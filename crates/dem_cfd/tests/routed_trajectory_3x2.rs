//! Live decomposed SOIL/FIELD trajectory over the routed 3-DEM/2-CFD seam.

#![cfg(feature = "mpi-routing")]

use std::any::TypeId;
use std::collections::{BTreeMap, BTreeSet};
use std::process::Command;

use cfd_boundary::{BoundaryPlugin, BoundaryRegistry};
use cfd_ibm::coupling::{self, InterphaseForces, ParticleKinematics, ParticleSet};
use cfd_solver::{FluxPlugin, IntegratorPlugin, SolverConfig, SolverPlugin};
use cfd_state::CfdState;
use dem_cfd::config::GasCfg;
use dem_cfd::routing::{decode_particles, reduce_forces, route_forces, route_particles};
use dem_cfd::routing::{RoutedForce, RoutedParticle};
use dem_cfd::seam::{build_cfd_base, build_soil_bed, point_particle_exchange};
use dem_cfd::seam::{FluidForces, SeamCtx};
use field_core::{FieldRegistry, FvMesh, MeshScheduleSet, PartitionDirectory};
use field_core::{UniformMesh, UniformMeshConfig};
use grass_app::App;
use grass_mpi::{CommBackend, MpiCommBackend};
use grass_multi::{CoupledPairRunner, CouplingEpoch, RoleLaunch};
use soil_core::Atom;

const CHILD_ENV: &str = "DEM_CFD_ROUTED_TRAJECTORY_CHILD";
const STEPS: u64 = 12;
const DT: f64 = 2.0e-5;
const RADIUS: f64 = 0.01;
const DENSITY: f64 = 2_500.0;
const CONFIG: &str = r#"
    [topology]
    mode = "split"
    [[topology.role]]
    name = "dem"
    ranks = 3
    [[topology.role]]
    name = "cfd"
    ranks = 2
"#;

fn mesh_config() -> UniformMeshConfig {
    UniformMeshConfig {
        nx: 10,
        ny: 2,
        nz: 2,
        ng: 1,
        bounds_lo: [0.0; 3],
        bounds_hi: [1.0; 3],
        y_edges: None,
        z_edges: None,
    }
}

fn directory() -> PartitionDirectory {
    PartitionDirectory::from_uniform_config(&mesh_config(), [2, 1, 1])
}

fn local_initial(rank: i32) -> Vec<(u64, [f64; 3], [f64; 3])> {
    match rank {
        0 => vec![(100, [0.25, 0.5, 0.5], [0.12, 0.0, 0.0])],
        // One DEM partition genuinely communicates with both FIELD owners.
        1 => vec![
            (110, [0.40, 0.5, 0.5], [0.08, 0.0, 0.0]),
            (111, [0.60, 0.5, 0.5], [0.06, 0.0, 0.0]),
        ],
        2 => vec![(120, [0.75, 0.5, 0.5], [0.10, 0.0, 0.0])],
        _ => unreachable!(),
    }
}

fn borrow_mut<T: 'static>(app: &App, f: impl FnOnce(&mut T)) {
    let cell = app.resource_cell(TypeId::of::<T>()).expect("resource");
    f(cell
        .borrow_mut()
        .downcast_mut::<T>()
        .expect("resource type"));
}

fn build_dem(rank: i32) -> App {
    let initial = local_initial(rank);
    let positions: Vec<_> = initial.iter().map(|(_, x, _)| *x).collect();
    let mut app = build_soil_bed(&positions, RADIUS, DENSITY, 0.0, DT);
    borrow_mut::<Atom>(&app, |atoms| {
        for (i, (id, _, velocity)) in initial.iter().enumerate() {
            atoms.tag[i] = *id as u32;
            atoms.vel[i] = velocity.map(|x| x as _);
        }
    });
    app.prepare();
    app
}

fn build_cfd() -> App {
    let gas = GasCfg {
        rho: 1.2,
        p: 101_325.0,
        mu: 1.8e-5,
    };
    let mut app = build_cfd_base(
        &gas,
        mesh_config(),
        SeamCtx {
            mu: gas.mu,
            rho: gas.rho,
            eps: 1.0,
            g: [0.0; 3],
            dt: DT,
            mode: Default::default(),
        },
    );
    app.add_plugins(BoundaryPlugin::<UniformMesh>::new(
        BoundaryRegistry::default(),
    ))
    .add_plugins(FluxPlugin::<UniformMesh>::hllc())
    .add_plugins(IntegratorPlugin::euler())
    .add_plugins(SolverPlugin::<UniformMesh>::new(SolverConfig {
        fixed_dt: Some(DT),
        ..SolverConfig::default()
    }));
    app.add_update_system(point_particle_exchange, MeshScheduleSet::Output);
    app.prepare();
    app
}

fn exported_particles(app: &App, dem_owner: i32) -> Vec<RoutedParticle> {
    let atoms = app.get_resource_ref::<Atom>().expect("DEM Atom");
    (0..atoms.nlocal as usize)
        .map(|i| RoutedParticle {
            id: atoms.tag[i] as u64,
            dem_owner,
            center: atoms.pos[i].map(|x| x as f64),
            velocity: atoms.vel[i].map(|x| x as f64),
            radius: RADIUS,
        })
        .collect()
}

fn set_particles(app: &App, particles: &[RoutedParticle]) {
    borrow_mut::<ParticleSet>(app, |set| {
        set.particles = particles
            .iter()
            .map(|p| ParticleKinematics {
                center: p.center,
                velocity: p.velocity,
                radius: p.radius,
            })
            .collect();
    });
}

fn fluid_forces(app: &App, particles: &[RoutedParticle]) -> Vec<RoutedForce> {
    let forces = app
        .get_resource_ref::<InterphaseForces>()
        .expect("FIELD forces");
    assert_eq!(forces.force.len(), particles.len());
    particles
        .iter()
        .zip(&forces.force)
        .map(|(particle, force)| RoutedForce {
            id: particle.id,
            dem_owner: particle.dem_owner,
            force: *force,
        })
        .collect()
}

fn apply_forces(app: &App, incoming: &[grass_multi::ReceivedPayload], prime_verlet: bool) -> f64 {
    let forces = reduce_forces(incoming).expect("reduce returned FIELD forces");
    let atoms = app.get_resource_ref::<Atom>().expect("DEM Atom");
    assert_eq!(
        forces.len(),
        atoms.nlocal as usize,
        "one force per local particle"
    );
    let expected: BTreeSet<_> = (0..atoms.nlocal as usize)
        .map(|i| atoms.tag[i] as u64)
        .collect();
    let actual: BTreeSet<_> = forces.iter().map(|force| force.id).collect();
    assert_eq!(actual, expected, "exact stable-ID force return set");
    let ordered: Vec<_> = (0..atoms.nlocal as usize)
        .map(|i| {
            forces
                .iter()
                .find(|force| force.id == atoms.tag[i] as u64)
                .expect("force returned to current DEM owner")
                .force
        })
        .collect();
    drop(atoms);
    borrow_mut::<FluidForces>(app, |fluid| {
        fluid.f.clone_from(&ordered);
    });
    if prime_verlet {
        borrow_mut::<Atom>(app, |atoms| {
            for (atom_force, fluid_force) in atoms.force.iter_mut().zip(&ordered) {
                *atom_force = fluid_force.map(|value| value as _);
            }
        });
    }
    ordered
        .iter()
        .flat_map(|force| force.iter())
        .map(|value| value.abs())
        .sum()
}

/// Correct the already-deposited rectangle-rule reaction to the same temporal
/// impulse used by velocity Verlet: F0 on the primed first step, then
/// 0.5*(F_previous + F_current). This makes the cross-role conservation check a
/// property of the actual coupled discretization, not a tolerance accident.
fn match_verlet_reaction(
    app: &App,
    particles: &[RoutedParticle],
    forces: &[RoutedForce],
    previous: &mut BTreeMap<u64, [f64; 3]>,
) -> [f64; 3] {
    let mut particle_impulse = [0.0; 3];
    let corrections: Vec<[f64; 3]> = forces
        .iter()
        .map(|force| {
            let old = previous.get(&force.id).copied().unwrap_or(force.force);
            for axis in 0..3 {
                particle_impulse[axis] += 0.5 * (old[axis] + force.force[axis]) * DT;
            }
            std::array::from_fn(|axis| 0.5 * (old[axis] - force.force[axis]))
        })
        .collect();
    let kinematics: Vec<_> = particles
        .iter()
        .map(|p| ParticleKinematics {
            center: p.center,
            velocity: p.velocity,
            radius: p.radius,
        })
        .collect();
    let mesh = app.get_resource_ref::<UniformMesh>().expect("FIELD mesh");
    let registry = app
        .get_resource_ref::<FieldRegistry>()
        .expect("FIELD registry");
    let mut state = registry.expect_mut::<CfdState>("CfdState");
    coupling::apply_momentum_sink(&*mesh, &mut state, &kinematics, &corrections, DT);
    for force in forces {
        previous.insert(force.id, force.force);
    }
    particle_impulse
}

fn particle_momentum(app: &App) -> [f64; 3] {
    let atoms = app.get_resource_ref::<Atom>().expect("DEM Atom");
    let mut momentum = [0.0; 3];
    for i in 0..atoms.nlocal as usize {
        for (axis, total) in momentum.iter_mut().enumerate() {
            *total += atoms.mass[i] as f64 * atoms.vel[i][axis] as f64;
        }
    }
    momentum
}

fn gas_momentum(app: &App) -> [f64; 3] {
    let mesh = app.get_resource_ref::<UniformMesh>().expect("FIELD mesh");
    let registry = app
        .get_resource_ref::<FieldRegistry>()
        .expect("FIELD registry");
    let state = registry.expect::<CfdState>("CfdState");
    let mut momentum = [0.0; 3];
    for cell in 0..mesh.n_cells_total() {
        if mesh.is_local_cell(cell) {
            let volume = mesh.cell_volume(cell);
            momentum[0] += state.u[cell].rho_u * volume;
            momentum[1] += state.u[cell].rho_v * volume;
            momentum[2] += state.u[cell].rho_w * volume;
        }
    }
    momentum
}

fn global_vector(comm: &dyn CommBackend, local: [f64; 3]) -> [f64; 3] {
    local.map(|value| comm.all_reduce_sum_f64(value))
}

fn run_role(launch: RoleLaunch) {
    let role = launch.role().to_owned();
    let exchange = launch.into_routed_exchange();
    let (rank, size) = exchange.role_position();
    let solver_comm = MpiCommBackend::new(grass_mpi::get_mpi_world());
    let mut app = if role == "dem" {
        assert_eq!(size, 3);
        build_dem(rank)
    } else {
        assert_eq!(role, "cfd");
        assert_eq!(size, 2);
        build_cfd()
    };
    assert_eq!(solver_comm.rank(), rank);
    assert_eq!(solver_comm.size(), size);

    let initial_particle_momentum =
        (role == "dem").then(|| global_vector(&solver_comm, particle_momentum(&app)));
    let initial_gas_momentum =
        (role == "cfd").then(|| global_vector(&solver_comm, gas_momentum(&app)));
    let initial_dem_velocity = if role == "dem" {
        exported_particles(&app, rank)
            .into_iter()
            .map(|particle| (particle.id, particle.velocity))
            .collect::<BTreeMap<_, _>>()
    } else {
        BTreeMap::new()
    };
    let mut previous_cfd_force = BTreeMap::new();
    let mut cfd_particle_impulse = [0.0; 3];
    let mut exchanged_force_l1 = 0.0;

    for step in 0..STEPS {
        let outgoing = if role == "dem" {
            route_particles(&directory(), &exported_particles(&app, rank))
                .expect("route live SOIL particles")
        } else {
            Vec::new()
        };
        let incoming = exchange
            .exchange(CouplingEpoch(2 * step), &outgoing)
            .expect("route SOIL state to FIELD owners");

        let outgoing = if role == "cfd" {
            let particles = decode_particles(&incoming).expect("decode live SOIL particles");
            assert!(particles
                .iter()
                .all(|p| directory().owner_rank(p.center) == Some(rank)));
            set_particles(&app, &particles);
            app.run();
            let forces = fluid_forces(&app, &particles);
            exchanged_force_l1 += forces
                .iter()
                .flat_map(|force| force.force)
                .map(f64::abs)
                .sum::<f64>();
            let impulse = match_verlet_reaction(&app, &particles, &forces, &mut previous_cfd_force);
            for axis in 0..3 {
                cfd_particle_impulse[axis] += impulse[axis];
            }
            route_forces(&forces)
        } else {
            assert!(incoming.is_empty());
            Vec::new()
        };
        let incoming = exchange
            .exchange(CouplingEpoch(2 * step + 1), &outgoing)
            .expect("return live FIELD forces to SOIL owners");
        if role == "dem" {
            exchanged_force_l1 += apply_forces(&app, &incoming, step == 0);
            app.run();
            for particle in exported_particles(&app, rank) {
                assert!(particle.center.into_iter().all(f64::is_finite));
                assert!(particle.velocity.into_iter().all(f64::is_finite));
            }
        } else {
            assert!(incoming.is_empty());
        }
    }

    let delta: [f64; 3] = if role == "dem" {
        let final_momentum = global_vector(&solver_comm, particle_momentum(&app));
        let initial = initial_particle_momentum.expect("DEM initial momentum");
        std::array::from_fn(|axis| final_momentum[axis] - initial[axis])
    } else {
        let final_gas = global_vector(&solver_comm, gas_momentum(&app));
        let initial_gas = initial_gas_momentum.expect("CFD initial momentum");
        assert!(
            (0..3).any(|axis| final_gas[axis] != initial_gas[axis]),
            "the CFD solver trajectory must evolve its conserved state"
        );
        // The Euler domain momentum also contains boundary fluxes. The exact
        // cross-role invariant is the equal-and-opposite coupling source.
        global_vector(&solver_comm, cfd_particle_impulse).map(|impulse| -impulse)
    };
    let global_force_l1 = solver_comm.all_reduce_sum_f64(exchanged_force_l1);
    assert!(
        global_force_l1 > 0.0,
        "coupling must exchange a nonzero force"
    );
    if role == "dem" {
        let changed = exported_particles(&app, rank)
            .into_iter()
            .any(|particle| particle.velocity != initial_dem_velocity[&particle.id]);
        assert!(
            solver_comm.all_reduce_sum_f64(f64::from(changed)) > 0.0,
            "at least one SOIL particle must evolve"
        );
    }
    // Cross-role diagnostic exchange: every peer root sees both role-global
    // momentum changes without leaking solver collectives across communicators.
    let diagnostic = if rank == 0 {
        let mut bytes = Vec::with_capacity(24);
        for value in delta {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        vec![grass_multi::RoutedPayload::new(
            0,
            grass_multi::EntityId(9_999),
            bytes,
        )]
    } else {
        Vec::new()
    };
    let peer = exchange
        .exchange(CouplingEpoch(2 * STEPS), &diagnostic)
        .expect("exchange role-global momentum diagnostics");
    if rank == 0 {
        assert_eq!(peer.len(), 1);
        let peer_delta: [f64; 3] = std::array::from_fn(|axis| {
            f64::from_le_bytes(peer[0].payload[8 * axis..8 * axis + 8].try_into().unwrap())
        });
        for axis in 0..3 {
            let residual = delta[axis] + peer_delta[axis];
            let scale = delta[axis].abs().max(peer_delta[axis].abs());
            assert!(
                scale > 0.0 || axis != 0,
                "x momentum exchange must be nonzero"
            );
            assert!(
                residual.abs() <= 1e-9 * scale.max(f64::MIN_POSITIVE),
                "cross-role momentum residual on axis {axis}: {residual:e} (scale {scale:e})"
            );
        }
    } else {
        assert!(peer.is_empty());
    }
    app.run_cleanup();
}

#[test]
fn routed_3x2_real_soil_field_trajectory_conserves_momentum() {
    if std::env::var_os(CHILD_ENV).is_some() {
        CoupledPairRunner::from_source(CONFIG)
            .and_then(|runner| runner.run(run_role, run_role))
            .expect("run real 3x2 SOIL-FIELD trajectory");
        return;
    }
    if Command::new("mpirun").arg("--version").output().is_err() {
        eprintln!("SKIP real 3x2 SOIL-FIELD trajectory: `mpirun` not found");
        return;
    }
    let executable = std::env::current_exe().expect("locate routed trajectory test binary");
    let status = Command::new("mpirun")
        .args(["--oversubscribe", "-np", "5"])
        .arg(executable)
        .args([
            "--exact",
            "routed_3x2_real_soil_field_trajectory_conserves_momentum",
            "--nocapture",
            "--test-threads=1",
        ])
        .env(CHILD_ENV, "1")
        .env("OMPI_MCA_btl", "self,vader")
        .status()
        .expect("spawn real 3x2 SOIL-FIELD trajectory");
    assert!(status.success(), "real 3x2 trajectory failed: {status}");
}
