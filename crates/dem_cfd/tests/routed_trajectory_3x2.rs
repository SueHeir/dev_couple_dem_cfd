//! Live decomposed SOIL/FIELD trajectory over the routed 3-DEM/2-CFD seam.

#![cfg(feature = "mpi-routing")]

use std::any::TypeId;
use std::process::Command;

use cfd_ibm::coupling::{InterphaseForces, ParticleKinematics, ParticleSet};
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
        0 => vec![(100, [0.10, 0.5, 0.5], [0.12, 0.0, 0.0])],
        // One DEM partition genuinely communicates with both FIELD owners.
        1 => vec![
            (110, [0.40, 0.5, 0.5], [0.08, 0.0, 0.0]),
            (111, [0.60, 0.5, 0.5], [-0.08, 0.0, 0.0]),
        ],
        2 => vec![(120, [0.90, 0.5, 0.5], [-0.12, 0.0, 0.0])],
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

fn apply_forces(app: &App, incoming: &[grass_multi::ReceivedPayload]) {
    let forces = reduce_forces(incoming).expect("reduce returned FIELD forces");
    borrow_mut::<FluidForces>(app, |fluid| {
        let atoms = app.get_resource_ref::<Atom>().expect("DEM Atom");
        fluid.f = (0..atoms.nlocal as usize)
            .map(|i| {
                forces
                    .iter()
                    .find(|force| force.id == atoms.tag[i] as u64)
                    .expect("force returned to current DEM owner")
                    .force
            })
            .collect();
    });
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

    let initial = if role == "dem" {
        global_vector(&solver_comm, particle_momentum(&app))
    } else {
        global_vector(&solver_comm, gas_momentum(&app))
    };

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
            route_forces(&fluid_forces(&app, &particles))
        } else {
            assert!(incoming.is_empty());
            Vec::new()
        };
        let incoming = exchange
            .exchange(CouplingEpoch(2 * step + 1), &outgoing)
            .expect("return live FIELD forces to SOIL owners");
        if role == "dem" {
            apply_forces(&app, &incoming);
            app.run();
            for particle in exported_particles(&app, rank) {
                assert!(particle.center.into_iter().all(f64::is_finite));
                assert!(particle.velocity.into_iter().all(f64::is_finite));
            }
        } else {
            assert!(incoming.is_empty());
        }
    }

    let final_momentum = if role == "dem" {
        global_vector(&solver_comm, particle_momentum(&app))
    } else {
        global_vector(&solver_comm, gas_momentum(&app))
    };
    let delta: [f64; 3] = std::array::from_fn(|axis| final_momentum[axis] - initial[axis]);

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
            assert!(
                residual.abs() <= 5e-12 * delta[axis].abs().max(peer_delta[axis].abs()).max(1.0),
                "cross-role momentum residual on axis {axis}: {residual:e}"
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
