//! Live decomposed SOIL/FIELD trajectory over the routed 3-DEM/2-CFD seam,
//! with particles that genuinely cross both CFD partition boundaries and DEM
//! ownership boundaries mid-trajectory.
//!
//! The mesh is split `[2, 1, 1]`, so the CFD partition boundary sits at
//! `x = 0.5`. The coupling package additionally owns a 3-way DEM ownership
//! decomposition in `x` (thresholds [`DEM_SPLIT`]); a particle that drifts past
//! one of those thresholds is physically migrated to the new DEM owner over the
//! DEM *role* communicator (real point-to-point MPI, isolated from both the
//! coupling exchange and the CFD role), carrying its stable ID and full
//! kinematic + force state.
//!
//! The initial conditions launch four grains so that, over the trajectory:
//!   * id 110 crosses the CFD partition boundary (owner 0 → 1) *and* later a DEM
//!     ownership boundary (rank 1 → 2) — a particle that crosses both seams;
//!   * id 111 crosses a DEM ownership boundary only (rank 1 → 2);
//!   * id 100 crosses a DEM ownership boundary only (rank 0 → 1);
//!   * id 120 crosses nothing (a stationary control).
//!
//! Everything the original conservation gate checked still holds exactly:
//! stable-ID force return, the trapezoidal (velocity-Verlet) temporal impulse
//! matched equal-and-opposite across roles to 1e-9, a nonzero exchanged force,
//! and communicator isolation. The crossing correctness rests on carrying the
//! previous fluid force with the particle ([`RoutedParticle::prev_force`]) so the
//! receiving CFD owner can form the trapezoid on the first step it sees a
//! migrated particle instead of degrading to a rectangle rule.

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
use soil_core::{Accum, Atom, Real};

const CHILD_ENV: &str = "DEM_CFD_ROUTED_TRAJECTORY_CHILD";
const STEPS: u64 = 320;
const DT: f64 = 2.0e-5;
const RADIUS: f64 = 0.01;
const DENSITY: f64 = 2_500.0;

/// DEM ownership thresholds in `x`: rank 0 owns `x < 0.35`, rank 1 owns
/// `0.35 ≤ x < 0.55`, rank 2 owns `x ≥ 0.55`. Deliberately independent of the
/// FIELD `[2,1,1]` CFD partition boundary at `x = 0.5` so a particle can cross a
/// CFD partition boundary and a DEM ownership boundary at different times.
const DEM_SPLIT: [f64; 2] = [0.35, 0.55];

/// Flat migration record: tag, x, y, z, vx, vy, vz, mass, radius, fx, fy, fz.
const MIG_RECORD: usize = 12;

/// Stable IDs seeded across the DEM role (used for the migration conservation
/// checks).
const INITIAL_IDS: [u64; 4] = [100, 110, 111, 120];

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

/// DEM role-rank that owns a particle at `center`, from its world `x`.
fn dem_owner_of(center: [f64; 3]) -> i32 {
    let x = center[0];
    if x < DEM_SPLIT[0] {
        0
    } else if x < DEM_SPLIT[1] {
        1
    } else {
        2
    }
}

fn local_initial(rank: i32) -> Vec<(u64, [f64; 3], [f64; 3])> {
    match rank {
        // Crosses a DEM ownership boundary (rank 0 → 1) but stays CFD owner 0.
        0 => vec![(100, [0.30, 0.5, 0.5], [20.0, 0.0, 0.0])],
        1 => vec![
            // Crosses BOTH the CFD partition boundary (owner 0 → 1) and, later,
            // a DEM ownership boundary (rank 1 → 2).
            (110, [0.40, 0.5, 0.5], [30.0, 0.0, 0.0]),
            // Crosses a DEM ownership boundary only (rank 1 → 2); already CFD
            // owner 1.
            (111, [0.52, 0.5, 0.5], [30.0, 0.0, 0.0]),
        ],
        // Stationary control: crosses nothing.
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

/// Export the live local bed as routed particles. `prev_force` carries the fluid
/// force this atom left the previous coupling step with (`atoms.force`), so the
/// trapezoidal reaction survives a CFD-owner or DEM-owner change. Zero before the
/// first primed step.
fn exported_particles(app: &App, dem_owner: i32) -> Vec<RoutedParticle> {
    let atoms = app.get_resource_ref::<Atom>().expect("DEM Atom");
    (0..atoms.nlocal as usize)
        .map(|i| RoutedParticle {
            id: atoms.tag[i] as u64,
            dem_owner,
            center: atoms.pos[i].map(|x| x as f64),
            velocity: atoms.vel[i].map(|x| x as f64),
            radius: RADIUS,
            prev_force: atoms.force[i].map(|x| x as f64),
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

/// Deposit the same temporal impulse used by velocity Verlet: `F0` on the primed
/// first step, then `0.5·(F_prev + F_now)` afterwards. `F_prev` is taken from the
/// particle's carried [`RoutedParticle::prev_force`], which travels with the
/// particle across CFD-partition and DEM-ownership boundaries, so the cross-role
/// conservation check is a property of the discretization even when the CFD owner
/// changes mid-trajectory (a fresh owner has no local force history).
fn match_verlet_reaction(
    app: &App,
    particles: &[RoutedParticle],
    forces: &[RoutedForce],
    first_step: bool,
) -> [f64; 3] {
    let mut particle_impulse = [0.0; 3];
    let corrections: Vec<[f64; 3]> = forces
        .iter()
        .zip(particles)
        .map(|(force, particle)| {
            let old = if first_step {
                force.force
            } else {
                particle.prev_force
            };
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

/// All-gather variable-length `f64` buffers around the role ring. Each rank sends
/// its own `local` to every peer across the rounds and collects each peer's
/// buffer; deadlock-free (immediate send + probing receive). Used for global
/// initial-state maps and the final stable-ID set.
fn ring_all_gather(comm: &MpiCommBackend, local: &[f64]) -> Vec<f64> {
    let rank = comm.rank() as usize;
    let size = comm.size() as usize;
    let mut collected: Vec<Vec<f64>> = vec![Vec::new(); size];
    collected[rank] = local.to_vec();
    for dist in 1..size {
        let dst = (rank + dist) % size;
        let src = (rank + size - dist) % size;
        collected[src] = comm.sendrecv_f64(dst as i32, local, src as i32);
    }
    collected.concat()
}

/// Migrate every local atom that has drifted out of this rank's DEM ownership
/// slab to its new owner over the DEM role communicator, and receive any atoms
/// migrating in. Real point-to-point MPI on the role communicator — isolated
/// from the coupling exchange (a duplicated communicator) and from the CFD role.
/// Stable ID and full state (position, velocity, mass, radius, and the fluid
/// force carried for the next trapezoid) travel with the atom. Returns the number
/// of atoms this rank sent away.
fn migrate_dem_atoms(app: &App, comm: &MpiCommBackend) -> u64 {
    let rank = comm.rank() as usize;
    let size = comm.size() as usize;

    let mut keep: Vec<f64> = Vec::new();
    let mut send: Vec<Vec<f64>> = vec![Vec::new(); size];
    {
        let atoms = app.get_resource_ref::<Atom>().expect("DEM Atom");
        for i in 0..atoms.nlocal as usize {
            let center = [
                atoms.pos[i][0] as f64,
                atoms.pos[i][1] as f64,
                atoms.pos[i][2] as f64,
            ];
            let record = [
                atoms.tag[i] as f64,
                center[0],
                center[1],
                center[2],
                atoms.vel[i][0] as f64,
                atoms.vel[i][1] as f64,
                atoms.vel[i][2] as f64,
                atoms.mass[i] as f64,
                atoms.cutoff_radius[i] as f64,
                atoms.force[i][0] as f64,
                atoms.force[i][1] as f64,
                atoms.force[i][2] as f64,
            ];
            let owner = dem_owner_of(center) as usize;
            if owner == rank {
                keep.extend_from_slice(&record);
            } else {
                send[owner].extend_from_slice(&record);
            }
        }
    }
    let migrated_out = (send.iter().map(Vec::len).sum::<usize>() / MIG_RECORD) as u64;

    let mut received: Vec<f64> = Vec::new();
    for dist in 1..size {
        let dst = (rank + dist) % size;
        let src = (rank + size - dist) % size;
        let got = comm.sendrecv_f64(dst as i32, &send[dst], src as i32);
        received.extend_from_slice(&got);
    }

    borrow_mut::<Atom>(app, |atoms| rebuild_local_atoms(atoms, &keep, &received));
    migrated_out
}

/// Rebuild the local `Atom` store from the kept and freshly received flat
/// records. The bed uses only core per-atom fields (velocity Verlet + the seam
/// force), so no `AtomDataRegistry` mirroring is required.
fn rebuild_local_atoms(atoms: &mut Atom, keep: &[f64], received: &[f64]) {
    atoms.tag.clear();
    atoms.atom_type.clear();
    atoms.origin_index.clear();
    atoms.pos.clear();
    atoms.vel.clear();
    atoms.force.clear();
    atoms.cutoff_radius.clear();
    atoms.mass.clear();
    atoms.inv_mass.clear();
    atoms.image.clear();
    atoms.is_ghost.clear();

    for record in keep
        .chunks_exact(MIG_RECORD)
        .chain(received.chunks_exact(MIG_RECORD))
    {
        let mass = record[7];
        atoms.tag.push(record[0] as u32);
        atoms.atom_type.push(0);
        atoms.origin_index.push(0);
        atoms
            .pos
            .push([record[1] as Real, record[2] as Real, record[3] as Real]);
        atoms
            .vel
            .push([record[4] as Real, record[5] as Real, record[6] as Real]);
        atoms
            .force
            .push([record[9] as Accum, record[10] as Accum, record[11] as Accum]);
        atoms.cutoff_radius.push(record[8] as Real);
        atoms.mass.push(mass as Real);
        atoms.inv_mass.push((1.0 / mass) as Real);
        atoms.image.push([0, 0, 0]);
        atoms.is_ghost.push(false);
    }

    let n = atoms.tag.len() as u32;
    atoms.nlocal = n;
    atoms.natoms = n as u64;
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

    // Initial DEM ownership must be self-consistent, and the global initial
    // velocity map lets the evolution check survive migration (a rank may end
    // holding an ID it did not start with).
    let global_initial_vel: BTreeMap<u64, [f64; 3]> = if role == "dem" {
        for particle in exported_particles(&app, rank) {
            assert_eq!(
                dem_owner_of(particle.center),
                rank,
                "initial DEM ownership of {} must match its seeding rank",
                particle.id
            );
        }
        let mut local = Vec::new();
        for particle in exported_particles(&app, rank) {
            local.extend_from_slice(&[
                particle.id as f64,
                particle.velocity[0],
                particle.velocity[1],
                particle.velocity[2],
            ]);
        }
        ring_all_gather(&solver_comm, &local)
            .chunks_exact(4)
            .map(|c| (c[0] as u64, [c[1], c[2], c[3]]))
            .collect()
    } else {
        BTreeMap::new()
    };

    let initial_particle_momentum =
        (role == "dem").then(|| global_vector(&solver_comm, particle_momentum(&app)));
    let initial_gas_momentum =
        (role == "cfd").then(|| global_vector(&solver_comm, gas_momentum(&app)));
    let mut cfd_particle_impulse = [0.0; 3];
    let mut exchanged_force_l1 = 0.0;

    // DEM-side observation of live routing: per-ID last CFD owner seen while the
    // particle was local here, and counters for the two kinds of crossing.
    let mut cfd_owner_seen: BTreeMap<u64, i32> = BTreeMap::new();
    let mut cfd_crossings = 0.0;
    let mut dem_migrations = 0.0;

    for step in 0..STEPS {
        let exported = if role == "dem" {
            exported_particles(&app, rank)
        } else {
            Vec::new()
        };
        if role == "dem" {
            for particle in &exported {
                let cfd_owner = directory()
                    .owner_rank(particle.center)
                    .expect("live particle left the FIELD domain");
                if let Some(previous) = cfd_owner_seen.insert(particle.id, cfd_owner) {
                    if previous != cfd_owner {
                        cfd_crossings += 1.0;
                    }
                }
            }
        }
        let outgoing = if role == "dem" {
            route_particles(&directory(), &exported).expect("route live SOIL particles")
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
            let impulse = match_verlet_reaction(&app, &particles, &forces, step == 0);
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
            // Migrate any atom that drifted out of this rank's ownership slab.
            dem_migrations += migrate_dem_atoms(&app, &solver_comm) as f64;
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
        // Both kinds of boundary crossing must actually have happened.
        let global_cfd_crossings = solver_comm.all_reduce_sum_f64(cfd_crossings);
        let global_migrations = solver_comm.all_reduce_sum_f64(dem_migrations);
        assert!(
            global_cfd_crossings >= 1.0,
            "at least one particle must cross the CFD partition boundary (observed {global_cfd_crossings})"
        );
        assert!(
            global_migrations >= 1.0,
            "at least one particle must cross a DEM ownership boundary (observed {global_migrations})"
        );
        if rank == 0 {
            eprintln!(
                "[routed-trajectory] CFD-partition crossings observed: {global_cfd_crossings}; \
                 DEM ownership migrations: {global_migrations}"
            );
        }

        // Migration must conserve the global particle count and the stable-ID
        // set exactly.
        let local_count = app.get_resource_ref::<Atom>().expect("DEM Atom").nlocal as f64;
        let global_count = solver_comm.all_reduce_sum_f64(local_count) as u64;
        assert_eq!(
            global_count,
            INITIAL_IDS.len() as u64,
            "DEM migration must conserve the global particle count"
        );
        let local_ids: Vec<f64> = exported_particles(&app, rank)
            .iter()
            .map(|p| p.id as f64)
            .collect();
        let global_ids: BTreeSet<u64> = ring_all_gather(&solver_comm, &local_ids)
            .iter()
            .map(|&value| value as u64)
            .collect();
        assert_eq!(
            global_ids,
            BTreeSet::from(INITIAL_IDS),
            "stable IDs must be conserved across DEM migration"
        );

        // At least one particle's velocity must have evolved.
        let changed = exported_particles(&app, rank).into_iter().any(|particle| {
            global_initial_vel
                .get(&particle.id)
                .map_or(true, |initial| particle.velocity != *initial)
        });
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

    // Wrap the launch in `timeout` when available so a hung MPI rank cannot run
    // unbounded in CI; the trajectory itself completes in a couple of seconds.
    let use_timeout = Command::new("timeout")
        .arg("--version")
        .output()
        .map_or(false, |o| o.status.success());
    let mut command = if use_timeout {
        let mut command = Command::new("timeout");
        command.args(["--signal=KILL", "180", "mpirun"]);
        command
    } else {
        Command::new("mpirun")
    };
    let status = command
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
