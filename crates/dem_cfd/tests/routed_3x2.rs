//! Live 3-DEM/2-CFD manufactured conservative routing validation.

#![cfg(feature = "mpi-routing")]

use dem_cfd::routing::{
    decode_particles, reduce_forces, route_forces, route_particles, RoutedForce, RoutedParticle,
};
use field_core::{PartitionDirectory, UniformMeshConfig};
use grass_multi::{CoupledPairRunner, CouplingEpoch, RoleLaunch};
use std::process::Command;

const CHILD_ENV: &str = "DEM_CFD_ROUTED_3X2_CHILD";
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

fn directory() -> PartitionDirectory {
    PartitionDirectory::from_uniform_config(
        &UniformMeshConfig {
            nx: 10,
            ny: 2,
            nz: 2,
            ng: 1,
            bounds_lo: [0.0; 3],
            bounds_hi: [1.0; 3],
            y_edges: None,
            z_edges: None,
        },
        [2, 1, 1],
    )
}

fn particle(id: u64, dem_owner: i32, x: f64, vx: f64) -> RoutedParticle {
    RoutedParticle {
        id,
        dem_owner,
        center: [x, 0.5, 0.5],
        velocity: [vx, 0.0, 0.0],
        radius: 0.01,
    }
}

fn local_particles(rank: i32) -> Vec<RoutedParticle> {
    match rank {
        0 => vec![particle(100, 0, 0.10, 0.01)],
        // This DEM partition straddles the CFD split and must route both ways.
        1 => vec![particle(110, 1, 0.40, 0.02), particle(111, 1, 0.60, 0.03)],
        2 => vec![particle(120, 2, 0.90, 0.04)],
        _ => unreachable!(),
    }
}

fn manufactured_drag(particle: RoutedParticle) -> [f64; 3] {
    let volume = 4.0 / 3.0 * std::f64::consts::PI * particle.radius.powi(3);
    let gas_velocity = [particle.center[0], -0.25, 0.125];
    let coefficient = 2.5;
    [
        coefficient * volume * (gas_velocity[0] - particle.velocity[0]),
        coefficient * volume * (gas_velocity[1] - particle.velocity[1]),
        coefficient * volume * (gas_velocity[2] - particle.velocity[2]),
    ]
}

fn assert_force(actual: [f64; 3], expected: [f64; 3]) {
    for axis in 0..3 {
        assert!(
            (actual[axis] - expected[axis]).abs() < 1e-18,
            "force axis {axis}: {} != {}",
            actual[axis],
            expected[axis]
        );
    }
}

fn run_role(launch: RoleLaunch) {
    let role = launch.role().to_owned();
    let exchange = launch.into_routed_exchange();
    let (rank, size) = exchange.role_position();
    let particles = if role == "dem" {
        assert_eq!(size, 3);
        local_particles(rank)
    } else {
        assert_eq!(role, "cfd");
        assert_eq!(size, 2);
        Vec::new()
    };

    let outgoing_particles = if role == "dem" {
        route_particles(&directory(), &particles).expect("route DEM particles")
    } else {
        Vec::new()
    };
    let incoming_particles = exchange
        .exchange(CouplingEpoch(0), &outgoing_particles)
        .expect("exchange routed particles");

    let outgoing_forces = if role == "cfd" {
        let owned = decode_particles(&incoming_particles).expect("decode routed particles");
        for particle in &owned {
            assert_eq!(directory().owner_rank(particle.center), Some(rank));
        }
        let mut particle_force = [0.0; 3];
        let mut gas_reaction = [0.0; 3];
        let forces: Vec<RoutedForce> = owned
            .into_iter()
            .map(|particle| {
                let force = manufactured_drag(particle);
                for axis in 0..3 {
                    particle_force[axis] += force[axis];
                    gas_reaction[axis] -= force[axis];
                }
                RoutedForce {
                    id: particle.id,
                    dem_owner: particle.dem_owner,
                    force,
                }
            })
            .collect();
        for axis in 0..3 {
            assert_eq!(particle_force[axis] + gas_reaction[axis], 0.0);
        }
        route_forces(&forces)
    } else {
        assert!(incoming_particles.is_empty());
        Vec::new()
    };

    let incoming_forces = exchange
        .exchange(CouplingEpoch(1), &outgoing_forces)
        .expect("return routed forces");
    if role == "dem" {
        let returned = reduce_forces(&incoming_forces).expect("reduce returned forces");
        assert_eq!(returned.len(), particles.len());
        for (actual, particle) in returned.iter().zip(&particles) {
            assert_eq!(actual.id, particle.id);
            assert_eq!(actual.dem_owner, rank);
            assert_force(actual.force, manufactured_drag(*particle));
        }
    } else {
        assert!(incoming_forces.is_empty());
    }
}

#[test]
fn routed_3x2_particle_force_map_is_conservative_and_owner_correct() {
    if std::env::var_os(CHILD_ENV).is_some() {
        CoupledPairRunner::from_source(CONFIG)
            .and_then(|runner| runner.run(run_role, run_role))
            .expect("run 3x2 routed DEM-CFD validation");
        return;
    }
    if Command::new("mpirun").arg("--version").output().is_err() {
        eprintln!("SKIP 3x2 routed DEM-CFD validation: `mpirun` not found");
        return;
    }
    let executable = std::env::current_exe().expect("locate routed test binary");
    let status = Command::new("mpirun")
        .args(["--oversubscribe", "-np", "5"])
        .arg(executable)
        .args([
            "--exact",
            "routed_3x2_particle_force_map_is_conservative_and_owner_correct",
            "--nocapture",
            "--test-threads=1",
        ])
        .env(CHILD_ENV, "1")
        .env("OMPI_MCA_btl", "self,vader")
        .status()
        .expect("spawn 3x2 routed DEM-CFD validation");
    assert!(status.success(), "3x2 routed validation failed: {status}");
}
