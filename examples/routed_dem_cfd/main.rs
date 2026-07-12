//! The same binary and TOML run a local DEM-CFD route or a 3+2 MPI split.

use dem_cfd::routing::{
    decode_particles, reduce_forces, route_forces, route_particles, RoutedForce, RoutedParticle,
};
use field_core::{PartitionDirectory, UniformMeshConfig};
use grass_multi::{CoupledPairRunner, CouplingEpoch, PairRun, RoleLaunch};

const DEFAULT_CONFIG: &str = include_str!("config.toml");

fn directory(partitions: i32) -> PartitionDirectory {
    PartitionDirectory::from_uniform_config(
        &UniformMeshConfig {
            nx: 12,
            ny: 2,
            nz: 2,
            ng: 1,
            bounds_lo: [0.0; 3],
            bounds_hi: [1.0; 3],
            y_edges: None,
            z_edges: None,
        },
        [partitions, 1, 1],
    )
}

fn particle(rank: i32, size: i32) -> RoutedParticle {
    RoutedParticle {
        id: 100 + rank as u64,
        dem_owner: rank,
        center: [(rank as f64 + 0.5) / size as f64, 0.5, 0.5],
        velocity: [0.01 * (rank + 1) as f64, 0.0, 0.0],
        radius: 0.01,
        prev_force: [0.0; 3],
    }
}

fn drag(p: RoutedParticle) -> [f64; 3] {
    let coefficient = 2.5;
    let fluid_velocity = [p.center[0], -0.25, 0.125];
    std::array::from_fn(|axis| coefficient * (fluid_velocity[axis] - p.velocity[axis]))
}

fn run_role(launch: RoleLaunch) -> usize {
    let role = launch.role().to_owned();
    let exchange = launch.into_routed_exchange();
    let (rank, role_size) = exchange.role_position();
    let local_particles = if role == "dem" {
        vec![particle(rank, role_size)]
    } else {
        Vec::new()
    };
    let particle_routes = if role == "dem" {
        route_particles(&directory(exchange.peer_size()), &local_particles)
            .expect("route particles to FIELD owners")
    } else {
        Vec::new()
    };
    let incoming = exchange
        .exchange(CouplingEpoch(0), &particle_routes)
        .expect("exchange particle records");
    let force_routes = if role == "cfd" {
        let owned = decode_particles(&incoming).expect("decode particles");
        let owners = directory(role_size);
        let forces = owned
            .into_iter()
            .map(|p| {
                assert_eq!(owners.owner_rank(p.center), Some(rank));
                RoutedForce {
                    id: p.id,
                    dem_owner: p.dem_owner,
                    force: drag(p),
                }
            })
            .collect::<Vec<_>>();
        route_forces(&forces)
    } else {
        Vec::new()
    };
    let returned = exchange
        .exchange(CouplingEpoch(1), &force_routes)
        .expect("return stable-ID forces");
    if role == "dem" {
        let forces = reduce_forces(&returned).expect("reduce returned forces");
        assert_eq!(forces.len(), local_particles.len());
        assert_eq!(forces[0].force, drag(local_particles[0]));
        forces.len()
    } else {
        incoming.len()
    }
}

fn main() {
    let run = CoupledPairRunner::from_cli_or(DEFAULT_CONFIG)
        .and_then(|runner| runner.run(run_role, run_role))
        .unwrap_or_else(|error| panic!("run routed DEM-CFD example: {error}"));
    match run {
        PairRun::Local { first, second } => {
            println!("LOCAL dem_records={first} cfd_records={second}");
            println!("PASS same-binary local routed coupling");
        }
        PairRun::Split { role, result } => {
            println!("MPI role={role} local_records={result}");
        }
    }
}
