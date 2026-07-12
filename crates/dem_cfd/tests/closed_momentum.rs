//! Closed-domain check of the gas's realized equal-and-opposite impulse.

use cfd_ibm::coupling::{apply_momentum_sink, ParticleKinematics};
use cfd_state::{CfdState, ConsVar};
use field_core::{FvMesh, UniformMesh, UniformMeshConfig};

fn momentum(mesh: &UniformMesh, state: &CfdState) -> [f64; 3] {
    let mut total = [0.0; 3];
    for cell in 0..mesh.n_cells_total() {
        if mesh.is_local_cell(cell) {
            let volume = mesh.cell_volume(cell);
            total[0] += state.u[cell].rho_u * volume;
            total[1] += state.u[cell].rho_v * volume;
            total[2] += state.u[cell].rho_w * volume;
        }
    }
    total
}

#[test]
fn closed_gas_absorbs_the_realized_opposite_particle_impulse() {
    let mesh = UniformMesh::from_config(&UniformMeshConfig {
        nx: 4,
        ny: 3,
        nz: 2,
        ng: 1,
        ..Default::default()
    });
    let n = mesh.n_cells_total();
    let mut state = CfdState {
        u: vec![ConsVar::new(1.2, 0.0, 0.0, 0.0, 2.5); n],
        rhs: vec![ConsVar::default(); n],
        u0: vec![ConsVar::default(); n],
    };
    let particles = [
        ParticleKinematics { center: [0.2, 0.4, 0.4], ..Default::default() },
        ParticleKinematics { center: [0.8, 0.6, 0.6], ..Default::default() },
    ];
    let forces = [[2.0, -3.0, 0.5], [-0.25, 1.0, 4.0]];
    let dt = 0.0125;
    let before = momentum(&mesh, &state);
    apply_momentum_sink(&mesh, &mut state, &particles, &forces, dt);
    let after = momentum(&mesh, &state);

    for axis in 0..3 {
        let particle_impulse = forces.iter().map(|force| force[axis]).sum::<f64>() * dt;
        let realized_gas_impulse = after[axis] - before[axis];
        assert!(
            (realized_gas_impulse + particle_impulse).abs() < 1e-14,
            "axis {axis}: gas {realized_gas_impulse} particle {particle_impulse}"
        );
    }
}
