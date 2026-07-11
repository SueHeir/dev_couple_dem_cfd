//! Bed-scale mesh plumbing shared by every unresolved DEM↔CFD case: containment
//! binning of the packing's solid volume onto the (coarse) gas mesh, FCC packing
//! construction, imposing the interstitial gas velocity, and the two-way momentum
//! sink + conservation check. The seam's own interpolation locator is tuned for a
//! single SUB-CELL particle and mis-bins at bed scale, so bed-scale containment
//! binning lives here rather than in the seam crate.

use cfd_ibm::coupling::{self, ParticleKinematics};
use cfd_state::CfdState;
use field_core::{FvMesh, StructuredMesh, UniformMesh, Vec3};

/// Interior cell-center coordinates along each axis (uniform, separable grid) plus
/// the ghost width.
pub fn axis_centers(mesh: &UniformMesh) -> ([Vec<f64>; 3], usize) {
    let [ni, nj, nk] = mesh.dims();
    let ng = mesh.n_ghost();
    let xc = (0..ni)
        .map(|i| mesh.cell_centroid(mesh.idx_raw(i + ng, ng, ng))[0])
        .collect();
    let yc = (0..nj)
        .map(|j| mesh.cell_centroid(mesh.idx_raw(ng, j + ng, ng))[1])
        .collect();
    let zc = (0..nk)
        .map(|k| mesh.cell_centroid(mesh.idx_raw(ng, ng, k + ng))[2])
        .collect();
    ([xc, yc, zc], ng)
}

/// Nearest cell-center index along one axis (uniform spacing).
#[inline]
pub fn nearest_center(cs: &[f64], v: f64) -> usize {
    if cs.len() < 2 {
        return 0;
    }
    let dx = cs[1] - cs[0];
    (((v - cs[0]) / dx).round() as isize).clamp(0, cs.len() as isize - 1) as usize
}

/// Raw cell index (with ghosts) of the cell containing point `p`.
pub fn containing_cell(mesh: &UniformMesh, centers: &[Vec<f64>; 3], ng: usize, p: Vec3) -> usize {
    let i = nearest_center(&centers[0], p[0]);
    let j = nearest_center(&centers[1], p[1]);
    let k = nearest_center(&centers[2], p[2]);
    mesh.idx_raw(i + ng, j + ng, k + ng)
}

/// Per-cell void fraction `ε = 1 − ΣV_p/V_cell` by containment deposition, plus the
/// per-particle containing-cell index (reused to drive the drag). In the unresolved
/// regime (gas cell ≫ particle) each cell holds many spheres, so this field
/// converges to the bed porosity and is the field that DRIVES the drag.
pub fn deposit_bed_void_fraction(
    mesh: &UniformMesh,
    particles: &[ParticleKinematics],
) -> (Vec<f64>, Vec<usize>) {
    let (centers, ng) = axis_centers(mesh);
    let total = mesh.n_cells_total();
    let mut solid = vec![0.0f64; total];
    let mut cell_of_particle = Vec::with_capacity(particles.len());
    for p in particles {
        let c = containing_cell(mesh, &centers, ng, p.center);
        solid[c] += p.volume();
        cell_of_particle.push(c);
    }
    let mut eps = vec![1.0f64; total];
    for c in 0..total {
        let v = mesh.cell_volume(c);
        if v > 0.0 {
            eps[c] = (1.0 - solid[c] / v).clamp(1e-6, 1.0);
        }
    }
    (eps, cell_of_particle)
}

/// Impose a uniform interstitial gas velocity `u_g` (world axes) in every interior
/// cell, leaving ρ untouched. `u_g = U/ε_bed` is the pore-scale velocity the drag
/// closure sees.
pub fn impose_interstitial_velocity(mesh: &UniformMesh, state: &mut CfdState, u_g: [f64; 3]) {
    for c in 0..mesh.n_cells_total() {
        if !mesh.is_local_cell(c) {
            continue;
        }
        let rho = state.u[c].rho;
        state.u[c].rho_u = rho * u_g[0];
        state.u[c].rho_v = rho * u_g[1];
        state.u[c].rho_w = rho * u_g[2];
    }
}

/// Apply the two-way momentum sink (reaction of the DRAG part of the interphase
/// force) to the gas and return the conservation error `‖Δp_gas − (−ΣF_drag·dt)‖ /
/// ‖ΣF_drag·dt‖`. `drag_on_particle` is the per-particle drag (world axes).
pub fn momentum_sink_and_check(
    mesh: &UniformMesh,
    state: &mut CfdState,
    particles: &[ParticleKinematics],
    drag_on_particle: &[[f64; 3]],
    dt: f64,
) -> f64 {
    let mut m0 = [0.0f64; 3];
    for c in 0..mesh.n_cells_total() {
        if mesh.is_local_cell(c) {
            let v = mesh.cell_volume(c);
            m0[0] += state.u[c].rho_u * v;
            m0[1] += state.u[c].rho_v * v;
            m0[2] += state.u[c].rho_w * v;
        }
    }
    coupling::apply_momentum_sink(mesh, state, particles, drag_on_particle, dt);
    let mut m1 = [0.0f64; 3];
    for c in 0..mesh.n_cells_total() {
        if mesh.is_local_cell(c) {
            let v = mesh.cell_volume(c);
            m1[0] += state.u[c].rho_u * v;
            m1[1] += state.u[c].rho_v * v;
            m1[2] += state.u[c].rho_w * v;
        }
    }
    let (mut dn, mut sc) = (0.0f64, 0.0f64);
    for k in 0..3 {
        let dm = m1[k] - m0[k];
        let imp = -drag_on_particle.iter().map(|f| f[k]).sum::<f64>() * dt;
        dn += (dm - imp) * (dm - imp);
        sc += imp * imp;
    }
    dn.sqrt() / sc.sqrt().max(1e-30)
}

/// FCC sphere centers filling `[0,Lx]×[0,Ly]×[0,Lz]`: `nc` conventional cells (4
/// spheres each), lattice constant `a`. Returns `(positions, bounds)`.
pub fn fcc_packing(nc: [usize; 3], a: f64) -> (Vec<[f64; 3]>, [f64; 3]) {
    let s = 0.25;
    let basis = [
        [s, s, s],
        [s, 0.5 + s, 0.5 + s],
        [0.5 + s, s, 0.5 + s],
        [0.5 + s, 0.5 + s, s],
    ];
    let mut pos = Vec::with_capacity(4 * nc[0] * nc[1] * nc[2]);
    for i in 0..nc[0] {
        for j in 0..nc[1] {
            for k in 0..nc[2] {
                for b in &basis {
                    pos.push([
                        (i as f64 + b[0]) * a,
                        (j as f64 + b[1]) * a,
                        (k as f64 + b[2]) * a,
                    ]);
                }
            }
        }
    }
    let bounds = [nc[0] as f64 * a, nc[1] as f64 * a, nc[2] as f64 * a];
    (pos, bounds)
}

/// FCC lattice constant that yields solid volume fraction `phi` for bead diameter
/// `d`: `a = d (2π / (3φ))^{1/3}`.
pub fn fcc_lattice_constant(d: f64, phi: f64) -> f64 {
    d * (2.0 * std::f64::consts::PI / (3.0 * phi)).cbrt()
}
