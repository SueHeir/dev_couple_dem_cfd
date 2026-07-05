//! Declarative case blocks every unresolved DEM↔CFD example deserializes from its
//! TOML. Kept here so the examples don't each redeclare an identical `[gas]`,
//! `[particle]`, `[mesh]`, `[packing]` struct.

use field_core::UniformMeshConfig;
use serde::Deserialize;

/// `[gas]` — carrier-phase state and transport.
#[derive(Deserialize, Default, Clone, Copy)]
pub struct GasCfg {
    pub rho: f64,
    pub p: f64,
    /// Dynamic viscosity μ [Pa·s].
    pub mu: f64,
}

/// `[particle]` — the DEM grain the drag closure sees.
#[derive(Deserialize, Default, Clone, Copy)]
pub struct ParticleCfg {
    /// Bead diameter d [m].
    pub diameter: f64,
    pub density: f64,
}

/// `[gravity]`.
#[derive(Deserialize, Default, Clone, Copy)]
pub struct GravityCfg {
    pub gz: f64,
}

/// `[packing]` — an FCC bed (4 spheres per conventional cell). The lattice
/// constant is chosen elsewhere so the solid fraction equals `solid_fraction`
/// (bed porosity ε = 1 − solid_fraction); `solid_fraction` must be ≤ 0.7405 (FCC max).
#[derive(Deserialize, Default, Clone, Copy)]
pub struct PackingCfg {
    pub ncx: usize,
    pub ncy: usize,
    pub ncz: usize,
    /// Target solid volume fraction φ = 1 − ε.
    pub solid_fraction: f64,
}

/// `[mesh]` — the gas mesh, **coarser than the particles** in the unresolved
/// regime. `nx*` must divide the packing's conventional-cell counts so the domain
/// tiles evenly.
#[derive(Deserialize, Default, Clone, Copy)]
pub struct MeshCfg {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    #[serde(default = "default_ng")]
    pub ng: usize,
}

/// Default ghost-cell width for a gas mesh.
pub fn default_ng() -> usize {
    2
}

impl MeshCfg {
    /// Build a `[0, bounds_hi]` uniform mesh config from this block.
    pub fn to_uniform(&self, bounds_hi: [f64; 3]) -> UniformMeshConfig {
        UniformMeshConfig {
            nx: self.nx,
            ny: self.ny,
            nz: self.nz,
            ng: self.ng,
            bounds_lo: [0.0, 0.0, 0.0],
            bounds_hi,
            y_edges: None,
            z_edges: None,
        }
    }
}
