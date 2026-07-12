//! # dem_cfd — reusable unresolved DEM↔CFD coupling on GRASS
//!
//! The pieces every *unresolved* (Euler–Lagrange, void-fraction) DEM↔CFD
//! simulation repeats, lifted out of the examples so a case `main` reads as its
//! physics, not its plumbing:
//!
//! - [`config`] — the declarative `[gas]`/`[particle]`/`[mesh]`/`[packing]` blocks.
//! - [`drag`] — the void-fraction β closures (MacDonald/Ergun) + the [`drag::SeamMode`].
//! - [`reference`] — published references (Wen & Yu, Ergun, Archimedes) kept out of
//!   the measured path so a gate can't quietly check itself.
//! - [`bed`] — deposit the packing onto the coarse gas mesh, impose the interstitial
//!   flow, the two-way momentum sink + conservation check, FCC packing.
//! - [`seam`] — the `grass_multi` scaffold: seam resources, the CFD sub-App base,
//!   the dynamic two-way schedule ([`seam::couple_two_way`]) + its systems, and
//!   driver-side accessors.
//!
//! What stays in a case: its **force model** (the seam system — drag-only vs
//! drag+∇P+buoyancy), its **topology** if non-standard (a static packed bed is a
//! two-phase export-once schedule, not the dynamic four-phase one), and its
//! **validation** tolerances. Fully *resolved* IBM couplings (a body meshed into
//! the gas, e.g. an immersed fiber) are a different pattern and do not use this crate.

pub mod bed;
pub mod config;
pub mod drag;
pub mod reference;
#[cfg(feature = "mpi-routing")]
pub mod routing;
pub mod seam;

pub use seam::DemCfdCouplingPlugin;

/// Common imports for a case `main`.
pub mod prelude {
    pub use crate::bed::{
        axis_centers, containing_cell, deposit_bed_void_fraction, fcc_lattice_constant,
        fcc_packing, impose_interstitial_velocity, momentum_sink_and_check, nearest_center,
    };
    pub use crate::config::{default_ng, GasCfg, GravityCfg, MeshCfg, PackingCfg, ParticleCfg};
    pub use crate::drag::{beta_for, ergun_beta, macdonald_beta, SeamMode};
    pub use crate::reference::{
        archimedes, ergun_dp_per_length, modified_reynolds, u_mf_balance, u_mf_wen_yu,
    };
    pub use crate::seam::{
        bed_force, build_cfd_base, build_soil_bed, couple_two_way, export_kinematics, import_force,
        import_force_to_dem, import_force_typed, read_subapp_resource, set_seam_mode,
        set_superficial, with_subapp_resource, BodyAccel, CfdNs, CouplePhase, DemCfdCouplingPlugin,
        DemNs, FluidForces, ParticleSpec, SeamCtx, Superficial, R_GAS,
    };
}
