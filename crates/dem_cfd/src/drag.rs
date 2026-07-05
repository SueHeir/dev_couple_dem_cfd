//! Void-fraction (Euler–Lagrange) interphase drag closures β, and the seam mode
//! that selects between them / injects the negative-control faults. These are the
//! packed-bed drag laws every unresolved DEM↔CFD bed reuses verbatim.

/// MacDonald et al. (1979) interphase coefficient β — the "Ergun revisited" 180/1.8
/// re-fit. The INDEPENDENT measured closure (shares no constant with Wen & Yu 1966).
pub fn macdonald_beta(eps: f64, rho_f: f64, mu: f64, d: f64, rel_speed: f64) -> f64 {
    let eps = eps.clamp(1e-6, 1.0);
    let om = 1.0 - eps;
    180.0 * om * om * mu / (eps * d * d) + 1.8 * om * rho_f * rel_speed / d
}

/// Ergun (1952) β (150/1.75) — same functional form with the original constants.
/// Reported only for the exact-Ergun bracket / harness-sanity baseline, never the
/// measured validation closure (using it against an Ergun reference is the tautology
/// the negative controls guard against).
pub fn ergun_beta(eps: f64, rho_f: f64, mu: f64, d: f64, rel_speed: f64) -> f64 {
    let eps = eps.clamp(1e-6, 1.0);
    let om = 1.0 - eps;
    150.0 * om * om * mu / (eps * d * d) + 1.75 * om * rho_f * rel_speed / d
}

/// Which β closure the seam assembles, and whether to inject a fault (negative
/// controls). Set by the parent per pass. The superset of what the packed-bed and
/// fluidized-bed cases each need — a case simply leaves the flags it does not use
/// at their default.
#[derive(Clone, Copy)]
pub struct SeamMode {
    /// `true` → MacDonald(1979) measured closure; `false` → Ergun(1952) constants
    /// (the reported bracket / baseline, never the pass gate).
    pub macdonald: bool,
    /// Negative control A (fluidized only): drop the ∇P (pressure-gradient buoyancy)
    /// force — a real CFD–DEM mistake that shifts U_mf by ~1/ε.
    pub omit_pressure_grad: bool,
    /// Negative control B: the ε²-instead-of-ε³ reduction bug (scale the assembled
    /// force by 1/ε).
    pub corrupt_eps_power: bool,
}

impl Default for SeamMode {
    fn default() -> Self {
        Self { macdonald: true, omit_pressure_grad: false, corrupt_eps_power: false }
    }
}

/// Select the β closure named by `mode`.
pub fn beta_for(mode: SeamMode, eps: f64, rho_f: f64, mu: f64, d: f64, rel_speed: f64) -> f64 {
    if mode.macdonald {
        macdonald_beta(eps, rho_f, mu, d, rel_speed)
    } else {
        ergun_beta(eps, rho_f, mu, d, rel_speed)
    }
}
