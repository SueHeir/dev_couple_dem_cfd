//! Published/analytical references for unresolved DEM↔CFD beds. These are the
//! independent truths the seam is validated against — kept in the library so a case
//! never quietly re-derives the number it is supposed to be checked by.

/// Archimedes number `Ar = ρ_f (ρ_p − ρ_f) g d³ / μ²`.
pub fn archimedes(rho_f: f64, rho_p: f64, g: f64, d: f64, mu: f64) -> f64 {
    rho_f * (rho_p - rho_f) * g.abs() * d.powi(3) / (mu * mu)
}

/// Wen & Yu (1966) minimum fluidization velocity:
/// `Re_mf = sqrt(33.7² + 0.0408 Ar) − 33.7`, then `U_mf = Re_mf μ / (ρ_f d)`.
pub fn u_mf_wen_yu(rho_f: f64, rho_p: f64, g: f64, d: f64, mu: f64) -> f64 {
    let ar = archimedes(rho_f, rho_p, g, d, mu);
    let re_mf = (33.7f64 * 33.7 + 0.0408 * ar).sqrt() - 33.7;
    re_mf * mu / (rho_f * d)
}

/// Superficial velocity at which a packed-bed pressure drop (`c1` viscous, `c2`
/// inertial constants, porosity `eps`) equals the buoyant weight per unit length
/// `(1−ε)(ρ_p−ρ_f)g` — the incipient-fluidization criterion. Closed form of the
/// Ergun/MacDonald quadratic `a_inert U² + a_visc U − target = 0`. Used for the
/// analytic Ergun/MacDonald brackets reported alongside a SEAM-measured U_mf.
#[allow(clippy::too_many_arguments)]
pub fn u_mf_balance(
    c1: f64,
    c2: f64,
    eps: f64,
    rho_f: f64,
    rho_p: f64,
    g: f64,
    d: f64,
    mu: f64,
) -> f64 {
    let om = 1.0 - eps;
    let e3 = eps.powi(3);
    let a_visc = c1 * om / e3 * mu / (d * d); // × U
    let a_inert = c2 / e3 * rho_f / d; // × U²
    let target = (rho_p - rho_f) * g.abs(); // (dP/L)/(1−ε) at balance
    (-a_visc + (a_visc * a_visc + 4.0 * a_inert * target).sqrt()) / (2.0 * a_inert)
}

/// Ergun (1952) pressure drop per unit length for a superficial velocity.
pub fn ergun_dp_per_length(eps: f64, mu: f64, rho: f64, d: f64, u_superficial: f64) -> f64 {
    let om = 1.0 - eps;
    let e3 = eps * eps * eps;
    let viscous = 150.0 * om * om / e3 * mu * u_superficial / (d * d);
    let inertial = 1.75 * om / e3 * rho * u_superficial * u_superficial / d;
    viscous + inertial
}

/// Modified particle Reynolds number `Re_p = ρ U d / (μ (1−ε))` — the standard
/// packed-bed Reynolds that places a sweep on the Ergun viscous↔inertial map.
pub fn modified_reynolds(rho: f64, u: f64, d: f64, mu: f64, eps: f64) -> f64 {
    rho * u * d / (mu * (1.0 - eps))
}
