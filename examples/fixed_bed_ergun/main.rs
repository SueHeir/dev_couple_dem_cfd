//! **Fixed-bed pressure drop vs the Ergun (1952) correlation** — the dense-regime
//! validation of the *unresolved* (point-particle) DEM–CFD seam, one rung up the
//! coupling ladder from the single-sphere settling case (`settling_sphere`).
//!
//! A static packing of spheres fills a column; gas moves up through it at a
//! prescribed **superficial velocity** `U`. At steady state the streamwise
//! pressure gradient balances the interphase drag the packing exerts on the gas,
//! and the classic packed-bed result is the **Ergun (1952)** correlation
//!
//! ```text
//!   dP/L = 150 (1−ε)²/ε³ · μ U / d²   +   1.75 (1−ε)/ε³ · ρ_f U² / d
//!          └──────── viscous (Blake–Kozeny) ───────┘   └──── inertial (Burke–Plummer) ────┘
//! ```
//!
//! with bed porosity `ε`, particle diameter `d`, superficial velocity `U`.
//!
//! ## Why the *measured* drag is an INDEPENDENT closure (no self-comparison)
//!
//! An earlier version of this example computed the "measured" `dP/L` from the
//! seam's own `beta_gidaspow` closure and compared it to the Ergun formula. That
//! is an **algebraic tautology**: the dense branch of `beta_gidaspow` *is* Ergun
//! (the very same 150 / 1.75 constants), so the measured `dP/L = |ΣF|/(V_bed ε)`
//! reduces IDENTICALLY to `ergun_dp_per_length` for any `ε, U, d` — the relative
//! error is machine-epsilon at every sweep point and the tolerance can never fail.
//! Such a check validates nothing.
//!
//! This version breaks that self-comparison. The **measured** per-particle drag is
//! assembled through the seam from an **independent** packed-bed closure —
//! **MacDonald et al. (1979)**, the "Ergun equation revisited" re-fit of the
//! Blake–Kozeny / Burke–Plummer constants to a much larger data set, which
//! recommends `180` (viscous) and `1.8` (inertial) in place of Ergun's `150` and
//! `1.75` (see [`dem_cfd::drag::macdonald_beta`]). Those constants come from data
//! Ergun never saw, so the measured value shares **no constant** with the Ergun
//! *reference*; the relative error is a genuine, non-zero, Reynolds-dependent spread
//! (~5–20 %, the documented inter-correlation difference — largest in the
//! viscous-dominated regime where `180/150 ≈ 1.20`, shrinking toward the inertial
//! end where `1.8/1.75 ≈ 1.03`).
//!
//! ### What that makes falsifiable
//!
//! The pass gate is `|dP/L_meas − dP/L_Ergun| / dP/L_Ergun ≤ tol` with `tol` set to
//! bracket the *documented* MacDonald-vs-Ergun spread. This is **not** a test of
//! the drag constants (Ergun's `150/1.75` and MacDonald's `180/1.8` are two
//! literature fits; neither is "the" truth). It tests that the unresolved seam's
//! **assembly** — deposit a discrete packing to a per-cell void fraction, read the
//! interstitial slip `u_g = U/ε`, evaluate a per-particle drag, sum it, and reduce
//! `Σ_i F_i → dP/L` — reconstructs the macroscopic packed-bed pressure-drop *law*
//! to within that spread. Any of the mistakes a CFD–DEM coupling routinely makes
//! breaks it by **factors far exceeding** the tolerance and the check FAILS:
//!
//!   * a wrong void-fraction power (`ε²` vs `ε³`) → error `~1/ε ≈ 2.5×`;
//!   * superficial `U` mistaken for interstitial `U/ε` in the slip → error `~ε`;
//!   * per-particle → bulk summation wrong (e.g. `V_bed` or the `1−ε` factor) → order-unity error.
//!
//! To prove that (and that the tolerance is not tuned to pass), this example also
//! runs a **negative control**: it re-assembles the SAME sweep with a deliberately
//! corrupted seam (the `ε²`-instead-of-`ε³` bug above) and asserts it would FAIL.
//!
//! ## The deposited void fraction is load-bearing
//!
//! The per-particle drag is driven by a per-cell void fraction **deposited from the
//! actual packing** ([`dem_cfd::bed::deposit_bed_void_fraction`]), not an analytic
//! `ε_geo`. For that field to be usable the gas mesh must be **coarser than the
//! particles** — the defining requirement of the *unresolved* regime — so each cell
//! contains *many* sub-cell spheres and its deposited void fraction converges to the
//! bed porosity. The `[mesh]` section sets a gas grid several particle diameters per
//! cell (`d/Δx ≪ 1`); the example gates that the deposited per-cell field is uniform
//! to within `tol_deposit_cell` of `ε_bed` (a real deposition-fidelity check — a
//! mis-binning deposit fails it), and then *uses* that field to evaluate the drag.
//! The binning uses a **containment** locator (in `dem_cfd::bed`): the seam's own
//! `coupling::deposit_void_fraction` uses an interpolation locator tuned for a
//! single SUB-CELL particle and mis-bins at bed scale. `ε_bed` itself is the exact
//! geometric porosity of the FCC packing, used only for the Ergun reference and the
//! superficial↔interstitial conversion.
//!
//! ## Momentum conservation (kept, but NOT the validation)
//!
//! The two-way momentum sink is still exercised and its conservation checked
//! (`|ΔP_gas + ΣF·dt| / |ΣF·dt| ~ 1e-13`). That is a legitimate property of the
//! seam but it is **not** Ergun-specific — it holds for any drag whatsoever — so it
//! is reported and gated as a sanity check, never counted as the pressure-drop
//! validation.
//!
//! This example imposes the superficial flow (interstitial `u_g = U/ε`) rather than
//! running a full driven compressible flow-through solve past the packing: a
//! porosity-weighted (volume-averaged) momentum solve is the resolved-track story,
//! out of scope for the unresolved seam (cf. the `Cd(Re)` deferral in
//! `cfd_ibm::coupling`). The falsifiability here comes from the *independent* drag
//! reference, not from a flow solve.
//!
//! Everything case-specific is declarative TOML from `argv[1]`:
//!
//! ```text
//! cargo run --release --example fixed_bed_ergun -- \
//!     examples/fixed_bed_ergun/config.toml
//! ```
//!
//! The shared unresolved DEM↔CFD machinery — config blocks, the MacDonald/Ergun β
//! closures + [`SeamMode`], the Ergun references, void-fraction deposition, the
//! interstitial-flow / momentum-sink plumbing, FCC packing, and the seam resources —
//! lives in the `dem_cfd` crate. This case keeps its STATIC topology (a two-phase
//! export-once schedule; the packing never moves) and its drag-only force model.
//!
//! References:
//! * S. Ergun, "Fluid flow through packed columns", *Chem. Eng. Prog.* 48(2):89–94 (1952).
//! * I.F. MacDonald, M.S. El-Sayed, K. Mow, F.A.L. Dullien, "Flow through Porous
//!   Media — the Ergun Equation Revisited", *Ind. Eng. Chem. Fundam.* 18(3):199–208 (1979).

use cfd_eos::{Eos, EosResource};
use cfd_ibm::coupling::{self, drag_force_from_beta, InterphaseForces, ParticleKinematics, ParticleSet};
use cfd_state::CfdState;
use field_core::{FieldRegistry, FvMesh, MeshScheduleSet, UniformMesh, UniformMeshConfig, Vec3};
use grass_app::prelude::*;
use grass_io::Config;
use grass_multi::{tick_subapp, Multi, MultiAppExt};
use grass_scheduler::prelude::*;
use grass_scheduler::{Res, ResMut};
use serde::Deserialize;
use soil_core::Atom;
use soil_verlet::VelocityVerletPlugin;

use dem_cfd::prelude::*;

// ─── Case-specific config (the shared blocks come from dem_cfd::config) ───────

#[derive(Deserialize, Default)]
struct FlowCfg {
    /// Superficial velocities U [m/s] to sweep.
    superficial: Vec<f64>,
    /// Coupling timestep for the momentum-sink integral / conservation check.
    dt: f64,
}

#[derive(Deserialize, Default)]
struct ValidationCfg {
    /// |dP/L_measured − dP/L_Ergun| / dP/L_Ergun, per sweep point. Set to bracket
    /// the documented MacDonald(1979)-vs-Ergun(1952) inter-correlation spread.
    tol_rel_ergun: f64,
    /// Regime gate: packed-bed (dense) branch requires ε ≤ this.
    eps_max: f64,
    /// The sweep must reach DOWN to at least this modified Re_p (viscous-dominated end).
    rep_viscous_below: f64,
    /// …and UP to at least this modified Re_p (inertia-comparable end).
    rep_inertial_above: f64,
    /// Two-way momentum-exchange conservation tolerance (sink vs −ΣF_drag·dt).
    tol_momentum: f64,
    /// Max |ε_cell(deposit) − ε_bed| / ε_bed over interior cells (deposition fidelity).
    tol_deposit_cell: f64,
}

/// Result of one fixed-bed evaluation, read back by the parent.
#[derive(Clone, Copy, Default)]
struct BedResult {
    /// Total interphase drag on the packing, Σ_i F_drag,i (world axes).
    f_total: Vec3,
    /// Max |ε_cell − ε_bed| / ε_bed over interior cells (deposition fidelity).
    eps_cell_err: f64,
    eps_min: f64,
    eps_max: f64,
    /// Momentum-conservation error of the two-way sink, |ΔP_gas + ΣF·dt|/|ΣF·dt|.
    mom_err: f64,
}

// ─── FIELD sub-App: impose the flow, run the seam, report the drag ───────────
//
// Plumbing (impose the interstitial velocity, deposit the void fraction, the
// momentum sink + conservation check) is from `dem_cfd::bed`; the case-specific
// part is the STATIC (u_p = 0) drag-only force model + its deposition diagnostics.

fn fixed_bed_seam_system(
    mesh: Res<UniformMesh>,
    reg: Res<FieldRegistry>,
    eos: Res<EosResource>,
    ctx: Res<SeamCtx>,
    sup: Res<Superficial>,
    pset: Res<ParticleSet>,
    mut forces: ResMut<InterphaseForces>,
    mut result: ResMut<BedResult>,
) {
    let eos: &dyn Eos = &*eos.0;
    let mut state = reg.expect_mut::<CfdState>("CfdState not registered");
    let parts = &pset.particles;
    forces.reset(parts.len());
    let mode = ctx.mode;

    // Interstitial gas velocity u_g = U/ε_bed, imposed uniformly (ρ untouched).
    let eps_bed = ctx.eps;
    let inv_eps = 1.0 / eps_bed;
    let u_g = [sup.u[0] * inv_eps, sup.u[1] * inv_eps, sup.u[2] * inv_eps];
    impose_interstitial_velocity(&mesh, &mut state, u_g);

    // Deposit the packing's solid volume onto the mesh; report worst interior-cell
    // deviation from ε_bed as the deposition-fidelity diagnostic.
    let (eps_field, cell_of_particle) = deposit_bed_void_fraction(&mesh, parts);
    let (mut e_min, mut e_max, mut e_err) = (f64::INFINITY, 0.0f64, 0.0f64);
    for (c, &e) in eps_field.iter().enumerate() {
        if !mesh.is_local_cell(c) || e >= 1.0 - 1e-9 {
            continue;
        }
        e_min = e_min.min(e);
        e_max = e_max.max(e);
        e_err = e_err.max((e - eps_bed).abs() / eps_bed);
    }
    result.eps_min = if e_min.is_finite() { e_min } else { 1.0 };
    result.eps_max = e_max;
    result.eps_cell_err = e_err;

    // Interphase drag on each particle through the selected closure, driven by the
    // DEPOSITED per-cell void fraction (load-bearing), summed.
    let mut drag_on_particle = vec![[0.0f64; 3]; parts.len()];
    let mut f_total = [0.0f64; 3];
    for (i, p) in parts.iter().enumerate() {
        let u_gas =
            coupling::sample_gas_velocity(&*mesh, &state, eos, p.center).unwrap_or([0.0; 3]);
        let rho_f = coupling::sample_gas_density(&*mesh, &state, p.center).unwrap_or(ctx.rho);
        let eps = eps_field[cell_of_particle[i]];

        // Static packing: u_p = 0, so the slip IS the local interstitial velocity.
        let rel = u_gas;
        let rel_speed = (rel[0] * rel[0] + rel[1] * rel[1] + rel[2] * rel[2]).sqrt();
        let d = p.diameter();
        let beta = beta_for(mode, eps, rho_f, ctx.mu, d, rel_speed);
        let drag = drag_force_from_beta(beta, p.volume(), eps, rel);

        forces.force[i] = drag;
        drag_on_particle[i] = drag;
        for k in 0..3 {
            f_total[k] += drag[k];
        }
    }
    // Negative control: corrupt the reduction with the ε²-instead-of-ε³ bug.
    if mode.corrupt_eps_power {
        for k in 0..3 {
            f_total[k] /= eps_bed;
        }
    }
    result.f_total = f_total;

    // Two-way momentum sink + conservation check (NOT the pressure-drop validation).
    result.mom_err = momentum_sink_and_check(&mesh, &mut state, parts, &drag_on_particle, ctx.dt);
}

fn build_cfd(gas: &GasCfg, mesh_cfg: UniformMeshConfig, eps: f64, dt: f64) -> App {
    let ctx = SeamCtx {
        mu: gas.mu,
        rho: gas.rho,
        eps,
        g: [0.0, 0.0, 0.0], // static bed: no hydrostatic buoyancy term
        dt,
        mode: SeamMode::default(),
    };
    let mut app = build_cfd_base(gas, mesh_cfg, ctx);
    app.add_resource(BedResult::default());
    app.add_update_system(fixed_bed_seam_system, MeshScheduleSet::Output);
    app
}

// ─── SOIL sub-App: the STATIC packing (never integrated) ─────────────────────

fn build_soil(positions: &[[f64; 3]], radius: f64, density: f64) -> App {
    let mut atoms = Atom::new();
    atoms.dt = 1.0; // never integrated; particles are held fixed
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
    for (tag, pos) in positions.iter().enumerate() {
        atoms.push_test_atom(tag as u32, *pos, radius, mass);
    }
    atoms.nlocal = positions.len() as u32;
    atoms.natoms = positions.len() as u64;

    let mut app = App::new();
    app.add_resource(atoms);
    // Present for parity with the settling case; the parent never ticks SOIL, so
    // the packing stays static.
    app.add_plugins(VelocityVerletPlugin::new());
    app
}

/// Run the full superficial-velocity sweep once for a given seam mode, returning
/// `(dP/L_meas per U, worst deposition-cell error, worst momentum error)`.
fn run_sweep(parent: &mut App, us: &[f64], v_bed: f64, eps_bed: f64, mode: SeamMode) -> (Vec<f64>, f64, f64) {
    set_seam_mode(parent, mode);
    let mut dpdl = Vec::with_capacity(us.len());
    let mut worst_dep = 0.0f64;
    let mut worst_mom = 0.0f64;
    for &u in us {
        set_superficial(parent, [0.0, 0.0, u]);
        parent.run();
        let res = read_result(parent);
        let f_mag =
            (res.f_total[0].powi(2) + res.f_total[1].powi(2) + res.f_total[2].powi(2)).sqrt();
        // Steady two-fluid balance: ε dP/dx = −β u_g ⇒ dP/L = |ΣF_drag| / (V_bed ε).
        dpdl.push(f_mag / (v_bed * eps_bed));
        worst_dep = worst_dep.max(res.eps_cell_err);
        worst_mom = worst_mom.max(res.mom_err);
    }
    (dpdl, worst_dep, worst_mom)
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: fixed_bed_ergun <case.toml>");
    let toml_src =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let cfg = Config::from_str(&toml_src);

    let gas: GasCfg = cfg.section("gas");
    let pc: ParticleCfg = cfg.section("particle");
    let pack: PackingCfg = cfg.section("packing");
    let meshc: MeshCfg = cfg.section("mesh");
    let flow: FlowCfg = cfg.section("flow");
    let valid: ValidationCfg = cfg.section("validation");

    let d = pc.diameter;
    let radius = 0.5 * d;
    let v_p = 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);

    // FCC lattice constant giving the requested solid fraction φ.
    assert!(
        pack.solid_fraction <= 0.7405,
        "solid_fraction {} exceeds the FCC maximum 0.7405 (spheres would overlap)",
        pack.solid_fraction
    );
    // The gas mesh must tile the packing evenly and be COARSER than the particles.
    assert!(
        pack.ncx % meshc.nx == 0 && pack.ncy % meshc.ny == 0 && pack.ncz % meshc.nz == 0,
        "mesh {}x{}x{} must divide packing {}x{}x{} evenly",
        meshc.nx, meshc.ny, meshc.nz, pack.ncx, pack.ncy, pack.ncz
    );
    let a = fcc_lattice_constant(d, pack.solid_fraction);
    let (positions, bounds) = fcc_packing([pack.ncx, pack.ncy, pack.ncz], a);
    let n = positions.len();
    let v_bed = bounds[0] * bounds[1] * bounds[2];
    let eps_geo = 1.0 - n as f64 * v_p / v_bed;
    let nn_gap = a / std::f64::consts::SQRT_2 - d; // nearest-neighbour surface gap

    let mesh_cfg = meshc.to_uniform(bounds);
    let n_cells = (meshc.nx * meshc.ny * meshc.nz) as f64;

    let soil = build_soil(&positions, radius, pc.density);
    let cfd = build_cfd(&gas, mesh_cfg, eps_geo, flow.dt);

    let mut parent = App::new();
    parent.add_subapp("soil", soil);
    parent.add_subapp("cfd", cfd);
    parent.add_update_system(export_kinematics, Phase::Export);
    parent.add_update_system(tick_subapp("cfd", 1), Phase::TickCfd);
    parent.prepare();

    // Cell size for the point-particle regime report (d/Δx should be ≪ 1).
    let dx = (v_bed / n_cells).cbrt();
    let parts_per_cell = n as f64 / n_cells;

    println!("# Fixed-bed pressure drop vs Ergun (1952) — unresolved DEM-CFD seam");
    println!("# MEASURED drag: INDEPENDENT MacDonald et al. (1979) closure (180/1.8), assembled through the seam");
    println!("# REFERENCE:     Ergun (1952) correlation (150/1.75)");
    println!(
        "# packing: FCC {}x{}x{} cells, N = {} spheres   d = {:.3e} m   a = {:.3e} m",
        pack.ncx, pack.ncy, pack.ncz, n, d, a
    );
    println!(
        "# gas mesh: {}x{}x{} cells   dx ~ {:.3e} m   d/dx = {:.3} (≪1 ⇒ unresolved)   ~{:.0} spheres/cell",
        meshc.nx, meshc.ny, meshc.nz, dx, d / dx, parts_per_cell
    );
    println!(
        "# column: {:.3e} x {:.3e} x {:.3e} m   V_bed = {:.3e} m^3   nn surface gap = {:.3e} m (>0 ⇒ non-overlapping)",
        bounds[0], bounds[1], bounds[2], v_bed, nn_gap
    );
    println!(
        "# porosity eps_geo = {:.4}   (solid frac phi = {:.4})   gas: rho = {} mu = {:.3e}",
        eps_geo,
        1.0 - eps_geo,
        gas.rho,
        gas.mu
    );
    println!("#");
    println!(
        "#     U [m/s]     Re_p      dP/L_meas [Pa/m]   dP/L_Ergun [Pa/m]   rel.err   visc:inert"
    );

    // ── The validation: MEASURED = independent MacDonald closure through the seam.
    let mode = SeamMode { macdonald: true, corrupt_eps_power: false, ..SeamMode::default() };
    let (dpdl_meas, worst_dep, worst_mom) = run_sweep(&mut parent, &flow.superficial, v_bed, eps_geo, mode);

    let mut all_ok = true;
    let mut rep_min = f64::INFINITY;
    let mut rep_max = 0.0f64;
    let mut worst_err = 0.0f64;
    let mut best_err = f64::INFINITY;
    for (idx, &u) in flow.superficial.iter().enumerate() {
        let dpdl_ergun = ergun_dp_per_length(eps_geo, gas.mu, gas.rho, d, u);
        let rel = (dpdl_meas[idx] - dpdl_ergun).abs() / dpdl_ergun;
        let rep = modified_reynolds(gas.rho, u, d, gas.mu, eps_geo);
        let om = 1.0 - eps_geo;
        let e3 = eps_geo.powi(3);
        let visc = 150.0 * om * om / e3 * gas.mu * u / (d * d);
        let inert = 1.75 * om / e3 * gas.rho * u * u / d;

        rep_min = rep_min.min(rep);
        rep_max = rep_max.max(rep);
        worst_err = worst_err.max(rel);
        best_err = best_err.min(rel);
        all_ok &= rel <= valid.tol_rel_ergun;

        println!(
            "  {u:>10.4}  {rep:>7.2}   {:>16.4}   {dpdl_ergun:>16.4}   {:>6.3}%   {:>6.2}:1",
            dpdl_meas[idx],
            100.0 * rel,
            visc / inert.max(1e-30)
        );
    }

    // ── Negative control: corrupt the seam reduction (ε² instead of ε³) and show
    // it blows past the tolerance — proof the pass is genuinely capable of failing
    // and the tolerance is not tuned to admit a broken seam.
    let corrupt = SeamMode { macdonald: true, corrupt_eps_power: true, ..SeamMode::default() };
    let (dpdl_bad, _, _) = run_sweep(&mut parent, &flow.superficial, v_bed, eps_geo, corrupt);
    let mut worst_bad = 0.0f64;
    for (idx, &u) in flow.superficial.iter().enumerate() {
        let dpdl_ergun = ergun_dp_per_length(eps_geo, gas.mu, gas.rho, d, u);
        worst_bad = worst_bad.max((dpdl_bad[idx] - dpdl_ergun).abs() / dpdl_ergun);
    }
    let neg_control_fails = worst_bad > valid.tol_rel_ergun;

    println!("#");
    println!("# ── result ─────────────────────────────────────────────");
    println!(
        "# deposition fidelity: worst |eps_cell - eps_bed|/eps_bed = {:.2}% (tol {:.1}%)  [the deposited field DRIVES the drag]",
        100.0 * worst_dep,
        100.0 * valid.tol_deposit_cell
    );
    println!(
        "# Re_p sweep span: [{:.2}, {:.2}]   (need ≤ {} viscous end AND ≥ {} inertial end)",
        rep_min, rep_max, valid.rep_viscous_below, valid.rep_inertial_above
    );
    println!(
        "# Ergun rel.err spread over sweep: [{:.3}%, {:.3}%]  (tol {:.1}%)  — MacDonald-vs-Ergun, non-zero ⇒ not a tautology",
        100.0 * best_err,
        100.0 * worst_err,
        100.0 * valid.tol_rel_ergun
    );
    println!(
        "# negative control (corrupted eps-power seam): worst rel.err {:.1}%  ⇒ {} (must exceed tol {:.1}%)",
        100.0 * worst_bad,
        if neg_control_fails { "FAILS as required" } else { "DID NOT FAIL — gate is vacuous!" },
        100.0 * valid.tol_rel_ergun
    );
    println!(
        "# momentum conservation err: {worst_mom:.2e}  (tol {:.0e})  [sanity only — not the Ergun validation]",
        valid.tol_momentum
    );

    let pass_regime = eps_geo <= valid.eps_max;
    let pass_span = rep_min <= valid.rep_viscous_below && rep_max >= valid.rep_inertial_above;
    let pass_dep = worst_dep <= valid.tol_deposit_cell;
    let pass_mom = worst_mom <= valid.tol_momentum;
    let pass_ergun = all_ok;
    // The pass is only meaningful if the rel.err is genuinely non-zero (independent
    // reference, not a self-comparison) AND the corrupted seam would fail.
    let pass_nontrivial = worst_err > 1e-4 && neg_control_fails;

    if pass_ergun && pass_regime && pass_span && pass_dep && pass_mom && pass_nontrivial {
        println!(
            "VALIDATION: PASS  (independent MacDonald vs Ergun: rel.err∈[{:.2}%,{:.2}%]<={:.1}% over Re_p∈[{:.1},{:.1}]; eps={:.3}<={}; deposit cell {:.1}%; neg-control fails at {:.0}%; mom {:.1e})",
            100.0 * best_err,
            100.0 * worst_err,
            100.0 * valid.tol_rel_ergun,
            rep_min,
            rep_max,
            eps_geo,
            valid.eps_max,
            100.0 * worst_dep,
            100.0 * worst_bad,
            worst_mom
        );
    } else {
        println!(
            "VALIDATION: FAIL  (ergun_ok={pass_ergun} regime_ok={pass_regime} span_ok={pass_span} dep_ok={pass_dep} mom_ok={pass_mom} nontrivial_ok={pass_nontrivial})"
        );
        std::process::exit(1);
    }
}

// ─── Parent coupling schedule (STATIC: export the packing once, tick the gas) ─

#[derive(Debug, Clone, Copy)]
enum Phase {
    Export,
    TickCfd,
}
impl ScheduleSet for Phase {
    fn to_index(&self) -> u32 {
        match self {
            Self::Export => 0,
            Self::TickCfd => 1,
        }
    }
    fn name(&self) -> &'static str {
        match self {
            Self::Export => "Export",
            Self::TickCfd => "TickCfd",
        }
    }
}

/// SOIL→FIELD half of the seam: hand the static packing's kinematics across (once).
fn export_kinematics(world: Multi) {
    let atoms = world.expect_read::<Atom>("soil");
    let n = atoms.nlocal as usize;
    let mut set = world.expect_write::<ParticleSet>("cfd");
    if set.particles.len() == n {
        return; // static packing — export once
    }
    let radius = atoms.cutoff_radius.first().copied().unwrap_or(0.0) as f64;
    set.particles.clear();
    for i in 0..n {
        set.particles.push(ParticleKinematics {
            center: [
                atoms.pos[i][0] as f64,
                atoms.pos[i][1] as f64,
                atoms.pos[i][2] as f64,
            ],
            velocity: [0.0; 3],
            radius,
        });
    }
}

/// (Σ drag, deposition err, momentum err, …) read back from the `cfd` sub-App.
fn read_result(parent: &App) -> BedResult {
    read_subapp_resource::<BedResult>(parent, "cfd")
}
