//! **Minimum fluidization velocity `U_mf` via the unresolved DEM–CFD seam** — the
//! *dynamic* two-way-coupling capstone of the coupling ladder, one rung above the
//! static fixed-bed pressure drop (`fixed_bed_ergun`).
//!
//! A packing of spheres (the SOIL/Lagrangian bed) is immersed in an upward gas
//! stream (the FIELD/Eulerian gas). Gravity pulls the bed down; the interphase
//! force pushes it up. The **minimum fluidization velocity** `U_mf` is the
//! superficial gas velocity at which the total upward fluid force on the bed just
//! balances its net (buoyant) weight — below it the bed stays packed on the
//! distributor, above it the bed lifts and fluidizes. The classic result equates
//! the Ergun pressure drop to the buoyant weight per unit bed length,
//!
//! ```text
//!   (dP/L)|_{U_mf}  =  (1 − ε)(ρ_p − ρ_f) g,
//! ```
//!
//! and the reference used here is the **Wen & Yu (1966)** closed-form correlation
//!
//! ```text
//!   Re_mf = sqrt(33.7² + 0.0408 Ar) − 33.7,     Ar = ρ_f (ρ_p − ρ_f) g d³ / μ²,
//!   U_mf  = Re_mf · μ / (ρ_f d),
//! ```
//!
//! (`Ar` = Archimedes number). This is the standard CFD–DEM fluidization
//! validation, cf. Goniva et al. (2012) for CFDEM.
//!
//! ## Why this is DYNAMIC (and what it adds over `fixed_bed_ergun`)
//!
//! The fixed-bed rung imposes a **static** packing and checks a pressure drop. Here
//! the bed particles are **free to move**: each is a `soil::Atom` integrated by
//! `soil_verlet` under gravity + the fluid force handed across the seam. The whole
//! coupled system is stepped, and the bed's own **centre-of-mass acceleration**
//! `a_z(U)` is *measured from the integration*, not computed on paper. That curve
//! is the fluidization signature:
//!
//!   * `U < U_mf` → `a_z < 0`: net force is downward, the bed compacts onto the
//!     distributor and stays a packed fixed bed;
//!   * `U > U_mf` → `a_z > 0`: net force is upward, the bed lifts — it fluidizes;
//!   * the zero-crossing `a_z(U_mf) = 0` is incipient fluidization.
//!
//! `U_mf` is extracted two ways that must agree: (i) a bisection on the **measured**
//! net bed force read back through the live seam at rest, and (ii) the zero-crossing
//! of the **integrated** `a_z(U)` sweep. Both are then compared to Wen & Yu.
//!
//! ## The total fluid force on a bed particle — drag AND the ∇P (buoyancy) force
//!
//! A particle in a bed feels two fluid forces, and getting `U_mf` right requires
//! **both**:
//!
//!   1. **Interphase drag**, `F_drag = β V_p /(1−ε) · (u_g − u_p)` (the same
//!      [`coupling::drag_force_from_beta`] the settling and fixed-bed cases use).
//!      Summed over the packing this is only `ε · (dP/L) V_bed`.
//!   2. **The pressure-gradient (generalized buoyancy) force**, `−V_p ∇P`. In the
//!      imposed-flow (porous-medium) model the gas-phase momentum balance fixes the
//!      streamwise gradient at `ε ∇P = −β (u_g − u_p)`, so the force on each
//!      particle is `+V_p β (u_g − u_p)/ε`. Summed this is the remaining
//!      `(1−ε) · (dP/L) V_bed`. This is the standard unresolved-CFD–DEM "pressure
//!      gradient force" (e.g. Kafui et al. 2002; Zhou et al. 2010); it is the
//!      bed-scale generalization of the tiny `−ρ_f V g` hydrostatic buoyancy the
//!      settling case carries (recovered as the `∇P → ρ_f g`, no-flow limit).
//!
//! Only `F_drag + F_∇P = (dP/L) V_bed` balances the buoyant weight at the correct
//! velocity. **Dropping the ∇P force is a real, commonly-made CFD–DEM mistake**; it
//! shifts `U_mf` by a factor `~1/ε` (here ≈ +80 %), so this example runs it as a
//! **negative control** and asserts the Wen-&-Yu gate then FAILS. A second negative
//! control injects the fixed-bed's `ε²`-instead-of-`ε³` reduction bug (≈ −53 %).
//!
//! ## Non-tautology (independent measured closure)
//!
//! Exactly as in `fixed_bed_ergun`, the **measured** drag is assembled from the
//! *independent* **MacDonald et al. (1979)** packed-bed closure (`180 / 1.8`), which
//! shares **no constant** with the Ergun-based Wen & Yu (1966) reference
//! (`150 / 1.75` → `33.7 / 0.0408`). The resulting `U_mf` therefore differs from
//! Wen & Yu by a genuine, non-zero spread (≈ 4.5 % here) — the documented
//! inter-correlation difference, well inside Wen & Yu's own ~34 % reported scatter —
//! not by construction. For orientation the example also reports the crossover
//! obtained with the *exact* Ergun (1952) constants, which brackets Wen & Yu from
//! the other side (≈ +4 %).
//!
//! The gas flow is **imposed** (interstitial `u_g = U/ε`), not solved — a porosity-
//! weighted momentum solve is the resolved-track story, out of scope for the
//! unresolved seam (cf. the `Cd(Re)` deferral in `cfd_ibm::coupling`); the
//! falsifiability here comes from the independent drag reference and the dynamic
//! force balance, not from a flow solve. Everything case-specific is declarative
//! TOML from `argv[1]`:
//!
//! ```text
//! cargo run --release --example fluidized_bed_umf -- \
//!     examples/fluidized_bed_umf/config.toml
//! ```
//!
//! The shared unresolved DEM↔CFD machinery — the `[gas]`/`[particle]`/`[mesh]`/
//! `[packing]` config blocks, the MacDonald/Ergun β closures and [`SeamMode`], the
//! Wen&Yu/Archimedes references, void-fraction deposition, the interstitial-flow /
//! momentum-sink plumbing, and the `grass_multi` two-way seam scaffold — lives in
//! the `dem_cfd` crate. What remains here is this case's *force model* (drag + ∇P +
//! buoyancy), its dynamic `U_mf` measurement, and its validation gates.
//!
//! References:
//! * C.Y. Wen & Y.H. Yu, "A generalized method for predicting the minimum
//!   fluidization velocity", *AIChE J.* 12(3):610–612 (1966).
//! * S. Ergun, "Fluid flow through packed columns", *Chem. Eng. Prog.* 48(2):89 (1952).
//! * I.F. MacDonald, M.S. El-Sayed, K. Mow, F.A.L. Dullien, "Flow through Porous
//!   Media — the Ergun Equation Revisited", *Ind. Eng. Chem. Fundam.* 18(3):199 (1979).
//! * C. Goniva et al., "Influence of rolling friction on single spout fluidized bed
//!   simulation", *Particuology* 10(5):582–591 (2012) — CFDEM CFD–DEM validation.

use cfd_eos::{Eos, EosResource};
use cfd_ibm::coupling::{self, drag_force_from_beta, InterphaseForces, ParticleSet};
use cfd_state::CfdState;
use field_core::{FieldRegistry, FvMesh, MeshScheduleSet, UniformMesh, Vec3};
use grass_app::prelude::*;
use grass_io::Config;
use grass_scheduler::{Res, ResMut};
use serde::Deserialize;
use soil_core::Atom;

use dem_cfd::prelude::*;

// ─── Case-specific config (the shared blocks come from dem_cfd::config) ───────

#[derive(Deserialize, Default)]
struct FlowCfg {
    /// Superficial velocities U [m/s] to sweep (must bracket U_mf: some below, some above).
    superficial: Vec<f64>,
    /// Macro coupling interval to compare across strategies.
    macro_dt: f64,
    /// Number of macro intervals integrated from rest to measure a_z(U).
    macro_steps: usize,
}

#[derive(Clone, Deserialize)]
#[serde(rename_all = "snake_case")]
enum StrategyKind {
    FixedExplicit,
    FixedSubcycle,
    ResidualGated,
}

impl Default for StrategyKind {
    fn default() -> Self {
        Self::FixedExplicit
    }
}

#[derive(Clone, Deserialize, Default)]
struct StrategyCfg {
    name: String,
    kind: StrategyKind,
    /// Fixed substep count for fixed strategies.
    substeps: Option<usize>,
    /// Candidate substep counts for residual-gated coupling.
    candidates: Option<Vec<usize>>,
    /// Accept when max |a_z(measured)-a_z(force)|/g is below this value.
    residual_tol: Option<f64>,
}

#[derive(Deserialize, Default)]
struct ValidationCfg {
    /// |U_mf_meas − U_mf_WenYu| / U_mf_WenYu. Brackets the documented MacDonald(1979)-
    /// vs-Wen&Yu(1966) closure spread; well inside Wen&Yu's own ~34% scatter.
    tol_rel_umf: f64,
    /// The measured rel.err must be at least this (non-tautology floor: an
    /// independent closure, not a self-comparison, so the error is genuinely > 0).
    umf_err_floor: f64,
    /// |a_z_measured(integrated) − a_z_force(net seam force / M_bed)| / g, per sweep
    /// point — the two-way handoff: the integrator must deliver the seam force to
    /// the freely-moving DEM bed.
    tol_handoff: f64,
    /// Regime gate: packed-bed (dense Ergun/MacDonald) branch requires ε ≤ this.
    eps_max: f64,
    /// Two-way momentum-exchange conservation tolerance (sink vs −ΣF_drag·dt).
    tol_momentum: f64,
    /// Worst |ε_cell(deposit) − ε_bed| / ε_bed over interior cells (deposition
    /// fidelity — the deposited void fraction DRIVES the drag).
    tol_deposit_cell: f64,
    /// The residual-gated strategy must reduce worst residual by at least this factor
    /// versus the one-shot explicit strategy.
    min_residual_improvement: f64,
}

/// Read back by the parent each coupled step.
#[derive(Clone, Copy, Default)]
struct BedResult {
    /// Σ drag force on the packing (world axes) — used for the momentum check.
    f_drag_total: Vec3,
    /// Σ TOTAL fluid force (drag + ∇P buoyancy + hydrostatic) — the force that
    /// balances the bed weight at U_mf.
    f_fluid_total: Vec3,
    eps_cell_err: f64,
    mom_err: f64,
}

// ─── FIELD sub-App: impose the flow, run the seam, deliver per-particle force ─
//
// The plumbing (impose the interstitial velocity, deposit the void fraction, the
// momentum sink + conservation check) comes from `dem_cfd::bed`; what is written
// out longhand here is this case's FORCE MODEL — drag + the ∇P (pressure-gradient
// buoyancy) force + hydrostatic buoyancy — and the two negative-control faults.

#[allow(clippy::too_many_arguments)]
fn fluidized_seam_system(
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

    // Impose the interstitial gas velocity u_g = U/ε_bed in every interior cell.
    let inv_eps = 1.0 / ctx.eps;
    let u_g = [sup.u[0] * inv_eps, sup.u[1] * inv_eps, sup.u[2] * inv_eps];
    impose_interstitial_velocity(&mesh, &mut state, u_g);

    // Deposit the packing's void fraction (drives the drag) and report fidelity.
    let (eps_field, cell_of_particle) = deposit_bed_void_fraction(&mesh, parts);
    let mut e_err = 0.0f64;
    for (c, &e) in eps_field.iter().enumerate() {
        if !mesh.is_local_cell(c) || e >= 1.0 - 1e-9 {
            continue;
        }
        e_err = e_err.max((e - ctx.eps).abs() / ctx.eps);
    }
    result.eps_cell_err = e_err;

    let mut drag_on_particle = vec![[0.0f64; 3]; parts.len()];
    let mut f_drag = [0.0f64; 3];
    let mut f_fluid = [0.0f64; 3];
    for (i, p) in parts.iter().enumerate() {
        let u_gas = coupling::sample_gas_velocity(&*mesh, &state, eos, p.center).unwrap_or(u_g);
        let rho_f = coupling::sample_gas_density(&*mesh, &state, p.center).unwrap_or(ctx.rho);
        let eps = eps_field[cell_of_particle[i]];

        // Slip = local interstitial gas velocity − particle velocity.
        let rel = [
            u_gas[0] - p.velocity[0],
            u_gas[1] - p.velocity[1],
            u_gas[2] - p.velocity[2],
        ];
        let rel_speed = (rel[0] * rel[0] + rel[1] * rel[1] + rel[2] * rel[2]).sqrt();
        let d = p.diameter();
        let beta = beta_for(mode, eps, rho_f, ctx.mu, d, rel_speed);
        let v_p = p.volume();

        // (1) Interphase drag F = β V_p/(1−ε) · u_rel  (Σ = ε (dP/L) V_bed).
        let drag = drag_force_from_beta(beta, v_p, eps, rel);

        // (2) Pressure-gradient (generalized buoyancy) force −V_p ∇P. The gas-phase
        //     momentum balance in the imposed-flow model fixes ε ∇P = −β u_rel, so
        //     the force on the particle is +V_p β u_rel / ε  (Σ = (1−ε)(dP/L) V_bed).
        //     Dropped by the omit_pressure_grad negative control.
        let pg_coeff = if mode.omit_pressure_grad {
            0.0
        } else {
            v_p * beta / eps
        };
        let f_pg = [pg_coeff * rel[0], pg_coeff * rel[1], pg_coeff * rel[2]];

        // (3) Hydrostatic buoyancy −ρ_f V g (tiny in gas; the ∇P→ρ_f g no-flow limit).
        let buoy = [
            -rho_f * v_p * ctx.g[0],
            -rho_f * v_p * ctx.g[1],
            -rho_f * v_p * ctx.g[2],
        ];

        let mut total = [
            drag[0] + f_pg[0] + buoy[0],
            drag[1] + f_pg[1] + buoy[1],
            drag[2] + f_pg[2] + buoy[2],
        ];
        // Negative control B: ε²-instead-of-ε³ reduction bug (scale total by 1/ε).
        if mode.corrupt_eps_power {
            for k in 0..3 {
                total[k] /= ctx.eps;
            }
        }

        forces.force[i] = total;
        drag_on_particle[i] = drag;
        for k in 0..3 {
            f_drag[k] += drag[k];
            f_fluid[k] += total[k];
        }
    }
    result.f_drag_total = f_drag;
    result.f_fluid_total = f_fluid;

    // Two-way momentum sink (reaction of the DRAG part) + conservation check.
    result.mom_err = momentum_sink_and_check(&mesh, &mut state, parts, &drag_on_particle, ctx.dt);
}

fn build_cfd(
    gas: &GasCfg,
    mesh_cfg: field_core::UniformMeshConfig,
    eps: f64,
    gz: f64,
    dt: f64,
) -> App {
    let ctx = SeamCtx {
        mu: gas.mu,
        rho: gas.rho,
        eps,
        g: [0.0, 0.0, gz],
        dt,
        mode: SeamMode::default(),
    };
    let mut app = build_cfd_base(gas, mesh_cfg, ctx);
    app.add_resource(BedResult::default());
    app.add_update_system(fluidized_seam_system, MeshScheduleSet::Output);
    app
}

#[derive(Clone, Copy)]
struct Sample {
    u: f64,
    a_meas: f64,
    a_force: f64,
    residual: f64,
    dep: f64,
    mom: f64,
    substeps: usize,
    accepted: bool,
}

struct StrategyOutcome {
    name: String,
    kind: &'static str,
    samples: Vec<Sample>,
    u_mf_dyn: Option<f64>,
    worst_residual: f64,
    worst_dep: f64,
    worst_mom: f64,
    max_substeps: usize,
    accepted_all: bool,
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: adaptive_umf_strategy <case.toml>");
    let toml_src =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let cfg = Config::from_str(&toml_src);

    let gas: GasCfg = cfg.section("gas");
    let pc: ParticleCfg = cfg.section("particle");
    let pack: PackingCfg = cfg.section("packing");
    let meshc: MeshCfg = cfg.section("mesh");
    let grav: GravityCfg = cfg.section("gravity");
    let flow: FlowCfg = cfg.section("flow");
    let valid: ValidationCfg = cfg.section("validation");
    let mut strategies: Vec<StrategyCfg> = cfg.parse_array("strategy");
    if strategies.is_empty() {
        strategies = vec![
            StrategyCfg {
                name: "fixed_explicit".into(),
                kind: StrategyKind::FixedExplicit,
                substeps: Some(1),
                candidates: None,
                residual_tol: None,
            },
            StrategyCfg {
                name: "subcycled_4".into(),
                kind: StrategyKind::FixedSubcycle,
                substeps: Some(4),
                candidates: None,
                residual_tol: None,
            },
        ];
    }

    let d = pc.diameter;
    let radius = 0.5 * d;
    let v_p = 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
    let g = grav.gz.abs();

    assert!(
        pack.solid_fraction <= 0.7405,
        "solid_fraction {} exceeds FCC max 0.7405",
        pack.solid_fraction
    );
    assert!(
        pack.ncx % meshc.nx == 0 && pack.ncy % meshc.ny == 0 && pack.ncz % meshc.nz == 0,
        "mesh {}x{}x{} must divide packing {}x{}x{} evenly",
        meshc.nx,
        meshc.ny,
        meshc.nz,
        pack.ncx,
        pack.ncy,
        pack.ncz
    );
    let a = fcc_lattice_constant(d, pack.solid_fraction);
    let (positions, bounds) = fcc_packing([pack.ncx, pack.ncy, pack.ncz], a);
    let n = positions.len();
    let v_bed = bounds[0] * bounds[1] * bounds[2];
    let eps = 1.0 - n as f64 * v_p / v_bed;
    let nn_gap = a / std::f64::consts::SQRT_2 - d;

    let mass = pc.density * v_p;
    let m_bed = mass * n as f64;
    let w_buoy = (pc.density - gas.rho) * (1.0 - eps) * v_bed * g;

    // Analytic references.
    let ar = archimedes(gas.rho, pc.density, g, d, gas.mu);
    let u_wy = u_mf_wen_yu(gas.rho, pc.density, g, d, gas.mu);
    let u_erg = u_mf_balance(150.0, 1.75, eps, gas.rho, pc.density, g, d, gas.mu);
    let u_mac = u_mf_balance(180.0, 1.8, eps, gas.rho, pc.density, g, d, gas.mu);
    let re_wy = gas.rho * u_wy * d / gas.mu;

    let n_cells = (meshc.nx * meshc.ny * meshc.nz) as f64;
    let dx = (v_bed / n_cells).cbrt();

    println!("# Adaptive coupling strategy near U_mf — unresolved DEM-CFD seam");
    println!("# MEASURED drag: INDEPENDENT MacDonald et al. (1979) closure (180/1.8), assembled through the seam");
    println!(
        "# REFERENCE:     Wen & Yu (1966) correlation  Re_mf = sqrt(33.7^2 + 0.0408 Ar) - 33.7"
    );
    println!(
        "# packing: FCC {}x{}x{} cells, N = {} spheres   d = {:.3e} m   a = {:.3e} m   nn gap = {:.2e} m",
        pack.ncx, pack.ncy, pack.ncz, n, d, a, nn_gap
    );
    println!(
        "# gas mesh: {}x{}x{} cells   dx ~ {:.3e} m   d/dx = {:.3} (<<1 => unresolved)   ~{:.0} spheres/cell",
        meshc.nx, meshc.ny, meshc.nz, dx, d / dx, n as f64 / n_cells
    );
    println!(
        "# porosity eps = {:.4}   rho_p = {}   rho_f = {}   mu = {:.3e}   g = {}",
        eps, pc.density, gas.rho, gas.mu, g
    );
    println!(
        "# M_bed = {:.4e} kg   W_buoy = (1-eps)(rho_p-rho_f)V_bed g = {:.4e} N",
        m_bed, w_buoy
    );
    println!("# Ar = {ar:.1}");
    println!("# U_mf (Wen&Yu 1966, REFERENCE) = {u_wy:.5} m/s   (Re_mf = {re_wy:.3})",);
    println!(
        "# U_mf (Ergun 1952 exact balance, bracket) = {u_erg:.5} m/s  ({:+.2}% vs Wen&Yu)",
        100.0 * (u_erg / u_wy - 1.0)
    );
    println!(
        "# macro_dt = {:.3e} s   macro_steps = {}",
        flow.macro_dt, flow.macro_steps
    );
    println!("# U_mf (MacDonald 1979 exact balance)       = {u_mac:.5} m/s  ({:+.2}% vs Wen&Yu)  <- the seam should measure this",
        100.0 * (u_mac / u_wy - 1.0));
    println!("#");

    // ── (1) Measure U_mf through the LIVE seam: bisection on the net bed force
    // read back at rest (MacDonald closure, full physics). f_net(U) = ΣF_fluid_z − W_full.
    let mode = SeamMode::default();
    // Bracket then bisect.
    let (mut lo, mut hi) = (1e-4, 10.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if rest_net_force(
            &positions,
            radius,
            pc.density,
            &gas,
            &meshc,
            bounds,
            eps,
            grav.gz,
            flow.macro_dt,
            mid,
            mode,
            mass,
            g,
            n,
        ) < 0.0
        {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let u_mf_seam = 0.5 * (lo + hi);
    let rel_umf = (u_mf_seam - u_wy).abs() / u_wy;

    // ── (2) Dynamic strategy sweep: same live seam, different schedule cadence.
    println!("# CSV strategy,kind,U,Re,a_meas,a_force,residual_g,substeps,accepted");
    let mut outcomes = Vec::new();
    for scfg in &strategies {
        let out = run_strategy(
            scfg, &flow, &positions, radius, pc.density, &gas, &meshc, bounds, eps, grav.gz, mode,
            mass, m_bed, g, n,
        );
        for s in &out.samples {
            let re = gas.rho * s.u * d / gas.mu;
            println!(
                "CSV,{},{},{:.6},{:.6},{:.8},{:.8},{:.8},{},{}",
                out.name,
                out.kind,
                s.u,
                re,
                s.a_meas,
                s.a_force,
                s.residual,
                s.substeps,
                s.accepted
            );
        }
        outcomes.push(out);
    }

    // ── (3) Negative controls (RUN, not asserted on paper): each shifts the seam
    // U_mf far outside the Wen&Yu tolerance.
    let u_mf_nopg = bisect_umf(
        &positions,
        radius,
        pc.density,
        &gas,
        &meshc,
        bounds,
        eps,
        grav.gz,
        flow.macro_dt,
        SeamMode {
            omit_pressure_grad: true,
            ..mode
        },
        mass,
        g,
        n,
    );
    let u_mf_epsbug = bisect_umf(
        &positions,
        radius,
        pc.density,
        &gas,
        &meshc,
        bounds,
        eps,
        grav.gz,
        flow.macro_dt,
        SeamMode {
            corrupt_eps_power: true,
            ..mode
        },
        mass,
        g,
        n,
    );
    let err_nopg = (u_mf_nopg - u_wy).abs() / u_wy;
    let err_epsbug = (u_mf_epsbug - u_wy).abs() / u_wy;
    let neg_ok = err_nopg > valid.tol_rel_umf && err_epsbug > valid.tol_rel_umf;

    println!("#");
    println!("# ── result ─────────────────────────────────────────────");
    println!("# U_mf MEASURED (seam, MacDonald via live net-force bisection): {u_mf_seam:.5} m/s",);
    println!(
        "# U_mf Wen&Yu (1966) REFERENCE:                                 {u_wy:.5} m/s   rel.err {:.2}%  (tol {:.1}%)",
        100.0 * rel_umf, 100.0 * valid.tol_rel_umf
    );
    for out in &outcomes {
        match out.u_mf_dyn {
            Some(u) => println!(
                "# strategy {:>16}: U_mf_dyn {u:.5} m/s ({:.2}% vs Wen&Yu), worst residual {:.4} g, max substeps {}, accepted_all={}",
                out.name,
                100.0 * (u - u_wy).abs() / u_wy,
                out.worst_residual,
                out.max_substeps,
                out.accepted_all,
            ),
            None => println!(
                "# strategy {:>16}: no dynamic zero crossing, worst residual {:.4} g",
                out.name, out.worst_residual
            ),
        }
    }
    println!(
        "# non-tautology: measured rel.err {:.2}% is genuinely > {:.1}% floor (independent MacDonald closure, not self-comparison)",
        100.0 * rel_umf, 100.0 * valid.umf_err_floor
    );
    println!(
        "# negative controls: omit-∇P U_mf {u_mf_nopg:.4} ({:+.1}%)  eps-power-bug U_mf {u_mf_epsbug:.4} ({:+.1}%)  => {} (must exceed tol {:.1}%)",
        100.0 * (u_mf_nopg / u_wy - 1.0),
        100.0 * (u_mf_epsbug / u_wy - 1.0),
        if neg_ok { "both FAIL as required" } else { "a control DID NOT FAIL — gate vacuous!" },
        100.0 * valid.tol_rel_umf
    );
    let fixed_residual = outcomes
        .iter()
        .find(|o| {
            matches!(
                strategies
                    .iter()
                    .find(|s| s.name == o.name)
                    .map(|s| &s.kind),
                Some(StrategyKind::FixedExplicit)
            )
        })
        .map(|o| o.worst_residual)
        .unwrap_or(f64::INFINITY);
    let adaptive_residual = outcomes
        .iter()
        .find(|o| {
            matches!(
                strategies
                    .iter()
                    .find(|s| s.name == o.name)
                    .map(|s| &s.kind),
                Some(StrategyKind::ResidualGated)
            )
        })
        .map(|o| o.worst_residual)
        .unwrap_or(f64::INFINITY);
    let residual_improvement = fixed_residual / adaptive_residual.max(1e-30);
    println!(
        "# strategy residual improvement fixed/adaptive = {residual_improvement:.2}x (min {:.2}x)",
        valid.min_residual_improvement
    );
    println!(
        "# deposition fidelity: worst |eps_cell − eps_bed|/eps_bed = {:.2}%  (tol {:.1}%)",
        100.0 * outcomes.iter().map(|o| o.worst_dep).fold(0.0, f64::max),
        100.0 * valid.tol_deposit_cell
    );
    println!(
        "# momentum conservation err: {:.2e}  (tol {:.0e})  [sanity only]",
        outcomes.iter().map(|o| o.worst_mom).fold(0.0, f64::max),
        valid.tol_momentum
    );

    let pass_umf = rel_umf <= valid.tol_rel_umf;
    let pass_nontrivial = rel_umf > valid.umf_err_floor && neg_ok;
    let pass_dyn = outcomes.iter().all(|o| {
        let sign_ok = o.samples.first().map(|x| x.a_meas < 0.0).unwrap_or(false)
            && o.samples.last().map(|x| x.a_meas > 0.0).unwrap_or(false)
            && o.samples.windows(2).all(|w| w[1].a_meas > w[0].a_meas);
        let dyn_matches = o
            .u_mf_dyn
            .map(|u| (u - u_mf_seam).abs() / u_mf_seam < 0.05)
            .unwrap_or(false);
        sign_ok && dyn_matches
    });
    let pass_handoff = adaptive_residual <= valid.tol_handoff
        && residual_improvement >= valid.min_residual_improvement
        && outcomes
            .iter()
            .filter(|o| o.kind == "residual_gated")
            .all(|o| o.accepted_all);
    let pass_regime = eps <= valid.eps_max;
    let pass_dep = outcomes
        .iter()
        .all(|o| o.worst_dep <= valid.tol_deposit_cell);
    let pass_mom = outcomes.iter().all(|o| o.worst_mom <= valid.tol_momentum);

    if pass_umf
        && pass_nontrivial
        && pass_dyn
        && pass_handoff
        && pass_regime
        && pass_dep
        && pass_mom
    {
        println!(
            "VALIDATION: PASS  (U_mf seam {u_mf_seam:.4} vs Wen&Yu {u_wy:.4}, {:.2}%<={:.1}%; all dynamic onsets match; adaptive residual {:.3}<={:.3} and improves {:.1}x; neg-controls fail at {:+.0}%/{:+.0}%; eps={eps:.3}<={}; dep {:.1}%; mom {:.1e})",
            100.0 * rel_umf,
            100.0 * valid.tol_rel_umf,
            adaptive_residual,
            valid.tol_handoff,
            residual_improvement,
            100.0 * (u_mf_nopg / u_wy - 1.0),
            100.0 * (u_mf_epsbug / u_wy - 1.0),
            valid.eps_max,
            100.0 * outcomes.iter().map(|o| o.worst_dep).fold(0.0, f64::max),
            outcomes.iter().map(|o| o.worst_mom).fold(0.0, f64::max),
        );
    } else {
        println!(
            "VALIDATION: FAIL  (umf_ok={pass_umf} nontrivial_ok={pass_nontrivial} dynamic_ok={pass_dyn} handoff_ok={pass_handoff} regime_ok={pass_regime} dep_ok={pass_dep} mom_ok={pass_mom})"
        );
        std::process::exit(1);
    }
}

fn run_strategy(
    scfg: &StrategyCfg,
    flow: &FlowCfg,
    positions: &[[f64; 3]],
    radius: f64,
    density: f64,
    gas: &GasCfg,
    meshc: &MeshCfg,
    bounds: [f64; 3],
    eps: f64,
    gz: f64,
    mode: SeamMode,
    mass: f64,
    m_bed: f64,
    g: f64,
    n: usize,
) -> StrategyOutcome {
    let kind = match scfg.kind {
        StrategyKind::FixedExplicit => "fixed_explicit",
        StrategyKind::FixedSubcycle => "fixed_subcycle",
        StrategyKind::ResidualGated => "residual_gated",
    };
    let mut samples = Vec::new();
    for &u in &flow.superficial {
        let sample = match scfg.kind {
            StrategyKind::FixedExplicit => measure_strategy_sample(
                positions,
                radius,
                density,
                gas,
                meshc,
                bounds,
                eps,
                gz,
                mode,
                mass,
                m_bed,
                g,
                n,
                u,
                flow.macro_dt,
                flow.macro_steps,
                scfg.substeps.unwrap_or(1).max(1),
            ),
            StrategyKind::FixedSubcycle => measure_strategy_sample(
                positions,
                radius,
                density,
                gas,
                meshc,
                bounds,
                eps,
                gz,
                mode,
                mass,
                m_bed,
                g,
                n,
                u,
                flow.macro_dt,
                flow.macro_steps,
                scfg.substeps.unwrap_or(4).max(1),
            ),
            StrategyKind::ResidualGated => {
                let tol = scfg.residual_tol.unwrap_or(0.01);
                let candidates = scfg.candidates.clone().unwrap_or_else(|| vec![1, 2, 4, 8]);
                let mut last = None;
                let mut chosen = None;
                for sub in candidates {
                    let s = measure_strategy_sample(
                        positions,
                        radius,
                        density,
                        gas,
                        meshc,
                        bounds,
                        eps,
                        gz,
                        mode,
                        mass,
                        m_bed,
                        g,
                        n,
                        u,
                        flow.macro_dt,
                        flow.macro_steps,
                        sub.max(1),
                    );
                    last = Some(s);
                    if s.residual <= tol {
                        chosen = Some(Sample {
                            accepted: true,
                            ..s
                        });
                        break;
                    }
                }
                chosen.unwrap_or_else(|| Sample {
                    accepted: false,
                    ..last.expect("residual-gated strategy requires candidates")
                })
            }
        };
        samples.push(sample);
    }

    let sweep: Vec<(f64, f64)> = samples.iter().map(|s| (s.u, s.a_meas)).collect();
    StrategyOutcome {
        name: scfg.name.clone(),
        kind,
        u_mf_dyn: zero_crossing(&sweep),
        worst_residual: samples.iter().map(|s| s.residual).fold(0.0, f64::max),
        worst_dep: samples.iter().map(|s| s.dep).fold(0.0, f64::max),
        worst_mom: samples.iter().map(|s| s.mom).fold(0.0, f64::max),
        max_substeps: samples.iter().map(|s| s.substeps).max().unwrap_or(0),
        accepted_all: samples.iter().all(|s| s.accepted),
        samples,
    }
}

#[allow(clippy::too_many_arguments)]
fn measure_strategy_sample(
    positions: &[[f64; 3]],
    radius: f64,
    density: f64,
    gas: &GasCfg,
    meshc: &MeshCfg,
    bounds: [f64; 3],
    eps: f64,
    gz: f64,
    mode: SeamMode,
    mass: f64,
    m_bed: f64,
    g: f64,
    n: usize,
    u: f64,
    macro_dt: f64,
    macro_steps: usize,
    substeps: usize,
) -> Sample {
    let dt = macro_dt / substeps as f64;
    let net = rest_net_force(
        positions, radius, density, gas, meshc, bounds, eps, gz, dt, u, mode, mass, g, n,
    );
    let mut parent = make_parent(positions, radius, density, gas, meshc, bounds, eps, gz, dt);
    set_seam_mode(&mut parent, mode);
    set_superficial(&mut parent, [0.0, 0.0, u]);
    for _ in 0..(macro_steps * substeps) {
        parent.run();
    }
    let vz = bed_com_vz(&parent);
    let a_meas = vz / (macro_steps as f64 * macro_dt);
    let a_force = net / m_bed;
    let (_fluid, _drag, dep, mom) = read_result(&parent);
    Sample {
        u,
        a_meas,
        a_force,
        residual: (a_meas - a_force).abs() / g,
        dep,
        mom,
        substeps,
        accepted: true,
    }
}

#[allow(clippy::too_many_arguments)]
fn make_parent(
    positions: &[[f64; 3]],
    radius: f64,
    density: f64,
    gas: &GasCfg,
    meshc: &MeshCfg,
    bounds: [f64; 3],
    eps: f64,
    gz: f64,
    dt: f64,
) -> App {
    let soil = build_soil_bed(positions, radius, density, gz, dt);
    let cfd = build_cfd(gas, meshc.to_uniform(bounds), eps, gz, dt);
    couple_two_way(soil, cfd, radius)
}

#[allow(clippy::too_many_arguments)]
fn rest_net_force(
    positions: &[[f64; 3]],
    radius: f64,
    density: f64,
    gas: &GasCfg,
    meshc: &MeshCfg,
    bounds: [f64; 3],
    eps: f64,
    gz: f64,
    dt: f64,
    u: f64,
    mode: SeamMode,
    mass: f64,
    g: f64,
    n: usize,
) -> f64 {
    let mut parent = make_parent(positions, radius, density, gas, meshc, bounds, eps, gz, dt);
    set_seam_mode(&mut parent, mode);
    set_superficial(&mut parent, [0.0, 0.0, u]);
    export_and_seam(&mut parent);
    let (fluid, ..) = read_result(&parent);
    fluid[2] - mass * g * n as f64
}

/// Step the coupled system once from the current (rest) bed state and record the
/// seam force. `BedResult` is written by the FIELD seam in the `TickCfd` phase —
/// i.e. from the bed kinematics exported at the START of the step (rest) — before
/// `TickSoil` advances the particles, so the recorded force is the rest-state force
/// regardless of the trailing integration (which the caller discards via `reset_bed`).
fn export_and_seam(parent: &mut App) {
    parent.run();
}

/// (Σ total fluid force, Σ drag force, deposition err, momentum err).
fn read_result(parent: &App) -> (Vec3, Vec3, f64, f64) {
    let r = read_subapp_resource::<BedResult>(parent, "cfd");
    (r.f_fluid_total, r.f_drag_total, r.eps_cell_err, r.mom_err)
}

/// Bed centre-of-mass vertical velocity (uniform-mass bed → mean v_z).
fn bed_com_vz(parent: &App) -> f64 {
    let mut sum = 0.0f64;
    let mut n = 0usize;
    with_subapp_resource::<Atom>(parent, "soil", |atoms| {
        n = atoms.nlocal as usize;
        for i in 0..n {
            sum += atoms.vel[i][2] as f64;
        }
    });
    sum / n.max(1) as f64
}

/// Bisect the superficial velocity where the seam's net bed force (from rest)
/// crosses zero, for a given (possibly corrupted) seam mode.
fn bisect_umf(
    positions: &[[f64; 3]],
    radius: f64,
    density: f64,
    gas: &GasCfg,
    meshc: &MeshCfg,
    bounds: [f64; 3],
    eps: f64,
    gz: f64,
    dt: f64,
    mode: SeamMode,
    mass: f64,
    g: f64,
    n: usize,
) -> f64 {
    let (mut lo, mut hi) = (1e-4, 20.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if rest_net_force(
            positions, radius, density, gas, meshc, bounds, eps, gz, dt, mid, mode, mass, g, n,
        ) < 0.0
        {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Linear-interpolated zero crossing of a monotone-increasing (U, a_z) sweep.
fn zero_crossing(sweep: &[(f64, f64)]) -> Option<f64> {
    for w in sweep.windows(2) {
        let (u0, a0) = w[0];
        let (u1, a1) = w[1];
        if a0 <= 0.0 && a1 > 0.0 {
            return Some(u0 + (u1 - u0) * (-a0) / (a1 - a0));
        }
    }
    None
}
