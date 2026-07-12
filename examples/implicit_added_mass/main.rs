//! Added-mass-style strong coupling on the real DEM-CFD seam.
//!
//! This is a deliberately small dense-regime validation case: one SOIL particle
//! and one FIELD sub-App exchange interface velocity and fluid load through the
//! existing `seam.rs` resources. The FIELD force model is the linearized
//! added-mass interface map
//!
//! ```text
//!   v_tilde = a v + b,  with a = -3
//! ```
//!
//! Plain explicit partitioning (`couple_two_way`, omega=1) therefore oscillates
//! and diverges. The identical seam phases converge when the parent driver wraps
//! them in `grass_multi::converge_outer_iter` with Aitken relaxation.

use cfd_ibm::coupling::{InterphaseForces, ParticleSet};
use field_core::{MeshScheduleSet, UniformMeshConfig};
use grass_app::prelude::*;
use grass_io::Config;
use grass_multi::{
    converge_outer_iter, tick_subapp, Multi, MultiAppExt, OuterIteration, Relaxation, SubApps,
};
use grass_scheduler::prelude::*;
use grass_scheduler::{Res, ResMut};
use serde::Deserialize;
use soil_core::Atom;

use dem_cfd::prelude::*;

#[derive(Deserialize, Default)]
struct ParticleCfg {
    radius: f64,
    density: f64,
    z0: f64,
}

#[derive(Deserialize, Default)]
struct CouplingCfg {
    dt: f64,
    map_slope: f64,
    map_intercept: f64,
    explicit_steps: usize,
    implicit_max_iters: u32,
    implicit_tol: f64,
    aitken_omega0: f64,
}

#[derive(Deserialize, Default)]
struct ValidationCfg {
    explicit_min_growth: f64,
    implicit_solution_tol: f64,
    implicit_residual_tol: f64,
    max_implicit_iters: u32,
}

#[derive(Clone, Copy)]
struct AddedMassMap {
    slope: f64,
    intercept: f64,
    mass: f64,
    dt: f64,
}

#[derive(Clone, Copy, Default)]
struct InterfaceTrace {
    v_in: f64,
    v_tilde: f64,
}

fn added_mass_force_system(
    map: Res<AddedMassMap>,
    pset: Res<ParticleSet>,
    mut forces: ResMut<InterphaseForces>,
    mut trace: ResMut<InterfaceTrace>,
) {
    forces.reset(pset.particles.len());
    let Some(p) = pset.particles.first() else {
        return;
    };
    let v = p.velocity[2];
    let v_tilde = map.slope * v + map.intercept;
    let fz = map.mass * (v_tilde - v) / map.dt;
    forces.force[0] = [0.0, 0.0, fz];
    *trace = InterfaceTrace { v_in: v, v_tilde };
}

fn build_cfd(gas: &GasCfg, mesh_cfg: UniformMeshConfig, map: AddedMassMap) -> App {
    let ctx = SeamCtx {
        mu: gas.mu,
        rho: gas.rho,
        eps: 1.0,
        g: [0.0, 0.0, 0.0],
        dt: map.dt,
        mode: SeamMode::default(),
    };
    let mut app = build_cfd_base(gas, mesh_cfg, ctx);
    app.add_resource(map);
    app.add_resource(InterfaceTrace::default());
    app.add_update_system(added_mass_force_system, MeshScheduleSet::Output);
    app
}

#[derive(Clone, Copy, Debug, ScheduleSet)]
enum ImplicitPhase {
    Export,
    TickCfd,
    Import,
    TickSoil,
    Converge,
}

fn coupled_parent(
    gas: &GasCfg,
    mesh_cfg: UniformMeshConfig,
    particle: &ParticleCfg,
    c: &CouplingCfg,
    implicit: bool,
) -> App {
    let mass = particle.density * 4.0 / 3.0 * std::f64::consts::PI * particle.radius.powi(3);
    let map = AddedMassMap {
        slope: c.map_slope,
        intercept: c.map_intercept,
        mass,
        dt: c.dt,
    };
    let soil = build_soil_bed(
        &[[0.5, 0.5, particle.z0]],
        particle.radius,
        particle.density,
        0.0,
        c.dt,
    );
    let cfd = build_cfd(gas, mesh_cfg, map);
    if !implicit {
        return couple_two_way(soil, cfd, particle.radius);
    }

    let mut parent = App::new();
    parent.add_subapp("soil", soil);
    parent.add_subapp("cfd", cfd);
    parent.add_resource(ParticleSpec {
        radius: particle.radius,
    });
    parent.add_resource(OuterIteration::new(
        vec![0.0],
        Relaxation::Aitken {
            omega0: c.aitken_omega0,
        },
        c.implicit_tol,
        c.implicit_max_iters,
    ));
    parent.add_update_system(export_kinematics, ImplicitPhase::Export);
    parent.add_update_system(tick_subapp("cfd", 1), ImplicitPhase::TickCfd);
    parent.add_update_system(import_force, ImplicitPhase::Import);
    parent.add_update_system(tick_subapp("soil", 1), ImplicitPhase::TickSoil);
    parent.add_update_system(
        converge_outer_iter(
            |w: &Multi| vec![w.expect_read::<Atom>("soil").vel[0][2] as f64],
            |w: &Multi, x: &[f64]| {
                let mut atoms = w.expect_write::<Atom>("soil");
                atoms.vel[0][2] = x[0] as _;
                atoms.force[0] = [0.0, 0.0, 0.0];
            },
        ),
        ImplicitPhase::Converge,
    );
    parent.prepare();
    parent
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: implicit_added_mass <case.toml>");
    let toml_src =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let cfg = Config::from_str(&toml_src);

    let gas: GasCfg = cfg.section("gas");
    let particle: ParticleCfg = cfg.section("particle");
    let mesh_cfg: UniformMeshConfig = cfg.section("grid");
    let coupling: CouplingCfg = cfg.section("coupling");
    let valid: ValidationCfg = cfg.section("validation");

    let x_star = coupling.map_intercept / (1.0 - coupling.map_slope);

    println!("# Added-mass strong-coupling seam validation");
    println!(
        "# interface map: v_tilde = a v + b, a = {:.3}, b = {:.3}; fixed point v* = {:.12}",
        coupling.map_slope, coupling.map_intercept, x_star
    );
    println!("# explicit path: dem_cfd::seam::couple_two_way (omega = 1 / plain Picard)");
    println!("# implicit path: same Export/TickCfd/Import/TickSoil seam phases + grass_multi::converge_outer_iter(Aitken)");
    println!("#");
    println!("# explicit trace");
    println!("# step  v_in          v_tilde      v_after      residual");

    let mut explicit = coupled_parent(&gas, mesh_cfg.clone(), &particle, &coupling, false);
    let mut explicit_rows = Vec::new();
    let mut first_res = 0.0;
    let mut last_res = 0.0;
    for step in 0..coupling.explicit_steps {
        explicit.run();
        let tr = read_subapp_resource::<InterfaceTrace>(&explicit, "cfd");
        let v_after = read_particle_vz(&explicit);
        let res = (tr.v_tilde - tr.v_in).abs();
        if step == 0 {
            first_res = res;
        }
        last_res = res;
        explicit_rows.push((step, tr.v_in, tr.v_tilde, v_after, res));
        println!(
            "{step:>6}  {:>12.6}  {:>12.6}  {:>12.6}  {:>12.6}",
            tr.v_in, tr.v_tilde, v_after, res
        );
    }
    cleanup(&mut explicit);

    println!("#");
    println!("# implicit Aitken trace");
    println!("# iter  v_in          v_tilde      v_relaxed    residual     omega");
    let mut implicit = coupled_parent(&gas, mesh_cfg, &particle, &coupling, true);
    implicit.start();
    let (implicit_iters, implicit_res, implicit_converged, implicit_omega) = {
        let it = implicit.get_resource_ref::<OuterIteration>().unwrap();
        (it.iters(), it.residual_norm(), it.converged(), it.omega())
    };
    let tr = read_subapp_resource::<InterfaceTrace>(&implicit, "cfd");
    let v_final = read_particle_vz(&implicit);
    println!(
        "{:>6}  {:>12.6}  {:>12.6}  {:>12.6}  {:>12.3e}  {:>8.4}",
        implicit_iters, tr.v_in, tr.v_tilde, v_final, implicit_res, implicit_omega
    );
    cleanup(&mut implicit);

    let growth = last_res / first_res.max(1e-30);
    let explicit_diverged = growth >= valid.explicit_min_growth && last_res > first_res;
    let implicit_solution_err = (v_final - x_star).abs();
    let implicit_ok = implicit_converged
        && implicit_solution_err <= valid.implicit_solution_tol
        && implicit_res <= valid.implicit_residual_tol
        && implicit_iters <= valid.max_implicit_iters;

    println!("#");
    println!("# result");
    println!(
        "# explicit residual growth: {:.3e} / {:.3e} = {:.3e} (required >= {:.3e})",
        last_res, first_res, growth, valid.explicit_min_growth
    );
    println!(
        "# implicit v_final: {:.12} vs analytic fixed point {:.12}; abs err {:.3e}",
        v_final, x_star, implicit_solution_err
    );
    println!(
        "# implicit residual: {:.3e}; iters {}; converged {}",
        implicit_res, implicit_iters, implicit_converged
    );

    if explicit_diverged && implicit_ok {
        println!(
            "VALIDATION: PASS  (explicit diverged {:.2e}x; Aitken converged in {} iters to {:.12} with residual {:.1e})",
            growth, implicit_iters, v_final, implicit_res
        );
    } else {
        println!(
            "VALIDATION: FAIL  (explicit_diverged={explicit_diverged} implicit_ok={implicit_ok})"
        );
        std::process::exit(1);
    }

    write_csv(
        &explicit_rows,
        v_final,
        x_star,
        implicit_iters,
        implicit_res,
    );
}

fn read_particle_vz(parent: &App) -> f64 {
    let subs = parent.get_resource_ref::<SubApps>().unwrap();
    let soil = subs.find("soil").unwrap();
    let cell = soil
        .resource_cell(std::any::TypeId::of::<Atom>())
        .unwrap()
        .borrow();
    cell.downcast_ref::<Atom>().unwrap().vel[0][2] as f64
}

fn cleanup(parent: &mut App) {
    if let Some(cell) = parent.get_mut_resource(std::any::TypeId::of::<SubApps>()) {
        cell.borrow_mut()
            .downcast_mut::<SubApps>()
            .unwrap()
            .cleanup_all();
    }
}

fn write_csv(
    rows: &[(usize, f64, f64, f64, f64)],
    v_final: f64,
    x_star: f64,
    iters: u32,
    res: f64,
) {
    let mut out = String::from("kind,step,v_in,v_tilde,v_after,residual\n");
    for (step, v_in, v_tilde, v_after, residual) in rows {
        out.push_str(&format!(
            "explicit,{step},{v_in:.12e},{v_tilde:.12e},{v_after:.12e},{residual:.12e}\n"
        ));
    }
    out.push_str(&format!(
        "implicit,{iters},{v_final:.12e},{x_star:.12e},{v_final:.12e},{res:.12e}\n"
    ));
    std::fs::create_dir_all("examples/implicit_added_mass/plots").unwrap();
    std::fs::write("examples/implicit_added_mass/plots/trace.csv", out).unwrap();
}
