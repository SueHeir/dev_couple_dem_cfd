//! Concrete solver setup for the short lesson in `main.rs`.
//!
//! This file is intentionally the second page of the lesson: `main.rs` teaches
//! composition first; here the curious reader can inspect the actual
//! Wen–Yu/Gidaspow point-particle force exchange.

use std::any::TypeId;

use cfd_eos::{Eos, EosResource};
use cfd_ibm::coupling::{
    self, beta_gidaspow, cd_schiller_naumann, deposit_void_fraction, particle_reynolds,
    InterphaseForces, ParticleSet,
};
use cfd_state::CfdState;
use dem_cfd::config::GasCfg;
use dem_cfd::drag::SeamMode;
use dem_cfd::seam::{build_cfd_base, build_soil_bed, SeamCtx};
use field_core::{FieldRegistry, MeshScheduleSet, UniformMesh, UniformMeshConfig};
use grass_app::prelude::*;
use grass_multi::SubApps;
use grass_scheduler::{Res, ResMut};
use soil_core::Atom;

pub const RADIUS: f64 = 5.0e-4;
const DT: f64 = 1.0e-4;
const GAS_DENSITY: f64 = 1.2;
const GAS_VISCOSITY: f64 = 1.8e-5;
const GRAVITY: f64 = -9.81;

pub fn dem() -> App {
    build_soil_bed(&[[0.005, 0.005, 0.015]], RADIUS, 2_500.0, GRAVITY, DT)
}

pub fn cfd() -> App {
    let gas = GasCfg {
        rho: GAS_DENSITY,
        mu: GAS_VISCOSITY,
        p: 101_325.0,
    };
    let mesh = UniformMeshConfig {
        nx: 4,
        ny: 4,
        nz: 8,
        ng: 2,
        bounds_lo: [0.0, 0.0, 0.0],
        bounds_hi: [0.01, 0.01, 0.02],
        y_edges: None,
        z_edges: None,
    };
    let ctx = SeamCtx {
        mu: GAS_VISCOSITY,
        rho: GAS_DENSITY,
        eps: 1.0,
        g: [0.0, 0.0, GRAVITY],
        dt: DT,
        mode: SeamMode::default(),
    };
    let mut app = build_cfd_base(&gas, mesh, ctx);
    app.add_update_system(point_particle_drag, MeshScheduleSet::Output);
    app
}

#[allow(clippy::too_many_arguments)]
fn point_particle_drag(
    mesh: Res<UniformMesh>,
    reg: Res<FieldRegistry>,
    eos: Res<EosResource>,
    ctx: Res<SeamCtx>,
    particles: Res<ParticleSet>,
    mut forces: ResMut<InterphaseForces>,
) {
    let eos: &dyn Eos = &*eos.0;
    let mut state = reg.expect_mut::<CfdState>("CfdState not registered");
    forces.reset(particles.particles.len());
    let void_fraction = deposit_void_fraction(&*mesh, &particles.particles, 1e-3);
    let mut drag_reaction = vec![[0.0; 3]; particles.particles.len()];

    for (i, particle) in particles.particles.iter().enumerate() {
        let gas_velocity = coupling::sample_gas_velocity(&*mesh, &state, eos, particle.center)
            .unwrap_or([0.0; 3]);
        let rho = coupling::sample_gas_density(&*mesh, &state, particle.center).unwrap_or(ctx.rho);
        let eps = coupling::sample_void_fraction(&*mesh, &void_fraction, particle.center)
            .unwrap_or(1.0);
        let slip = [
            gas_velocity[0] - particle.velocity[0],
            gas_velocity[1] - particle.velocity[1],
            gas_velocity[2] - particle.velocity[2],
        ];
        let speed = slip.iter().map(|x| x * x).sum::<f64>().sqrt();
        let reynolds = particle_reynolds(rho, speed, particle.diameter(), ctx.mu) * eps;
        let beta = beta_gidaspow(
            eps,
            rho,
            ctx.mu,
            particle.diameter(),
            speed,
            cd_schiller_naumann(reynolds.max(1e-12)),
        );
        let drag = coupling::drag_force_from_beta(beta, particle.volume(), eps, slip);
        let buoyancy = [
            -rho * particle.volume() * ctx.g[0],
            -rho * particle.volume() * ctx.g[1],
            -rho * particle.volume() * ctx.g[2],
        ];
        forces.force[i] = [
            drag[0] + buoyancy[0],
            drag[1] + buoyancy[1],
            drag[2] + buoyancy[2],
        ];
        drag_reaction[i] = drag;
    }

    coupling::apply_momentum_sink(
        &*mesh,
        &mut state,
        &particles.particles,
        &drag_reaction,
        ctx.dt,
    );
}

pub fn particle_height(parent: &App) -> f64 {
    let subapps = parent
        .get_resource_ref::<SubApps>()
        .expect("coupling sub-Apps");
    let atoms = subapps
        .find("dem")
        .expect("dem sub-App")
        .resource_cell(TypeId::of::<Atom>())
        .expect("DEM particles")
        .borrow();
    atoms.downcast_ref::<Atom>().expect("Atom resource").pos[0][2] as f64
}
