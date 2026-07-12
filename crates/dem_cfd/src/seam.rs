//! The `grass_multi` seam scaffold: the resources that cross the SOIL↔FIELD
//! namespace boundary, the CFD sub-App base builder, the dynamic two-way coupling
//! systems + parent schedule, and small accessors for poking sub-App resources
//! from the driver. A case supplies only its own force model (the seam *system*)
//! and, if its topology differs (e.g. a static packed bed), its own schedule.

use std::any::TypeId;

use cfd_eos::{Eos, EosResource, IdealGas, Viscosity};
use cfd_ibm::coupling::{
    self, beta_gidaspow, cd_schiller_naumann, deposit_void_fraction, particle_reynolds,
    InterphaseForces, ParticleKinematics, ParticleSet,
};
use cfd_solver::{CfdStatePlugin, IdealGasPlugin};
use cfd_state::{CfdState, PrimVar};
use field_core::{
    FieldDefaultPlugins, FieldRegistry, MeshScheduleSet, UniformMesh, UniformMeshConfig, Vec3,
};
use grass_app::prelude::*;
use grass_multi::{
    namespace, tick_n_times, tick_subapp, Multi, MultiAppExt, MultiRes, MultiResMut,
    OuterIterStopPlugin, SubApps,
};
use grass_scheduler::prelude::*;
use grass_scheduler::{Res, ResMut};
use soil_core::{Accum, Atom, ParticleSimScheduleSet};
use soil_verlet::VelocityVerletPlugin;

use crate::config::GasCfg;
use crate::drag::SeamMode;

/// Matches `IdealGas::air()`.
pub const R_GAS: f64 = 287.058;

// ─── Seam resources shared across the namespace boundary ─────────────────────

/// Per-particle total fluid force (drag + optional ∇P buoyancy + hydrostatic),
/// FIELD→SOIL.
#[derive(Clone, Debug, Default)]
pub struct FluidForces {
    pub f: Vec<Vec3>,
}

/// Body acceleration (gravity) applied to the bed atoms on the SOIL side.
#[derive(Clone, Copy)]
pub struct BodyAccel {
    pub g: Vec3,
}

/// Imposed superficial velocity for the current pass (world axes), set by the
/// driver; the FIELD side turns it into the interstitial gas velocity `U/ε`.
#[derive(Clone, Copy, Default)]
pub struct Superficial {
    pub u: Vec3,
}

/// Radius handed to the exported particle kinematics.
#[derive(Clone, Copy)]
pub struct ParticleSpec {
    pub radius: f64,
}

/// Static seam context on the FIELD sub-App: gas transport, bed porosity ε (for
/// the superficial↔interstitial conversion `u_g = U/ε`), gravity (hydrostatic
/// buoyancy), the coupling timestep, and the current [`SeamMode`]. Bundled so a
/// seam system stays within the scheduler's parameter-count limit.
#[derive(Clone, Copy)]
pub struct SeamCtx {
    pub mu: f64,
    pub rho: f64,
    /// Bed porosity ε = 1 − Σ V_p / V_bed.
    pub eps: f64,
    /// Gravity (for the hydrostatic buoyancy term); `[0,0,0]` for a static bed.
    pub g: Vec3,
    pub dt: f64,
    pub mode: SeamMode,
}

namespace!(pub DemNs = "dem");
namespace!(pub CfdNs = "cfd");

// ─── FIELD sub-App base ──────────────────────────────────────────────────────

/// Build the CFD sub-App with the standard plugins (uniform mesh, ideal-gas state,
/// constant viscosity) and the seam resources every case reads/writes
/// (`SeamCtx`, `Superficial`, `ParticleSet`, `InterphaseForces`). The caller adds
/// its own seam *system* (the force model) and result resource, e.g.
/// `app.add_resource(MyResult::default()); app.add_update_system(my_seam, MeshScheduleSet::Output);`.
pub fn build_cfd_base(gas: &GasCfg, mesh_cfg: UniformMeshConfig, ctx: SeamCtx) -> App {
    let (rho, p) = (gas.rho, gas.p);
    let t = p / (rho * R_GAS);
    let init = move |_x: Vec3| {
        let eos = IdealGas::air();
        eos.prim_to_cons(&PrimVar::new(rho, 0.0, 0.0, 0.0, p, t))
    };
    let mut app = App::new();
    app.add_plugins(FieldDefaultPlugins { mesh: mesh_cfg })
        .add_plugins(CfdStatePlugin::new(init))
        .add_plugins(IdealGasPlugin);
    app.add_resource(Viscosity::Constant(gas.mu));
    app.add_resource(ctx);
    app.add_resource(Superficial::default());
    app.add_resource(ParticleSet::default());
    app.add_resource(InterphaseForces::default());
    app
}

// ─── SOIL sub-App: a freely-moving (integrated) bed ──────────────────────────

/// `Force` phase: force = m·gravity + the per-particle fluid load from the seam.
// Casts bridge soil's precision-generic `Real` to the CFD side's fixed `f64`;
// they are no-ops only under the `precision-double` build this coupling pins.
#[allow(clippy::unnecessary_cast)]
pub fn bed_force(mut atoms: ResMut<Atom>, ff: Res<FluidForces>, body: Res<BodyAccel>) {
    let n = atoms.nlocal as usize;
    for i in 0..n {
        let m = atoms.mass[i] as f64;
        let f = ff.f.get(i).copied().unwrap_or([0.0; 3]);
        atoms.force[i] = [
            (m * body.g[0] + f[0]) as _,
            (m * body.g[1] + f[1]) as _,
            (m * body.g[2] + f[2]) as _,
        ];
    }
}

/// Add the imported fluid force to whatever forces the standalone DEM solver
/// already computes (gravity, contacts, bonds, constraints, ...).
pub fn add_fluid_force(mut atoms: ResMut<Atom>, ff: Res<FluidForces>) {
    for i in 0..atoms.nlocal as usize {
        let f = ff.f.get(i).copied().unwrap_or([0.0; 3]);
        atoms.force[i][0] += f[0] as Accum;
        atoms.force[i][1] += f[1] as Accum;
        atoms.force[i][2] += f[2] as Accum;
    }
}

/// Point-particle Wen–Yu/Gidaspow exchange installed into the CFD sub-App by
/// [`DemCfdCouplingPlugin`]. It deposits void fraction, samples local gas state,
/// computes drag and buoyancy, and applies the equal-and-opposite gas momentum
/// sink.
#[allow(clippy::too_many_arguments)]
pub fn point_particle_exchange(
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
        let gas_velocity =
            coupling::sample_gas_velocity(&*mesh, &state, eos, particle.center).unwrap_or([0.0; 3]);
        let rho = coupling::sample_gas_density(&*mesh, &state, particle.center).unwrap_or(ctx.rho);
        let eps =
            coupling::sample_void_fraction(&*mesh, &void_fraction, particle.center).unwrap_or(1.0);
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

/// Build an integrated (velocity-Verlet) SOIL bed of equal spheres at `positions`.
pub fn build_soil_bed(positions: &[[f64; 3]], radius: f64, density: f64, gz: f64, dt: f64) -> App {
    let mut atoms = Atom::new();
    atoms.dt = dt;
    let mass = density * 4.0 / 3.0 * std::f64::consts::PI * radius.powi(3);
    for (tag, pos) in positions.iter().enumerate() {
        atoms.push_test_atom(tag as u32, *pos, radius, mass);
    }
    atoms.nlocal = positions.len() as u32;
    atoms.natoms = positions.len() as u64;

    let mut app = App::new();
    app.add_resource(atoms);
    app.add_resource(FluidForces::default());
    app.add_resource(BodyAccel { g: [0.0, 0.0, gz] });
    app.add_update_system(bed_force, ParticleSimScheduleSet::Force);
    app.add_plugins(VelocityVerletPlugin::new());
    app
}

// ─── Dynamic two-way coupling schedule (moving bed) ──────────────────────────

/// The parent phases for a two-way coupled moving bed.
#[derive(Debug, Clone, Copy)]
pub enum CouplePhase {
    Export,
    TickCfd,
    Import,
    TickSoil,
    Check,
}
impl ScheduleSet for CouplePhase {
    fn to_index(&self) -> u32 {
        match self {
            Self::Export => 0,
            Self::TickCfd => 1,
            Self::Import => 2,
            Self::TickSoil => 3,
            Self::Check => 4,
        }
    }
    fn name(&self) -> &'static str {
        match self {
            Self::Export => "Export",
            Self::TickCfd => "TickCfd",
            Self::Import => "Import",
            Self::TickSoil => "TickSoil",
            Self::Check => "Check",
        }
    }
}

/// SOIL→FIELD: hand the (moving) bed kinematics across each step.
// Casts bridge soil's precision-generic `Real` to the CFD side's fixed `f64`;
// they are no-ops only under the `precision-double` build this coupling pins.
#[allow(clippy::unnecessary_cast)]
pub fn export_kinematics(world: Multi, spec: Res<ParticleSpec>) {
    let atoms = world.expect_read::<Atom>("dem");
    let n = atoms.nlocal as usize;
    let mut set = world.expect_write::<ParticleSet>("cfd");
    set.particles.clear();
    for i in 0..n {
        set.particles.push(ParticleKinematics {
            center: [
                atoms.pos[i][0] as f64,
                atoms.pos[i][1] as f64,
                atoms.pos[i][2] as f64,
            ],
            velocity: [
                atoms.vel[i][0] as f64,
                atoms.vel[i][1] as f64,
                atoms.vel[i][2] as f64,
            ],
            radius: spec.radius,
        });
    }
}

/// FIELD→SOIL: copy the per-particle total fluid force back to the bed.
pub fn import_force(world: Multi) {
    let forces = world.expect_read::<InterphaseForces>("cfd");
    let v = forces.force.clone();
    drop(forces);
    world.expect_write::<FluidForces>("dem").f = v;
}

/// Typed parent-owned FIELD→SOIL force handoff.
///
/// The explicit mesh schedule ends with [`MeshScheduleSet::Output`], so running
/// this system in the parent's [`CouplePhase::Import`] phase immediately after
/// the typed CFD tick observes the same completed force result as a child
/// adapter in `Output`. Keeping the handoff on the parent makes the coupling
/// order explicit and avoids cross-App access from inside a child scheduler.
pub fn import_force_typed(
    forces: MultiRes<InterphaseForces, CfdNs>,
    mut fluid_forces: MultiResMut<FluidForces, DemNs>,
) {
    fluid_forces.f.clone_from(&forces.force);
}

/// The standard dynamic unresolved DEM↔CFD coupling loop.
///
/// Add independently configured DEM and CFD solvers as the `"dem"` and `"cfd"`
/// sub-Apps, then add this plugin to the parent. Neither solver needs coupling
/// code. Before either is prepared, the plugin installs the seam resources and
/// systems into both Apps. It then owns particle export, Wen–Yu/Gidaspow force
/// exchange, the equal-and-opposite gas momentum sink, fluid-force addition on
/// the DEM side, solver advancement, termination, and cleanup.
pub struct DemCfdCouplingPlugin {
    /// Radius exported with every DEM particle.
    pub particle_radius: f64,
    /// Number of coupled outer steps executed by [`App::start`].
    pub steps: u32,
    /// Carrier density used if a particle lies outside the sampled mesh.
    pub gas_density: f64,
    /// Dynamic viscosity used by the point-particle closure.
    pub gas_viscosity: f64,
    /// Body acceleration used for generalized buoyancy.
    pub gravity: Vec3,
    /// Coupling timestep used by the equal-and-opposite gas momentum sink.
    pub dt: f64,
}

impl DemCfdCouplingPlugin {
    /// Construct a coupling loop for equal-radius particles.
    pub fn new(
        particle_radius: f64,
        steps: u32,
        dt: f64,
        gas_density: f64,
        gas_viscosity: f64,
        gravity: Vec3,
    ) -> Self {
        assert!(particle_radius > 0.0, "particle radius must be positive");
        assert!(steps > 0, "coupling steps must be positive");
        assert!(dt > 0.0, "coupling timestep must be positive");
        assert!(gas_density > 0.0, "gas density must be positive");
        assert!(gas_viscosity > 0.0, "gas viscosity must be positive");
        Self {
            particle_radius,
            steps,
            gas_density,
            gas_viscosity,
            gravity,
            dt,
        }
    }

    /// Convenience constructor for the teaching case's room-temperature air
    /// (`rho = 1.2 kg/m³`, `mu = 1.8e-5 Pa·s`). Production cases should use
    /// [`new`](Self::new) with properties matching their CFD state.
    pub fn for_air(particle_radius: f64, steps: u32, dt: f64, gravity: Vec3) -> Self {
        Self::new(particle_radius, steps, dt, 1.2, 1.8e-5, gravity)
    }
}

impl Plugin for DemCfdCouplingPlugin {
    fn build(&self, app: &mut App) {
        let ctx = SeamCtx {
            mu: self.gas_viscosity,
            rho: self.gas_density,
            eps: 1.0,
            g: self.gravity,
            dt: self.dt,
            mode: SeamMode::default(),
        };
        app.configure_subapp("dem", |dem| {
            dem.add_resource(FluidForces::default());
            dem.add_update_system(add_fluid_force, ParticleSimScheduleSet::Force);
        });
        app.configure_subapp("cfd", |cfd| {
            cfd.add_resource(ctx);
            cfd.add_resource(ParticleSet::default());
            cfd.add_resource(InterphaseForces::default());
            cfd.add_update_system(point_particle_exchange, MeshScheduleSet::Output);
        });
        app.add_resource(ParticleSpec {
            radius: self.particle_radius,
        });
        app.add_update_system(export_kinematics, CouplePhase::Export);
        app.add_update_system(tick_n_times::<CfdNs>(1), CouplePhase::TickCfd);
        app.add_update_system(import_force_typed, CouplePhase::Import);
        app.add_update_system(tick_n_times::<DemNs>(1), CouplePhase::TickSoil);
        app.add_plugins(OuterIterStopPlugin {
            n_iters: self.steps,
            phase: CouplePhase::Check,
        });
        app.add_cleanup_with_app(|parent| {
            if let Some(cell) = parent.get_mut_resource(TypeId::of::<SubApps>()) {
                cell.borrow_mut()
                    .downcast_mut::<SubApps>()
                    .expect("SubApps resource type")
                    .cleanup_all();
            }
        });
    }
}

/// Assemble the parent App for a dynamic two-way coupling: mount the two sub-Apps
/// and wire the four-phase [`CouplePhase`] schedule (export → gas → import → bed).
/// Returns the prepared parent, ready to `update()`.
pub fn couple_two_way(soil: App, cfd: App, radius: f64) -> App {
    let mut parent = App::new();
    parent.add_subapp("dem", soil);
    parent.add_subapp("cfd", cfd);
    // Manual drivers own their loop and cleanup, so wire only the four phases
    // here. New self-driving examples should prefer DemCfdCouplingPlugin.
    parent.add_resource(ParticleSpec { radius });
    parent.add_update_system(export_kinematics, CouplePhase::Export);
    parent.add_update_system(tick_subapp("cfd", 1), CouplePhase::TickCfd);
    parent.add_update_system(import_force, CouplePhase::Import);
    parent.add_update_system(tick_subapp("dem", 1), CouplePhase::TickSoil);
    parent.prepare();
    parent
}

// ─── Driver-side access to a sub-App resource (outside any system) ────────────

/// Mutate a resource of type `T` on the named sub-App from the driver.
pub fn with_subapp_resource<T: 'static>(parent: &App, sub: &str, f: impl FnOnce(&mut T)) {
    let subs = parent.get_resource_ref::<SubApps>().unwrap();
    let participant = subs.find(sub).unwrap();
    let cell = participant.resource_cell(TypeId::of::<T>()).unwrap();
    f(cell.borrow_mut().downcast_mut::<T>().unwrap());
}

/// Read a `Copy` resource of type `T` from the named sub-App.
pub fn read_subapp_resource<T: Copy + 'static>(parent: &App, sub: &str) -> T {
    let subs = parent.get_resource_ref::<SubApps>().unwrap();
    let participant = subs.find(sub).unwrap();
    let cell = participant
        .resource_cell(TypeId::of::<T>())
        .unwrap()
        .borrow();
    *cell.downcast_ref::<T>().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parent_imports_force_after_typed_cfd_tick() {
        let expected = vec![[1.0, -2.0, 3.5], [4.0, 5.0, -6.0]];

        fn publish_in_output(mut forces: ResMut<InterphaseForces>) {
            forces.force = vec![[1.0, -2.0, 3.5], [4.0, 5.0, -6.0]];
        }

        let mut dem = App::new();
        dem.add_resource(FluidForces::default());

        let mut cfd = App::new();
        cfd.add_resource(InterphaseForces::default());
        cfd.add_update_system(publish_in_output, MeshScheduleSet::Output);

        let mut parent = App::new();
        parent.add_subapp_typed::<DemNs>(dem);
        parent.add_subapp_typed::<CfdNs>(cfd);
        parent.add_update_system(tick_n_times::<CfdNs>(1), CouplePhase::TickCfd);
        parent.add_update_system(import_force_typed, CouplePhase::Import);
        parent.prepare();
        parent.run();

        let subs = parent.get_resource_ref::<SubApps>().unwrap();
        let dem = subs.find("dem").unwrap();
        let forces = dem
            .resource_cell(TypeId::of::<FluidForces>())
            .unwrap()
            .borrow();
        assert_eq!(forces.downcast_ref::<FluidForces>().unwrap().f, expected);
    }
}

/// Set the imposed superficial velocity on the `cfd` sub-App.
pub fn set_superficial(parent: &App, u: Vec3) {
    with_subapp_resource::<Superficial>(parent, "cfd", |s| s.u = u);
}

/// Set the seam mode on the `cfd` sub-App.
pub fn set_seam_mode(parent: &App, mode: SeamMode) {
    with_subapp_resource::<SeamCtx>(parent, "cfd", |c| c.mode = mode);
}
