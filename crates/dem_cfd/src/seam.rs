//! The `grass_multi` seam scaffold: the resources that cross the SOIL↔FIELD
//! namespace boundary, the CFD sub-App base builder, the dynamic two-way coupling
//! systems + parent schedule, and small accessors for poking sub-App resources
//! from the driver. A case supplies only its own force model (the seam *system*)
//! and, if its topology differs (e.g. a static packed bed), its own schedule.

use std::any::TypeId;

use cfd_eos::{Eos, IdealGas, Viscosity};
use cfd_ibm::coupling::{InterphaseForces, ParticleKinematics, ParticleSet};
use cfd_solver::{CfdStatePlugin, IdealGasPlugin};
use cfd_state::PrimVar;
use field_core::{
    FieldDefaultPlugins, UniformMeshConfig, Vec3,
};
use grass_app::prelude::*;
use grass_multi::{tick_subapp, Multi, MultiAppExt, SubApps};
use grass_scheduler::prelude::*;
use grass_scheduler::{Res, ResMut};
use soil_core::{Atom, ParticleSimScheduleSet};
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

/// The four-phase schedule for a two-way coupled moving bed: export kinematics →
/// tick the gas (and its seam system) → import the fluid force → tick the bed.
#[derive(Debug, Clone, Copy)]
pub enum CouplePhase {
    Export,
    TickCfd,
    Import,
    TickSoil,
}
impl ScheduleSet for CouplePhase {
    fn to_index(&self) -> u32 {
        match self {
            Self::Export => 0,
            Self::TickCfd => 1,
            Self::Import => 2,
            Self::TickSoil => 3,
        }
    }
    fn name(&self) -> &'static str {
        match self {
            Self::Export => "Export",
            Self::TickCfd => "TickCfd",
            Self::Import => "Import",
            Self::TickSoil => "TickSoil",
        }
    }
}

/// SOIL→FIELD: hand the (moving) bed kinematics across each step.
pub fn export_kinematics(world: Multi, spec: Res<ParticleSpec>) {
    let atoms = world.expect_read::<Atom>("soil");
    let n = atoms.nlocal as usize;
    let mut set = world.expect_write::<ParticleSet>("cfd");
    set.particles.clear();
    for i in 0..n {
        set.particles.push(ParticleKinematics {
            center: [atoms.pos[i][0] as f64, atoms.pos[i][1] as f64, atoms.pos[i][2] as f64],
            velocity: [atoms.vel[i][0] as f64, atoms.vel[i][1] as f64, atoms.vel[i][2] as f64],
            radius: spec.radius,
        });
    }
}

/// FIELD→SOIL: copy the per-particle total fluid force back to the bed.
pub fn import_force(world: Multi) {
    let forces = world.expect_read::<InterphaseForces>("cfd");
    let v = forces.force.clone();
    drop(forces);
    world.expect_write::<FluidForces>("soil").f = v;
}

/// Assemble the parent App for a dynamic two-way coupling: mount the two sub-Apps
/// and wire the four-phase [`CouplePhase`] schedule (export → gas → import → bed).
/// Returns the prepared parent, ready to `update()`.
pub fn couple_two_way(soil: App, cfd: App, radius: f64) -> App {
    let mut parent = App::new();
    parent.add_subapp("soil", soil);
    parent.add_subapp("cfd", cfd);
    parent.add_resource(ParticleSpec { radius });
    parent.add_update_system(export_kinematics, CouplePhase::Export);
    parent.add_update_system(tick_subapp("cfd", 1), CouplePhase::TickCfd);
    parent.add_update_system(import_force, CouplePhase::Import);
    parent.add_update_system(tick_subapp("soil", 1), CouplePhase::TickSoil);
    parent.prepare();
    parent
}

// ─── Driver-side access to a sub-App resource (outside any system) ────────────

/// Mutate a resource of type `T` on the named sub-App from the driver.
pub fn with_subapp_resource<T: 'static>(parent: &App, sub: &str, f: impl FnOnce(&mut T)) {
    let subs = parent.get_resource_ref::<SubApps>().unwrap();
    let cell = subs.find(sub).unwrap().resource_cell(TypeId::of::<T>()).unwrap();
    f(cell.borrow_mut().downcast_mut::<T>().unwrap());
}

/// Read a `Copy` resource of type `T` from the named sub-App.
pub fn read_subapp_resource<T: Copy + 'static>(parent: &App, sub: &str) -> T {
    let subs = parent.get_resource_ref::<SubApps>().unwrap();
    let cell = subs.find(sub).unwrap().resource_cell(TypeId::of::<T>()).unwrap().borrow();
    *cell.downcast_ref::<T>().unwrap()
}

/// Set the imposed superficial velocity on the `cfd` sub-App.
pub fn set_superficial(parent: &App, u: Vec3) {
    with_subapp_resource::<Superficial>(parent, "cfd", |s| s.u = u);
}

/// Set the seam mode on the `cfd` sub-App.
pub fn set_seam_mode(parent: &App, mode: SeamMode) {
    with_subapp_resource::<SeamCtx>(parent, "cfd", |c| c.mode = mode);
}
