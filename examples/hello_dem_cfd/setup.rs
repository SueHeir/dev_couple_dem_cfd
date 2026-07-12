//! Concrete solver setup for the short lesson in `main.rs`.
//!
//! These are ordinary standalone solver Apps. There are no coupling resources,
//! exchange systems, or references to the other solver in this file.

use std::any::TypeId;

use cfd_eos::{Eos, IdealGas, Viscosity};
use cfd_solver::{CfdStatePlugin, IdealGasPlugin};
use cfd_state::PrimVar;
use field_core::{FieldDefaultPlugins, UniformMeshConfig, Vec3};
use grass_app::prelude::*;
use grass_multi::SubApps;
use grass_scheduler::{CurrentState, ResMut};
use soil_core::{Accum, Atom, AtomPlugin, CommState, ParticleSimScheduleSet};
use soil_verlet::VelocityVerletPlugin;

pub const RADIUS: f64 = 5.0e-4;
pub const DT: f64 = 1.0e-4;
pub const GAS_DENSITY: f64 = 1.2;
pub const GAS_VISCOSITY: f64 = 1.8e-5;
pub const GRAVITY: Vec3 = [0.0, 0.0, -9.81];

fn gravity(mut atoms: ResMut<Atom>) {
    for i in 0..atoms.nlocal as usize {
        atoms.force[i][0] += (atoms.mass[i] as f64 * GRAVITY[0]) as Accum;
        atoms.force[i][1] += (atoms.mass[i] as f64 * GRAVITY[1]) as Accum;
        atoms.force[i][2] += (atoms.mass[i] as f64 * GRAVITY[2]) as Accum;
    }
}

pub fn dem() -> App {
    let mut atoms = Atom::new();
    atoms.dt = DT;
    let mass = 2_500.0 * 4.0 / 3.0 * std::f64::consts::PI * RADIUS.powi(3);
    atoms.push_test_atom(0, [0.005, 0.005, 0.015], RADIUS, mass);
    atoms.nlocal = 1;
    atoms.natoms = 1;

    let mut app = App::new();
    app.add_plugins(AtomPlugin)
        .add_resource(CurrentState(CommState::FullRebuild))
        .add_resource(atoms)
        .add_update_system(gravity, ParticleSimScheduleSet::Force)
        .add_plugins(VelocityVerletPlugin::new());
    app
}

pub fn cfd() -> App {
    let pressure = 101_325.0;
    let temperature = pressure / (GAS_DENSITY * 287.058);
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
    let initial_state = move |_x: Vec3| {
        IdealGas::air().prim_to_cons(&PrimVar::new(
            GAS_DENSITY,
            0.0,
            0.0,
            0.0,
            pressure,
            temperature,
        ))
    };

    let mut app = App::new();
    app.add_plugins(FieldDefaultPlugins { mesh })
        .add_plugins(CfdStatePlugin::new(initial_state))
        .add_plugins(IdealGasPlugin)
        .add_resource(Viscosity::Constant(GAS_VISCOSITY));
    app
}

pub fn particle_height(parent: &App) -> f64 {
    let subapps = parent
        .get_resource_ref::<SubApps>()
        .expect("coupling sub-Apps");
    let dem = subapps.find("dem").expect("dem sub-App");
    let atoms = dem
        .resource_cell(TypeId::of::<Atom>())
        .expect("DEM particles")
        .borrow();
    atoms.downcast_ref::<Atom>().expect("Atom resource").pos[0][2] as f64
}
