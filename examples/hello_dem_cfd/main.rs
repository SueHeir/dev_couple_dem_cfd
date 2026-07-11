//! The smallest useful DEM↔CFD lesson: construct two ordinary GRASS Apps,
//! mount them, add their coupling, and run.

mod setup;

use dem_cfd::DemCfdCouplingPlugin;
use grass_app::prelude::*;
use grass_multi::MultiAppExt;
use setup::{DT, GRAVITY, RADIUS};

fn main() {
    let mut app = App::new();

    app.add_subapp("dem", setup::dem())
        .add_subapp("cfd", setup::cfd())
        .add_plugins(DemCfdCouplingPlugin::for_air(RADIUS, 200, DT, GRAVITY))
        .start();

    println!(
        "particle height after coupling: {:.6} m",
        setup::particle_height(&app)
    );
}
