#!/usr/bin/env bash
# dev_couple_dem_cfd validation harness — DEM<->CFD interphase-drag couplings.
#
#   ./validation/run.sh          # smoke: fast coupling gates (every CI push/PR)
#   ./validation/run.sh full     # smoke + heavier bed/fiber gates (scheduled / on-demand)
#
# Each gate checks the live drag seam against a literature reference with an
# independent closure + a negative control (see VALIDATION.md). Every runnable
# example lives under examples/<name>/ (no crate-embedded examples).
set -euo pipefail
cd "$(dirname "$0")/.."
MODE="${1:-smoke}"

echo "=== dev_couple_dem_cfd validation set (mode=$MODE) ==="

# --- Smoke gate: fast point-particle seam vs a closed-form reference. ---

# Stokes terminal velocity, low-Reynolds settling sphere through the drag seam.
cargo run --release --example settling_sphere -- examples/settling_sphere/config.toml

if [ "$MODE" = "full" ]; then
  # --- Full gates: heavier bed/fiber couplings. Scheduled / on-demand only. ---

  # Ergun (1952) packed-bed pressure drop across a Reynolds sweep (MacDonald control).
  cargo run --release --example fixed_bed_ergun -- examples/fixed_bed_ergun/config.toml

  # Wen & Yu (1966) minimum fluidization via live net-force bisection (neg controls).
  cargo run --release --example fluidized_bed_umf -- examples/fluidized_bed_umf/config.toml

  # Resolved DIRT bonded-clump fiber in the gas: Archimedes + Tirado slender-body drag.
  cargo run --release --example cfd_ibm_fiber -- examples/cfd_ibm_fiber/config.toml
fi

echo "=== all dev_couple_dem_cfd validation gates passed (mode=$MODE) ==="
