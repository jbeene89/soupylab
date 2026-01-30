#!/bin/bash
# Asteroid Impact Engine WASM Build Script
# Uses the generic build_wasm.sh script

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"$SCRIPT_DIR/build_wasm.sh" "asteroid_impact" "☄️ Building Asteroid Impact" "createAsteroidModule"
