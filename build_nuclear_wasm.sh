#!/bin/bash
# Nuclear War Engine WASM Build Script
# Uses the generic build_wasm.sh script

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"$SCRIPT_DIR/build_wasm.sh" "nuclear_war" "☢️ Building Nuclear War" "createNuclearModule"
