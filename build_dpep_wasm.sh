#!/bin/bash
# DPEP Engine WASM Build Script
# Uses the generic build_wasm.sh script with SIMD support

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"$SCRIPT_DIR/build_wasm.sh" "dpep" "🚀 Building DPEP" "createDPEPModule" "-s SIMD=1"
