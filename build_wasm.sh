#!/bin/bash
# Generic WASM Build Script
# Usage: ./build_wasm.sh <engine_name> <display_icon> <export_name> [extra_flags]
#
# Example: ./build_wasm.sh asteroid_impact "☄️ Asteroid Impact" "createAsteroidModule"
#

set -e

# Validate arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <engine_name> <display_icon> <export_name> [extra_flags]"
    echo ""
    echo "Arguments:"
    echo "  engine_name   - Name of the engine (e.g., asteroid_impact)"
    echo "  display_icon  - Display name with icon (e.g., '☄️ Asteroid Impact')"
    echo "  export_name   - ES6 export name (e.g., createAsteroidModule)"
    echo "  extra_flags   - Optional additional emcc flags (e.g., '-s SIMD=1')"
    echo ""
    echo "Example:"
    echo "  $0 asteroid_impact '☄️ Asteroid Impact' createAsteroidModule"
    exit 1
fi

ENGINE_NAME="$1"
DISPLAY_ICON="$2"
EXPORT_NAME="$3"
EXTRA_FLAGS="${4:-}"

echo "$DISPLAY_ICON Engine WebAssembly Module..."

# Activate Emscripten
if [ -d "/tmp/emsdk" ]; then
    source /tmp/emsdk/emsdk_env.sh
else
    echo "❌ Emscripten not found. Installing..."
    cd /tmp
    git clone --depth 1 https://github.com/emscripten-core/emsdk.git
    cd emsdk
    ./emsdk install latest
    ./emsdk activate latest
    source ./emsdk_env.sh
fi

# Build with optimizations
echo "📦 Compiling C++ to WebAssembly..."
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

emcc ${ENGINE_NAME}_engine.cpp -o ${ENGINE_NAME}_engine.js \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_ES6=1 \
    -s EXPORT_NAME="$EXPORT_NAME" \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s MAXIMUM_MEMORY=512MB \
    -lembind \
    -O3 \
    $EXTRA_FLAGS \
    -s ASSERTIONS=0 \
    -s MALLOC=emmalloc \
    --closure 1 \
    -s ENVIRONMENT=web

echo "✅ Build complete!"
echo "📁 Generated files:"
echo "   - ${ENGINE_NAME}_engine.js"
echo "   - ${ENGINE_NAME}_engine.wasm"
ls -lh ${ENGINE_NAME}_engine.js ${ENGINE_NAME}_engine.wasm 2>/dev/null || echo "   (Build in progress...)"
