#!/bin/bash
# Asteroid Impact Engine WASM Build Script

set -e

echo "☄️ Building Asteroid Impact Engine WebAssembly Module..."

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
cd /home/user/soupylab

emcc asteroid_impact_engine.cpp -o asteroid_impact_engine.js \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_ES6=1 \
    -s EXPORT_NAME="createAsteroidModule" \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s MAXIMUM_MEMORY=512MB \
    -lembind \
    -O3 \
    -s ASSERTIONS=0 \
    -s MALLOC=emmalloc \
    --closure 1 \
    -s ENVIRONMENT=web

echo "✅ Build complete!"
echo "📁 Generated files:"
echo "   - asteroid_impact_engine.js"
echo "   - asteroid_impact_engine.wasm"
ls -lh asteroid_impact_engine.js asteroid_impact_engine.wasm 2>/dev/null || echo "   (Build in progress...)"
