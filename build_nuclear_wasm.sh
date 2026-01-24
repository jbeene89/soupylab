#!/bin/bash
# Nuclear War Engine WASM Build Script

set -e

echo "☢️ Building Nuclear War Engine WebAssembly Module..."

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

emcc nuclear_war_engine.cpp -o nuclear_war_engine.js \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_ES6=1 \
    -s EXPORT_NAME="createNuclearModule" \
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
echo "   - nuclear_war_engine.js"
echo "   - nuclear_war_engine.wasm"
ls -lh nuclear_war_engine.js nuclear_war_engine.wasm 2>/dev/null || echo "   (Build in progress...)"
