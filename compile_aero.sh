#!/bin/bash
# WebAssembly Compilation Script for Aerodynamics Simulators
# Requires Emscripten SDK (emsdk)

echo "🚀 Compiling Aerodynamics Simulators to WebAssembly..."
echo ""

# Check if emcc is available
if ! command -v emcc &> /dev/null; then
    echo "❌ Error: emcc not found!"
    echo "Please install Emscripten SDK:"
    echo "  git clone https://github.com/emscripten-core/emsdk.git"
    echo "  cd emsdk"
    echo "  ./emsdk install latest"
    echo "  ./emsdk activate latest"
    echo "  source ./emsdk_env.sh"
    exit 1
fi

# Common compilation flags
COMMON_FLAGS="-s WASM=1 -s MODULARIZE=1 -s EXPORT_ES6=1 -s ALLOW_MEMORY_GROWTH=1 -lembind -O3 -s ENVIRONMENT=web"

echo "📦 Compiling Stratospheric Glider..."
emcc glider_engine.cpp -o glider_engine.js \
    -s EXPORT_NAME="createGliderModule" \
    $COMMON_FLAGS

if [ $? -eq 0 ]; then
    echo "✅ Glider compiled successfully!"
    ls -lh glider_engine.{js,wasm}
else
    echo "❌ Glider compilation failed!"
    exit 1
fi

echo ""
echo "📦 Compiling Hypersonic Flight..."
emcc hypersonic_engine.cpp -o hypersonic_engine.js \
    -s EXPORT_NAME="createHypersonicModule" \
    $COMMON_FLAGS

if [ $? -eq 0 ]; then
    echo "✅ Hypersonic compiled successfully!"
    ls -lh hypersonic_engine.{js,wasm}
else
    echo "❌ Hypersonic compilation failed!"
    exit 1
fi

echo ""
echo "📦 Compiling Vehicle Downforce..."
emcc vehicle_aero_engine.cpp -o vehicle_aero_engine.js \
    -s EXPORT_NAME="createVehicleModule" \
    $COMMON_FLAGS

if [ $? -eq 0 ]; then
    echo "✅ Vehicle compiled successfully!"
    ls -lh vehicle_aero_engine.{js,wasm}
else
    echo "❌ Vehicle compilation failed!"
    exit 1
fi

echo ""
echo "📦 Compiling Ground Effect Vehicle..."
emcc ground_effect_engine.cpp -o ground_effect_engine.js \
    -s EXPORT_NAME="createEkranoplanModule" \
    $COMMON_FLAGS

if [ $? -eq 0 ]; then
    echo "✅ Ground Effect compiled successfully!"
    ls -lh ground_effect_engine.{js,wasm}
else
    echo "❌ Ground Effect compilation failed!"
    exit 1
fi

echo ""
echo "📊 Compilation Summary:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Module                    | JS Size | WASM Size"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

if [ -f glider_engine.js ]; then
    JS_SIZE=$(ls -lh glider_engine.js | awk '{print $5}')
    WASM_SIZE=$(ls -lh glider_engine.wasm | awk '{print $5}')
    echo "Glider                    | $JS_SIZE    | $WASM_SIZE"
fi

if [ -f hypersonic_engine.js ]; then
    JS_SIZE=$(ls -lh hypersonic_engine.js | awk '{print $5}')
    WASM_SIZE=$(ls -lh hypersonic_engine.wasm | awk '{print $5}')
    echo "Hypersonic                | $JS_SIZE    | $WASM_SIZE"
fi

if [ -f vehicle_aero_engine.js ]; then
    JS_SIZE=$(ls -lh vehicle_aero_engine.js | awk '{print $5}')
    WASM_SIZE=$(ls -lh vehicle_aero_engine.wasm | awk '{print $5}')
    echo "Vehicle Downforce         | $JS_SIZE    | $WASM_SIZE"
fi

if [ -f ground_effect_engine.js ]; then
    JS_SIZE=$(ls -lh ground_effect_engine.js | awk '{print $5}')
    WASM_SIZE=$(ls -lh ground_effect_engine.wasm | awk '{print $5}')
    echo "Ground Effect             | $JS_SIZE    | $WASM_SIZE"
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "🎉 All simulators compiled successfully!"
echo ""
echo "🌐 To run:"
echo "  python3 -m http.server 8000"
echo "  Then open http://localhost:8000 in your browser"
echo ""
