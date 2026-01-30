#!/bin/bash
# WebAssembly Compilation Script for Aerodynamics Simulators
# Requires Emscripten SDK (emsdk)

echo "🚀 Compiling Aerodynamics Simulators to WebAssembly..."
echo ""

# Get the script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Check if the generic build script exists
if [ ! -f "$SCRIPT_DIR/build_wasm.sh" ]; then
    echo "❌ Error: build_wasm.sh not found!"
    exit 1
fi

# Array of engines to compile: engine_name, display_name, export_name
declare -a engines=(
    "glider:📦 Compiling Stratospheric Glider:createGliderModule"
    "hypersonic:📦 Compiling Hypersonic Flight:createHypersonicModule"
    "vehicle_aero:📦 Compiling Vehicle Downforce:createVehicleModule"
    "ground_effect:📦 Compiling Ground Effect Vehicle:createEkranoplanModule"
)

# Compile each engine
SUCCESS_COUNT=0
FAIL_COUNT=0
for engine_spec in "${engines[@]}"; do
    IFS=':' read -r engine_name display_name export_name <<< "$engine_spec"
    
    echo ""
    echo "$display_name..."
    if "$SCRIPT_DIR/build_wasm.sh" "$engine_name" "$display_name" "$export_name"; then
        echo "✅ ${engine_name} compiled successfully!"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "❌ ${engine_name} compilation failed!"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

echo ""
echo "📊 Compilation Summary:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Module                    | JS Size | WASM Size"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

for engine_spec in "${engines[@]}"; do
    IFS=':' read -r engine_name display_name export_name <<< "$engine_spec"
    
    if [ -f "${engine_name}_engine.js" ]; then
        JS_SIZE=$(ls -lh "${engine_name}_engine.js" | awk '{print $5}')
        WASM_SIZE=$(ls -lh "${engine_name}_engine.wasm" | awk '{print $5}')
        printf "%-25s | %-7s | %-9s\n" "${engine_name}" "$JS_SIZE" "$WASM_SIZE"
    fi
done

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

if [ $FAIL_COUNT -eq 0 ]; then
    echo "🎉 All $SUCCESS_COUNT simulators compiled successfully!"
    echo ""
    echo "🌐 To run:"
    echo "  python3 -m http.server 8000"
    echo "  Then open http://localhost:8000 in your browser"
    echo ""
    exit 0
else
    echo "⚠️  $SUCCESS_COUNT succeeded, $FAIL_COUNT failed"
    exit 1
fi
