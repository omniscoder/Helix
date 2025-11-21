#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN=${PYTHON_BIN:-python3}
BUILD_DIR=${BUILD_DIR:-build/native_cuda}
PYTHON_EXE=$($PYTHON_BIN -c 'import sys; print(sys.executable)')
PYBIND_INCLUDES=$($PYTHON_BIN -m pybind11 --includes 2>/dev/null || echo "")

cmake_args=(
  -DHELIX_ENGINE_ENABLE_CUDA=ON
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  -DPython_EXECUTABLE="$PYTHON_EXE"
)
if [ -n "$PYBIND_INCLUDES" ]; then
  cmake_args+=("-DCMAKE_CXX_FLAGS:STRING=$PYBIND_INCLUDES" "-DCMAKE_CUDA_FLAGS:STRING=$PYBIND_INCLUDES")
fi

cmake -S src/helix_engine -B "$BUILD_DIR" "${cmake_args[@]}"
cmake --build "$BUILD_DIR" --target helix_engine_py
EXT_SUFFIX=$($PYTHON_BIN - <<'PY'
import sysconfig
print(sysconfig.get_config_var('EXT_SUFFIX') or '')
PY
)
MODULE=$(ls "$BUILD_DIR"/*helix_engine_py*.so | head -n1)
cp "$MODULE" "src/helix_engine/_native${EXT_SUFFIX}"
echo "Native CUDA module copied to src/helix_engine/_native${EXT_SUFFIX}"
