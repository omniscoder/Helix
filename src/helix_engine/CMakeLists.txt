cmake_minimum_required(VERSION 3.15)
project(helix_engine LANGUAGES CXX)

option(HELIX_ENGINE_ENABLE_CUDA "Build CUDA kernels for Helix engine" OFF)
if(HELIX_ENGINE_ENABLE_CUDA)
    enable_language(CUDA)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(helix_engine STATIC
    src/physics.cpp
)
if(HELIX_ENGINE_ENABLE_CUDA)
    target_sources(helix_engine PRIVATE src/physics_cuda.cu)
    target_compile_definitions(helix_engine PUBLIC HELIX_ENGINE_HAS_CUDA=1)
endif()
target_include_directories(helix_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(helix_engine_py MODULE
    src/binding.cpp
)
target_link_libraries(helix_engine_py PRIVATE helix_engine)
# pybind11 will normally be discovered by the build environment; this file is a placeholder
# until the full native build system is wired up.
