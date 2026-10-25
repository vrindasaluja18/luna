cmake_minimum_required(VERSION 3.16)
project(luna VERSION 1.0.0 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# LUNA: CPU-only, 8 GB, no-CUDA SLAM + navigation stack.
# Hard rules enforced by this build:
#   * C++17, no GPU code paths, no ROS, no network at runtime
#   * all long-lived allocations go through luna::alloc pools
# ---------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(LUNA_ASAN "Build with AddressSanitizer/UBSan" OFF)
option(LUNA_TESTS "Build the test suite" ON)

set(LUNA_WARN -Wall -Wextra -Wno-unused-parameter)

add_library(luna_flags INTERFACE)
target_compile_options(luna_flags INTERFACE ${LUNA_WARN} -O3 -funroll-loops)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64|AMD64")
  target_compile_options(luna_flags INTERFACE -mssse3 -msse4.1)
endif()
if(LUNA_ASAN)
  target_compile_options(luna_flags INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(luna_flags INTERFACE -fsanitize=address,undefined)
endif()

# Eigen is the only external dependency: header-only, so no link step and no
# runtime shared-object coupling.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(Eigen3_FOUND)
  get_target_property(LUNA_EIGEN_INCLUDE Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
else()
  find_path(LUNA_EIGEN_INCLUDE Eigen/Core PATHS /usr/include/eigen3 /usr/local/include/eigen3)
endif()
if(NOT LUNA_EIGEN_INCLUDE)
  message(FATAL_ERROR "Eigen3 not found. Install it with: apt-get install libeigen3-dev")
endif()
message(STATUS "Eigen3 headers: ${LUNA_EIGEN_INCLUDE}")

add_library(luna_core
  src/core/se3.cpp
  src/core/se2.cpp
  src/core/qp.cpp
  src/core/sparse.cpp
  src/core/alloc.cpp
  src/core/stats.cpp)
target_include_directories(luna_core PUBLIC include ${LUNA_EIGEN_INCLUDE})
target_link_libraries(luna_core PUBLIC luna_flags)

add_library(luna_sim
  src/sim/world.cpp
  src/sim/simulator.cpp
  src/sim/vehicle.cpp
  src/sim/sensors.cpp
  src/sim/faults.cpp)
target_link_libraries(luna_sim PUBLIC luna_core)

add_library(luna_frontend
  src/frontend/fast.cpp
  src/frontend/brief.cpp
  src/frontend/klt.cpp
  src/frontend/stereo.cpp
  src/frontend/odometry.cpp)
target_link_libraries(luna_frontend PUBLIC luna_core luna_sim)

add_library(luna_backend
  src/backend/pose_graph.cpp
  src/backend/loop_close.cpp)
target_link_libraries(luna_backend PUBLIC luna_core luna_frontend luna_sim)

add_library(luna_mapping
  src/mapping/occupancy.cpp)
target_link_libraries(luna_mapping PUBLIC luna_core luna_frontend luna_sim)

add_library(luna_topology
  src/topology/morse.cpp)
target_link_libraries(luna_topology PUBLIC luna_core luna_mapping luna_sim)

add_library(luna_planning
  src/planning/planner.cpp)
target_link_libraries(luna_planning PUBLIC luna_core luna_topology luna_sim)

add_library(luna_control
  src/control/mpc.cpp)
target_link_libraries(luna_control PUBLIC luna_core luna_planning luna_sim)

add_library(luna_bench
  src/bench/trajectory.cpp)
target_link_libraries(luna_bench PUBLIC luna_core)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
add_library(luna_app
  src/app/estimator.cpp
  src/app/navigation.cpp)
target_link_libraries(luna_app PUBLIC luna_core luna_sim luna_frontend luna_backend luna_bench luna_mapping luna_topology luna_planning luna_control)
target_include_directories(luna_app PUBLIC include)

add_executable(luna src/app/main.cpp)
target_link_libraries(luna PRIVATE luna_app)

# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------
if(LUNA_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

install(TARGETS luna RUNTIME DESTINATION bin)
