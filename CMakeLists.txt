cmake_minimum_required(VERSION 3.16)
project(slamcore VERSION 1.0.0 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# SLAMCore-8G: CPU-only, 8 GB, no-CUDA SLAM + navigation stack.
# Hard rules enforced by this build:
#   * C++17, no GPU code paths, no ROS, no network at runtime
#   * all long-lived allocations go through sc::alloc pools
# ---------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(SLAMCORE_ASAN "Build with AddressSanitizer/UBSan" OFF)
option(SLAMCORE_TESTS "Build the test suite" ON)

set(SC_WARN -Wall -Wextra -Wno-unused-parameter)

add_library(slamcore_flags INTERFACE)
target_compile_options(slamcore_flags INTERFACE ${SC_WARN} -O3 -funroll-loops)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64|AMD64")
  target_compile_options(slamcore_flags INTERFACE -mssse3 -msse4.1)
endif()
if(SLAMCORE_ASAN)
  target_compile_options(slamcore_flags INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_options(slamcore_flags INTERFACE -fsanitize=address,undefined)
endif()

# Eigen is the only external dependency: header-only, so no link step and no
# runtime shared-object coupling.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(Eigen3_FOUND)
  get_target_property(SC_EIGEN_INCLUDE Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
else()
  find_path(SC_EIGEN_INCLUDE Eigen/Core PATHS /usr/include/eigen3 /usr/local/include/eigen3)
endif()
if(NOT SC_EIGEN_INCLUDE)
  message(FATAL_ERROR "Eigen3 not found. Install it with: apt-get install libeigen3-dev")
endif()
message(STATUS "Eigen3 headers: ${SC_EIGEN_INCLUDE}")

add_library(sc_core
  src/core/se3.cpp
  src/core/se2.cpp
  src/core/qp.cpp
  src/core/sparse.cpp
  src/core/alloc.cpp
  src/core/stats.cpp)
target_include_directories(sc_core PUBLIC include ${SC_EIGEN_INCLUDE})
target_link_libraries(sc_core PUBLIC slamcore_flags)

add_library(sc_sim
  src/sim/world.cpp
  src/sim/simulator.cpp
  src/sim/vehicle.cpp
  src/sim/sensors.cpp
  src/sim/faults.cpp)
target_link_libraries(sc_sim PUBLIC sc_core)

add_library(sc_frontend
  src/frontend/fast.cpp
  src/frontend/brief.cpp
  src/frontend/klt.cpp
  src/frontend/stereo.cpp
  src/frontend/odometry.cpp)
target_link_libraries(sc_frontend PUBLIC sc_core sc_sim)

add_library(sc_backend
  src/backend/pose_graph.cpp
  src/backend/loop_close.cpp)
target_link_libraries(sc_backend PUBLIC sc_core sc_frontend sc_sim)

add_library(sc_mapping
  src/mapping/occupancy.cpp)
target_link_libraries(sc_mapping PUBLIC sc_core sc_frontend sc_sim)

add_library(sc_topology
  src/topology/morse.cpp)
target_link_libraries(sc_topology PUBLIC sc_core sc_mapping sc_sim)

add_library(sc_planning
  src/planning/planner.cpp)
target_link_libraries(sc_planning PUBLIC sc_core sc_topology sc_sim)

add_library(sc_control
  src/control/mpc.cpp)
target_link_libraries(sc_control PUBLIC sc_core sc_planning sc_sim)

add_library(sc_bench
  src/bench/trajectory.cpp)
target_link_libraries(sc_bench PUBLIC sc_core)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
add_library(sc_app
  src/app/estimator.cpp
  src/app/navigation.cpp)
target_link_libraries(sc_app PUBLIC sc_core sc_sim sc_frontend sc_backend sc_bench sc_mapping sc_topology sc_planning sc_control)
target_include_directories(sc_app PUBLIC include)

add_executable(slamcore src/app/main.cpp)
target_link_libraries(slamcore PRIVATE sc_app)

# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------
if(SLAMCORE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

install(TARGETS slamcore RUNTIME DESTINATION bin)
