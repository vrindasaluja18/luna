# LUNA test suite.
#
# One executable per module, registered with CTest. Tests are deliberately
# analytic where possible: SE(3) identities, Morse test functions, and the
# simulator's own exact landmark projections. A test that only re-derives the
# implementation's own behaviour proves nothing, so the reference value is
# always computed independently of the code under test.

add_executable(test_core test_core.cpp)
target_link_libraries(test_core PRIVATE luna_core)

add_executable(test_sim test_sim.cpp)
target_link_libraries(test_sim PRIVATE luna_core luna_sim)

add_executable(test_frontend test_frontend.cpp)
target_link_libraries(test_frontend PRIVATE luna_core luna_sim luna_frontend)

add_executable(test_backend test_backend.cpp)
target_link_libraries(test_backend PRIVATE luna_core luna_sim luna_frontend luna_backend luna_bench)

add_executable(test_determinism test_determinism.cpp)
target_link_libraries(test_determinism PRIVATE luna_core luna_sim luna_frontend)

add_executable(test_estimator test_estimator.cpp)
target_link_libraries(test_estimator PRIVATE luna_app)

add_executable(test_navigation test_navigation.cpp)
target_link_libraries(test_navigation PRIVATE luna_app)

add_test(NAME core COMMAND test_core)
add_test(NAME sim COMMAND test_sim)
add_test(NAME frontend COMMAND test_frontend)
add_test(NAME backend COMMAND test_backend)
add_test(NAME determinism COMMAND test_determinism)
add_test(NAME estimator COMMAND test_estimator)
add_test(NAME navigation COMMAND test_navigation)
