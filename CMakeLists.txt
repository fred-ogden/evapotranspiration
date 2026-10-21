cmake_minimum_required(VERSION 3.12)

project(petbmi VERSION 1.0.0 DESCRIPTION "OWP PET BMI Module Shared Library")

# PET
set(PET_LIB_NAME_CMAKE petbmi)
set(PET_LIB_DESC_CMAKE "OWP PET BMI Module Shared Library")

# Make sure these are compiled with this directive
add_compile_definitions(BMI_ACTIVE)

set(PET_SOURCES
    src/bmi_pet.c
    src/pet.c
    src/noahmp_stomata.c
    src/ball_berry_state.c
    src/ball_berry_helpers.c
    src/ball_berry_method6_inputs.c
)

if(WIN32)
    add_library(petbmi ${PET_SOURCES})
else()
    add_library(petbmi SHARED ${PET_SOURCES})
endif()

target_include_directories(petbmi PRIVATE include)

# FLO 10/26: The shared library previously omitted libm, leaving C math
# functions unresolved when linking the BMI regression-test executable.
if(NOT WIN32)
    target_link_libraries(petbmi PRIVATE m)
endif()

set_target_properties(petbmi PROPERTIES VERSION ${PROJECT_VERSION})

set_target_properties(petbmi PROPERTIES PUBLIC_HEADER include/bmi_pet.h)

# Code requires minimum of C99 standard to compile
set_target_properties(petbmi PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

install(TARGETS petbmi
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

configure_file(petbmi.pc.in petbmi.pc @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/petbmi.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

# ---------------------------------------------------------------------------
# Tests (ctest). Run from the build directory with: ctest --output-on-failure
# ---------------------------------------------------------------------------
include(CTest)

if(BUILD_TESTING)
    # Ball-Berry vegetation-state (LAI, SAI, SHDFAC, ...) validation.
    add_executable(test_ball_berry_state
        test/test_ball_berry_state.c
        src/ball_berry_state.c
    )
    target_include_directories(test_ball_berry_state PRIVATE include)
    set_target_properties(test_ball_berry_state PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
    add_test(NAME ball_berry_state COMMAND test_ball_berry_state)

    # Method 6 through the BMI only: initialize, set_value (forcing and time), update,
    # get_value, finalize. Paths in the configuration are relative to the repository root.
    add_executable(test_ball_berry_method6 test/main_ball_berry_bmi.c)
    target_include_directories(test_ball_berry_method6 PRIVATE include)
    target_link_libraries(test_ball_berry_method6 PRIVATE petbmi m)
    set_target_properties(test_ball_berry_method6 PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
    add_test(NAME ball_berry_method6
        COMMAND test_ball_berry_method6 ${CMAKE_SOURCE_DIR}/configs/method_6_ball_berry.cfg)
    set_tests_properties(ball_berry_method6 PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()
