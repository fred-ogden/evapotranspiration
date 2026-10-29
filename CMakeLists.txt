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
    src/ball_berry_canopy_physics.c
    src/ball_berry_struct_copies.c
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
# Only when this is the top-level project (not when built inside ngen).
# ---------------------------------------------------------------------------
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    include(CTest)
endif()

if(BUILD_TESTING AND CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
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
        COMMAND test_ball_berry_method6 ${PROJECT_SOURCE_DIR}/configs/method_6_ball_berry.cfg)
    set_tests_properties(ball_berry_method6 PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

    # BMI unit test, one ctest per method configuration file.
    add_executable(test_pet_bmi_unit test/main_unit_test_bmi.c)
    target_include_directories(test_pet_bmi_unit PRIVATE include)
    target_link_libraries(test_pet_bmi_unit PRIVATE petbmi m)
    set_target_properties(test_pet_bmi_unit PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
    file(GLOB PET_METHOD_CONFIGURATION_FILES RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/configs/method_*.cfg)
    foreach(configuration_file ${PET_METHOD_CONFIGURATION_FILES})
        get_filename_component(configuration_name ${configuration_file} NAME_WE)
        add_test(NAME bmi_unit_${configuration_name} COMMAND test_pet_bmi_unit ${configuration_file})
        set_tests_properties(bmi_unit_${configuration_name} PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    endforeach()

    # FLO 10/24: revised canopy energy physics (canopy_energy_physics_option = 2): physics checks
    # through the BMI, and the BMI unit test with that configuration.
    add_executable(test_ball_berry_canopy_energy test/test_ball_berry_canopy_energy.c)
    target_include_directories(test_ball_berry_canopy_energy PRIVATE include)
    target_link_libraries(test_ball_berry_canopy_energy PRIVATE petbmi m)
    set_target_properties(test_ball_berry_canopy_energy PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
    add_test(NAME ball_berry_revised_canopy_energy
        COMMAND test_ball_berry_canopy_energy configs/ball_berry_revised_canopy_energy.cfg)
    set_tests_properties(ball_berry_revised_canopy_energy PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    add_test(NAME bmi_unit_ball_berry_revised_canopy_energy
        COMMAND test_pet_bmi_unit configs/ball_berry_revised_canopy_energy.cfg)

    # FLO 10/26: options for the decoupled Noah-MP stack (solar geometry, RSURF with snow, DVEG = 1,
    # buried canopy, soil albedo from soil moisture), and the BMI unit test with that configuration.
    add_executable(test_ball_berry_noahmp_options test/test_ball_berry_noahmp_options.c)
    target_include_directories(test_ball_berry_noahmp_options PRIVATE include)
    target_link_libraries(test_ball_berry_noahmp_options PRIVATE petbmi m)
    set_target_properties(test_ball_berry_noahmp_options PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
    add_test(NAME ball_berry_noahmp_options
        COMMAND test_ball_berry_noahmp_options configs/ball_berry_noahmp_cat-18560.cfg)
    set_tests_properties(ball_berry_noahmp_options PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    add_test(NAME bmi_unit_ball_berry_noahmp_cat-18560
        COMMAND test_pet_bmi_unit configs/ball_berry_noahmp_cat-18560.cfg)
    set_tests_properties(bmi_unit_ball_berry_noahmp_cat-18560 PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    set_tests_properties(bmi_unit_ball_berry_revised_canopy_energy PROPERTIES WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

    # Method 6 against the stored NWM Noah-MP (Fortran r8) results; needs bash and python3.
    if(NOT WIN32)
        add_test(NAME ball_berry_method6_vs_stored_noahmp_reference
            COMMAND bash ${PROJECT_SOURCE_DIR}/noahmp_comparison/check_C_against_stored_reference.sh)
    endif()
endif()
