#!/bin/bash
# ==========================================================================================
# check_C_against_stored_reference.sh
#
# Regression test of PET method 6 against NWM Noah-MP without the NWM source: builds and runs
# the C exact solver (VEGE_FLUX loop1 replica) through the BMI with the committed
# configs/method_6_ball_berry.cfg, then compares its solution summary with the stored Fortran
# r8 (64-bit) results in noahmp_comparison/reference. Fails if any compared variable differs
# from r8 by more than ${maximum_relative_difference_vs_r8} (relative).
#
# The stored results hold only for the committed configuration and forcing; after changing
# either, rerun noahmp_comparison/run_ball_berry_C_and_Fortran.sh with the NWM source and
# copy its output/Fortran/noahmp_canopy_* files to noahmp_comparison/reference.
# Run from anywhere.  -h prints this usage.
# ==========================================================================================

# ---------------- directories and files ----------------
comparison_directory="$(cd "$(dirname "$0")" && pwd)"
project_root_directory="$(cd "${comparison_directory}/.." && pwd)"
configuration_file="configs/method_6_ball_berry.cfg"                       # relative to the repository root
c_source_directory="${project_root_directory}/src"
c_include_directory="${project_root_directory}/include"
c_driver_source_file="${project_root_directory}/test/main_ball_berry_bmi.c"
c_build_directory="${comparison_directory}/build/C"
c_output_directory="${comparison_directory}/output/C"
check_configuration_file="noahmp_comparison/output/C/configuration_exact.txt" # relative to the repository root
maximum_relative_difference_vs_r8="1.0e-12"

if [ "${1:-}" = "-h" ]; then sed -n '3,15p' "$0" | sed 's/^# \{0,1\}//'; exit 0; fi
stop_with_error() { echo "ERROR: $1" >&2; exit 1; }

cd "${project_root_directory}" || stop_with_error "cannot cd to ${project_root_directory}"
CC="${CC:-gcc}"
mkdir -p "${c_build_directory}" "${c_output_directory}"
rm -f "${c_output_directory}"/*.csv "${c_output_directory}"/*.txt "${comparison_directory}/output/effective_configuration.txt"

echo "=== C: building with ${CC}"
${CC} -std=c99 -O2 -DBMI_ACTIVE -I"${c_include_directory}" \
    "${c_source_directory}/bmi_pet.c" "${c_source_directory}/pet.c" "${c_source_directory}/noahmp_stomata.c" \
    "${c_source_directory}/ball_berry_state.c" "${c_source_directory}/ball_berry_helpers.c" \
    "${c_source_directory}/ball_berry_method6_inputs.c" "${c_driver_source_file}" \
    -lm -o "${c_build_directory}/test_ball_berry_method6" > "${c_build_directory}/build.log" 2>&1 \
    || stop_with_error "C build failed, see ${c_build_directory}/build.log"

echo "=== C: running the exact solver with ${configuration_file}"
{
    sed -e '/^[[:space:]]*bb_noahmp_exact_iteration[[:space:]]*=/d' \
        -e '/^[[:space:]]*bb_diagnostic_summary_file[[:space:]]*=/d' \
        -e '/^[[:space:]]*bb_iteration_trace_file[[:space:]]*=/d' "${configuration_file}"
    echo ""
    echo "# ---------------- set by check_C_against_stored_reference.sh ----------------"
    echo "bb_noahmp_exact_iteration = 1"
    echo "bb_diagnostic_summary_file = ./noahmp_comparison/output/C/c_summary_exact.csv"
    echo "bb_iteration_trace_file = ./noahmp_comparison/output/C/c_canopy_iteration_trace_exact.csv"
} > "${check_configuration_file}"
"${c_build_directory}/test_ball_berry_method6" "${check_configuration_file}" > "${c_output_directory}/c_bmi_outputs_exact.txt" \
    || stop_with_error "C exact run failed"

echo "=== comparing with the stored Fortran results in noahmp_comparison/reference"
python3 "${comparison_directory}/compare_c_to_noahmp_reference.py" \
    --maximum-relative-difference-vs-r8 "${maximum_relative_difference_vs_r8}" \
    "${c_output_directory}/c_summary_exact.csv" \
    "${comparison_directory}/reference/noahmp_canopy_reference_summary_r4.csv" \
    "${comparison_directory}/reference/noahmp_canopy_reference_summary_r8.csv" \
    "${c_output_directory}/c_canopy_iteration_trace_exact.csv" \
    "${comparison_directory}/reference/noahmp_canopy_iteration_trace_r4.csv"
