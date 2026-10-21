#!/bin/bash
# ==========================================================================================
# compare_C_to_Fortran.sh
#
# Compares the C (PET BMI method 6) results in noahmp_comparison/output/C with the Fortran
# NWM Noah-MP reference results (r4 = 32-bit REAL as NWM, r8 = 64-bit): those of the last run
# in noahmp_comparison/output/Fortran, or with -r the stored ones in noahmp_comparison/reference.
# For the C exact solver it also compares the per-iteration traces.
# Produce the results first with noahmp_comparison/run_ball_berry_C_and_Fortran.sh, or the C
# results alone with noahmp_comparison/check_C_against_stored_reference.sh.
# ==========================================================================================

# ---------------- directories and files (edit here) ----------------
comparison_directory="$(cd "$(dirname "$0")" && pwd)"
output_directory="${comparison_directory}/output"
c_output_directory="${output_directory}/C"
fortran_output_directory="${output_directory}/Fortran"
stored_reference_directory="${comparison_directory}/reference"
comparison_helper="${comparison_directory}/compare_c_to_noahmp_reference.py"
fortran_summary_r4="${fortran_output_directory}/noahmp_canopy_reference_summary_r4.csv"
fortran_summary_r8="${fortran_output_directory}/noahmp_canopy_reference_summary_r8.csv"
fortran_trace_r4="${fortran_output_directory}/noahmp_canopy_iteration_trace_r4.csv"
c_trace_exact="${c_output_directory}/c_canopy_iteration_trace_exact.csv"
script_name="$(basename "$0")"

print_usage() {
    cat <<EOF
usage: ./noahmp_comparison/${script_name} [-r] [exact|nested|both]

Compares C method 6 results with the Fortran NWM Noah-MP reference (r4 and r8).
  -r       use the stored Fortran results in noahmp_comparison/reference (made with the
           committed configs/method_6_ball_berry.cfg) instead of those of the last run
  (none)   compare every C solver mode that has results in output/C
  exact    compare the C exact solver results only (summary and iteration trace)
  nested   compare the C nested solver results only
  both     compare both; it is an error if either has no results

Produce results first with:  ./noahmp_comparison/run_ball_berry_C_and_Fortran.sh [exact|nested|both]
EOF
}

stop_with_error() {
    echo "ERROR: $1" >&2
    exit 1
}

if [ "${1:-}" = "-r" ]; then fortran_output_directory="${stored_reference_directory}"; shift; fi
if [ $# -gt 1 ]; then print_usage; exit 1; fi
requested_mode="${1:-}"
case "${requested_mode}" in
    -h|--help) print_usage; exit 0 ;;
    "")        c_modes_to_compare=""
               for c_mode in exact nested; do
                   [ -f "${c_output_directory}/c_summary_${c_mode}.csv" ] && c_modes_to_compare="${c_modes_to_compare} ${c_mode}"
               done ;;
    exact|nested) c_modes_to_compare="${requested_mode}" ;;
    both)      c_modes_to_compare="exact nested" ;;
    *)         echo "unrecognized argument '${requested_mode}'" >&2; print_usage; exit 1 ;;
esac

fortran_summary_r4="${fortran_output_directory}/noahmp_canopy_reference_summary_r4.csv"
fortran_summary_r8="${fortran_output_directory}/noahmp_canopy_reference_summary_r8.csv"
fortran_trace_r4="${fortran_output_directory}/noahmp_canopy_iteration_trace_r4.csv"
echo "Fortran results: ${fortran_output_directory}"
[ -f "${fortran_summary_r4}" ] && [ -f "${fortran_summary_r8}" ] \
    || stop_with_error "no Fortran results in ${fortran_output_directory}; run ./noahmp_comparison/run_ball_berry_C_and_Fortran.sh first"
[ -n "${c_modes_to_compare}" ] \
    || stop_with_error "no C results in noahmp_comparison/output/C; run ./noahmp_comparison/run_ball_berry_C_and_Fortran.sh first"

effective_configuration_file="${output_directory}/effective_configuration.txt"
if [ -f "${effective_configuration_file}" ]; then
    echo "configuration used by the last run: noahmp_comparison/output/effective_configuration.txt"
    override_lines="$(grep "OVERRIDDEN on the command line" "${effective_configuration_file}" | sed 's/[[:space:]]*#.*//')"
    if [ -n "${override_lines}" ]; then echo "overrides:"; echo "${override_lines}" | sed 's/^/  /'; else echo "overrides: none"; fi
    echo ""
fi
for c_mode in ${c_modes_to_compare}; do
    c_summary_file="${c_output_directory}/c_summary_${c_mode}.csv"
    [ -f "${c_summary_file}" ] || stop_with_error "no C ${c_mode} results; run ./noahmp_comparison/run_ball_berry_C_and_Fortran.sh ${c_mode}"
    echo "================ C ${c_mode} solver vs Fortran NWM Noah-MP"
    if [ "${c_mode}" = "exact" ] && [ -f "${c_trace_exact}" ] && [ -f "${fortran_trace_r4}" ]; then
        python3 "${comparison_helper}" "${c_summary_file}" "${fortran_summary_r4}" "${fortran_summary_r8}" \
                "${c_trace_exact}" "${fortran_trace_r4}" || stop_with_error "comparison failed"
    else
        python3 "${comparison_helper}" "${c_summary_file}" "${fortran_summary_r4}" "${fortran_summary_r8}" \
            || stop_with_error "comparison failed"
    fi
    echo ""
done

# reminder of any C solver mode not compared
for c_mode in exact nested; do
    case " ${c_modes_to_compare} " in
        *" ${c_mode} "*) ;;
        *) if [ -f "${c_output_directory}/c_summary_${c_mode}.csv" ]; then
               echo "REMINDER: C ${c_mode} results exist but were not compared:  ./noahmp_comparison/${script_name} ${c_mode}"
           else
               echo "REMINDER: no C ${c_mode} results yet. Run:  ./noahmp_comparison/run_ball_berry_C_and_Fortran.sh ${c_mode}   then  ./noahmp_comparison/${script_name} ${c_mode}"
           fi ;;
    esac
done
