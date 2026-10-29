#!/bin/bash
# ==========================================================================================
# run_ball_berry_C_and_Fortran.sh
#
# Builds and runs the C (PET BMI method 6) and Fortran (NWM Noah-MP reference) Ball-Berry
# canopy transpiration codes, both reading the same configuration and forcing files.
# Results go to noahmp_comparison/output/C and noahmp_comparison/output/Fortran; compare them
# with noahmp_comparison/compare_C_to_Fortran.sh.
# Optional -s key=value arguments override configuration values for this run only: the
# configuration file is not edited; both codes read noahmp_comparison/output/effective_configuration.txt.
#
# The NWM Noah-MP source module_sf_noahmplsm.F is NOT in this repository. Set its path in
# noahmp_source_file below, or with the environment variable NOAHMP_SOURCE_FILE.
# Without it, check_C_against_stored_reference.sh compares the C code with the stored
# Fortran results in noahmp_comparison/reference instead.
#
# Compilers: CC (default gcc) and FC (default gfortran), e.g.  FC=flang ./run_ball_berry_C_and_Fortran.sh
# ==========================================================================================

# ---------------- directories and files (edit here) ----------------
comparison_directory="$(cd "$(dirname "$0")" && pwd)"
project_root_directory="$(cd "${comparison_directory}/.." && pwd)"   # repository root; the scripts run here
noahmp_source_file="${NOAHMP_SOURCE_FILE:-}"   # NWM module_sf_noahmplsm.F, not in this repository
configuration_file="configs/method_6_ball_berry.cfg"                              # relative to the repository root
effective_configuration_file="noahmp_comparison/output/effective_configuration.txt"   # configuration with -s overrides applied
c_source_directory="${project_root_directory}/src"
c_include_directory="${project_root_directory}/include"
c_driver_source_file="${project_root_directory}/test/main_ball_berry_bmi.c"
c_build_directory="${comparison_directory}/build/C"
fortran_source_directory="${comparison_directory}/Fortran"
fortran_build_directory="${comparison_directory}/build/Fortran"
c_output_directory="${comparison_directory}/output/C"
fortran_output_directory="${comparison_directory}/output/Fortran"
script_name="noahmp_comparison/$(basename "$0")"

print_usage() {
    cat <<EOF
usage: ./${script_name} [exact|nested|both] [-s key=value ...]
       ./${script_name} -l      list the configuration keys that -s can override
       ./${script_name} -h      this message

Builds and runs the C and Fortran Ball-Berry codes with
  configuration: ${configuration_file}
  forcing:       the forcing_csv_file named in the configuration

C solver mode (the Fortran is always NWM Noah-MP):
  (none)   use bb_noahmp_exact_iteration from the configuration file
  exact    Noah-MP VEGE_FLUX loop1 replica, this run only
  nested   nested Picard/Newton solver, this run only
  both     run the C code in both modes

Overrides (this run only; the configuration file is not changed):
  -s key=value   replace the value of an existing configuration key; repeat as needed.
                 Both codes read the result, written to ${effective_configuration_file}.
  examples:  ./${script_name} both -s root_stress_factor_0_to_1=0.5
             ./${script_name} exact -s snow_cover_fraction=0.4 -s snow_depth_m=6.0

NWM Noah-MP source (not in this repository): noahmp_source_file at the top of this script,
           or the environment variable NOAHMP_SOURCE_FILE=/path/to/module_sf_noahmplsm.F
Compilers: CC (default gcc), FC (default gfortran), e.g.  FC=flang ./${script_name} both
Outputs:   noahmp_comparison/output/C and output/Fortran.
Compare:   ./noahmp_comparison/compare_C_to_Fortran.sh

To list the configuration keys you can override with -s, run:  ./${script_name} -l
EOF
}

# Print every key of the configuration file with its current value and comment, under the
# section headings of the file. The list is read from the file itself, so it is always current.
# Keys the script sets for each C solver mode are marked: they cannot usefully be overridden.
print_override_list() {
    echo "Configuration keys that -s key=value can override, from ${configuration_file}:"
    echo "(current value and comment; a value set with -s applies to one run only)"
    awk '
        /^#[[:space:]]*-+[[:space:]].*[[:space:]]-+[[:space:]]*$/ {
            heading = $0
            gsub(/^#[[:space:]]*-+[[:space:]]*/, "", heading)
            gsub(/[[:space:]]*-+[[:space:]]*$/, "", heading)
            printf "\n%s\n", heading
            next
        }
        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
        /=/ {
            line = $0
            comment = ""
            comment_position = index(line, "#")
            if (comment_position > 0) {
                comment = substr(line, comment_position + 1)
                line = substr(line, 1, comment_position - 1)
                gsub(/^[[:space:]]+/, "", comment)
            }
            equals_position = index(line, "=")
            key = substr(line, 1, equals_position - 1)
            value = substr(line, equals_position + 1)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", key)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            if (key == "bb_noahmp_exact_iteration" || key == "bb_diagnostic_summary_file" || key == "bb_iteration_trace_file")
                comment = "SET BY THE RUN SCRIPT for each C solver mode; use exact|nested|both instead"
            if (comment != "") printf "  %-50s = %-14s  %s\n", key, value, comment
            else               printf "  %-50s = %s\n", key, value
        }
    ' "${configuration_file}"
    echo ""
    echo "example:  ./${script_name} both -s root_stress_factor_0_to_1=0.5 -s snow_cover_fraction=0.4"
}

stop_with_error() {
    echo "ERROR: $1" >&2
    exit 1
}

# ---------------- command line ----------------
requested_mode=""
configuration_overrides=()
while [ $# -gt 0 ]; do
    case "$1" in
        -h|--help) print_usage; exit 0 ;;
        -l|--list)
            cd "${project_root_directory}" || stop_with_error "cannot cd to ${project_root_directory}"
            [ -f "${configuration_file}" ] || stop_with_error "configuration file ${configuration_file} not found"
            print_override_list; exit 0 ;;
        exact|nested|both)
            [ -z "${requested_mode}" ] || { echo "more than one solver mode given" >&2; print_usage; exit 1; }
            requested_mode="$1" ;;
        -s)
            [ $# -ge 2 ] || { echo "-s needs key=value" >&2; print_usage; exit 1; }
            case "$2" in
                *=*) configuration_overrides+=("$2") ;;
                *) echo "-s needs key=value, got '$2'" >&2; print_usage; exit 1 ;;
            esac
            shift ;;
        *) echo "unrecognized argument '$1'" >&2; print_usage; exit 1 ;;
    esac
    shift
done

cd "${project_root_directory}" || stop_with_error "cannot cd to ${project_root_directory}"
[ -f "${configuration_file}" ] || stop_with_error "configuration file ${configuration_file} not found"

# ---------------- effective configuration: the configuration file with -s overrides ----------------
mkdir -p "$(dirname "${effective_configuration_file}")"
cp "${configuration_file}" "${effective_configuration_file}" || stop_with_error "cannot write ${effective_configuration_file}"
for override in "${configuration_overrides[@]}"; do
    override_key="${override%%=*}"
    override_value="${override#*=}"
    grep -q "^[[:space:]]*${override_key}[[:space:]]*=" "${effective_configuration_file}" \
        || stop_with_error "-s ${override}: key '${override_key}' is not in ${configuration_file}"
    sed -i "s|^[[:space:]]*${override_key}[[:space:]]*=.*|${override_key} = ${override_value}   # OVERRIDDEN on the command line (-s)|" "${effective_configuration_file}"
done

# mode set in the configuration file
configured_exact_flag="$(sed -n 's/^[[:space:]]*bb_noahmp_exact_iteration[[:space:]]*=[[:space:]]*\([0-9]*\).*/\1/p' "${effective_configuration_file}" | tail -n 1)"
if [ "${configured_exact_flag}" = "1" ]; then configured_mode="exact"; else configured_mode="nested"; fi

case "${requested_mode}" in
    "")   c_modes_to_run="${configured_mode}"; mode_source="configuration file" ;;
    both) c_modes_to_run="exact nested"; mode_source="command line" ;;
    *)    c_modes_to_run="${requested_mode}"; mode_source="command line" ;;
esac

CC="${CC:-gcc}"
FC="${FC:-gfortran}"
case "$(basename "${FC}")" in
    gfortran*) fortran_free_form_flags="-ffree-form -ffree-line-length-none" ;;
    *)         fortran_free_form_flags="-ffree-form" ;;
esac
command -v "${CC}" > /dev/null || stop_with_error "C compiler ${CC} not found"
command -v "${FC}" > /dev/null || stop_with_error "Fortran compiler ${FC} not found (set FC)"
[ -n "${noahmp_source_file}" ] || stop_with_error "set noahmp_source_file in ${script_name} or NOAHMP_SOURCE_FILE to the NWM module_sf_noahmplsm.F (without it, use noahmp_comparison/check_C_against_stored_reference.sh)"
[ -f "${noahmp_source_file}" ] || stop_with_error "NWM Noah-MP source ${noahmp_source_file} not found"

mkdir -p "${c_build_directory}" "${fortran_build_directory}" "${c_output_directory}" "${fortran_output_directory}"
# every result in output/ comes from this run: remove results of earlier runs
rm -f "${c_output_directory}"/*.csv "${c_output_directory}"/*.txt "${fortran_output_directory}"/*.csv

# ---------------- Fortran NWM Noah-MP reference, r4 and r8 ----------------
echo "=== Fortran: patching NWM module_sf_noahmplsm.F (HARNESS PATCH changes only)"
python3 "${fortran_source_directory}/make_harness_module.py" \
    "${noahmp_source_file}" \
    "${fortran_build_directory}/module_sf_noahmplsm_harness.F" || stop_with_error "patching failed"

for precision_tag in r4 r8; do
    precision_build_directory="${fortran_build_directory}/${precision_tag}"
    mkdir -p "${precision_build_directory}"
    if [ "${precision_tag}" = "r8" ]; then precision_flags="-fdefault-real-8"; else precision_flags=""; fi
    echo "=== Fortran: building ${precision_tag} with ${FC}"
    (
        cd "${precision_build_directory}" || exit 1
        ${FC} -cpp ${fortran_free_form_flags} -O0 ${precision_flags} -c "${fortran_build_directory}/module_sf_noahmplsm_harness.F" -o module_sf_noahmplsm_harness.o || exit 1
        for fortran_source_file in harness_text_utilities noahmp_reference_setup; do
            ${FC} -O0 ${precision_flags} -c "${fortran_source_directory}/${fortran_source_file}.f90" -o "${fortran_source_file}.o" || exit 1
        done
        for driver_name in noahmp_stomata_sweep_driver noahmp_canopy_reference_driver; do
            ${FC} -O0 ${precision_flags} "${fortran_source_directory}/${driver_name}.f90" module_sf_noahmplsm_harness.o harness_text_utilities.o noahmp_reference_setup.o -o "${driver_name}_${precision_tag}" || exit 1
        done
    ) > "${precision_build_directory}/build.log" 2>&1 || stop_with_error "Fortran ${precision_tag} build failed, see ${precision_build_directory}/build.log"
    echo "=== Fortran: running ${precision_tag}"
    "${precision_build_directory}/noahmp_stomata_sweep_driver_${precision_tag}" "${effective_configuration_file}" "${fortran_output_directory}" || stop_with_error "Fortran ${precision_tag} stomata sweep failed"
    "${precision_build_directory}/noahmp_canopy_reference_driver_${precision_tag}" "${effective_configuration_file}" "${fortran_output_directory}" || stop_with_error "Fortran ${precision_tag} canopy run failed"
done

# ---------------- C, PET BMI method 6 ----------------
echo "=== C: building with ${CC}"
${CC} -std=c99 -O2 -DBMI_ACTIVE -I"${c_include_directory}" \
    "${c_source_directory}/bmi_pet.c" "${c_source_directory}/pet.c" "${c_source_directory}/noahmp_stomata.c" \
    "${c_source_directory}/ball_berry_state.c" "${c_source_directory}/ball_berry_helpers.c" \
    "${c_source_directory}/ball_berry_canopy_physics.c" "${c_source_directory}/ball_berry_struct_copies.c" \
    "${c_source_directory}/ball_berry_method6_inputs.c" "${c_driver_source_file}" \
    -lm -o "${c_build_directory}/test_ball_berry_method6" > "${c_build_directory}/build.log" 2>&1 \
    || stop_with_error "C build failed, see ${c_build_directory}/build.log"

# Each C solver mode reads the effective configuration with three keys set for that mode:
# the solver selection and the names of its summary and trace files (any existing lines for
# these keys are removed and the three are appended).
for c_mode in ${c_modes_to_run}; do
    echo "=== C: running ${c_mode} solver (BMI driver)"
    if [ "${c_mode}" = "exact" ]; then exact_flag=1; else exact_flag=0; fi
    c_mode_configuration_file="noahmp_comparison/output/C/configuration_${c_mode}.txt"
    {
        sed -e '/^[[:space:]]*bb_noahmp_exact_iteration[[:space:]]*=/d' \
            -e '/^[[:space:]]*bb_diagnostic_summary_file[[:space:]]*=/d' \
            -e '/^[[:space:]]*bb_iteration_trace_file[[:space:]]*=/d' "${effective_configuration_file}"
        echo ""
        echo "# ---------------- set by ${script_name} for the ${c_mode} solver ----------------"
        echo "bb_noahmp_exact_iteration = ${exact_flag}"
        echo "bb_diagnostic_summary_file = ./noahmp_comparison/output/C/c_summary_${c_mode}.csv"
        echo "bb_iteration_trace_file = ./noahmp_comparison/output/C/c_canopy_iteration_trace_${c_mode}.csv"
    } > "${c_mode_configuration_file}"
    "${c_build_directory}/test_ball_berry_method6" "${c_mode_configuration_file}" > "${c_output_directory}/c_bmi_outputs_${c_mode}.txt" \
        || stop_with_error "C ${c_mode} run failed"
    grep "canopy_transpiration_flux_m_per_s" "${c_output_directory}/c_bmi_outputs_${c_mode}.txt"
done

# ---------------- summary and reminder ----------------
echo ""
echo "Done. C solver mode(s) run: ${c_modes_to_run} (from the ${mode_source}). Fortran r4 and r8 run."
echo "Outputs: noahmp_comparison/output/C and output/Fortran. Configuration used: ${effective_configuration_file}"
if [ ${#configuration_overrides[@]} -gt 0 ]; then
    echo "Overrides applied: ${configuration_overrides[*]}"
else
    echo "Overrides applied: none"
fi
modes_not_run=""
for c_mode in exact nested; do
    case " ${c_modes_to_run} " in *" ${c_mode} "*) ;; *) modes_not_run="${modes_not_run} ${c_mode}" ;; esac
done
if [ -n "${modes_not_run}" ]; then
    for c_mode in ${modes_not_run}; do
        echo "REMINDER: the C ${c_mode} solver was not run. To run it:  ./${script_name} ${c_mode}   (or ./${script_name} both)"
    done
fi
echo "To compare C with Fortran:  ./noahmp_comparison/compare_C_to_Fortran.sh"
