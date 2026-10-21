#!/bin/bash
# ==========================================================================================
# make_and_run_read_forcings.sh
#
# Standalone test, read-forcing mode: PET reads its own forcing file. Methods 1-5 only;
# method 6 takes its forcing through the BMI (see test/make_and_run_bmi_unit_test.sh).
#
# The method configuration files set forcing_file = BMI. For this test each one is copied
# to ${standalone_configuration_directory} with three keys changed:
#   forcing_file   = the configuration's forcing_csv_file (PET reads the forcing itself)
#   run_unit_tests = 1 (finalize prints the results next to their expected values)
#   verbose        = 1
# Run from the repository root.  -h prints this usage.
# ==========================================================================================
set -e

# ---------------- directories and files ----------------
configuration_directory="./configs"
standalone_configuration_directory="./test/output/standalone_configs"
source_directory="./src"
include_directory="./include"
executable="./run_bmi_forcings_read"
method_configuration_files="method_1_radiation_balance.cfg method_2_aerodynamic_method.cfg method_3_combination_method.cfg method_4_priestley_taylor_method.cfg method_5_penman_monteith_method.cfg"
catchment_configuration_file="${configuration_directory}/pet_config_cat_67.txt"   # method 5, catchment 67 AORC forcing

PET_SOURCES="${source_directory}/pet.c ${source_directory}/bmi_pet.c ${source_directory}/noahmp_stomata.c ${source_directory}/ball_berry_state.c ${source_directory}/ball_berry_helpers.c ${source_directory}/ball_berry_method6_inputs.c"

if [ "${1:-}" = "-h" ]; then sed -n '3,13p' "$0" | sed 's/^# \{0,1\}//'; exit 0; fi

set -x
gcc -Werror -Og -I${include_directory} ${source_directory}/main_read_forcing.c ${PET_SOURCES} -lm -o ${executable}
set +x

mkdir -p "${standalone_configuration_directory}"
for method_configuration_file in ${method_configuration_files}; do
    source_configuration="${configuration_directory}/${method_configuration_file}"
    standalone_configuration="${standalone_configuration_directory}/${method_configuration_file}"
    forcing_csv_file="$(sed -n 's/^[[:space:]]*forcing_csv_file[[:space:]]*=[[:space:]]*\([^#[:space:]]*\).*/\1/p' "${source_configuration}")"
    sed -e "s|^[[:space:]]*forcing_file[[:space:]]*=.*|forcing_file = ${forcing_csv_file}   # set by make_and_run_read_forcings.sh|" \
        -e "s|^[[:space:]]*run_unit_tests[[:space:]]*=.*|run_unit_tests = 1   # set by make_and_run_read_forcings.sh|" \
        -e "s|^[[:space:]]*verbose[[:space:]]*=.*|verbose = 1   # set by make_and_run_read_forcings.sh|" \
        "${source_configuration}" > "${standalone_configuration}"
    echo "=== ${method_configuration_file} (standalone copy ${standalone_configuration})"
    ${executable} "${standalone_configuration}"
done

echo "=== ${catchment_configuration_file}"
${executable} "${catchment_configuration_file}"
echo
echo "Read-forcing test done: methods 1-5 and catchment 67. Method 6 is BMI only: test/make_and_run_bmi_unit_test.sh"
