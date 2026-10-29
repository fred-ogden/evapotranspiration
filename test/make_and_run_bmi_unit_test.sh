#!/bin/bash
# ==========================================================================================
# make_and_run_bmi_unit_test.sh
#
# BMI unit test of the PET module for all six methods. Builds run_pet_bmi_test and runs it
# once per method configuration file. Each run calls every BMI function, passes the first row
# of the configuration's forcing_csv_file through set_value, and checks the results (see
# main_unit_test_bmi.c). Runs from the test directory; returns non-zero if any method fails.
# Usage: ./make_and_run_bmi_unit_test.sh [-h] [method_number ...]   (default: all six)
# ==========================================================================================
set -e
cd "$(dirname "$0")"

# ---------------- directories and files (relative to test/) ----------------
repository_root_directory=".."
source_directory="${repository_root_directory}/src"
include_directory="${repository_root_directory}/include"
configuration_directory="${repository_root_directory}/configs"
executable="./run_pet_bmi_test"

PET_SOURCES="${source_directory}/bmi_pet.c ${source_directory}/pet.c ${source_directory}/noahmp_stomata.c ${source_directory}/ball_berry_state.c ${source_directory}/ball_berry_helpers.c ${source_directory}/ball_berry_canopy_physics.c ${source_directory}/ball_berry_struct_copies.c ${source_directory}/ball_berry_method6_inputs.c"

if [ "${1:-}" = "-h" ]; then sed -n '3,10p' "$0" | sed 's/^# \{0,1\}//'; exit 0; fi
methods_to_run="${*:-1 2 3 4 5 6}"

set -x
gcc -Werror -Og -I${include_directory} ./main_unit_test_bmi.c ${PET_SOURCES} -lm -o ${executable}
set +x

# configuration paths inside the files are relative to the repository root
cd "${repository_root_directory}"
failed_methods=""
for method_number in ${methods_to_run}; do
    configuration_file="$(ls configs/method_${method_number}_*.cfg)"
    echo "=== method ${method_number}: ${configuration_file}"
    ./test/${executable#./} "${configuration_file}" || failed_methods="${failed_methods} ${method_number}"
done

echo
if [ -n "${failed_methods}" ]; then echo "BMI unit test FAILED for method(s):${failed_methods}"; exit 1; fi
echo "BMI unit test passed for method(s): ${methods_to_run}"
not_run=""
for method_number in 1 2 3 4 5 6; do case " ${methods_to_run} " in *" ${method_number} "*) ;; *) not_run="${not_run} ${method_number}";; esac; done
[ -n "${not_run}" ] && echo "Not run this time: method(s)${not_run}"
exit 0
