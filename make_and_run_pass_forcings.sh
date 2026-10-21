#!/bin/bash
# ==========================================================================================
# make_and_run_pass_forcings.sh
#
# Standalone test, pass-forcing mode: the AORC BMI module (extern/forcing_code, a git
# submodule) reads the catchment 67 forcing and the driver passes it to PET through the BMI,
# 100 hourly steps. Runs all six methods; for method 6 the driver also passes the calendar
# time inputs. Run "git submodule update --init" once first.
#
# Each method configuration file is copied to ${pass_configuration_directory} with the site of
# catchment 67 (latitude, longitude, elevation) and its number of time steps.
# Run from the repository root.  Usage: ./make_and_run_pass_forcings.sh [-h] [method_number ...]
# ==========================================================================================
set -e

# ---------------- directories and files ----------------
configuration_directory="./configs"
pass_configuration_directory="./test/output/pass_forcing_configs"
source_directory="./src"
include_directory="./include"
aorc_directory="./extern/forcing_code"
aorc_configuration_file="${configuration_directory}/aorc_config_cat_67.txt"
executable="./run_bmi_forcings_pass"

# ---------------- catchment 67 site ----------------
catchment_latitude_degrees="37.25"
catchment_longitude_degrees="-97.5554"
catchment_site_elevation_m="303.33"
catchment_num_timesteps="720"

PET_SOURCES="${source_directory}/pet.c ${source_directory}/bmi_pet.c ${source_directory}/noahmp_stomata.c ${source_directory}/ball_berry_state.c ${source_directory}/ball_berry_helpers.c ${source_directory}/ball_berry_method6_inputs.c"
AORC_SOURCES="${aorc_directory}/src/aorc.c ${aorc_directory}/src/bmi_aorc.c"

if [ "${1:-}" = "-h" ]; then sed -n '3,13p' "$0" | sed 's/^# \{0,1\}//'; exit 0; fi
methods_to_run="${*:-1 2 3 4 5 6}"
[ -f "${aorc_directory}/src/aorc.c" ] || { echo "ERROR: ${aorc_directory} is empty: run git submodule update --init"; exit 1; }

set -x
gcc -I${include_directory} ${source_directory}/main_pass_forcing.c ${PET_SOURCES} ${AORC_SOURCES} -lm -o ${executable}
set +x

mkdir -p "${pass_configuration_directory}"
for method_number in ${methods_to_run}; do
    method_configuration="$(ls ${configuration_directory}/method_${method_number}_*.cfg)"
    pass_configuration="${pass_configuration_directory}/$(basename "${method_configuration}")"
    sed -e "s|^[[:space:]]*latitude_degrees[[:space:]]*=.*|latitude_degrees = ${catchment_latitude_degrees}   # set by make_and_run_pass_forcings.sh|" \
        -e "s|^[[:space:]]*longitude_degrees[[:space:]]*=.*|longitude_degrees = ${catchment_longitude_degrees}   # set by make_and_run_pass_forcings.sh|" \
        -e "s|^[[:space:]]*site_elevation_m[[:space:]]*=.*|site_elevation_m = ${catchment_site_elevation_m}   # set by make_and_run_pass_forcings.sh|" \
        -e "s|^[[:space:]]*num_timesteps[[:space:]]*=.*|num_timesteps = ${catchment_num_timesteps}   # set by make_and_run_pass_forcings.sh|" \
        "${method_configuration}" > "${pass_configuration}"
    echo "=== method ${method_number}: ${pass_configuration} with ${aorc_configuration_file}"
    ${executable} "${pass_configuration}" "${aorc_configuration_file}"
done
echo
echo "Pass-forcing test done for method(s): ${methods_to_run}"
