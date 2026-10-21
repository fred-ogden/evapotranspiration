#!/bin/bash
# Standalone test, pass-forcing mode: the AORC BMI module passes forcing to PET through the BMI.
set -ex

PET_SOURCES="./src/pet.c ./src/bmi_pet.c ./src/noahmp_stomata.c ./src/ball_berry_state.c ./src/ball_berry_helpers.c ./src/ball_berry_method6_inputs.c"
AORC_SOURCES="./extern/forcing_code/src/aorc.c ./extern/forcing_code/src/bmi_aorc.c"

gcc -I./include ./src/main_pass_forcing.c ${PET_SOURCES} ${AORC_SOURCES} -lm -o run_bmi_forcings_pass
./run_bmi_forcings_pass ./configs/pet_config_bmi.txt ./configs/aorc_config_cat_67.txt 
