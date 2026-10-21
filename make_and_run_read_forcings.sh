#!/bin/bash
# Standalone test, read-forcing mode: PET reads its own forcing file (methods 1-5).
set -ex

PET_SOURCES="./src/pet.c ./src/bmi_pet.c ./src/noahmp_stomata.c ./src/ball_berry_state.c ./src/ball_berry_helpers.c ./src/ball_berry_method6_inputs.c"

gcc -Werror -Og -I./include ./src/main_read_forcing.c ${PET_SOURCES} -lm -o run_bmi_forcings_read
./run_bmi_forcings_read ./configs/pet_config_unit_test1.txt 
./run_bmi_forcings_read ./configs/pet_config_unit_test2.txt 
./run_bmi_forcings_read ./configs/pet_config_unit_test3.txt 
./run_bmi_forcings_read ./configs/pet_config_unit_test4.txt 
./run_bmi_forcings_read ./configs/pet_config_unit_test5.txt 
./run_bmi_forcings_read ./configs/pet_config_cat_67.txt
