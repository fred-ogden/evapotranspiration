#!/bin/bash

set -ex

gcc -Werror -Og -I../include ./main_unit_test_bmi.c ../src/bmi_pet.c ../src/pet.c ../src/noahmp_stomata.c ../src/ball_berry_state.c ../src/ball_berry_helpers.c ../src/ball_berry_method6_inputs.c -lm -o run_pet_bmi_test
./run_pet_bmi_test ../configs/pet_config_bmi_unit_test.txt
#./run_pet_bmi_test ../configs/pet_config_cat_67.txt
