#!/bin/sh
# FLO 10/26: Run isolated BMI interface and negative tests without AORC CSVs.
set -eu
cd "$(dirname "$0")"
CC=${CC:-gcc}
"$CC" -std=c99 -O0 -g -Wall -Wno-unused-parameter -Iinclude \
    test/bmi_regression.c src/bmi_pet.c src/pet.c src/noahmp_stomata.c src/ball_berry_state.c src/ball_berry_helpers.c src/ball_berry_method6_inputs.c -lm -o test/run_bmi_regression
./test/run_bmi_regression configs/method_5_penman_monteith_method.cfg
