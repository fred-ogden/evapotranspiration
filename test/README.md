# PET BMI tests

## BMI unit test, all six methods

```
cd test && ./make_and_run_bmi_unit_test.sh          # all six methods
cd test && ./make_and_run_bmi_unit_test.sh 2 6      # methods 2 and 6 only
cd test && ./make_and_run_bmi_unit_test.sh -h       # usage
```

The script builds `run_pet_bmi_test` from [main_unit_test_bmi.c](./main_unit_test_bmi.c) and runs
it once per method configuration file, `../configs/method_N_*.cfg`. Every configuration file sets
`forcing_file = BMI`; the test reads the first data row of the file named by `forcing_csv_file`
and passes it to PET with `set_value`, as a model engine would.

Each run calls every BMI function and checks what it returns:

| Group | Checks |
|---|---|
| control | initialize, update, update_until, finalize |
| model information | component name, input and output item counts and names |
| variable information | grid, type, itemsize, location, units, nbytes of every variable; unknown names rejected |
| grid information | rank, size, type of grid 0 |
| time information | start, end, step, units; update advances the current time by one step |
| getters | get_value, get_value_at_indices, get_value_ptr agree and are finite for every output |
| setters | every input accepts its current value; index 1, outputs and unknown names are rejected |
| results | methods 1-5 give the unit test PET values (also printed by `pet_unit_tests()`) |
| method 6 validation | out-of-range and NaN inputs (LAI, SHDFAC, root_stress_factor_0_to_1, time_month, ...) are rejected and leave the value unchanged |

Methods 1-5 expose 7 inputs (the AORC forcing) and 1 output (`water_potential_evaporation_flux`).
Method 6 exposes 31 inputs (forcing, vegetation state, coupling, time) and 9 outputs.

## Other tests

| Test | Command (from the repository root) |
|---|---|
| BMI regression and negative tests (method 5) | `./test_bmi.sh` |
| standalone, PET reads its forcing (methods 1-5) | `./make_and_run_read_forcings.sh` |
| standalone, AORC BMI passes the forcing | `./make_and_run_pass_forcings.sh` (needs `git submodule update --init`) |
| everything through CMake | `cmake -S . -B build && cmake --build build && (cd build && ctest --output-on-failure)` |
