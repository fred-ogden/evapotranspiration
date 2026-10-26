# PET BMI final interface fixes (FLO 10/26)

Files: src/bmi_pet.c, test/bmi_regression.c, test_bmi.sh.

- Correct int timestep printf format (%d instead of %g).
- Reject null names/pointers in indexed getter/setter before strcmp or dereference.
- Reject unknown or invalid variable item counts rather than using num_timesteps as a buffer length.
- Add independent, forcing-file-free negative regression tests.
- Leave the upstream AORC implementation, PET physics, and integration script unchanged.

Run from repository root: ./test_bmi.sh

The configured finite get_end_time behavior is unchanged pending a decision about
whether BMI-driven PET should advertise a finite horizon or an unbounded one.
The legacy test/main_unit_test_bmi.c has separate problems (attempts to write
read-only output and leaks allocations); the new regression test is independent.
