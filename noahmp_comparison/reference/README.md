# Stored NWM Noah-MP reference results

Fortran results of the NWM Noah-MP harness (`noahmp_comparison/Fortran`) for the committed
`configs/method_6_ball_berry.cfg` and `forcing/ball_berry_forcing.csv`, one time step
(2016-07-01 12:00:00 UTC, USGS type 15 Mixed Forest, 8 m canopy, no snow, dry canopy,
root_stress_factor_0_to_1 = 1).

| File | Contents |
|---|---|
| `noahmp_canopy_reference_summary_r8.csv` | solution summary, 64-bit REAL (`-fdefault-real-8`) |
| `noahmp_canopy_reference_summary_r4.csv` | solution summary, 32-bit REAL, as NWM runs |
| `noahmp_canopy_iteration_trace_r8.csv`, `..._r4.csv` | VEGE_FLUX loop1, one row per iteration |

Made with:

```
NOAHMP_SOURCE_FILE=/path/to/module_sf_noahmplsm.F FC=flang ./noahmp_comparison/run_ball_berry_C_and_Fortran.sh both
cp noahmp_comparison/output/Fortran/noahmp_canopy_* noahmp_comparison/reference/
```

- NWM `module_sf_noahmplsm.F`: 10198 lines, sha256 `1e84c5867b3b7b5b031283e2eef76722bd9067020697a2cc271b2f8631421205`
  (not in this repository)
- Fortran compiler: LLVM flang 20.1.8, `-O0`. gfortran gives the same results.

`noahmp_comparison/check_C_against_stored_reference.sh` compares the C exact solver with the r8
summary and fails above a relative difference of 1e-12; at this commit the largest is 2.9e-15.
Regenerate these files whenever the configuration or the forcing changes.
