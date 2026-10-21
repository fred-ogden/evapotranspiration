# PET method 6: Ball-Berry canopy transpiration (Noah-MP physics)

PET method 6 is a C99 port of the Noah-MP Ball-Berry canopy energy balance and transpiration
calculation, built into the PET BMI module as `pet_method = 6`. It is checked against the
Noah-MP code that runs in the current National Water Model (`module_sf_noahmplsm.F`) with a
Fortran harness in `noahmp_comparison/`. Both codes read the same configuration file,
`configs/method_6_ball_berry.cfg`, and the same forcing file, `forcing/ball_berry_forcing.csv`.

It is the first of the modules intended to replace Noah-OWP-Modular (NOM) in NextGen without
NOM's soil moisture (Richards equation) solver: this Ball-Berry module, a canopy interception,
evaporation, sublimation and melt module (CIESM), a snow ablation module, soil freeze-thaw (SFT)
with the soil skin temperature, and a bare soil evaporation module, coupled to a hydrologic
model such as CFE3.1.

## Layout

```
configs/method_6_ball_berry.cfg         the configuration, read by both codes
forcing/ball_berry_forcing.csv          AORC-format forcing (one time step)
src/ball_berry_helpers.c                the canopy solvers (exact VEGE_FLUX replica, nested solver)
src/ball_berry_method6_inputs.c         configuration parsing and per-step input assembly
src/noahmp_stomata.c, ball_berry_state.c   STOMATA, LAI/SAI/SHDFAC validation
src/bmi_pet.c, src/pet.c                the PET BMI module, all six methods
test/main_ball_berry_bmi.c              BMI-only method 6 driver (ctest ball_berry_method6)
test/main_unit_test_bmi.c               BMI unit test, all six methods
noahmp_comparison/
  run_ball_berry_C_and_Fortran.sh       builds and runs both codes (needs the NWM source)
  compare_C_to_Fortran.sh               compares output/C with output/Fortran, or with -r reference/
  check_C_against_stored_reference.sh   C exact solver vs the stored r8 results (no NWM source needed)
  compare_c_to_noahmp_reference.py      helper called by both compare scripts
  Fortran/                              patch script and Fortran drivers (no NWM source)
  reference/                            stored Fortran r4 and r8 results for the committed configuration
  build/, output/                       created by the scripts, not in git
```

The NWM `module_sf_noahmplsm.F` is not in this repository. Give its path in `noahmp_source_file`
at the top of `run_ball_berry_C_and_Fortran.sh` or with the environment variable
`NOAHMP_SOURCE_FILE`.

## Requirements

bash, a C99 compiler (default `gcc`), python3, and for the Fortran reference a Fortran compiler
(default `gfortran`; LLVM `flang` also works) and the NWM source.

## Running

Run from the repository root.

```
./noahmp_comparison/check_C_against_stored_reference.sh   # C vs stored Noah-MP r8, pass/fail (CI runs this)
NOAHMP_SOURCE_FILE=/path/to/module_sf_noahmplsm.F ./noahmp_comparison/run_ball_berry_C_and_Fortran.sh both
./noahmp_comparison/run_ball_berry_C_and_Fortran.sh            # C solver mode from the configuration file
./noahmp_comparison/run_ball_berry_C_and_Fortran.sh exact      # override for this run: Noah-MP exact iteration
./noahmp_comparison/run_ball_berry_C_and_Fortran.sh nested     # override for this run: nested solver
./noahmp_comparison/run_ball_berry_C_and_Fortran.sh both -s root_stress_factor_0_to_1=0.5    # override a key
./noahmp_comparison/run_ball_berry_C_and_Fortran.sh -h         # usage
./noahmp_comparison/run_ball_berry_C_and_Fortran.sh -l         # list the keys -s can override, with current values
FC=flang ./noahmp_comparison/run_ball_berry_C_and_Fortran.sh both   # choose compilers with CC and FC
./noahmp_comparison/compare_C_to_Fortran.sh             # compare every C mode that has results
./noahmp_comparison/compare_C_to_Fortran.sh -r exact    # compare with the stored reference; -h for usage
```

In tcsh, set the source path with `env NOAHMP_SOURCE_FILE=/path/to/module_sf_noahmplsm.F ./noahmp_comparison/run_ball_berry_C_and_Fortran.sh both`,
or edit `noahmp_source_file` at the top of the script once.

`-s key=value` replaces the value of an existing configuration key for that run only, and can
be repeated. `-l` lists every key in the configuration file with its current value and comment,
under the file's section headings; it reads the file each time, so it is always current. Three
keys are set by the run script for each C solver mode (`bb_noahmp_exact_iteration`,
`bb_diagnostic_summary_file`, `bb_iteration_trace_file`); use the exact|nested|both argument
instead of overriding them. Forcing values are in the forcing CSV, not the configuration, so
they are changed by editing that file. The configuration file is never edited: the run script writes
`noahmp_comparison/output/effective_configuration.txt` (marking each overridden line) and both
codes read it. A key that is not in the configuration file is an error.

Each run first deletes the results of the previous run from `noahmp_comparison/output/C` and
`output/Fortran`, so everything there comes from one run with one configuration. The run script
always rebuilds and reruns the Fortran reference in single (r4, 32-bit REAL as NWM) and double
(r8, `-fdefault-real-8`) precision. At the end, both scripts print the overrides used and a
reminder of any C solver mode that was not run or not compared.

The C test driver can also be run directly, from the repository root:
`build/test_ball_berry_method6 configs/method_6_ball_berry.cfg` (CMake build). It uses only BMI
functions (see "Model life cycle" below). The run script gives each C solver mode its own copy of
the effective configuration, `noahmp_comparison/output/C/configuration_<mode>.txt`, in which it
sets `bb_noahmp_exact_iteration` and the names of the summary and trace files for that mode.

## Outputs

| File | Contents |
|---|---|
| noahmp_comparison/output/C/c_summary_exact.csv, c_summary_nested.csv | C solution, Noah-MP names, one value per row (written by the model, key bb_diagnostic_summary_file) |
| noahmp_comparison/output/C/c_bmi_outputs_exact.txt, c_bmi_outputs_nested.txt | the nine method 6 BMI outputs, as the driver received them with get_value |
| noahmp_comparison/output/C/configuration_exact.txt, configuration_nested.txt | the configuration each C solver mode read |
| noahmp_comparison/output/C/c_canopy_iteration_trace_exact.csv | one row per canopy iteration (exact solver only) |
| noahmp_comparison/output/Fortran/noahmp_canopy_reference_summary_r4.csv, _r8.csv | Fortran solution: name, value, units, description |
| noahmp_comparison/output/Fortran/noahmp_canopy_iteration_trace_r4.csv, _r8.csv | one row per VEGE_FLUX canopy (loop1) iteration |
| noahmp_comparison/output/Fortran/noahmp_stomata_sweep_r4.csv, _r8.csv | STOMATA alone over a grid of APAR, TV, EA, RB; every input on each row |

## The configuration file

`configs/method_6_ball_berry.cfg` has the layout shared by all six method files (control, time,
site, then the method's parameters). Each code ignores keys it does not use, so one edit changes
both codes. In this file values must be plain numbers (no `[units]` annotations, because the
Fortran reads it too) except `forcing_file`, `forcing_csv_file`, `bb_iteration_trace_file` and
`bb_diagnostic_summary_file`. Relative paths are relative to the repository root. The two
diagnostic file keys are commented out: without them method 6 writes no files.

| Keys | Read by |
|---|---|
| verbose, pet_method, run_unit_tests, yes_aorc, forcing_file (= BMI) | C (PET BMI) |
| bb_noahmp_exact_iteration, bb_noahmp_compatibility, bb_iteration_trace_file, bb_diagnostic_summary_file | C |
| forcing_csv_file, num_timesteps | the drivers: the C test driver reads the forcing and passes it through BMI; the Fortran driver reads the first row |
| bb_maximum_stability_iterations, bb_maximum_canopy_iterations, bb_stability_relative_tolerance, bb_canopy_temperature_tolerance_k, bb_maximum_temperature_step_k | C, nested solver only |
| time_step_size_s, latitude_degrees (sign only), and all physics keys | both |
| coupling inputs (root_stress_factor_0_to_1 through canopy_ice_mm) | both: C initial values of the BMI inputs; Fortran imposed values |
| usgs_vegetation_type, tree_crown_radius_m, ground_surface_relative_humidity | Fortran only (see the comments in the file) |

Physics keys use long descriptive names with units; the Noah-MP variable name is in brackets
in the comment on each line. Vegetation values are the USGS type 15 (Mixed Forest) column of
NOM_MPTABLE.TBL, except canopy_top_height_m = 8.0 (table 16.0) and canopy_bottom_height_m = 5.0
(table 10.0).

`use_monthly_lai_sai_table = 1` interpolates LAI and SAI from the monthly table to the forcing
date exactly as Noah-MP PHENOLOGY does (DVEG = 4). With 0, `leaf_area_index` and
`stem_area_index` are used.

Both codes run one time step: the Fortran driver reads the first forcing row only.

## C solver modes

| Mode | Key value | What it does | Cost per call (gcc -O2) |
|---|---|---|---|
| exact | bb_noahmp_exact_iteration = 1 | line-by-line replica of VEGE_FLUX loop1 (`noahmp_vege_flux_exact_solve`) | about 4 us |
| nested | bb_noahmp_exact_iteration = 0 | outer Picard iteration on the resistances around a fully converged inner Newton solve of TV (`noahmp_stability_coupled_solve`) | about 150 us |

The exact solver reproduces NWM Noah-MP to machine precision against the double-precision (r8)
Fortran reference (largest relative difference about 5e-15), and to about 1e-5 or better against
the single-precision (r4, as NWM) reference, with the same number of canopy iterations. The nested solver converges further than
Noah-MP does, so it differs from NWM by up to about 4e-4 in the fluxes and 4e-6 in TV. Use exact
unless there is a reason to want the more fully converged answer.

`bb_noahmp_compatibility = 1` makes the nested solver freeze the stomatal resistances at the
first iteration, as Noah-MP does. The exact solver always does this.

## Model life cycle (method 6)

Method 6 follows the BMI life cycle a model engine expects:

- **initialize**: reads the configuration file once. All parameters and constants that are not
  BMI variables go into a per-instance parameter struct (`BallBerryMethod6Parameters`, including
  the 12-month LAI and SAI tables); the physiology goes into the STOMATA parameter struct; the
  BMI input variables get their initial values. The parsed text is freed before initialize
  returns, and nothing is stored in static variables, so several instances (catchments) can run
  in one process independently. Method 6 requires `forcing_file = BMI`.
- **set_value**: the driver sets the seven AORC forcing variables and the six time variables
  before every update, and the coupling inputs when it has them. Out-of-range values are
  rejected and leave the stored value unchanged.
- **update**: builds the canopy input from the current BMI values and the parameter struct,
  solves, and sets the outputs. It reads no files (except writing the optional summary and trace).
- **get_value**: the outputs listed below.
- **finalize**: frees everything the model allocated. Checked with valgrind: no leaks, for both
  solver modes and for a failed initialize.

## Coupling design

### Root stress (unique to this module)

Noah-MP's soil moisture factor BTRAN is not the ratio of actual to potential transpiration. It
enters only STOMATA, where it scales the maximum carboxylation rate,
VCMX = VCMX25/F2(TC)*FNF*BTRAN*F1(AVCMX,TC). Transpiration then responds through photosynthesis,
the Ball-Berry conductance, the cuticular limit set by BP, and the leaf energy balance. For the
test case:

| BTRAN | 1.00 | 0.75 | 0.50 | 0.25 | 0.10 | 0.00 |
|---|---|---|---|---|---|---|
| TR/TR(BTRAN=1) | 1.000 | 0.925 | 0.761 | 0.420 | 0.191 | 0.025 |

So the coupling is:

1. A hydrologic model coupled to this module must expose its root stress factor, 0 to 1, on its
   BMI. This module reads it as `root_stress_factor_0_to_1` (Noah-MP BTRAN).
2. This module applies it inside STOMATA and computes the transpiration consistent with its leaf
   energy balance, `canopy_transpiration_flux_m_per_s`. That is the actual transpiration.
3. The hydrologic model removes that amount from its root zone, distributed by its own method
   (for CFE3.1: from the wettest disc, or one quarter from each disc), and does NOT apply its
   stress factor to it again.

The stress factor from the previous time step is what Noah-MP itself uses: in NOAHMP_SFLX,
ENERGY (BTRAN and VEGE_FLUX) runs before WATER, so BTRAN is computed from the soil moisture at
the start of the step. If no stress factor is connected, it stays at its configured value (1),
and the output is the unstressed (potential) transpiration; the hydrologic model may then apply
its own stress, at the cost of a leaf energy balance that assumed no stress.

The other PET methods (1 to 5) return potential evapotranspiration, and the receiving model
applies its own root stress as before.

The meaning of the stress factor depends on the hydrologic model that supplies it. Noah-MP's own
form (OPT_BTR = 1) is the root-disc thickness-weighted sum of
GX = (SH2O - SMCWLT)/(SMCREF - SMCWLT), clipped to [0,1].

### Order of execution in an ngen realization

1. This module (Ball-Berry): uses the previous step's ground skin temperature (SFT), surface
   resistance (bare soil evaporation) and root stress factor (hydrologic model), and the current
   snow and canopy water state.
2. Soil freeze-thaw (skin temperature) and bare soil evaporation: use this step's canopy air
   temperature, canopy air vapor pressure, under-canopy resistance and absorbed radiation.
3. The hydrologic model: removes this step's transpiration.

This is the same sequential splitting Noah-MP already uses inside VEGE_FLUX: loop1 solves the
canopy with TG fixed at its start-of-step value, then loop2 solves TG with TAH and EAH fixed.
When snow is deeper than 0.05 m and TG is above freezing, Noah-MP (OPT_STC = 3) then resets TG to
(1-FSNO)*TG + FSNO*TFRZ; that belongs to the skin temperature module, not to this one.

### Canopy evaporation

Transpiration TR and the evaporation of intercepted water EVC come from the same leaf
temperature solve in VEGE_FLUX, so this module computes both. CIESM supplies the wetted
fraction and the intercepted liquid and ice, receives EVC, and owns the canopy water storage
mass balance. Until CIESM exists, these inputs keep their configured values (a dry canopy).

### Tiles

This module solves only the vegetated tile. The bare soil fraction (1 - FVEG - impervious
fraction in CFE3.1) is handled elsewhere: SFT runs both tiles, and bare soil evaporation runs
once per tile (with TAH, EAH and RAHG from this module for the vegetated tile, and the
reference-height air and a bare-ground surface-layer resistance for the bare tile). The
bare-ground resistance (SFCDIF1 over bare soil in Noah-MP BARE_FLUX) has no home yet.

Radiation follows NWM (OPT_RAD = 3, gap fraction 1 - FVEG), so SAG is a single grid-mean value
that Noah-MP applies to both the under-canopy and the bare ground.

## BMI variables

All variables are `double`, one value, grid 0, location node. Long descriptive names are used at
the BMI; inside the Noah-MP science code the Noah-MP names are kept for provenance, and the
mapping is made at the start of the method 6 input builder (`bb6_build_canopy_input`)
and in `run_pet` for the outputs.

### Inputs

The time inputs are UTC at the start of the time step (the forcing time stamp). They are
doubles holding whole numbers (`time_second` may be fractional); set_value rejects anything else,
and update checks that the day exists in the month. Method 6 uses them only for the monthly
LAI/SAI table (`use_monthly_lai_sai_table = 1`); otherwise they are ignored. They are NaN until
set, and an update that needs them before they are set fails with a message.

| BMI name | Noah-MP | Units | Supplied by |
|---|---|---|---|
| land_surface_radiation~incoming~longwave__energy_flux | LWDN | W m-2 | forcing (AORC DLWRF) |
| land_surface_air__pressure | SFCPRS | Pa | forcing (AORC PRES) |
| atmosphere_air_water~vapor__specific_humidity | QAIR | kg kg-1 | forcing (AORC SPFH) |
| land_surface_radiation~incoming~shortwave__energy_flux | SOLDN | W m-2 | forcing (AORC DSWRF) |
| land_surface_air__temperature | SFCTMP | K | forcing (AORC TMP) |
| land_surface_wind__x_component_of_velocity | UU | m s-1 | forcing (AORC UGRD) |
| land_surface_wind__y_component_of_velocity | VV | m s-1 | forcing (AORC VGRD) |
| LAI | LAI | 1 | vegetation model or data assimilation; monthly table if configured |
| SAI | SAI | 1 | as LAI |
| SHDFAC | FVEG | 1 | vegetation model or data assimilation |
| CWP | CWPVT | 1 | configuration, calibration |
| MP | MP | 1 | configuration, calibration |
| VCMX25 | VCMX25 | umol CO2 m-2 s-1 | configuration, calibration |
| root_stress_factor_0_to_1 | BTRAN | 1 | hydrologic model, previous step |
| ground_skin_temperature_K | TG | K | soil freeze-thaw, previous step |
| surface_resistance_s_per_m | RSURF | s m-1 | bare soil evaporation, previous step |
| snow_cover_fraction | FSNO | 1 | snow module |
| snow_depth_m | SNOWH | m | snow module |
| snow_albedo_direct_visible | ALBSND(1) | 1 | snow module |
| snow_albedo_direct_nir | ALBSND(2) | 1 | snow module |
| snow_albedo_diffuse_visible | ALBSNI(1) | 1 | snow module |
| snow_albedo_diffuse_nir | ALBSNI(2) | 1 | snow module |
| canopy_wetted_fraction | FWET | 1 | CIESM |
| canopy_liquid_water_mm | CANLIQ | mm | CIESM |
| canopy_ice_mm | CANICE | mm | CIESM |
| time_year | | year | driver |
| time_month | | month (1 to 12) | driver |
| time_day_of_month | | d (1 to 31) | driver |
| time_hour | | h (0 to 23) | driver |
| time_minute | | min (0 to 59) | driver |
| time_second | | s (0 to under 60) | driver |

### Outputs

| BMI name | Noah-MP | Units | Used by |
|---|---|---|---|
| water_potential_evaporation_flux | | m s-1 | PET methods 1 to 5: PET. Method 6: canopy transpiration, set to zero if negative (kept for compatibility) |
| canopy_transpiration_flux_m_per_s | ETRAN = TR/(LATHEAV*rho_water) | m s-1 | hydrologic model |
| canopy_evaporation_heat_flux_W_per_m2 | EVC | W m-2 | CIESM |
| ground_shortwave_absorbed_below_canopy_W_per_m2 | SAG | W m-2 | SFT |
| ground_longwave_absorbed_below_canopy_W_per_m2 | EMG*((1-EMV)*LWDN + EMV*SB*TV^4) | W m-2 | SFT |
| canopy_air_temperature_K | TAH | K | SFT, bare soil evaporation |
| canopy_air_vapor_pressure_Pa | EAH | Pa | SFT, bare soil evaporation |
| below_canopy_aerodynamic_resistance_s_per_m | RAHG (= RAWG) | s m-1 | SFT, bare soil evaporation |
| vegetation_temperature_K | TV | K | CIESM, snow module |

Conventions:

- TR and EVC are per unit grid-cell area: Noah-MP includes the FVEG weight in both. Positive
  transpiration and evaporation are water leaving the surface; the sign is passed through as
  Noah-MP computes it.
- The absorbed longwave is the absorbed part of Noah-MP's IRG for the under-canopy ground, at the
  final TV. The ground's own emission, EMG*SB*TG^4, belongs to SFT. Ground emissivity is
  EMG = EG*(1-FSNO) + FSNO.
- Outputs are NaN until the first update.

### How the inputs are treated

- Initial values of the vegetation and coupling inputs come from the configuration file during
  BMI Initialize, from keys with the same names (and `leaf_area_index`, `stem_area_index`, `green_vegetation_fraction`,
  `canopy_wind_parameter`, `conductance_photosynthesis_slope`,
  `maximum_carboxylation_rate_25c_umol_per_m2_per_s` for the six vegetation variables).
- Every update uses the current BMI values, so `set_value` after Initialize takes effect.
- `ground_skin_temperature_K` starts as NaN. Until it is set, TG = air temperature +
  `ground_temperature_offset_from_air_k`, the cold start used by the test.
- `set_value` rejects values out of range and leaves the stored value unchanged: fractions and
  the root stress factor in [0,1], albedos in (0,1), depths, water amounts and RSURF
  non-negative, TG between 150 and 400 K. LAI, SAI and SHDFAC are checked by
  `ball_berry_state_set`.
- With `use_monthly_lai_sai_table = 1`, LAI and SAI are interpolated from the table every step
  and written back to the BMI variables, overriding any set_value. Set the key to 0 when another
  model supplies them.
- With snow, the provider of RSURF supplies the snow-weighted value Noah-MP uses with OPT_RSF = 4,
  RSURF = 1/(FSNO/RSURF_SNOW + (1-FSNO)/RSURF_soil).

### Snow, as in Noah-MP

- Burial (PHENOLOGY): the exposed indices ELAI and ESAI are LAI and SAI times (1 - FB), with
  FB from SNOWH, HVB and HVT (a different form when HVT <= 1 m); small values are set to zero.
- ENERGY: ZPD is raised to SNOWH when the snow is deeper than 0.65*HVT, and ZLVL is raised if the
  snow surface reaches it. Z0MG blends the bare soil (0.002 m) and snow (Z0SNO = 0.002 m)
  roughness by FSNO.
- GROUNDALB: ground albedo = soil albedo*(1-FSNO) + snow albedo*FSNO, for direct and diffuse,
  visible and near infrared. The snow-free soil albedos are imposed by the configuration.

### Why LAI, SAI, SHDFAC, CWP, MP and VCMX25 are exposed

- LAI, SAI and SHDFAC are canopy state that changes through the season and can come from another
  model or from data assimilation instead of the static monthly table.
- MP and VCMX25 control stomatal conductance and photosynthesis and differ between vegetation
  types in MPTABLE (MP is 6 or 9, VCMX25 is 40 to 80), so they are the physiology parameters for
  parameter estimation. The other STOMATA parameters (KC25, AKC, KO25, AKO, AVCMX, QE25, C3PSN)
  have the same value for every vegetated type.
- CWP sets the within-canopy wind attenuation, and through it RB and RAHG; it also varies by
  vegetation type and is empirical.

## Differences from a full NWM NOAHMP_SFLX call (both codes)

1. The snow-free soil albedos are imposed from the configuration instead of computed from soil
   color and soil moisture by GROUNDALB; snow albedos are inputs instead of computed from snow age.
2. The ground skin temperature TG is an input. The Fortran VEGE_FLUX ground loop (loop2)
   evaluates IRG, SHG and EVG once at that TG, and GH is the residual SAG - IRG - SHG - EVG; soil
   temperature and thermal conductivity do not enter.
3. BTRAN and RSURF are inputs instead of computed from Noah-MP's soil moisture.
4. FSNO, SNOWH, FWET, CANLIQ and CANICE are inputs instead of computed by the snow and canopy
   water routines.
5. Cold start: TV = TAH = TMP and EAH = the air vapor pressure. Because STOMATA is called only on
   canopy iteration 1, RSSUN and RSSHA depend directly on these initial values.
6. ZLVL = HVT + ZREF with ZREF = 2 m, so the forcing is 10 m above the ground. NWM with a ZLVL
   namelist value of 10 m would place it at HVT + 10 m = 18 m.
7. COSZ is fixed by the configuration.

## Noah-MP behavior the C code reproduces

- STOMATA is called once, on canopy iteration 1, with the incoming TV, ESAT(TV), EAH and the
  iteration 1 RB, and the supplied BTRAN. RS is not updated as TV converges.
- Canopy loop: at most NITERC = 20 iterations. Once ITER >= 5 and |DTV| <= 0.01 K, one more
  iteration is taken and the loop exits. STOMATA Ci loop: NITER = 3, no tolerance.
- Each canopy iteration: SFCDIF1 with H from the previous iteration, RAGRB with HG and TAH from the
  previous iteration, ESAT at TV, one linearized Newton step for TV.
  RAHC = RAWC = MAX(1, 1/(CH*UR)). EVC is limited to CANLIQ*LATHEAV/DT (CANICE below freezing).
- Z0H = Z0M and Z0HG = Z0MG. UR = MAX(wind speed, 1 m/s). VAIE = MIN(6, ELAI + ESAI).
- ESAT argument TDC(T) = MIN(50, MAX(-50, T - TFRZ)) with TFRZ = 273.16 K.
- Constants: TFRZ = 273.16, HVAP = 2.5104E06, HSUB = 2.8440E06, CPAIR = 1004.64, RAIR = 287.04,
  VKC = 0.40, GRAV = 9.80616, SB = 5.67E-08, MPE = 1.E-6.
- Radiation: OPT_RAD = 3 (NWM), gap fraction = 1 - FVEG; HVB and RC are not used by the radiation.

## The Fortran reference harness

The NWM `module_sf_noahmplsm.F` is used as supplied, unmodified, and is not in this repository.
The run script runs `noahmp_comparison/Fortran/make_harness_module.py`, which writes
`noahmp_comparison/build/Fortran/module_sf_noahmplsm_harness.F` with
five changes, each marked HARNESS PATCH: RADIATION, VEGE_FLUX, STOMATA and ESAT made public;
driver control variables added; imposed soil and snow albedos blended by FSNO in GROUNDALB; a
per-iteration trace in VEGE_FLUX loop1; and the prescribed-TG option in loop2. The patch script
stops if any change does not apply exactly once. The Fortran summary reports both the prescribed
TG and `TG_AFTER_VEGE_FLUX`, the value after Noah-MP's snow reset.

## Results

Mixed Forest, HVT = 8 m, forcing 2016-07-01 12:00, C exact solver driven through the BMI.
Largest relative difference over all compared variables against the double-precision (r8)
Fortran reference:

| Case | Overrides or forcing change | Iterations (C, Fortran) | Largest relative difference vs r8 |
|---|---|---|---|
| baseline (285.3 K, relative humidity about 99.8 percent) | none | 7, 7 | 2.9e-15 |
| root stress | root_stress_factor_0_to_1=0.5 | 7, 7 | 5.3e-15 |
| snow and wet canopy | snow_cover_fraction=0.4 snow_depth_m=6.0 canopy_wetted_fraction=0.3 canopy_liquid_water_mm=0.02 | 9, 9 | 2.8e-15 |
| warm and drier | TMP = 296.3 K in the forcing (relative humidity about 50 percent) | 7, 7 | 0 |

Against the single-precision (r4) reference the differences are about 1e-5 or less, set by
32-bit roundoff (largest: EVC 5.3e-5 in the snow case, where the CANLIQ limit is active).

The snow case buries one third of the canopy (SNOWH 6 m against HVB 5 m and HVT 8 m) and raises
ZPD to the snow surface; it is a code test, not a realistic July state. The C test driver parses
the forcing in double precision. With the forcing parsed in single precision instead, the
BMI-driven C code reproduces the earlier file-driven results digit for digit, which confirms the
move to BMI forcing changed nothing else. The gfortran and LLVM flang builds of the reference give
identical results.

The forcing is nearly saturated (relative humidity about 99.8 percent), so it barely exercises
the vapor pressure deficit term in STOMATA. A drier, warmer case is the next test to add.

## Running more than one time step

`make_and_run_pass_forcings.sh` runs method 6 for 100 hourly AORC steps, with the calendar time
inputs set by the driver. Until the two items below are done, those steps are a test of the
coupling, not a Noah-MP simulation:

- `cosine_solar_zenith_angle` is a fixed configuration value. It should be computed each step from
  the time inputs and the site latitude and longitude.
- The canopy state is not carried between steps (below).

## Known issues, not yet changed

- Nested solver diagnostics: RS is frozen at iteration 1, but PSNSUN and PSNSHA are re-evaluated
  at the final TV, so the reported photosynthesis differs from Noah-MP. Transpiration is not
  affected. The compare script therefore reports PSNSUN as the largest difference for the
  nested mode.
- The canopy state (TV, TAH, EAH) is not carried from one time step to the next: every update is
  a cold start (TV = TAH = air temperature, EAH = air vapor pressure), or uses the configured
  initial values. Noah-MP carries them as prognostic state, and because STOMATA uses the incoming
  TV and EAH, multi-step runs will differ from Noah-MP until they are carried. This is the next
  change needed before multi-step runs.
- `parse_aorc_line_pet` does not subtract 1 from tm_mon and uses local-time mktime, so the stored
  forcing time is off by a month. Only methods 1 to 5 (file forcing) use it; method 6 gets its
  time through the BMI. The AORC BMI module (`extern/forcing_code`, `aorc.c`) has the same bug, so
  `main_pass_forcing.c` computes the calendar time from the first time stamp of the forcing file
  and the AORC BMI current time instead of using AORC's stored time.
