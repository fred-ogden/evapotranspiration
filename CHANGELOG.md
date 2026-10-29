All notable changes to this project will be documented in this file.
We follow the [Semantic Versioning 2.0.0](http://semver.org/) format.


## Unreleased (stack-phase1, dcfe-noahmp-stack) - 2026-10-28

### Changed

- PET method 6, the exact solver (bb_noahmp_exact_iteration = 1): one step is a preparation
  (`bb6_prepare_step_inputs`: the checks, solar geometry, LAI/SAI, NaN conventions, messages) and
  the physics on resolved inputs (`ball_berry_method6_physics_step`, new
  `src/ball_berry_canopy_physics.c`), in the form Tapenade differentiates: one return per function,
  no printf, isnan, ternaries, break, NAN, whole-structure assignment or memset (generated
  field-by-field copies, `src/ball_berry_struct_copies.c`), and every chain of multiplications and
  divisions written one operation per statement. The exact solver records its iteration trace and
  the caller writes it. EVC = MIN(latent limit, EVC) is a comparison (fmin of +0 and -0 may return
  either zero). The nested solver is unchanged and uses the same canopy input.
- Bitwise the same as before: the 13 ctests and their output files, a replay of the
  dcfe-noahmp-stack inputs of the six 01123000 divides over 2012-10-01 to 2019-07-12 (356,544
  divide-hours, every coupling output), the CIESM coupled ctests and the stack's core path test.

### Added

- BMI input `RSURF_SNOW` (s m-1): the snow surface resistance of the OPT_RSF = 4 RSURF blend, settable as a
  model_params value for parameter estimation; an independent of the Tapenade head.
- `tapenade/canopy_energy/`: Tapenade 3.16 tangent and reverse of one step (32 independents,
  14 dependents) with p/t/a tests on the stack's inputs; `tapenade/BBGeneralLib` (hypot).

## Unreleased (feature/ball-berry-canopy-state) - 2026-10-24

### Added

- PET method 6: the canopy state TV, TAH, EAH is carried between steps
  (`carry_canopy_state_between_steps = 1`), as Noah-MP carries it; exposed as the BMI inputs
  `vegetation_temperature_start_of_step_K`, `canopy_air_temperature_start_of_step_K`,
  `canopy_air_vapor_pressure_start_of_step_Pa` (hot start, or the CIESM Noah-MP reset TV).
- PET method 6: `canopy_energy_physics_option = 2` (nested solver): canopy biomass heat capacity,
  enthalpy-form fusion of intercepted water in the canopy energy balance, transpiration ramped to
  zero below TFRZ, radiation at the start-of-step TV, the EVC water limit fixed by LATHEAV, and the
  NWM outer-iteration rule. New outputs `canopy_freezing_rate_mm_per_s` (to CIESM canopy_storage),
  `canopy_fusion_heat_flux_W_per_m2`, `canopy_heat_storage_flux_W_per_m2`.
  `configs/ball_berry_revised_canopy_energy.cfg`; ctests `ball_berry_revised_canopy_energy` and
  `bmi_unit_ball_berry_revised_canopy_energy`.

### Changed

- Option 1 (default) is unchanged except that the state is now carried: the exact solver still
  matches NWM VEGE_FLUX (stored reference, r8, 2.9e-15) on the first step.

## Unreleased (feature/ciesm-coupling) - 2026-10-24

### Added

- PET method 6: BMI input `precipitation_advected_heat_to_canopy_W_per_m2` (Noah-MP PAHV, from
  CIESM canopy_interception), added to the canopy energy balance; it was fixed at 0. Configuration
  key of the same name (initial value) in `configs/method_6_ball_berry.cfg`; the Fortran reference
  driver reads it too. With PAHV = 15 W m-2 the exact solver matches NWM VEGE_FLUX (r8) to 2.7e-15.
- PET method 6: BMI output `canopy_latent_heat_J_per_kg` (LATHEAV), so CIESM canopy_storage
  converts EVC to water with the latent heat Ball-Berry used.

### Fixed

- `noahmp_comparison/Fortran/make_harness_module.py` matches the `private ::` declaration lines
  with a regular expression; copies of the NWM file differ in trailing blanks on those lines.

## x.y.z - YYYY-MM-DD

### Added

- Lorem ipsum dolor sit amet

### Deprecated

- Nothing.

### Removed

- Nothing.

### Fixed

- Nothing.
