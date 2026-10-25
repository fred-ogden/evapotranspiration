All notable changes to this project will be documented in this file.
We follow the [Semantic Versioning 2.0.0](http://semver.org/) format.


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
