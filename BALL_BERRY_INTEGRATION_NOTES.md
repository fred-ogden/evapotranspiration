# Ball-Berry method 6 integration (review build)

This archive updates the source supplied by the user; it does not modify Git history.

## Changes
- The method 6 driver prepares the canopy state from the current AORC forcing record before exactly one BMI update, checks errors, and prints potential transpiration in m/s. The historical two-update behavior of methods 1–5 is retained.
- The method 6 `bb_stability_input_ready` flag is consumed on an attempted solve, so stale forcing cannot silently be reused.
- Test canopy geometry, optical, physiological, solar, and numerical settings are configurable through `bb_*` keys. The six pre-existing BMI inputs (`LAI`, `SAI`, `SHDFAC`, `CWP`, `MP`, `VCMX25`) remain independent.
- Conflicting hardcoded canopy geometry assignments are removed from the test adapter. Existing Ball-Berry physics helpers are not modified.
- Added a synthetic AORC test record and a CTest end-to-end test.

## Tests
On the isolated source snapshot: CMake build successful, 3/3 CTests passed; standalone GCC driver compiles. Test output with synthetic forcing: potential transpiration approximately 3.035708973599e-08 m/s. This is a software smoke test, not physical validation.

## Important limitations
- `configs/pet_config_unit_test6.txt` uses `test/ball_berry_forcing.csv`, a **synthetic** AORC record. It is not observed forcing.
- Test solar zenith, direct/diffuse partition, ground temperature offset, and optical/physiological assumptions remain fixed **per configuration**, not derived dynamically from date/time or vegetation datasets.
- The test sets a 10 m canopy and 10 m wind reference (matching the AORC wind observation height). This is an unvalidated canopy-top reference-height approximation; no 10-to-20 m wind extrapolation is performed. Review the aerodynamic treatment before production.
- `bb_*` test parameters are distinct from the production BMI canopy state; the test adapter is only linked into the standalone driver and end-to-end test.
- Existing methods 1–5 were preserved in source but not exercised against their forcing data, which were not included in the uploaded archive.
- The original project file `petbmi.pc.in` was not in the uploaded archive. The isolated build used a temporary local stand-in; the returned archive does not replace your existing `petbmi.pc.in`.
- The configuration parser retains legacy behavior of ignoring unknown keys. Further validation and parser unit tests are recommended before production.

## Suggested installation
Inspect the changed files before copying them into the working repository. The archive includes the original `src`, `include`, `configs`, `test`, `CMakeLists.txt`, and forcing-driver script, plus these notes. Do not overwrite unrelated files or commit until your local tests pass.
