# PET BMI stabilization: review notes and proposed PR

## Scope

This patch is a targeted reliability pass over the supplied `src/` and `include/` snapshot. It is **not** a completed BMI conformance or PET numerical validation effort. No Ball–Berry/Noah-OWP-Modular code is included.

## Changes

- Correct file-forcing specific humidity mapping (previously read precipitation).
- Correct BMI input metadata from relative saturation to specific humidity (`kg kg-1`). **Breaking interface-name change:** coupled clients must update their variable mapping.
- Zero-initialize allocated model state and free the model, owned forcing-file string, and forcing arrays at finalization. A finalized `Bmi` cannot be reused without registering again.
- Propagate PET calculation failure through BMI `Update`; return errors instead of terminating the process for two input-parameter checks.
- Reject invalid PET method IDs; preserve configured wind measurement height between calls.
- Report scalar grid rank 0, not 1; restrict indexed access to scalar index 0; disallow setting output variables.
- Use elapsed model time for `Get_current_time`, with simulation start 0 and configured finite end time. This is **simulation-relative time**, not Unix epoch time.
- Check forcing record count, forcing array allocations, and some malformed configuration lines.
- Reject backward `Update_until` calls and propagate update failures. Fractional-step behavior remains provisional and needs explicit tests.
- Enable GNU/POSIX extensions in source for existing `strdup`/`strsep` usage; this is **not strict ISO C99 portability**.

## Known limitations / follow-up work before merging

1. **Run full executable tests:** only source syntax checks were performed on this snapshot. No Makefile, test fixtures, or representative config/forcing files were supplied.
2. Configuration parsing still uses unchecked numeric `strtod` conversions and requires comprehensive validation of parameter ranges and missing required keys.
3. `Update_until` fractional steps should be redesigned/tested with forcing-file and BMI-driven modes; file-mode fluxes are instantaneous and forcing records are indexed by update count.
4. Additional `exit()` calls remain in utility/solar routines in `pet.c`. They should be converted to error returns in a separate focused change.
5. Some physics inputs remain hard-coded in `pet_setup`; they are not changed to avoid altering reference calculations without tests.
6. Time convention and the `num_timesteps == 1` behavior have changed; review coupling assumptions.
7. Reinitialization without finalization, parser error cleanup, and all allocation failure paths need tests.
8. Confirm BMI standard names with the target BMI vocabulary and update coupling configuration mappings.
9. Validate all five PET methods against reference forcing and establish before/after numerical regression tests.

## Suggested PR title

`Stabilize PET BMI lifecycle, forcing metadata, and error handling`

## Suggested validation commands

```sh
gcc -std=c99 -Wall -Wextra -fsyntax-only -Iinclude src/bmi_pet.c src/pet.c
```

Review the source diff and run the repository's normal build/test suite before creating a pull request. This snapshot contains no Git history, so this archive is not itself a PR branch.

## Follow-up correction after five-method unit-test report (2026-10-19)

The initial review patch regressed aerodynamic and combination PET and caused the
Penman–Monteith method to terminate because the wind-speed forcing was normalized
to 2 m while aerodynamic formulas still used the original configured measurement
height (typically 10 m). This follow-up changes aerodynamic/combination logarithms
and Penman–Monteith aerodynamic resistance to use 2 m, consistent with the
normalized `pet_forcing.wind_speed_m_per_s`, while preserving the original
measurement height for subsequent forcing records. It also moves the null check
in `run_pet` before dereferencing the model pointer.

The source compiles with `gcc -std=c99 -Werror -Og -Iinclude` for the read-forcing
driver. Numerical unit tests have NOT been rerun because the supplied archive does
not include `configs/` or forcing data. Re-run all five original unit tests and
compare against the unmodified baseline before submitting the PR. This is a
correction to the patch, not yet a validated numerical fix.
