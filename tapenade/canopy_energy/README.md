# Tapenade differentiation of one PET method 6 (Ball-Berry) step

One step of the Ball-Berry canopy energy balance with the VEGE_FLUX loop1 replica, differentiated by Tapenade 3.16
in tangent and reverse mode, in the configuration of the dcfe-noahmp-stack for the 01123000 basin (cat-18560):
- the exact solver (`bb_noahmp_exact_iteration = 1`), with `canopy_energy_physics_option = 1` (Noah-MP);
- the canopy state carried between steps;
- solar geometry from the time inputs, and the monthly LAI/SAI table;
- the coupling inputs from the other modules: BTRAN from CFE3.1, TG from sft-C (fidelity 1), RSURF, the snow
  fraction, depth and albedos from NoahMP_snow_ablation, and the canopy water from CIESM.

This is the fifth and last module core of the dcfe-noahmp-stack (its DESIGN.md 7.2, ladder step 5), after sft-C
(fidelity 1), CIESM, NoahMP_snow_ablation and CFE3.1.

Run, from anywhere (on a machine with Tapenade):

    tapenade/canopy_energy/run_tapenade_ball_berry.csh

Tapenade location: `$TAPENADE_HOME`, or `/user2/ogden/opt/tapenade/tapenade_3.16`. Without Tapenade the script
builds the primal-only test (experiment p, and t as a central-difference self-consistency check).

## Files

- `ball_berry_head.h`, `.c`: the head, `ball_berry_head()`, over `ball_berry_method6_physics_step()`. Its
  independents and dependents are flat arrays, so the generated signature is predictable. The canopy state at the
  start of the step is among the independents and the state at the end among the dependents, so the head chains
  over steps.
  - 33 independents:
    - the canopy state at the start of the step: TV, TAH, EAH;
    - the forcing: air temperature, pressure, specific humidity, shortwave and longwave down, the two wind
      components;
    - the coupling inputs: BTRAN, TG, RSURF, FSNO, SNOWH, the four snow albedos, FWET, CANLIQ, CANICE, PAHV, SMC(1);
    - the parameters: MP, VCMX25, BP, QE25, CWP, SHDFAC, Z0M, DLEAF, and RSURF_SNOW (FLO 10/28).
  - 14 dependents: the coupling outputs (ETRAN, EVC, SAG, the ground longwave absorbed, TAH, EAH, RAHG, TV,
    LATHEAV, the option 2 fusion and storage heat, RSURF after the snow blend, FVEG) and the PET. Left out: the
    freezing rate (NaN in the BMI with the Noah-MP physics) and COSZ (set by the time inputs alone).
  - The rest comes in passive: the parameters and the step inputs the preparation resolved
    (`bb6_prepare_step_inputs()`). They are copied into local structures field by field.
- `../BBGeneralLib`: `hypot` with its derivative, for Tapenade (`-ext`).
- Two adapters, `ball_berry_{tangent,adjoint}_adapter.c`, with `ball_berry_adapters.h`. They are the only files
  that see generated code. The run script reads each generated prototype and tells the adapter which passive
  arguments got derivative slots, and of which type; zeros are passed in them.
- `ball_berry_ad_test.c`: the tests (below).
- `data/`:
  - the divide's Ball-Berry configuration for cat-18560;
  - the inputs dcfe-noahmp-stack gave Ball-Berry over the whole realization period, 2012-10-01 to 2019-07-12
    (`build/stack_driver --core-path --dump-slot-inputs`, the canopy_energy slot, 59,424 hours, %.17g).
- `run_tapenade_ball_berry.csh`: stage, differentiate, repair, compile, run.
  - Generated code goes to `tapenade_work/canopy_energy/`, and results to `output/tapenade_canopy_energy/`. Git
    ignores both.
  - The generated code is also saved to `output/tapenade_canopy_energy/generated_code.tar.gz`.

## Tests

Ball-Berry is initialized through its BMI from the configuration, as the stack does. Each hour, the dumped inputs
are set into the model as the stack's adapter sets them. The preparation is run on a copy of the model, the head's
independents are packed from the step inputs it gives, and the model is stepped with
`ball_berry_method6_core_step()`, as the stack's adapter steps it.

- `p`: the head against the production step at every hour, bitwise, every dependent.
- `t`: the tangent against central differences at the test states, each independent with a step-size scan from
  1e-3 to 1e-7 of the larger of its value and a typical scale.
  - The test states are the first hour of each regime and every 5000th hour (20 states): night, day, stressed
    (BTRAN < 1), wet canopy, wet canopy with under 0.05 mm of water, an SFCDIF1 neutral reset (the stability
    parameter reset to 0 after two sign changes), snow on the ground, a frozen canopy (LATHEAV = HSUB), and loop1
    at 20 iterations (not converged to 0.01 K). No hour of the record has an unexposed canopy or the dry bare soil
    resistance.
  - Each step is compared over the dependents it resolves (relative change above 1e-9).
  - Directions are classified as in CFE3.1's test: at a bound of the domain (zero snow, wetted fraction or canopy
    water, BTRAN = 1, a fraction or albedo at 1; one-sided, not tested), non-finite tangent (a failure), below
    finite-difference resolution (tested only for agreement), kinks (reported apart, with whether the tangent
    equals one side's difference), smooth (pass below 1e-4).
  - The tangent code's primal must equal the head bitwise.
- `a`: the reverse Jacobian against the tangent Jacobian, and the dot-product test (pass below 1e-10).

## Tapenade-ready form (byte-identical)

One step is a preparation and the physics:
- `bb6_prepare_step_inputs()` (`src/ball_berry_method6_inputs.c`, not differentiated): the checks of the coupling
  inputs and the forcing, the solar geometry, the leaf and stem area index, the NaN conventions of the BMI inputs
  resolved into flags (a start-of-step canopy state, or TG, not yet supplied), and the messages.
- `ball_berry_method6_physics_step()` (`src/ball_berry_canopy_physics.c`, new): the canopy input
  (`bb6_canopy_input_physics()`: ATM, PHENOLOGY, GROUNDALB, the snow burial, the RSURF snow blend, the ENERGY
  geometry, the start-of-step state), the no-canopy branch or the VEGE_FLUX loop1 replica, and the coupling outputs.
  The functions of the replica (ALBEDO/TWOSTREAM/SURRAD, ESAT, SFCDIF1, RAGRB, the foliage step and the exact
  solver) moved there from `src/ball_berry_helpers.c`, which keeps the nested solver (not differentiated; it uses
  the same canopy input).

The code Tapenade parses (`src/ball_berry_canopy_physics.c`, `src/noahmp_stomata.c`,
`src/ball_berry_struct_copies.c`) follows the rules of the other module cores:
- one return per function (a status variable), no break (the loop1 exit is a flag), no ternaries;
- no printf (the exact solver records its iteration trace; the caller writes the file), no `isnan` or `isfinite`
  (`bb_is_finite()`), no NAN (the BMI's NaN freezing rate is a flag the caller turns into NaN);
- no whole-structure assignment or memset: generated field-by-field copies and zeroing
  (`tools/generate_ball_berry_struct_copies.py` writes `src/ball_berry_struct_copies.c`);
- nested fmin and fmax written as two statements;
- every chain of two or more multiplications and divisions (every division by a product, and every product with a
  quotient in parentheses) written one operation per statement, in C's left-to-right order.

Two changes go beyond the form:
- EVC = MIN(latent limit, EVC) is a comparison, `if (latent_limit < evc) evc = latent_limit`. When both are zero
  with opposite signs, `fmin()` may return either, and gcc 13 swapped its arguments when the function changed; the
  comparison keeps EVC in a tie, as the code before did.
- TWOSTREAM's `copysign(1.0e-6, sigma)` is a comparison: the same except at sigma = -0.0.

The 13 ctests and their output files, a replay of the stack's inputs for the six 01123000 divides over the whole
period (356,544 divide-hours, every coupling output, %a), the CIESM coupled ctests and the dcfe-noahmp-stack core
path test (56 files) are bitwise the same as before the change.

## Results, Tapenade 3.16, FLO 10/28, all PASS (here, and on nordin with the same numbers)

- p: 59,424 hours, bitwise.
- t: worst tangent error away from kinks 7.4e-7 (BP on the first day), over 640 point-independent pairs; no kinks;
  the tangent code's primal equals the head bitwise; every tangent is finite.
  - 166 directions with no resolvable effect, for example: TAH at the start of the step (loop1 does not use it:
    RAGRB's stability on iteration 1 is neutral, and the iterations use the updated TAH); SMC(1) above 0.275,
    where GROUNDALB's increment is 0; the snow albedos and depth without snow (the canopy bottom is at 11.5 m);
    MP, VCMX25 and QE25 at night, where STOMATA's resistance is the cuticular one, cf/BP.
  - 140 directions at a bound of the domain, mostly BTRAN = 1 and zero snow or canopy water.
- a: reverse against tangent Jacobian 2.8e-12 at worst; dot products 4.4e-14 at worst.

The first run failed t only on the tangent code's primal, which differed from the head in the last bits: Tapenade
regrouped the chains of multiplications and divisions (for example STOMATA's `0.5*kc/ko*o2*0.21`). Writing the
chains one operation per statement fixed it.

## RSURF_SNOW added, FLO 10/28, all PASS (Tapenade 3.16)

RSURF_SNOW (`surface_resistance_snow_s_per_m`, the snow resistance of the OPT_RSF = 4 blend of RSURF) is among the
parameters estimated with the shared NOM-named set, with VCMX25, MP and CWP (Fred, 10/28). It is an independent of
the head and a BMI input (`RSURF_SNOW`, s/m, positive), so the stack can set it as a model_params value. The BMI now
has 37 inputs; the ctests' outputs differ only in that count and in the BMI unit test's check count (983, was 972).
- p: 59,424 hours, bitwise.
- t: worst tangent error away from kinks 7.4e-7, over 660 point-independent pairs; no kinks. RSURF_SNOW has an
  effect only with snow on the ground (5 of the 20 states; 1e-10 to 1e-8 there).
- a: reverse against tangent Jacobian 2.8e-12 at worst; dot products 5.2e-14 at worst.

FLO 10/28: first version.
