PET BMI final wind-height review (FLO 10/26)

1. Added missing FLO 10/26 explanations to the forcing-file SPFH correction,
   the preservation of configured wind measurement height, and the effective
   2 m wind reference in the aerodynamic, combination, and resistance formulas.
2. Verified run_pet normalizes incoming wind from the configured measurement
   height to 2 m without changing the configured height. The aerodynamic and
   combination methods and Penman-Monteith aerodynamic resistance use the
   normalized 2 m wind. Energy-balance and Priestley-Taylor do not use wind
   directly in their principal PET formula.
3. Extended test/bmi_regression.c to exercise all five methods with equivalent
   wind at 10 m and 2 m, repeated updates, and preservation of the configured
   height. Existing BMI invalid-argument and fractional-step tests remain.
4. This supplemental patch changes only comments in src/pet.c and three
   headers, plus the new regression coverage in test/bmi_regression.c.

Validation: ./test_bmi.sh PASS (compiler reports pre-existing warnings).
The external-forcing integration tests were not rerun in this environment.

Caveat: This test checks numerical consistency of the existing logarithmic
wind adjustment, not the meteorological validity of its displacement-height
parameterization for every vegetation or atmospheric stability regime.
