# PET BMI infrastructure fixes (FLO 10/26)

Modified: src/bmi_pet.c and src/main_pass_forcing.c.

- Reject fractional Update_until requests. The old implementation consumed a complete forcing index for a partial time increment. Full timestep updates are retained. Fractional stepping requires a separate forcing-interpolation design.
- Reject malformed non-finite numeric configuration values; preserve the legacy text-valued cloud_base_height_known and forcing_file entries.
- Validate positive timestep, remove duplicated albedo branch, fix floating-point printf format.
- Reject null BMI get/set arguments before strcmp/memcpy.
- Check coupling driver allocation, registration, initialization, update and transfer failures, with cleanup on exit.
- Preserve the AORC legacy humidity BMI name and PET specific-humidity BMI name.

Testing: compilation of both example drivers succeeded. Complete runtime regression was NOT possible because forcing/cat-67_2015.csv is absent from the uploaded repository. The existing read-forcing driver segfaults after missing-file initialization failures; this is outside the two-file patch. AORC finalize has an existing upstream memory-management defect (reallocates model data) that is not changed here.

Apply from repository root: tar -xzf pet_bmi_infrastructure_fixes_FLO_10_26.tar.gz
