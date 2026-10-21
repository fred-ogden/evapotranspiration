#!/usr/bin/env python3
# ------------------------------------------------------------------------------------------
# make_harness_module.py
#
# Purpose:
#   Create module_sf_noahmplsm_harness.F from the unmodified NWM Noah-MP source
#   module_sf_noahmplsm.F, with the smallest set of changes needed to drive the
#   radiation, vegetated canopy energy balance, and Ball-Berry stomata routines
#   from a standalone reference driver.  Every change is bracketed by
#   "HARNESS PATCH BEGIN" / "HARNESS PATCH END" comments in the output file.
#
# Changes made (and nothing else):
#   1. RADIATION, VEGE_FLUX, STOMATA, ESAT changed from private to public so the
#      driver can call them directly.
#   2. Module variables added after OPT_IMPERV:
#        HARNESS_FIX_GROUND_ALBEDO, HARNESS_SOIL_ALBEDO_DIRECT(2), HARNESS_SOIL_ALBEDO_DIFFUSE(2),
#        HARNESS_SNOW_ALBEDO_DIRECT(2), HARNESS_SNOW_ALBEDO_DIFFUSE(2)
#                                          -> impose snow-free soil and snow albedos
#        HARNESS_FIX_GROUND_TEMPERATURE    -> hold TG fixed (skip loop2 update)
#        HARNESS_TRACE_ENABLED, HARNESS_TRACE_UNIT -> per-iteration trace output
#   3. GROUNDALB: when HARNESS_FIX_GROUND_ALBEDO is true, ALBGRD and ALBGRI are
#      overwritten after the original computation by the imposed soil and snow albedos,
#      blended by FSNO exactly as GROUNDALB blends them:
#      ALBGRD = ALBSOD*(1-FSNO) + ALBSND*FSNO, ALBGRI = ALBSOI*(1-FSNO) + ALBSNI*FSNO.
#   4. VEGE_FLUX loop1: when HARNESS_TRACE_ENABLED, one CSV line per canopy
#      iteration is written to HARNESS_TRACE_UNIT (after TV is updated).
#   5. VEGE_FLUX loop2: when HARNESS_FIX_GROUND_TEMPERATURE is true, the ground
#      fluxes IRG, SHG, EVG are evaluated once at the prescribed TG (no Newton
#      update of TG) and GH is the residual SAG + PAHG - IRG - SHG - EVG, i.e.
#      the ground heat flux implied by the prescribed skin temperature.
#      When false, the original loop2 runs unchanged.
#
# Usage:
#   python3 make_harness_module.py module_sf_noahmplsm.F module_sf_noahmplsm_harness.F
# ------------------------------------------------------------------------------------------
import sys

if len(sys.argv) != 3:
    sys.exit("usage: make_harness_module.py input_noahmp_module.F output_harness_module.F")

input_path = sys.argv[1]
output_path = sys.argv[2]

with open(input_path, "r", encoding="ascii", errors="strict") as input_file:
    source_text = input_file.read()


def replace_exactly_once(text, old_text, new_text, description):
    """Replace old_text with new_text, requiring exactly one occurrence."""
    occurrence_count = text.count(old_text)
    if occurrence_count != 1:
        sys.exit("patch failed (%s): expected 1 occurrence, found %d" % (description, occurrence_count))
    return text.replace(old_text, new_text)


# ---- change 1: make the routines the driver calls public ----------------------------------
for routine_line in ["  private ::       RADIATION\n",
                     "  private ::       VEGE_FLUX\n",
                     "  private ::               STOMATA                  \n",
                     "  private ::               ESAT\n"]:
    public_line = routine_line.replace("private", "public ")
    source_text = replace_exactly_once(source_text, routine_line,
                                       "! HARNESS PATCH: was private\n" + public_line,
                                       "public " + routine_line.strip())

# ---- change 2: harness control variables ----------------------------------------------------
old_text = "  INTEGER :: OPT_IMPERV  ! options for imperviousness infiltration adjustment\n"
new_text = old_text + (
    "! HARNESS PATCH BEGIN: control variables for the standalone reference driver\n"
    "  LOGICAL :: HARNESS_FIX_GROUND_ALBEDO       = .FALSE. ! true: impose soil and snow albedos below\n"
    "  REAL    :: HARNESS_SOIL_ALBEDO_DIRECT(2)    = 0.     ! imposed snow-free soil albedo, direct (1=vis,2=nir)\n"
    "  REAL    :: HARNESS_SOIL_ALBEDO_DIFFUSE(2)   = 0.     ! imposed snow-free soil albedo, diffuse (1=vis,2=nir)\n"
    "  REAL    :: HARNESS_SNOW_ALBEDO_DIRECT(2)    = 0.     ! imposed snow albedo, direct (1=vis,2=nir)\n"
    "  REAL    :: HARNESS_SNOW_ALBEDO_DIFFUSE(2)   = 0.     ! imposed snow albedo, diffuse (1=vis,2=nir)\n"
    "  LOGICAL :: HARNESS_FIX_GROUND_TEMPERATURE  = .FALSE. ! true: TG prescribed, loop2 does not update TG\n"
    "  LOGICAL :: HARNESS_TRACE_ENABLED           = .FALSE. ! true: write per-iteration canopy trace\n"
    "  INTEGER :: HARNESS_TRACE_UNIT              = 0       ! Fortran unit for the canopy trace\n"
    "! HARNESS PATCH END\n")
source_text = replace_exactly_once(source_text, old_text, new_text, "control variables")

# ---- change 3: ground albedo override at end of GROUNDALB -----------------------------------
old_text = ("        ALBGRD(IB) = ALBSOD*(1.-FSNO) + ALBSND(IB)*FSNO\n"
            "        ALBGRI(IB) = ALBSOI*(1.-FSNO) + ALBSNI(IB)*FSNO\n"
            "  END DO\n")
new_text = old_text + (
    "! HARNESS PATCH BEGIN: impose soil and snow albedos supplied by the driver, blended by FSNO\n"
    "  IF (HARNESS_FIX_GROUND_ALBEDO) THEN\n"
    "     ALBGRD(1:2) = HARNESS_SOIL_ALBEDO_DIRECT(1:2)*(1.-FSNO)  + HARNESS_SNOW_ALBEDO_DIRECT(1:2)*FSNO\n"
    "     ALBGRI(1:2) = HARNESS_SOIL_ALBEDO_DIFFUSE(1:2)*(1.-FSNO) + HARNESS_SNOW_ALBEDO_DIFFUSE(1:2)*FSNO\n"
    "  END IF\n"
    "! HARNESS PATCH END\n")
source_text = replace_exactly_once(source_text, old_text, new_text, "GROUNDALB override")

# ---- change 4: per-iteration trace in VEGE_FLUX loop1 ---------------------------------------
old_text = ("! update vegetation surface temperature\n"
            "        TV  = TV + DTV\n")
new_text = old_text + (
    "! HARNESS PATCH BEGIN: per-iteration canopy trace (values after the TV update)\n"
    "        IF (HARNESS_TRACE_ENABLED) THEN\n"
    "           WRITE(HARNESS_TRACE_UNIT,'(I4,30(\",\",ES16.8E3))') ITER, TV, DTV, TAH, EAH, &\n"
    "                RSSUN, RSSHA, PSNSUN, PSNSHA, RB, RAHC, RAWC, RAHG, RAWG, CM, CH, FV, MOZ, &\n"
    "                FM, FH, FHG, IRC, SHC, EVC, TR, CVH, CTW, CAH, CGH, ESTV, DESTV\n"
    "        END IF\n"
    "! HARNESS PATCH END\n")
source_text = replace_exactly_once(source_text, old_text, new_text, "loop1 trace")

# ---- change 5: optional prescribed TG in VEGE_FLUX loop2 ------------------------------------
old_text = ("     loop2: DO ITER = 1, NITERG\n"
            "\n"
            "        T = TDC(TG)\n")
new_text = (
    "! HARNESS PATCH BEGIN: prescribed ground skin temperature, no Newton update of TG\n"
    "     IF (HARNESS_FIX_GROUND_TEMPERATURE) THEN\n"
    "        T = TDC(TG)\n"
    "        CALL ESAT(T, ESATW, ESATI, DSATW, DSATI)\n"
    "        IF (T .GT. 0.) THEN\n"
    "            ESTG  = ESATW\n"
    "            DESTG = DSATW\n"
    "        ELSE\n"
    "            ESTG  = ESATI\n"
    "            DESTG = DSATI\n"
    "        END IF\n"
    "        IRG = CIR*TG**4 + AIR\n"
    "        SHG = CSH * (TG         - TAH         )\n"
    "        EVG = CEV * (ESTG*RHSUR - EAH         )\n"
    "        GH  = SAG + PAHG - IRG - SHG - EVG   ! residual: ground heat flux implied by TG\n"
    "     ELSE\n"
    "! HARNESS PATCH END (original loop2 follows, unchanged)\n"
    + old_text)
source_text = replace_exactly_once(source_text, old_text, new_text, "loop2 prescribed TG, part 1")

old_text = ("        TG  = TG  + DTG\n"
            "\n"
            "     END DO loop2\n")
new_text = old_text + (
    "! HARNESS PATCH BEGIN: close the prescribed-TG branch\n"
    "     END IF\n"
    "! HARNESS PATCH END\n")
source_text = replace_exactly_once(source_text, old_text, new_text, "loop2 prescribed TG, part 2")

with open(output_path, "w", encoding="ascii") as output_file:
    output_file.write(source_text)

print("wrote %s" % output_path)
