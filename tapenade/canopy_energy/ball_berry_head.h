/*
 * tapenade/canopy_energy/ball_berry_head.h - the head Tapenade differentiates: one PET method 6
 * (Ball-Berry canopy energy balance) step with the VEGE_FLUX loop1 replica.
 *
 * Author: Fred L. Ogden, Ph.D., P.E.
 *         NOAA/National Weather Service
 *
 * This software was developed by a U.S. Government employee as part of
 * official duties and is not subject to U.S. copyright protection under
 * 17 U.S.C. Section 105.
 *
 * Licensed under the Apache License, Version 2.0.
 * SPDX-License-Identifier: Apache-2.0
 */

/*
 * ball_berry_head() runs ball_berry_method6_physics_step() (src/ball_berry_canopy_physics.c) with every
 * active quantity in one array of independents and one of dependents, so the generated signatures are
 * predictable.  The canopy state at the start of the step (TV, TAH, EAH) is among the independents and
 * the state at the end among the dependents, so the head chains over steps.  Everything else comes in
 * passive: the parameters and the step inputs the preparation resolved (bb6_prepare_step_inputs()); the
 * independents are written over them.
 *
 * FLO 10/28: first version.  ASCII only.
 */
#ifndef BALL_BERRY_HEAD_H
#define BALL_BERRY_HEAD_H

#include "ball_berry_canopy_physics.h"

/* ---- independents */
#define BB_HEAD_IN_VEGETATION_TEMPERATURE_START      0   /* [TV]  K, start of step */
#define BB_HEAD_IN_CANOPY_AIR_TEMPERATURE_START      1   /* [TAH] K */
#define BB_HEAD_IN_CANOPY_AIR_VAPOR_PRESSURE_START   2   /* [EAH] Pa */
#define BB_HEAD_IN_AIR_TEMPERATURE                   3   /* [SFCTMP] K */
#define BB_HEAD_IN_SURFACE_PRESSURE                  4   /* [SFCPRS] Pa */
#define BB_HEAD_IN_SPECIFIC_HUMIDITY                 5   /* [Q2] kg/kg */
#define BB_HEAD_IN_SHORTWAVE_DOWN                    6   /* [SOLDN] W/m2 */
#define BB_HEAD_IN_LONGWAVE_DOWN                     7   /* [LWDN] W/m2 */
#define BB_HEAD_IN_WIND_X                            8   /* [UU] m/s */
#define BB_HEAD_IN_WIND_Y                            9   /* [VV] m/s */
#define BB_HEAD_IN_ROOT_STRESS_FACTOR               10   /* [BTRAN] 0 to 1 */
#define BB_HEAD_IN_GROUND_TEMPERATURE               11   /* [TG] K */
#define BB_HEAD_IN_SURFACE_RESISTANCE               12   /* [RSURF] s/m */
#define BB_HEAD_IN_SNOW_COVER_FRACTION              13   /* [FSNO] */
#define BB_HEAD_IN_SNOW_DEPTH                       14   /* [SNOWH] m */
#define BB_HEAD_IN_SNOW_ALBEDO_DIRECT_VISIBLE       15   /* [ALBSND(1)] */
#define BB_HEAD_IN_SNOW_ALBEDO_DIRECT_NIR           16   /* [ALBSND(2)] */
#define BB_HEAD_IN_SNOW_ALBEDO_DIFFUSE_VISIBLE      17   /* [ALBSNI(1)] */
#define BB_HEAD_IN_SNOW_ALBEDO_DIFFUSE_NIR          18   /* [ALBSNI(2)] */
#define BB_HEAD_IN_CANOPY_WETTED_FRACTION           19   /* [FWET] */
#define BB_HEAD_IN_CANOPY_LIQUID_WATER              20   /* [CANLIQ] mm */
#define BB_HEAD_IN_CANOPY_ICE                       21   /* [CANICE] mm */
#define BB_HEAD_IN_PRECIPITATION_ADVECTED_HEAT      22   /* [PAHV] W/m2 */
#define BB_HEAD_IN_SOIL_TOP_DISC_MOISTURE           23   /* [SMC(1)] m3/m3 */
#define BB_HEAD_IN_STOMATAL_SLOPE                   24   /* [MP] */
#define BB_HEAD_IN_MAXIMUM_CARBOXYLATION_RATE       25   /* [VCMX25] umol/m2/s */
#define BB_HEAD_IN_MINIMUM_LEAF_CONDUCTANCE         26   /* [BP] umol/m2/s */
#define BB_HEAD_IN_QUANTUM_EFFICIENCY               27   /* [QE25] */
#define BB_HEAD_IN_CANOPY_WIND_PARAMETER            28   /* [CWP] */
#define BB_HEAD_IN_GREEN_VEGETATION_FRACTION        29   /* [SHDFAC] */
#define BB_HEAD_IN_MOMENTUM_ROUGHNESS_LENGTH        30   /* [Z0MVT] m */
#define BB_HEAD_IN_LEAF_DIMENSION                   31   /* [DLEAF] m */
#define BB_HEAD_INDEPENDENT_COUNT                   32

/* ---- dependents: the coupling outputs (the freezing rate, NaN with the Noah-MP physics, and COSZ, set
 *      by the time inputs alone, left out) and the PET */
#define BB_HEAD_OUT_TRANSPIRATION_FLUX               0   /* [ETRAN] m/s */
#define BB_HEAD_OUT_CANOPY_EVAPORATION_HEAT_FLUX     1   /* [EVC] W/m2 */
#define BB_HEAD_OUT_GROUND_SHORTWAVE_ABSORBED        2   /* [SAG] W/m2 */
#define BB_HEAD_OUT_GROUND_LONGWAVE_ABSORBED         3   /* W/m2 */
#define BB_HEAD_OUT_CANOPY_AIR_TEMPERATURE           4   /* [TAH] K, end of step */
#define BB_HEAD_OUT_CANOPY_AIR_VAPOR_PRESSURE        5   /* [EAH] Pa */
#define BB_HEAD_OUT_BELOW_CANOPY_RESISTANCE          6   /* [RAHG] s/m */
#define BB_HEAD_OUT_VEGETATION_TEMPERATURE           7   /* [TV] K, end of step */
#define BB_HEAD_OUT_CANOPY_LATENT_HEAT               8   /* [LATHEAV] J/kg */
#define BB_HEAD_OUT_CANOPY_FUSION_HEAT_FLUX          9   /* W/m2, option 2 only */
#define BB_HEAD_OUT_CANOPY_HEAT_STORAGE_FLUX        10   /* W/m2, option 2 only */
#define BB_HEAD_OUT_SURFACE_RESISTANCE_WITH_SNOW    11   /* [RSURF] s/m */
#define BB_HEAD_OUT_VEGETATED_FRACTION              12   /* [FVEG] */
#define BB_HEAD_OUT_PET                             13   /* m/s */
#define BB_HEAD_DEPENDENT_COUNT                     14

void ball_berry_head(double independent[BB_HEAD_INDEPENDENT_COUNT], double dependent[BB_HEAD_DEPENDENT_COUNT],
                     const BallBerryMethod6Parameters *passive_parameters, const BallBerryStepInputs *passive_step);

#endif /* BALL_BERRY_HEAD_H */
