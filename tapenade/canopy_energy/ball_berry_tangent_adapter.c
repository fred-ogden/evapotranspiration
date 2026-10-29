/*
 * tapenade/canopy_energy/ball_berry_tangent_adapter.c - the only file that sees Tapenade's generated TANGENT of ball_berry_head().
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
 * It #includes the generated source, so the call is compiled against the real generated definition
 * and a signature mismatch is a compile error.  Generated signature (suffix d): each active array
 * followed by its derivative array; a passive argument copied into a structure with active fields can
 * get a derivative slot too (sft-C, CIESM, NoahMP_snow_ablation, CFE3.1):
 *   ball_berry_head_d(independent, independentd, dependent, dependentd,
 *       passive_parameters, [passive_parametersd,] passive_step, [passive_stepd])
 * The run script reads the generated prototype and defines BB_HEAD_D_<ARGUMENT>_SLOT_BY_POINTER and
 * BB_HEAD_D_<ARGUMENT>_SLOT_TYPE (<name>_diff, or the original type for a structure without
 * differentiable fields) for the slots present; zeros are passed in them.
 *
 * FLO 10/28: first version.  ASCII only.
 */
#include <string.h>
#include "ball_berry_head_d.c"
#include "ball_berry_adapters.h"

#ifndef BB_HEAD_D_PASSIVE_PARAMETERS_SLOT_TYPE
#define BB_HEAD_D_PASSIVE_PARAMETERS_SLOT_TYPE BallBerryMethod6Parameters_diff
#endif
#ifndef BB_HEAD_D_PASSIVE_STEP_SLOT_TYPE
#define BB_HEAD_D_PASSIVE_STEP_SLOT_TYPE BallBerryStepInputs_diff
#endif

void ball_berry_tangent(const double independent[], const double independent_tangent[], double dependent[],
                        double dependent_tangent[], const BallBerryMethod6Parameters *passive_parameters,
                        const BallBerryStepInputs *passive_step)
{
    double independent_copy[BB_HEAD_INDEPENDENT_COUNT], independent_tangent_copy[BB_HEAD_INDEPENDENT_COUNT];
    int index;
#ifdef BB_HEAD_D_PASSIVE_PARAMETERS_SLOT_BY_POINTER
    static BB_HEAD_D_PASSIVE_PARAMETERS_SLOT_TYPE passive_parameters_tangent;   /* passive: zero */
#endif
#ifdef BB_HEAD_D_PASSIVE_STEP_SLOT_BY_POINTER
    static BB_HEAD_D_PASSIVE_STEP_SLOT_TYPE passive_step_tangent;               /* passive: zero */
#endif
#ifdef BB_HEAD_D_PASSIVE_PARAMETERS_SLOT_BY_POINTER
    memset(&passive_parameters_tangent, 0, sizeof(passive_parameters_tangent));
#endif
#ifdef BB_HEAD_D_PASSIVE_STEP_SLOT_BY_POINTER
    memset(&passive_step_tangent, 0, sizeof(passive_step_tangent));
#endif
    for (index = 0; index < BB_HEAD_INDEPENDENT_COUNT; ++index) {
        independent_copy[index] = independent[index];
        independent_tangent_copy[index] = independent_tangent[index];
    }
    for (index = 0; index < BB_HEAD_DEPENDENT_COUNT; ++index) {
        dependent[index] = 0.0;
        dependent_tangent[index] = 0.0;
    }
    ball_berry_head_d(independent_copy, independent_tangent_copy, dependent, dependent_tangent,
                      passive_parameters,
#if defined(BB_HEAD_D_PASSIVE_PARAMETERS_SLOT_BY_POINTER)
                      &passive_parameters_tangent,
#endif
                      passive_step
#if defined(BB_HEAD_D_PASSIVE_STEP_SLOT_BY_POINTER)
                      , &passive_step_tangent
#endif
                      );
}
