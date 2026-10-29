/*
 * tapenade/canopy_energy/ball_berry_adjoint_adapter.c - the only file that sees Tapenade's generated REVERSE of ball_berry_head().
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
 * and a signature mismatch is a compile error.  Generated signature (suffix b): each active array
 * followed by its derivative array; a passive argument copied into a structure with active fields can
 * get a derivative slot too (sft-C, CIESM, NoahMP_snow_ablation, CFE3.1):
 *   ball_berry_head_b(independent, independentb, dependent, dependentb,
 *       passive_parameters, [passive_parametersb,] passive_step, [passive_stepb])
 * The run script reads the generated prototype and defines BB_HEAD_B_<ARGUMENT>_SLOT_BY_POINTER and
 * BB_HEAD_B_<ARGUMENT>_SLOT_TYPE (<name>_diff, or the original type for a structure without
 * differentiable fields) for the slots present; zeros are passed in them. Reverse mode needs ADFirstAidKit (adStack.c); independentb is incremented, so it is zeroed first.
 *
 * FLO 10/28: first version.  ASCII only.
 */
#include <string.h>
#include "ball_berry_head_b.c"
#include "ball_berry_adapters.h"

#ifndef BB_HEAD_B_PASSIVE_PARAMETERS_SLOT_TYPE
#define BB_HEAD_B_PASSIVE_PARAMETERS_SLOT_TYPE BallBerryMethod6Parameters_diff
#endif
#ifndef BB_HEAD_B_PASSIVE_STEP_SLOT_TYPE
#define BB_HEAD_B_PASSIVE_STEP_SLOT_TYPE BallBerryStepInputs_diff
#endif

void ball_berry_adjoint(const double independent[], const double dependent_adjoint[], double independent_adjoint[],
                        const BallBerryMethod6Parameters *passive_parameters, const BallBerryStepInputs *passive_step)
{
    double independent_copy[BB_HEAD_INDEPENDENT_COUNT], dependent[BB_HEAD_DEPENDENT_COUNT];
    double dependent_adjoint_copy[BB_HEAD_DEPENDENT_COUNT];
    int index;
#ifdef BB_HEAD_B_PASSIVE_PARAMETERS_SLOT_BY_POINTER
    static BB_HEAD_B_PASSIVE_PARAMETERS_SLOT_TYPE passive_parameters_adjoint;   /* passive: zero */
#endif
#ifdef BB_HEAD_B_PASSIVE_STEP_SLOT_BY_POINTER
    static BB_HEAD_B_PASSIVE_STEP_SLOT_TYPE passive_step_adjoint;               /* passive: zero */
#endif
#ifdef BB_HEAD_B_PASSIVE_PARAMETERS_SLOT_BY_POINTER
    memset(&passive_parameters_adjoint, 0, sizeof(passive_parameters_adjoint));
#endif
#ifdef BB_HEAD_B_PASSIVE_STEP_SLOT_BY_POINTER
    memset(&passive_step_adjoint, 0, sizeof(passive_step_adjoint));
#endif
    for (index = 0; index < BB_HEAD_INDEPENDENT_COUNT; ++index) {
        independent_copy[index] = independent[index];
        independent_adjoint[index] = 0.0;
    }
    for (index = 0; index < BB_HEAD_DEPENDENT_COUNT; ++index) {
        dependent[index] = 0.0;
        dependent_adjoint_copy[index] = dependent_adjoint[index];
    }
    ball_berry_head_b(independent_copy, independent_adjoint, dependent, dependent_adjoint_copy,
                      passive_parameters,
#if defined(BB_HEAD_B_PASSIVE_PARAMETERS_SLOT_BY_POINTER)
                      &passive_parameters_adjoint,
#endif
                      passive_step
#if defined(BB_HEAD_B_PASSIVE_STEP_SLOT_BY_POINTER)
                      , &passive_step_adjoint
#endif
                      );
}
