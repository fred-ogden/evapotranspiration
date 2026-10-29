/*
 * tapenade/canopy_energy/ball_berry_adapters.h - calls into the Tapenade-generated tangent and reverse of ball_berry_head().
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
 * The two *_adapter.c files are the only ones that see generated code.  The adapters copy their
 * arguments, so the generated code never writes into the caller's arrays.
 *
 * FLO 10/28: first version.  ASCII only.
 */
#ifndef BALL_BERRY_ADAPTERS_H
#define BALL_BERRY_ADAPTERS_H

#include "ball_berry_head.h"

/* tangent: dependent_tangent = J * independent_tangent; dependent: the tangent code's primal */
void ball_berry_tangent(const double independent[], const double independent_tangent[], double dependent[],
                        double dependent_tangent[], const BallBerryMethod6Parameters *passive_parameters,
                        const BallBerryStepInputs *passive_step);

/* reverse: independent_adjoint = J^T * dependent_adjoint (set, not accumulated) */
void ball_berry_adjoint(const double independent[], const double dependent_adjoint[], double independent_adjoint[],
                        const BallBerryMethod6Parameters *passive_parameters, const BallBerryStepInputs *passive_step);

#endif
