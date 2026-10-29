/*
 * ball_berry_canopy_physics.h - the physics of one PET method 6 (Ball-Berry canopy energy balance)
 * time step on resolved inputs: the part Tapenade differentiates.
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
 * FLO 10/28: one method 6 step is a preparation and the physics.
 *
 * The preparation (bb6_prepare_step_inputs(), ball_berry_method6_inputs.c, not differentiated) checks
 * the coupling inputs and the forcing, computes the solar geometry (COSZ from the time inputs and the
 * site), takes the leaf and stem area index from the monthly table or the BMI values, resolves the NaN
 * conventions of the BMI inputs into flags (a start-of-step canopy state or a ground temperature not yet
 * supplied), and prints the messages.  It fills a BallBerryStepInputs.
 *
 * The physics (ball_berry_canopy_physics.c, differentiated) builds the canopy input from the step inputs
 * and the parameters (bb6_canopy_input_physics(): ATM, PHENOLOGY, GROUNDALB, the snow burial, the surface
 * resistance with snow, the ENERGY geometry, the start-of-step canopy state), and, with the VEGE_FLUX
 * loop1 replica (bb_noahmp_exact_iteration = 1), solves the canopy energy balance and forms the coupling
 * outputs (ball_berry_method6_physics_step()).  The nested Picard/Newton solver
 * (noahmp_stability_coupled_solve(), ball_berry_helpers.c) is not differentiated; it uses the same
 * canopy input.
 *
 * The code Tapenade parses has one return per function, no printf, no isnan or isfinite (bb_is_finite()
 * instead), no ternaries, no break, no whole-structure assignment or memset (the generated field-by-field
 * copies of ball_berry_struct_copies.c instead), and no NAN: the freezing rate the BMI reports as NaN
 * with the Noah-MP physics is returned as 0 with a flag, and the preparation writes the NaN.
 *
 * Byte-identical to the code before 10/28.  ASCII only.
 */
#ifndef BALL_BERRY_CANOPY_PHYSICS_H
#define BALL_BERRY_CANOPY_PHYSICS_H

#include "noahmp_stomata.h"
#include "ball_berry_helpers.h"
#include "ball_berry_state.h"

#define NOAHMP_EXACT_MAXIMUM_CANOPY_ITERATIONS 20   /* NITERC of VEGE_FLUX loop1 */

/* One row of the exact solver's optional iteration trace (bb_iteration_trace_file): the values of the
 * CSV row the solver wrote before 10/28, recorded so the preparation can write the file. */
typedef struct {
    int iteration;
    double vegetation_temperature_k;
    double vegetation_temperature_change_k;
    double canopy_air_temperature_k;
    double canopy_air_vapor_pressure_pa;
    double sunlit_stomatal_resistance_s_per_m;
    double shaded_stomatal_resistance_s_per_m;
    double leaf_boundary_layer_resistance_s_per_m;
    double above_canopy_heat_resistance_s_per_m;
    double under_canopy_heat_resistance_s_per_m;
    double momentum_drag_coefficient;
    double heat_exchange_coefficient;
    double friction_velocity_m_per_s;
    double monin_obukhov_stability_parameter;
    double canopy_net_longwave_w_per_m2;
    double canopy_sensible_heat_w_per_m2;
    double transpiration_heat_w_per_m2;
} NoahmpVegeFluxExactTraceRecord;
typedef struct {
    int started;        /* 1: the solver reached its iterations (the header line is written) */
    int count;          /* rows recorded */
    NoahmpVegeFluxExactTraceRecord record[NOAHMP_EXACT_MAXIMUM_CANOPY_ITERATIONS];
} NoahmpVegeFluxExactTrace;

/* The inputs of one step, resolved by the preparation. */
typedef struct {
    double air_temperature_K;                    /* [SFCTMP] */
    double surface_pressure_Pa;                  /* [SFCPRS] */
    double specific_humidity_kg_per_kg;          /* [Q2] */
    double shortwave_down_W_per_m2;              /* [SOLDN] before the terrain correction */
    double longwave_down_W_per_m2;               /* [LWDN] */
    double wind_x_m_per_s;                       /* [UU] */
    double wind_y_m_per_s;                       /* [VV] */
    double cosine_solar_zenith_angle;            /* [COSZ] from the time inputs (solar_geometry_option = 1)
                                                    or the configured value */
    double cosine_solar_zenith_angle_horizontal; /* [COSZ_HORIZ] solar_geometry_option = 1 only */
    double leaf_area_index;                      /* [LAI] table or BMI value, before the PHENOLOGY checks */
    double stem_area_index;                      /* [SAI] */
    double green_vegetation_fraction;            /* [SHDFAC] BMI variable */
    double canopy_wind_parameter;                /* [CWP] BMI variable */
    NoahmpStomataParameters stomata;             /* STOMATA parameters, including the BMI variables MP and VCMX25 */
    BallBerryCouplingInputs coupling;            /* the BMI coupling inputs; NaN fields are flagged below */
    int ground_temperature_known;                /* 0: TG not yet supplied; TG = SFCTMP + offset */
    int vegetation_temperature_known;            /* 0: the cold start sets TV, TAH, EAH ... */
    int canopy_air_temperature_known;            /* ... unless the start-of-step value is known */
    int canopy_air_vapor_pressure_known;
    double time_step_s;                          /* [DT] */
} BallBerryStepInputs;

/* The results of one step. */
typedef struct {
    BallBerryCouplingOutputs coupling_outputs;   /* every field; see freezing_rate_defined */
    int freezing_rate_defined;                   /* 0: canopy_freezing_rate_mm_per_s is NaN in the BMI (the
                                                    Noah-MP physics: CIESM uses its own rate); the physics
                                                    returns 0 */
    int no_canopy_this_step;                     /* 1: ELAI + ESAI = 0, VEGE_FLUX skipped */
    double no_canopy_ground_shortwave_W_per_m2;  /* the no-canopy values (set only when no_canopy_this_step) */
    double no_canopy_ground_emissivity;
    double no_canopy_air_temperature_K;
    double no_canopy_air_vapor_pressure_Pa;
    double no_canopy_vegetation_temperature_K;
    double transpiration_m_per_s;                /* [ETRAN] */
    double pet_m_per_s;                          /* ETRAN, at least 0 */
} BallBerryStepResults;

/* The canopy input of the step (bb_stability_input) from the step inputs and the parameters; sets the
 * results' no-canopy flag and values, and the coupling outputs COSZ, FVEG and RSURF with snow.  The nested
 * solver controls are set by the caller.  Returns 0. */
int bb6_canopy_input_physics(const BallBerryMethod6Parameters *parameters, const BallBerryStepInputs *step,
                             NoahmpStabilityInput *canopy_input, BallBerryStepResults *results);

/* The VEGE_FLUX loop1 replica (see noahmp_vege_flux_exact_solve(), ball_berry_helpers.h), recording the
 * iteration trace in trace instead of writing it.  Returns 0 or -1. */
int noahmp_vege_flux_exact_solve_traced(const NoahmpStabilityInput *input, NoahmpStabilityOutput *output,
                                        NoahmpVegeFluxExactTrace *trace);

/* The trace rows of an exact solve, written to stream as the solver wrote them before 10/28 (nothing if
 * stream is NULL); ball_berry_helpers.c, not differentiated. */
void noahmp_vege_flux_exact_write_trace(FILE *stream, const NoahmpVegeFluxExactTrace *trace);

/* One method 6 step with the exact solver: the canopy input, then the no-canopy outputs or the exact
 * solution (in solution) and the coupling outputs.  Returns 0 or -1 (the solver failed). */
int ball_berry_method6_physics_step(const BallBerryMethod6Parameters *parameters, const BallBerryStepInputs *step,
                                    NoahmpStabilityInput *canopy_input, NoahmpStabilityOutput *solution,
                                    NoahmpVegeFluxExactTrace *trace, BallBerryStepResults *results);

/* 1 if x is finite (neither infinite nor NaN), else 0: isfinite() without the macro. */
int bb_is_finite(double x);

#endif /* BALL_BERRY_CANOPY_PHYSICS_H */
