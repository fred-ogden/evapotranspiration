#include "ball_berry_helpers.h"
#include "ball_berry_canopy_physics.h"
#include <math.h>
#include <stddef.h>


/* Original: noahmp_canopy_aerodynamics.c */
#include <stddef.h>

int noahmp_canopy_aerodynamics(const NoahmpCanopyAerodynamicsInputs *inputs,
                              NoahmpSfcdif1State *state,
                              NoahmpCanopyAerodynamicsOutputs *outputs)
{
    NoahmpSfcdif1State next_state;
    NoahmpCanopyAerodynamicsOutputs next_outputs;

    if (inputs == NULL || state == NULL || outputs == NULL) {
        return -1;
    }
    next_state = *state;
    if (noahmp_sfcdif1(&inputs->surface_layer, &next_state,
                       &next_outputs.surface_layer) != 0) {
        return -1;
    }
    if (noahmp_leaf_boundary_resistance(&inputs->leaf_boundary,
             &next_outputs.leaf_boundary_resistance_s_per_m) != 0) {
        return -1;
    }
    *state = next_state;
    *outputs = next_outputs;
    return 0;
}


/* Original: noahmp_canopy_iterative.c */
#include <math.h>
#include <stddef.h>

static int canopy_evaluate_with_resistance(const NoahmpCanopyIterativeInput *in,
                           double foliage_temperature_k,
                           double *residual_w_per_m2,
                           double *transpiration_m_per_s,
                           NoahmpRadiationStomataOutput *radiation_out,
                           NoahmpVegeFluxOutputs *foliage_out,
                           const double *frozen_rs_sun,
                           const double *frozen_rs_shade)
{
    NoahmpRadiationStomataInput radiation;
    NoahmpRadiationStomataOutput radiation_result;
    NoahmpVegeFluxInputs foliage;
    NoahmpVegeFluxOutputs step;
    double slope, ground_slope;
    if (in == NULL || residual_w_per_m2 == NULL ||
        transpiration_m_per_s == NULL || !isfinite(foliage_temperature_k)) return -1;
    radiation = in->radiation_stomata;
    foliage = in->foliage;
    foliage.foliage_temperature_k = foliage_temperature_k;
    if (noahmp_esat(foliage.ground_temperature_k,
                    &foliage.ground_saturation_vapor_pressure_pa,
                    &ground_slope) != 0) return -1;
    if (noahmp_esat(foliage_temperature_k,
                    &foliage.foliage_saturation_vapor_pressure_pa,
                    &slope) != 0) return -1;
    foliage.foliage_saturation_slope_pa_per_k = slope;
    radiation.radiation.tv_k = foliage_temperature_k;
    if (in->hold_radiation_at_start_temperature) {
        /* FLO 10/24: Noah-MP ENERGY calls RADIATION once, at the TV of the start of the step. Holding
         * it there also keeps the TWOSTREAM TV > TFRZ snow-optics switch out of the Newton residual. */
        radiation.radiation.tv_k = in->foliage.previous_foliage_temperature_k;
    }
    radiation.canopy_air.tv = foliage_temperature_k;
    radiation.canopy_air.ei = foliage.foliage_saturation_vapor_pressure_pa;
    radiation.canopy_air.rb = foliage.rb_s_per_m;
    if (noahmp_radiation_stomata(&radiation, &radiation_result) != 0) return -1;
    foliage.rs_sun_s_per_m = radiation_result.sunlit.rs;
    foliage.rs_shade_s_per_m = radiation_result.shaded.rs;
    if (frozen_rs_sun != NULL && frozen_rs_shade != NULL) {
        foliage.rs_sun_s_per_m = *frozen_rs_sun;
        foliage.rs_shade_s_per_m = *frozen_rs_shade;
        radiation_result.sunlit.rs = *frozen_rs_sun;
        radiation_result.shaded.rs = *frozen_rs_shade;
    }
    foliage.lai_sun = radiation_result.radiation.radiation_inputs.laisun;
    foliage.lai_shade = radiation_result.radiation.radiation_inputs.laisha;
    foliage.absorbed_shortwave_w_per_m2 = radiation_result.radiation.radiation.sav;
    if (noahmp_vege_flux_foliage_step(&foliage, &step) != 0) return -1;
    *residual_w_per_m2 = step.energy_balance_residual_w_per_m2;
    *transpiration_m_per_s = step.unlinearized_transpiration_w_per_m2 /
                             foliage.latent_heat_j_per_kg / 1000.0;
    if (!isfinite(*residual_w_per_m2) || !isfinite(*transpiration_m_per_s)) return -1;
    if (radiation_out != NULL) *radiation_out = radiation_result;
    if (foliage_out != NULL) *foliage_out = step;
    return 0;
}

int noahmp_canopy_evaluate(const NoahmpCanopyIterativeInput *in,
                           double foliage_temperature_k,
                           double *residual_w_per_m2,
                           double *transpiration_m_per_s,
                           NoahmpRadiationStomataOutput *radiation_out,
                           NoahmpVegeFluxOutputs *foliage_out)
{
    return canopy_evaluate_with_resistance(in, foliage_temperature_k,
               residual_w_per_m2, transpiration_m_per_s, radiation_out,
               foliage_out, NULL, NULL);
}

/* In reference mode, freeze the two stomatal resistances at the initial
 * foliage temperature. They remain constant for all Newton trial states.
 * This emulates the first-iteration stomatal update policy, not the full
 * Noah-MP VEGE_FLUX stability/ground-temperature solution. */
static int initial_frozen_resistance(const NoahmpCanopyIterativeInput *in,
                                     double *rs_sun, double *rs_shade)
{
    NoahmpRadiationStomataOutput radiation;
    double residual, transpiration;
    if (in->stomata_update_mode != NOAHMP_STOMATA_REFERENCE_FROZEN) return 0;
    if (noahmp_canopy_evaluate(in, in->foliage.foliage_temperature_k,
            &residual, &transpiration, &radiation, NULL) != 0) return -1;
    *rs_sun = radiation.sunlit.rs;
    *rs_shade = radiation.shaded.rs;
    return 0;
}

int noahmp_canopy_heat_sensitivity(const NoahmpCanopyIterativeInput *in,
                                   double tv,
                                   double *dtdq,
                                   double *dedq)
{
    double h = 1.0e-3, fp, fm, ep, em, derivative;
    double rs_sun = 0.0, rs_shade = 0.0;
    const double *sun = NULL, *shade = NULL;
    if (in == NULL || dtdq == NULL || dedq == NULL) return -1;
    if (in->stomata_update_mode == NOAHMP_STOMATA_REFERENCE_FROZEN) {
        if (in->use_prescribed_frozen_stomata) {
            rs_sun = in->prescribed_frozen_sunlit_resistance_s_per_m;
            rs_shade = in->prescribed_frozen_shaded_resistance_s_per_m;
            if (!isfinite(rs_sun) || !isfinite(rs_shade) ||
                rs_sun <= 0.0 || rs_shade <= 0.0) return -1;
        } else if (initial_frozen_resistance(in, &rs_sun, &rs_shade) != 0) return -1;
        sun = &rs_sun;
        shade = &rs_shade;
    }
    if (canopy_evaluate_with_resistance(in,tv+h,&fp,&ep,NULL,NULL,sun,shade) != 0 ||
        canopy_evaluate_with_resistance(in,tv-h,&fm,&em,NULL,NULL,sun,shade) != 0) return -1;
    derivative = (fp-fm)/(2.0*h);
    if (!isfinite(derivative) || fabs(derivative) < 1.0e-10) return -1;
    *dtdq = -1.0/derivative;
    *dedq = (ep-em)/(2.0*h)*(*dtdq);
    return 0;
}

int noahmp_canopy_iterative_solve(const NoahmpCanopyIterativeInput *in,
                                  NoahmpCanopyIterativeOutput *out)
{
    NoahmpCanopyIterativeOutput candidate = {0};
    double tv, residual, transpiration, derivative, fp, fm, ep, em;
    double delta, h = 1.0e-3;
    int iteration;
    double rs_sun = 0.0, rs_shade = 0.0;
    const double *sun = NULL, *shade = NULL;
    if (in == NULL || out == NULL) return -1;
    if (in->maximum_iterations < 1 || in->maximum_iterations > 100 ||
        !isfinite(in->temperature_tolerance_k) || in->temperature_tolerance_k <= 0.0 ||
        !isfinite(in->maximum_temperature_step_k) || in->maximum_temperature_step_k <= 0.0)
        return -1;
    if (in->stomata_update_mode != NOAHMP_STOMATA_REEVALUATE &&
        in->stomata_update_mode != NOAHMP_STOMATA_REFERENCE_FROZEN) return -1;
    if (in->stomata_update_mode == NOAHMP_STOMATA_REFERENCE_FROZEN) {
        if (in->use_prescribed_frozen_stomata) {
            rs_sun = in->prescribed_frozen_sunlit_resistance_s_per_m;
            rs_shade = in->prescribed_frozen_shaded_resistance_s_per_m;
            if (!isfinite(rs_sun) || !isfinite(rs_shade) ||
                rs_sun <= 0.0 || rs_shade <= 0.0) return -1;
        } else if (initial_frozen_resistance(in, &rs_sun, &rs_shade) != 0) return -1;
        sun = &rs_sun;
        shade = &rs_shade;
    }
    tv = in->foliage.foliage_temperature_k;
    for (iteration = 0; iteration < in->maximum_iterations; ++iteration) {
        if (canopy_evaluate_with_resistance(in,tv,&residual,&transpiration,NULL,NULL,sun,shade) != 0 ||
            canopy_evaluate_with_resistance(in,tv+h,&fp,&ep,NULL,NULL,sun,shade) != 0 ||
            canopy_evaluate_with_resistance(in,tv-h,&fm,&em,NULL,NULL,sun,shade) != 0) return -1;
        derivative = (fp-fm)/(2.0*h);
        if (!isfinite(derivative) || fabs(derivative) < 1.0e-10) return -1;
        delta = -residual/derivative;
        if (!isfinite(delta) || fabs(delta) > in->maximum_temperature_step_k) return -1;
        tv += delta;
        if (fabs(delta) <= in->temperature_tolerance_k) break;
    }
    if (iteration == in->maximum_iterations) return -1;
    if (canopy_evaluate_with_resistance(in,tv,&residual,&transpiration,
                               &candidate.radiation_stomata,
                               &candidate.foliage.foliage,sun,shade) != 0) return -1;
    if (sun != NULL) {
        candidate.frozen_sunlit_resistance_s_per_m = rs_sun;
        candidate.frozen_shaded_resistance_s_per_m = rs_shade;
    }
    candidate.foliage.foliage_temperature_k = tv;
    candidate.foliage.ground_skin_temperature_k = in->foliage.ground_temperature_k;
    candidate.foliage.transpiration_m_per_s = transpiration;
    candidate.iterations = iteration+1;
    candidate.final_temperature_increment_k = delta;
    candidate.energy_balance_residual_w_per_m2 = residual;
    *out = candidate;
    return 0;
}


/* Original: noahmp_leaf_boundary_resistance.c */
#include <math.h>
#include <stddef.h>

int noahmp_leaf_boundary_resistance(const NoahmpLeafBoundaryInputs *inputs,
                                   double *resistance_s_per_m)
{
    double cwpc;
    double attenuation_ratio;
    double resistance;

    if (inputs == NULL || resistance_s_per_m == NULL) {
        return -1;
    }
    if (!isfinite(inputs->cwp) || inputs->cwp < 0.0 ||
        !isfinite(inputs->vegetation_area_index) || inputs->vegetation_area_index < 0.0 ||
        !isfinite(inputs->canopy_height_m) || inputs->canopy_height_m < 0.0 ||
        !isfinite(inputs->fhg) || inputs->fhg < 0.0 ||
        !isfinite(inputs->leaf_dimension_m) || inputs->leaf_dimension_m <= 0.0 ||
        !isfinite(inputs->canopy_wind_speed_m_per_s) ||
        inputs->canopy_wind_speed_m_per_s <= 0.0) {
        return -1;
    }

    /* Noah-MP RAGRB CWP branch:
     * CWPC = sqrt(CWP * VAI * HCAN * FHG)
     * RB = max(50 * CWPC/(1-exp(-CWPC/2)) * sqrt(DLEAF/UC), 20)
     * FLO 10/26: evaluate the exact zero-CWPC limit, which is 2.
     */
    cwpc = sqrt(inputs->cwp * inputs->vegetation_area_index *
                inputs->canopy_height_m * inputs->fhg);
    if (cwpc == 0.0) {
        attenuation_ratio = 2.0;
    } else {
        attenuation_ratio = cwpc / (-expm1(-0.5 * cwpc));
    }
    resistance = fmax(50.0 * attenuation_ratio *
                      sqrt(inputs->leaf_dimension_m /
                           inputs->canopy_wind_speed_m_per_s), 20.0);
    if (!isfinite(resistance)) {
        return -1;
    }
    *resistance_s_per_m = resistance;
    return 0;
}


/* Original: noahmp_prescribed_tg.c */
#include <math.h>

int noahmp_prescribed_tg_solve(const NoahmpVegeFluxInputs *inputs,
                              const NoahmpPrescribedTgOptions *options,
                              NoahmpPrescribedTgOutputs *outputs)
{
    NoahmpVegeFluxInputs state;
    NoahmpVegeFluxOutputs step;
    NoahmpPrescribedTgOutputs result;
    int iteration;
    if (inputs == 0 || options == 0 || outputs == 0) return -1;
    if (options->iterations < 1 || options->iterations > 100 ||
        !isfinite(options->maximum_temperature_step_k) ||
        options->maximum_temperature_step_k <= 0.0 ||
        !isfinite(inputs->ground_temperature_k) ||
        !isfinite(inputs->foliage_temperature_k)) return -1;
    state = *inputs;
    if (noahmp_esat(state.ground_temperature_k,
                    &state.ground_saturation_vapor_pressure_pa,
                    &result.foliage.canopy_air_vapor_pressure_pa) != 0) return -1;
    for (iteration = 0; iteration < options->iterations; ++iteration) {
        if (noahmp_esat(state.foliage_temperature_k,
                        &state.foliage_saturation_vapor_pressure_pa,
                        &state.foliage_saturation_slope_pa_per_k) != 0) return -1;
        if (noahmp_vege_flux_foliage_step(&state, &step) != 0) return -1;
        if (!isfinite(step.foliage_temperature_increment_k) ||
            fabs(step.foliage_temperature_increment_k) > options->maximum_temperature_step_k) return -1;
        state.foliage_temperature_k = step.updated_foliage_temperature_k;
    }
    /* Re-evaluate fluxes at the final TV, without taking another temperature update. */
    if (noahmp_esat(state.foliage_temperature_k,
                    &state.foliage_saturation_vapor_pressure_pa,
                    &state.foliage_saturation_slope_pa_per_k) != 0) return -1;
    if (noahmp_vege_flux_foliage_step(&state, &step) != 0) return -1;
    if (!isfinite(step.transpiration_w_per_m2)) return -1;
    result.foliage = step;
    result.foliage_temperature_k = state.foliage_temperature_k;
    result.ground_skin_temperature_k = state.ground_temperature_k;
    result.transpiration_m_per_s =
        step.transpiration_w_per_m2 / state.latent_heat_j_per_kg / 1000.0;
    if (!isfinite(result.transpiration_m_per_s)) return -1;
    *outputs = result;
    return 0;
}


/* Original: noahmp_radiation_stomata.c */
#include <math.h>
#include <stddef.h>
/* Follows VEGE_FLUX sunlit/shaded STOMATA calls: shared canopy-air state,
 * PARSUN/PARSHA from SURRAD, LAISUN/LAISHA for canopy integration.
 * No aerodynamic or foliage energy-balance assumptions are introduced.
 */
int noahmp_radiation_stomata(const NoahmpRadiationStomataInput *in,
                            NoahmpRadiationStomataOutput *out)
{
    NoahmpRadiationStomataOutput candidate;
    NoahmpStomataInputs leaf;
    double laisun, laisha;
    if (in == NULL || out == NULL) return -1;
    if (noahmp_albedo_radiation(&in->radiation, &candidate.radiation) != 0)
        return -1;
    leaf = in->canopy_air;
    leaf.apar = candidate.radiation.radiation.parsun;
    if (noahmp_stomata(&in->physiology, &leaf,
                                &candidate.sunlit) != 0) return -1;
    leaf.apar = candidate.radiation.radiation.parsha;
    if (noahmp_stomata(&in->physiology, &leaf,
                                &candidate.shaded) != 0) return -1;
    laisun = candidate.radiation.radiation_inputs.laisun;
    laisha = candidate.radiation.radiation_inputs.laisha;
    candidate.photosynthesis_umol_m2_s = candidate.sunlit.psn * laisun +
                                          candidate.shaded.psn * laisha;
    if (!isfinite(candidate.photosynthesis_umol_m2_s)) return -1;
    *out = candidate;
    return 0;
}


/* Original: noahmp_stability_coupling.c */
#include <math.h>
#include <stddef.h>
#ifdef NOAHMP_TRACE_STABILITY
#include <stdio.h>
#endif

/* SFCDIF1 CH is a dimensionless bulk coefficient. Its corresponding
 * resistance is 1/(CH*UR), not 1/CH. */
int noahmp_stability_coupled_solve(const NoahmpStabilityInput *input,
                                  NoahmpStabilityOutput *output)
{
    NoahmpSfcdif1State sf_state = {0};
    NoahmpRagrbState ground_state = {0};
    NoahmpStabilityOutput candidate = {0};
    NoahmpCanopyIterativeInput canopy;
    NoahmpSfcdif1Inputs sf;
    NoahmpRagrbInputs rg;
    double sensible_heat = 0.0;
    double ground_sensible_heat = 0.0;
    double canopy_air_temperature;
    double previous_rahc = 0.0;
    double previous_rahg = 0.0;
    double previous_rb = 0.0;
    double rahc, change = 0.0;
    int i;
    double first_rs_sun = 0.0;
    double first_rs_shade = 0.0;
    double vegetation_temperature_at_iteration_start_k;     /* FLO 10/24: NWM outer rule */
    int noahmp_last_iteration_flag = 0;
    int maximum_outer_iterations;
    if (input == NULL || output == NULL ||
        input->maximum_stability_iterations < 2 ||
        input->maximum_stability_iterations > 100 ||
        !isfinite(input->resistance_relative_tolerance) ||
        input->resistance_relative_tolerance <= 0.0 ||
        input->surface_layer.wind_speed_m_per_s <= 0.0) return -1;
    canopy = input->canopy;
    sf = input->surface_layer;
    rg = input->ground_canopy;
    canopy_air_temperature = canopy.foliage.air_temperature_k;
    maximum_outer_iterations = input->maximum_stability_iterations;
    if (input->noahmp_outer_iteration_rule) maximum_outer_iterations = 20;      /* NITERC */
    for (i = 1; i <= maximum_outer_iterations; ++i) {
        vegetation_temperature_at_iteration_start_k = canopy.foliage.foliage_temperature_k;
        sf.iteration = i;
        sf.sensible_heat_w_per_m2 = sensible_heat;
        if (noahmp_sfcdif1(&sf, &sf_state, &candidate.surface_layer) != 0) return -1;
        if (candidate.surface_layer.ch <= 0.0) return -1;
        rahc = 1.0 / (candidate.surface_layer.ch * sf.wind_speed_m_per_s);
        rg.iteration = i;
        rg.friction_velocity_m_per_s = candidate.surface_layer.friction_velocity_m_per_s;
        rg.ground_sensible_heat_w_per_m2 = ground_sensible_heat;
        rg.canopy_air_temperature_k = canopy_air_temperature;
        if (noahmp_ragrb(&rg, &ground_state, &candidate.ground_canopy) != 0) return -1;
        canopy.foliage.rahc_s_per_m = rahc;
        canopy.foliage.rawc_s_per_m = rahc;
        canopy.foliage.rahg_s_per_m = candidate.ground_canopy.ground_heat_resistance_s_per_m;
        canopy.foliage.rawg_s_per_m = candidate.ground_canopy.ground_vapor_resistance_s_per_m;
        canopy.foliage.rb_s_per_m = candidate.ground_canopy.leaf_boundary_resistance_s_per_m;
        if (canopy.stomata_update_mode == NOAHMP_STOMATA_REFERENCE_FROZEN && i > 1) {
            canopy.use_prescribed_frozen_stomata = 1;
            canopy.prescribed_frozen_sunlit_resistance_s_per_m = first_rs_sun;
            canopy.prescribed_frozen_shaded_resistance_s_per_m = first_rs_shade;
        }
        if (noahmp_canopy_iterative_solve(&canopy, &candidate.canopy) != 0) return -1;
        if (canopy.stomata_update_mode == NOAHMP_STOMATA_REFERENCE_FROZEN && i == 1) {
            first_rs_sun = candidate.canopy.frozen_sunlit_resistance_s_per_m;
            first_rs_shade = candidate.canopy.frozen_shaded_resistance_s_per_m;
        }
#ifdef NOAHMP_TRACE_STABILITY
        if (canopy.stomata_update_mode == NOAHMP_STOMATA_REFERENCE_FROZEN) {
            fprintf(stderr, "TRACE %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                i, candidate.canopy.radiation_stomata.radiation.radiation.sav,
                candidate.canopy.radiation_stomata.radiation.radiation_inputs.laisun,
                candidate.canopy.radiation_stomata.radiation.radiation_inputs.laisha,
                candidate.canopy.frozen_sunlit_resistance_s_per_m,
                candidate.canopy.frozen_shaded_resistance_s_per_m,
                rahc, canopy.foliage.rahg_s_per_m, canopy.foliage.rb_s_per_m,
                candidate.canopy.foliage.foliage_temperature_k,
                candidate.canopy.foliage.transpiration_m_per_s,
                candidate.surface_layer.friction_velocity_m_per_s,
                candidate.ground_canopy.canopy_wind_speed_m_per_s,
                candidate.ground_canopy.ground_heat_resistance_s_per_m);
        }
#endif
        /* FLO 10/20: Noah-MP passes SFCDIF1 the canopy-air-to-reference flux
         * H = RHOAIR*CPAIR*(TAH-SFCTMP)/RAHC = rho*cp*[CVH*(TV-TAH) + CGH*(TG-TAH)].
         * SHC carries the FVEG grid-fraction weight, so divide it back out
         * before adding the under-canopy ground flux (was SHC + HG, which
         * understated H by (1-FVEG)*SHC and biased RAHC high). */
        sensible_heat = candidate.canopy.foliage.foliage.sensible_heat_w_per_m2 /
                        canopy.foliage.fveg;
        canopy_air_temperature = candidate.canopy.foliage.foliage.canopy_air_temperature_k;
        ground_sensible_heat = canopy.foliage.air_density_kg_per_m3 *
             canopy.foliage.cp_air_j_per_kg_k *
             (canopy.foliage.ground_temperature_k - canopy_air_temperature) /
             canopy.foliage.rahg_s_per_m;
        sensible_heat += ground_sensible_heat;
        if (input->noahmp_outer_iteration_rule) {
            /* FLO 10/24: NWM loop1 exit rule, on the change of the converged TV over this iteration */
            double temperature_change_k = candidate.canopy.foliage.foliage_temperature_k -
                                          vegetation_temperature_at_iteration_start_k;
            if (i > 1) {
                change = fmax(fabs(rahc-previous_rahc)/fmax(rahc, 1.0),
                              fabs(canopy.foliage.rahg_s_per_m-previous_rahg)/fmax(fabs(canopy.foliage.rahg_s_per_m),1.0));
            }
            if (noahmp_last_iteration_flag == 1) break;
            if (i >= 5 && fabs(temperature_change_k) <= 0.01) noahmp_last_iteration_flag = 1;
            if (i == maximum_outer_iterations) break;               /* accepted, as NWM */
        } else if (i > 1) {
            double c1 = fabs(rahc-previous_rahc)/fmax(rahc, 1.0);
            double c2 = fabs(canopy.foliage.rahg_s_per_m-previous_rahg)/
                        fmax(fabs(canopy.foliage.rahg_s_per_m),1.0);
            double c3 = fabs(canopy.foliage.rb_s_per_m-previous_rb)/
                        fmax(canopy.foliage.rb_s_per_m,1.0);
            change = fmax(c1,fmax(c2,c3));
            if (change < input->resistance_relative_tolerance) break;
        }
        previous_rahc = rahc;
        previous_rahg = canopy.foliage.rahg_s_per_m;
        previous_rb = canopy.foliage.rb_s_per_m;
        canopy.foliage.foliage_temperature_k = candidate.canopy.foliage.foliage_temperature_k;
    }
    if (!input->noahmp_outer_iteration_rule && i > input->maximum_stability_iterations) return -1;
    candidate.stability_iterations = i;
    candidate.relative_resistance_change = change;
    candidate.reference_canopy_resistance_s_per_m = rahc;
    *output = candidate;
    return 0;
}


/* FLO 10/28: the functions of the VEGE_FLUX loop1 replica (ALBEDO/TWOSTREAM/SURRAD, ESAT, SFCDIF1,
 * RAGRB, the foliage step and the exact solver) moved to ball_berry_canopy_physics.c, in the form
 * Tapenade differentiates; the nested solver above calls them there.  This is the exact solver with its
 * optional per-iteration trace written to input->iteration_trace_stream, as before 10/28. */
int noahmp_vege_flux_exact_solve(const NoahmpStabilityInput *input,
                                 NoahmpStabilityOutput *output)
{
    static NoahmpVegeFluxExactTrace trace;
    NoahmpStabilityOutput candidate;
    int status;
    if (input == NULL || output == NULL) return -1;
    status = noahmp_vege_flux_exact_solve_traced(input, &candidate, &trace);
    noahmp_vege_flux_exact_write_trace(input->iteration_trace_stream, &trace);
    if (status == 0) *output = candidate;
    return status;
}

/* The trace rows of an exact solve, as the solver wrote them before 10/28 (nothing if stream is NULL). */
void noahmp_vege_flux_exact_write_trace(FILE *stream, const NoahmpVegeFluxExactTrace *trace)
{
    int row;
    const NoahmpVegeFluxExactTraceRecord *r;
    if (stream == NULL || trace == NULL || trace->started == 0) return;
    fprintf(stream,
            "iteration,vegetation_temperature_k,vegetation_temperature_change_k,"
            "canopy_air_temperature_k,canopy_air_vapor_pressure_pa,"
            "sunlit_stomatal_resistance_s_per_m,shaded_stomatal_resistance_s_per_m,"
            "leaf_boundary_layer_resistance_s_per_m,above_canopy_heat_resistance_s_per_m,"
            "under_canopy_heat_resistance_s_per_m,momentum_drag_coefficient,"
            "heat_exchange_coefficient,friction_velocity_m_per_s,"
            "monin_obukhov_stability_parameter,canopy_net_longwave_w_per_m2,"
            "canopy_sensible_heat_w_per_m2,transpiration_heat_w_per_m2\n");
    for (row = 0; row < trace->count; ++row) {
        r = &trace->record[row];
        fprintf(stream,
                "%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                r->iteration, r->vegetation_temperature_k, r->vegetation_temperature_change_k,
                r->canopy_air_temperature_k, r->canopy_air_vapor_pressure_pa,
                r->sunlit_stomatal_resistance_s_per_m, r->shaded_stomatal_resistance_s_per_m,
                r->leaf_boundary_layer_resistance_s_per_m, r->above_canopy_heat_resistance_s_per_m,
                r->under_canopy_heat_resistance_s_per_m, r->momentum_drag_coefficient,
                r->heat_exchange_coefficient, r->friction_velocity_m_per_s,
                r->monin_obukhov_stability_parameter, r->canopy_net_longwave_w_per_m2,
                r->canopy_sensible_heat_w_per_m2, r->transpiration_heat_w_per_m2);
    }
}
