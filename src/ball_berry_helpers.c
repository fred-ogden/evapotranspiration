#include "ball_berry_helpers.h"
#include <math.h>
#include <stddef.h>


/* Original: noahmp_albedo.c */
#include <math.h>
#include <stddef.h>
#include <string.h>

int noahmp_albedo_radiation(const NoahmpAlbedoInput *in, NoahmpAlbedoOutput *out)
{
    const double mpe = 1.0e-6;
    NoahmpAlbedoOutput result;
    NoahmpTwostreamInput tw_in;
    NoahmpTwostreamOutput tw_out;
    NoahmpSurradInput *rad;
    double vai, wl, ws, rho[2], tau[2], ext, gdir;
    int ib, ic;
    if (in == NULL || out == NULL) return 1;
    if (!isfinite(in->cosz) || !isfinite(in->elai) || !isfinite(in->esai) ||
        !isfinite(in->fveg) || !isfinite(in->fwet) ||
        in->elai < 0.0 || in->esai < 0.0 || in->fveg < 0.0 ||
        in->fveg > 1.0 || in->fwet < 0.0 || in->fwet > 1.0) return 1;
    vai = in->elai + in->esai;
    /* The existing TWOSTREAM port requires positive VAI and daylight. */
    if (vai <= 0.0 || in->cosz <= 0.0 || in->cosz > 1.0) return 1;
    memset(&result, 0, sizeof(result));
    memset(&tw_in, 0, sizeof(tw_in));
    memset(&tw_out, 0, sizeof(tw_out));
    rad = &result.radiation_inputs;
    wl = in->elai / fmax(vai, mpe);
    ws = in->esai / fmax(vai, mpe);
    tw_in.cosz = in->cosz;
    tw_in.vai = vai;
    tw_in.fveg = in->fveg;
    tw_in.fwet = in->fwet;
    tw_in.temperature_k = in->tv_k;
    tw_in.xl = in->xl;
    tw_in.betads = in->betads;
    tw_in.betais = in->betais;
    for (ib = 0; ib < 2; ib++) {
        rho[ib] = fmax(in->rhol[ib] * wl + in->rhos[ib] * ws, mpe);
        tau[ib] = fmax(in->taul[ib] * wl + in->taus[ib] * ws, mpe);
        tw_in.rho[ib] = rho[ib];
        tw_in.tau[ib] = tau[ib];
        tw_in.albgrd[ib] = in->albgrd[ib];
        tw_in.albgri[ib] = in->albgri[ib];
        tw_in.omegas[ib] = in->omegas[ib];
        rad->albgrd[ib] = in->albgrd[ib];
        rad->albgri[ib] = in->albgri[ib];
        rad->solad[ib] = in->solad[ib];
        rad->solai[ib] = in->solai[ib];
        if (!isfinite(rho[ib]) || !isfinite(tau[ib]) ||
            !isfinite(in->solad[ib]) || !isfinite(in->solai[ib])) return 1;
    }
    for (ib = 0; ib < 2; ib++) {
        for (ic = 0; ic < 2; ic++) {
            tw_in.ib = ib;
            tw_in.ic = ic;
            if (noahmp_twostream(&tw_in, &tw_out) != 0) return 1;
            if (ic == 0) {
                rad->fabd[ib] = tw_out.fab[ib];
                rad->albd[ib] = tw_out.fre[ib];
                rad->ftdd[ib] = tw_out.ftd[ib];
                rad->ftid[ib] = tw_out.fti[ib];
                rad->frevd[ib] = tw_out.frev[ib];
                rad->fregd[ib] = tw_out.freg[ib];
            } else {
                rad->fabi[ib] = tw_out.fab[ib];
                rad->albi[ib] = tw_out.fre[ib];
                rad->ftii[ib] = tw_out.fti[ib];
                rad->frevi[ib] = tw_out.frev[ib];
                rad->fregi[ib] = tw_out.freg[ib];
            }
        }
    }
    gdir = tw_out.gdir;
    ext = gdir / in->cosz * sqrt(1.0 - rho[0] - tau[0]);
    if (!isfinite(ext)) return 1;
    result.fsun = (1.0 - exp(-ext * vai)) / fmax(ext * vai, mpe);
    if (result.fsun < 0.01) result.fsun = 0.0;
    rad->fsun = result.fsun;
    rad->fsha = 1.0 - result.fsun;
    rad->elai = in->elai;
    rad->vai = vai;
    rad->laisun = in->elai * result.fsun;
    rad->laisha = in->elai * (1.0 - result.fsun);
    rad->mpe = mpe;
    if (noahmp_surrad(rad, &result.radiation) != 0) return 1;
    *out = result;
    return 0;
}


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

/* ESAT polynomial from module_sf_noahmplsm.F, with water/ice switch at 0 C.
 * FLO 10/20: argument formed exactly as VEGE_FLUX does before calling ESAT:
 * TDC(T) = MIN(50, MAX(-50, T - TFRZ)) with TFRZ = 273.16 K (was T - 273.15,
 * unclamped; that gave a 0.066 percent high EI at 285 K). */
int noahmp_esat(double temperature_k, double *esat_pa, double *slope_pa_per_k)
{
    const double noahmp_freezing_point_k = 273.16;   /* TFRZ */
    double t, es, slope;
    if (esat_pa == 0 || slope_pa_per_k == 0 || !isfinite(temperature_k)) return -1;
    t = fmin(50.0, fmax(-50.0, temperature_k - noahmp_freezing_point_k));
    if (t > 0.0) {
        es = 100.0 * (6.107799961 + t * (4.436518521e-1 + t *
             (1.428945805e-2 + t * (2.650648471e-4 + t *
             (3.031240396e-6 + t * (2.034080948e-8 + t * 6.136820929e-11))))));
        slope = 100.0 * (4.438099984e-1 + t * (2.857002636e-2 + t *
             (7.938054040e-4 + t * (1.215215065e-5 + t *
             (1.036561403e-7 + t * (3.532421810e-10 - t * 7.090244804e-13))))));
    } else {
        es = 100.0 * (6.109177956 + t * (5.034698970e-1 + t *
             (1.886013408e-2 + t * (4.176223716e-4 + t *
             (5.824720280e-6 + t * (4.838803174e-8 + t * 1.838826904e-10))))));
        slope = 100.0 * (5.030305237e-1 + t * (3.773255020e-2 + t *
             (1.267995369e-3 + t * (2.477563108e-5 + t *
             (3.005693132e-7 + t * (2.158542548e-9 + t * 7.131097725e-12))))));
    }
    if (!isfinite(es) || !isfinite(slope) || es <= 0.0 || slope <= 0.0) return -1;
    *esat_pa = es;
    *slope_pa_per_k = slope;
    return 0;
}

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


/* Original: noahmp_ragrb.c */
#include <math.h>
#include <stddef.h>

/* Port of VEGE_FLUX UC assignment and RAGRB, OPT_SFC=1 path.
 * Noah-MP constants VKC = 0.40, GRAV = 9.80616, CPAIR = 1004.64 (module_sf_noahmplsm.F).
 * FLO 10/21: GRAV and CPAIR were 9.81 and 1004.6, which made RAHG differ from Noah-MP by
 * about 6e-6 relative from the second canopy iteration on (through the under-canopy
 * Monin-Obukhov length); with the Noah-MP values the exact solver matches the r8 reference.
 * Invalid or singular domains are rejected instead of generating NaNs.
 */
int noahmp_ragrb(const NoahmpRagrbInputs *in,
                 NoahmpRagrbState *state, NoahmpRagrbOutputs *out)
{
    const double vkc = 0.4;
    const double grav = 9.80616;     /* GRAV  (m/s2) */
    const double cpair = 1004.64;    /* CPAIR (J/kg/K) */
    double uc, tmp1, molg, mozg, fhgnew, fhg, cwpc, tmp2;
    double tmprah2, kh, rb, rahg, log_denom;
    NoahmpRagrbOutputs result;
    NoahmpRagrbState next;
    if (in == NULL || state == NULL || out == NULL) return -1;
    if (in->iteration < 1 || !isfinite(in->reference_wind_speed_m_per_s) ||
        !isfinite(in->reference_height_m) || !isfinite(in->momentum_roughness_m) ||
        !isfinite(in->displacement_height_m) || !isfinite(in->canopy_height_m) ||
        !isfinite(in->vegetation_area_index) || !isfinite(in->air_density_kg_per_m3) ||
        !isfinite(in->ground_sensible_heat_w_per_m2) ||
        !isfinite(in->canopy_air_temperature_k) ||
        !isfinite(in->ground_momentum_roughness_m) ||
        !isfinite(in->ground_heat_roughness_m) ||
        !isfinite(in->canopy_heat_roughness_m) ||
        !isfinite(in->friction_velocity_m_per_s) || !isfinite(in->cwp) ||
        !isfinite(in->leaf_dimension_m) || !isfinite(in->mpe) ||
        in->reference_wind_speed_m_per_s <= 0.0 ||
        in->reference_height_m <= in->momentum_roughness_m ||
        in->momentum_roughness_m <= 0.0 || in->canopy_height_m <= in->displacement_height_m ||
        in->ground_momentum_roughness_m <= 0.0 || in->ground_heat_roughness_m < 0.0 ||
        in->canopy_heat_roughness_m < 0.0 || in->vegetation_area_index <= 0.0 ||
        in->air_density_kg_per_m3 <= 0.0 || in->canopy_air_temperature_k <= 0.0 ||
        in->friction_velocity_m_per_s <= 0.0 || in->cwp <= 0.0 ||
        in->leaf_dimension_m <= 0.0 || in->mpe <= 0.0) return -1;
    log_denom = log(in->reference_height_m / in->momentum_roughness_m);
    if (log_denom == 0.0) return -1;
    /* The second (overwriting) UC expression in VEGE_FLUX is authoritative. */
    uc = in->reference_wind_speed_m_per_s *
         log((in->canopy_height_m - in->displacement_height_m +
              in->momentum_roughness_m) / in->momentum_roughness_m) / log_denom;
    if (!isfinite(uc) || uc <= 0.0) return -1;
    mozg = 0.0;
    if (in->iteration > 1) {
        tmp1 = vkc * (grav / in->canopy_air_temperature_k) *
               in->ground_sensible_heat_w_per_m2 /
               (in->air_density_kg_per_m3 * cpair);
        if (fabs(tmp1) <= in->mpe) tmp1 = in->mpe;
        molg = -pow(in->friction_velocity_m_per_s, 3.0) / tmp1;
        mozg = fmin((in->displacement_height_m -
                     in->ground_momentum_roughness_m) / molg, 1.0);
    }
    if (mozg < 0.0) fhgnew = pow(1.0 - 15.0 * mozg, -0.25);
    else fhgnew = 1.0 + 4.7 * mozg;
    if (in->iteration == 1) fhg = fhgnew;
    else fhg = 0.5 * (state->fhg + fhgnew);
    cwpc = sqrt(in->cwp * in->vegetation_area_index * in->canopy_height_m * fhg);
    if (!isfinite(cwpc) || cwpc <= 0.0) return -1;
    tmp1 = exp(-cwpc * in->ground_heat_roughness_m / in->canopy_height_m);
    tmp2 = exp(-cwpc * (in->canopy_heat_roughness_m +
                        in->displacement_height_m) / in->canopy_height_m);
    tmprah2 = in->canopy_height_m * exp(cwpc) / cwpc * (tmp1 - tmp2);
    kh = fmax(vkc * in->friction_velocity_m_per_s *
              (in->canopy_height_m - in->displacement_height_m), in->mpe);
    rahg = tmprah2 / kh;
    rb = cwpc * 50.0 / (1.0 - exp(-cwpc / 2.0));
    rb = fmax(rb * sqrt(in->leaf_dimension_m / uc), 20.0);
    if (!isfinite(rahg) || !isfinite(rb)) return -1;
    next.mozg = mozg;
    next.fhg = fhg;
    result.canopy_wind_speed_m_per_s = uc;
    result.ground_momentum_resistance_s_per_m = 0.0;
    result.ground_heat_resistance_s_per_m = rahg;
    result.ground_vapor_resistance_s_per_m = rahg;
    result.leaf_boundary_resistance_s_per_m = rb;
    *state = next;
    *out = result;
    return 0;
}


/* Original: noahmp_sfcdif1.c */
#include <math.h>
#include <stddef.h>

/* Constants from supplied module_sf_noahmplsm.F. */
#define VKC 0.40
#define GRAV 9.80616
#define CPAIR 1004.64

static void stability_functions(double moz, double *fm, double *fh)
{
    if (moz < 0.0) {
        double x = pow(1.0 - 16.0 * moz, 0.25);
        double a = log((1.0 + x * x) / 2.0);
        double b = log((1.0 + x) / 2.0);
        *fm = 2.0 * b + a - 2.0 * atan(x) + 1.5707963;
        *fh = 2.0 * a;
    } else {
        *fm = -5.0 * moz;
        *fh = *fm;
    }
}

int noahmp_sfcdif1(const NoahmpSfcdif1Inputs *in,
                   NoahmpSfcdif1State *state, NoahmpSfcdif1Outputs *out)
{
    double old_moz, moz2, mol, tmp1, tvir;
    double tmpcm, tmpch, tmpcm2, tmpch2;
    double fmnew, fhnew, fm2new, fh2new;
    double cmfm, chfh, cm2fm2, ch2fh2;
    if (in == NULL || state == NULL || out == NULL) return -1;
    if (in->iteration < 1 || in->reference_height_m <= in->displacement_height_m ||
        in->momentum_roughness_m <= 0.0 || in->heat_roughness_m <= 0.0 ||
        in->air_density_kg_per_m3 <= 0.0 || in->air_temperature_k <= 0.0 ||
        in->wind_speed_m_per_s < 0.0 || in->mpe <= 0.0) return -1;
    old_moz = state->moz;
    tmpcm = log((in->reference_height_m-in->displacement_height_m)/in->momentum_roughness_m);
    tmpch = log((in->reference_height_m-in->displacement_height_m)/in->heat_roughness_m);
    tmpcm2 = log((2.0+in->momentum_roughness_m)/in->momentum_roughness_m);
    tmpch2 = log((2.0+in->heat_roughness_m)/in->heat_roughness_m);
    if (in->iteration == 1) {
        state->friction_velocity_m_per_s = 0.0;
        state->moz = 0.0;
        moz2 = 0.0;
    } else {
        tvir = (1.0+0.61*in->specific_humidity_kg_per_kg)*in->air_temperature_k;
        tmp1 = VKC*(GRAV/tvir)*in->sensible_heat_w_per_m2/
               (in->air_density_kg_per_m3*CPAIR);
        if (fabs(tmp1) <= in->mpe) tmp1 = in->mpe;
        mol = -pow(state->friction_velocity_m_per_s,3.0)/tmp1;
        if (mol == 0.0) return -1;
        state->moz = fmin((in->reference_height_m-in->displacement_height_m)/mol,1.0);
        moz2 = fmin((2.0+in->heat_roughness_m)/mol,1.0);
    }
    if (old_moz*state->moz < 0.0) state->moz_sign_changes++;
    if (state->moz_sign_changes >= 2) {
        state->moz=0.0; state->fm=0.0; state->fh=0.0;
        moz2=0.0; state->fm2=0.0; state->fh2=0.0;
    }
    stability_functions(state->moz,&fmnew,&fhnew);
    stability_functions(moz2,&fm2new,&fh2new);
    if (in->iteration == 1) {
        state->fm=fmnew; state->fh=fhnew;
        state->fm2=fm2new; state->fh2=fh2new;
    } else {
        state->fm=0.5*(state->fm+fmnew);
        state->fh=0.5*(state->fh+fhnew);
        state->fm2=0.5*(state->fm2+fm2new);
        state->fh2=0.5*(state->fh2+fh2new);
    }
    state->fh=fmin(state->fh,0.9*tmpch);
    state->fm=fmin(state->fm,0.9*tmpcm);
    state->fh2=fmin(state->fh2,0.9*tmpch2);
    state->fm2=fmin(state->fm2,0.9*tmpcm2);
    cmfm=tmpcm-state->fm; chfh=tmpch-state->fh;
    cm2fm2=tmpcm2-state->fm2; ch2fh2=tmpch2-state->fh2;
    if (fabs(cmfm)<=in->mpe) cmfm=in->mpe;
    if (fabs(chfh)<=in->mpe) chfh=in->mpe;
    if (fabs(cm2fm2)<=in->mpe) cm2fm2=in->mpe;
    if (fabs(ch2fh2)<=in->mpe) ch2fh2=in->mpe;
    out->cm=VKC*VKC/(cmfm*cmfm);
    out->ch=VKC*VKC/(cmfm*chfh);
    out->ch2=VKC*VKC/(cm2fm2*ch2fh2);
    state->friction_velocity_m_per_s=in->wind_speed_m_per_s*sqrt(out->cm);
    out->friction_velocity_m_per_s=state->friction_velocity_m_per_s;
    out->ch2=VKC*state->friction_velocity_m_per_s/ch2fh2;
    if (!isfinite(out->cm) || !isfinite(out->ch) || !isfinite(out->ch2)) return -1;
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
    for (i = 1; i <= input->maximum_stability_iterations; ++i) {
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
        if (i > 1) {
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
    if (i > input->maximum_stability_iterations) return -1;
    candidate.stability_iterations = i;
    candidate.relative_resistance_change = change;
    candidate.reference_canopy_resistance_s_per_m = rahc;
    *output = candidate;
    return 0;
}


/* Original: noahmp_surrad.c */
#include <math.h>
int noahmp_surrad(const NoahmpSurradInput *in, NoahmpSurradOutput *out)
{
    NoahmpSurradOutput result = {0};
    double cad[2], cai[2], trd, tri, absorbed, laifra;
    int ib;
    if (in == 0 || out == 0 || !isfinite(in->mpe) || in->mpe <= 0.0) return 1;
    for (ib = 0; ib < 2; ++ib) {
        cad[ib] = in->solad[ib] * in->fabd[ib];
        cai[ib] = in->solai[ib] * in->fabi[ib];
        result.sav = result.sav + cad[ib] + cai[ib];
        result.fsa = result.fsa + cad[ib] + cai[ib];
        trd = in->solad[ib] * in->ftdd[ib];
        tri = in->solad[ib] * in->ftid[ib] + in->solai[ib] * in->ftii[ib];
        absorbed = trd * (1.0 - in->albgrd[ib]) + tri * (1.0 - in->albgri[ib]);
        result.sag = result.sag + absorbed;
        result.fsa = result.fsa + absorbed;
    }
    laifra = in->elai / fmax(in->vai, in->mpe);
    if (in->fsun > 0.0) {
        result.parsun = (cad[0] + in->fsun * cai[0]) * laifra / fmax(in->laisun, in->mpe);
        result.parsha = (in->fsha * cai[0]) * laifra / fmax(in->laisha, in->mpe);
    } else {
        result.parsun = 0.0;
        result.parsha = (cad[0] + cai[0]) * laifra / fmax(in->laisha, in->mpe);
    }
    result.fsr = in->albd[0] * in->solad[0] + in->albi[0] * in->solai[0]
               + in->albd[1] * in->solad[1] + in->albi[1] * in->solai[1];
    for (ib = 0; ib < 2; ++ib) {
        result.fsrv += in->frevd[ib] * in->solad[ib] + in->frevi[ib] * in->solai[ib];
        result.fsrg += in->fregd[ib] * in->solad[ib] + in->fregi[ib] * in->solai[ib];
    }
    *out = result;
    return 0;
}


/* Original: noahmp_twostream.c */
#include <math.h>
#include <stddef.h>

/* Algebra and evaluation sequence follow TWOSTREAM in
 * module_sf_noahmplsm(3).F, lines 3095-3341, for OPT_RAD=3.
 * Noah-MP REAL constants are represented here in double precision.
 */
int noahmp_twostream(const NoahmpTwostreamInput *in, NoahmpTwostreamOutput *out)
{
    int ib;
    double gap, kopen, coszi, chil, phi1, phi2, gdir, ext, avmu;
    double omegal, asu, betadl, betail, omega, betad, betai;
    double tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
    double b, c, d, f, h, sigma, p1, p2, p3, p4, s1, s2, u1, u2, u3;
    double d1, d2, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10;
    double ftds, ftis, fres, freveg, frebar, alb;
    NoahmpTwostreamOutput result;
    if (in == NULL || out == NULL || in->ib < 0 || in->ib > 1 ||
        (in->ic != 0 && in->ic != 1) || !isfinite(in->cosz) ||
        in->cosz <= 0.0 || in->cosz > 1.0 || !isfinite(in->vai) ||
        in->vai <= 0.0 || !isfinite(in->fveg) || in->fveg < 0.0 ||
        in->fveg > 1.0 || !isfinite(in->temperature_k) ||
        !isfinite(in->fwet) || in->fwet < 0.0 || in->fwet > 1.0 ||
        !isfinite(in->xl)) return 1;
    ib = in->ib;
    if (!isfinite(in->rho[ib]) || !isfinite(in->tau[ib]) ||
        !isfinite(in->albgrd[ib]) || !isfinite(in->albgri[ib]) ||
        in->albgrd[ib] <= 0.0 || in->albgri[ib] <= 0.0 ||
        in->albgrd[ib] >= 1.0 || in->albgri[ib] >= 1.0) return 1;
    result = *out;
    gap = 1.0 - in->fveg;
    kopen = 1.0 - in->fveg;
    coszi = fmax(0.001, in->cosz);
    chil = fmin(fmax(in->xl, -0.4), 0.6);
    if (fabs(chil) <= 0.01) chil = 0.01;
    phi1 = 0.5 - 0.633 * chil - 0.330 * chil * chil;
    phi2 = 0.877 * (1.0 - 2.0 * phi1);
    gdir = phi1 + phi2 * coszi;
    ext = gdir / coszi;
    avmu = (1.0 - phi1 / phi2 * log((phi1 + phi2) / phi1)) / phi2;
    omegal = in->rho[ib] + in->tau[ib];
    if (omegal <= 0.0 || !isfinite(avmu)) return 1;
    tmp0 = gdir + phi2 * coszi;
    tmp1 = phi1 * coszi;
    asu = 0.5 * omegal * gdir / tmp0 * (1.0 - tmp1 / tmp0 * log((tmp1 + tmp0) / tmp1));
    betadl = (1.0 + avmu * ext) / (omegal * avmu * ext) * asu;
    betail = 0.5 * (in->rho[ib] + in->tau[ib] +
             (in->rho[ib] - in->tau[ib]) * pow((1.0 + chil) / 2.0, 2.0)) / omegal;
    if (in->temperature_k > 273.16) {
        tmp0 = omegal;
        tmp1 = betadl;
        tmp2 = betail;
    } else {
        tmp0 = (1.0 - in->fwet) * omegal + in->fwet * in->omegas[ib];
        if (tmp0 == 0.0) return 1;
        tmp1 = ((1.0 - in->fwet) * omegal * betadl +
                in->fwet * in->omegas[ib] * in->betads) / tmp0;
        tmp2 = ((1.0 - in->fwet) * omegal * betail +
                in->fwet * in->omegas[ib] * in->betais) / tmp0;
    }
    omega = tmp0; betad = tmp1; betai = tmp2;
    b = 1.0 - omega + omega * betai;
    c = omega * betai;
    tmp0 = avmu * ext;
    d = tmp0 * omega * betad;
    f = tmp0 * omega * (1.0 - betad);
    tmp1 = b * b - c * c;
    if (tmp1 < 0.0) return 1;
    h = sqrt(tmp1) / avmu;
    sigma = tmp0 * tmp0 - tmp1;
    if (fabs(sigma) < 1.0e-6) sigma = copysign(1.0e-6, sigma);
    p1 = b + avmu * h;
    p2 = b - avmu * h;
    p3 = b + tmp0;
    p4 = b - tmp0;
    s1 = exp(-h * in->vai);
    s2 = exp(-ext * in->vai);
    if (s1 == 0.0) return 1;
    if (in->ic == 0) alb = in->albgrd[ib];
    else alb = in->albgri[ib];
    u1 = b - c / alb;
    u2 = b - c * alb;
    u3 = f + c * alb;
    tmp2 = u1 - avmu * h;
    tmp3 = u1 + avmu * h;
    d1 = p1 * tmp2 / s1 - p2 * tmp3 * s1;
    tmp4 = u2 + avmu * h;
    tmp5 = u2 - avmu * h;
    d2 = tmp4 / s1 - tmp5 * s1;
    if (d1 == 0.0 || d2 == 0.0) return 1;
    h1 = -d * p4 - c * f;
    tmp6 = d - h1 * p3 / sigma;
    tmp7 = (d - c - h1 / sigma * (u1 + tmp0)) * s2;
    h2 = (tmp6 * tmp2 / s1 - p2 * tmp7) / d1;
    h3 = -(tmp6 * tmp3 * s1 - p1 * tmp7) / d1;
    h4 = -f * p3 - c * d;
    tmp8 = h4 / sigma;
    tmp9 = (u3 - tmp8 * (u2 - tmp0)) * s2;
    h5 = -(tmp8 * tmp4 / s1 + tmp9) / d2;
    h6 = (tmp8 * tmp5 * s1 + tmp9) / d2;
    h7 = c * tmp2 / (d1 * s1);
    h8 = -c * tmp3 * s1 / d1;
    h9 = tmp4 / (d2 * s1);
    h10 = -tmp5 * s1 / d2;
    if (in->ic == 0) {
        ftds = s2 * (1.0 - gap) + gap;
        ftis = (h4 * s2 / sigma + h5 * s1 + h6 / s1) * (1.0 - gap);
        fres = (h1 / sigma + h2 + h3) * (1.0 - gap) + in->albgrd[ib] * gap;
        freveg = (h1 / sigma + h2 + h3) * (1.0 - gap);
        frebar = in->albgrd[ib] * gap;
    } else {
        ftds = 0.0;
        ftis = (h9 * s1 + h10 / s1) * (1.0 - kopen) + kopen;
        fres = (h7 + h8) * (1.0 - kopen) + in->albgri[ib] * kopen;
        freveg = fres;
        frebar = 0.0;
    }
    result.ftd[ib] = ftds;
    result.fti[ib] = ftis;
    result.fre[ib] = fres;
    result.frev[ib] = freveg;
    result.freg[ib] = frebar;
    result.fab[ib] = 1.0 - fres - (1.0 - in->albgrd[ib]) * ftds -
                     (1.0 - in->albgri[ib]) * ftis;
    result.gdir = gdir;
    if (!isfinite(result.fab[ib]) || !isfinite(result.fre[ib]) ||
        !isfinite(result.ftd[ib]) || !isfinite(result.fti[ib])) return 1;
    *out = result;
    return 0;
}


/* Original: noahmp_vege_flux.c */
#include <math.h>

int noahmp_vege_flux_foliage_step(const NoahmpVegeFluxInputs *in,
                                 NoahmpVegeFluxOutputs *out)
{
    double cah, cvh, cgh, cond, ata, bta, csh;
    double caw, cew, ctw, cgw, aea, bea, cev, ctr;
    double air_term, cir, residual, jacobian, dtv, tv, estv, tah, eah;
    double irc, shc, evc, tr, latent_limit;
    NoahmpVegeFluxOutputs result;
    if (in == 0 || out == 0) return -1;
    if (!(in->fveg > 0.0 && in->fveg <= 1.0 && in->vai >= 0.0 &&
          in->lai_sun >= 0.0 && in->lai_shade >= 0.0 &&
          in->wet_fraction >= 0.0 && in->wet_fraction <= 1.0 &&
          in->rb_s_per_m > 0.0 && in->rs_sun_s_per_m >= 0.0 &&
          in->rs_shade_s_per_m >= 0.0 && in->rahc_s_per_m > 0.0 &&
          in->rahg_s_per_m > 0.0 && in->rawc_s_per_m > 0.0 &&
          in->rawg_s_per_m > 0.0 && in->rsurf_s_per_m >= 0.0 &&
          in->air_density_kg_per_m3 > 0.0 && in->cp_air_j_per_kg_k > 0.0 &&
          in->psychrometric_constant_pa_per_k > 0.0 &&
          in->timestep_s > 0.0 && in->latent_heat_j_per_kg > 0.0 &&
          in->vegetation_emissivity >= 0.0 && in->vegetation_emissivity <= 1.0 &&
          in->ground_emissivity >= 0.0 && in->ground_emissivity <= 1.0 &&
          in->stefan_boltzmann_w_per_m2_k4 > 0.0)) return -1;
    tv = in->foliage_temperature_k;
    estv = in->foliage_saturation_vapor_pressure_pa;
    cah = 1.0 / in->rahc_s_per_m;
    cvh = 2.0 * fmin(6.0, in->vai) / in->rb_s_per_m;
    cgh = 1.0 / in->rahg_s_per_m;
    cond = cah + cvh + cgh;
    ata = (in->air_temperature_k * cah + in->ground_temperature_k * cgh) / cond;
    bta = cvh / cond;
    csh = (1.0 - bta) * in->air_density_kg_per_m3 * in->cp_air_j_per_kg_k * cvh;
    caw = 1.0 / in->rawc_s_per_m;
    cew = in->wet_fraction * fmin(6.0, in->vai) / in->rb_s_per_m;
    ctw = (1.0 - in->wet_fraction) *
          (fmin(6.0, in->lai_sun) / (in->rb_s_per_m + in->rs_sun_s_per_m) +
           fmin(6.0, in->lai_shade) / (in->rb_s_per_m + in->rs_shade_s_per_m));
    cgw = 1.0 / (in->rawg_s_per_m + in->rsurf_s_per_m);
    cond = caw + cew + ctw + cgw;
    aea = (in->air_vapor_pressure_pa * caw +
           in->ground_saturation_vapor_pressure_pa * cgw) / cond;
    bea = (cew + ctw) / cond;
    cev = (1.0 - bea) * cew * in->air_density_kg_per_m3 *
          in->cp_air_j_per_kg_k / in->psychrometric_constant_pa_per_k;
    ctr = (1.0 - bea) * ctw * in->air_density_kg_per_m3 *
          in->cp_air_j_per_kg_k / in->psychrometric_constant_pa_per_k;
    tah = ata + bta * tv;
    eah = aea + bea * estv;
    air_term = -in->vegetation_emissivity *
               (1.0 + (1.0 - in->vegetation_emissivity) *
               (1.0 - in->ground_emissivity)) * in->incoming_longwave_w_per_m2 -
               in->vegetation_emissivity * in->ground_emissivity *
               in->stefan_boltzmann_w_per_m2_k4 * pow(in->ground_temperature_k, 4.0);
    cir = (2.0 - in->vegetation_emissivity * (1.0 - in->ground_emissivity)) *
          in->vegetation_emissivity * in->stefan_boltzmann_w_per_m2_k4;
    irc = in->fveg * (air_term + cir * pow(tv, 4.0));
    shc = in->fveg * in->air_density_kg_per_m3 * in->cp_air_j_per_kg_k * cvh * (tv - tah);
    evc = in->fveg * in->air_density_kg_per_m3 * in->cp_air_j_per_kg_k * cew *
          (estv - eah) / in->psychrometric_constant_pa_per_k;
    tr = in->fveg * in->air_density_kg_per_m3 * in->cp_air_j_per_kg_k * ctw *
         (estv - eah) / in->psychrometric_constant_pa_per_k;
    if (tv > 273.16) latent_limit = in->canopy_liquid_mm * in->latent_heat_j_per_kg / in->timestep_s;
    else latent_limit = in->canopy_ice_mm * in->latent_heat_j_per_kg / in->timestep_s;
    evc = fmin(latent_limit, evc);
    residual = in->absorbed_shortwave_w_per_m2 - irc - shc - evc - tr + in->additional_heat_w_per_m2;
    jacobian = in->fveg * (4.0 * cir * pow(tv, 3.0) + csh +
               (cev + ctr) * in->foliage_saturation_slope_pa_per_k);
    if (!isfinite(jacobian) || fabs(jacobian) < 1e-15) return -1;
    dtv = residual / jacobian;
    result.energy_balance_residual_w_per_m2 = residual;
    result.unlinearized_transpiration_w_per_m2 = tr;
    result.net_longwave_w_per_m2 = irc + in->fveg * 4.0 * cir * pow(tv, 3.0) * dtv;
    result.sensible_heat_w_per_m2 = shc + in->fveg * csh * dtv;
    result.interception_evaporation_w_per_m2 = evc + in->fveg * cev *
                                                 in->foliage_saturation_slope_pa_per_k * dtv;
    result.transpiration_w_per_m2 = tr + in->fveg * ctr *
                                     in->foliage_saturation_slope_pa_per_k * dtv;
    result.canopy_air_temperature_k = tah;
    result.canopy_air_vapor_pressure_pa = eah;
    result.foliage_temperature_increment_k = dtv;
    result.updated_foliage_temperature_k = tv + dtv;
    result.transpiration_conductance_m_per_s = ctw;
    if (!isfinite(result.transpiration_w_per_m2) || !isfinite(dtv)) return -1;
    *out = result;
    return 0;
}


/* FLO 10/20: Noah-MP exact iteration mode (VEGE_FLUX loop1 replica).
 * See the description in ball_berry_helpers.h. Variable names in comments
 * are the Noah-MP names in module_sf_noahmplsm.F. */
#include <math.h>
#include <stddef.h>
#include <string.h>

int noahmp_vege_flux_exact_solve(const NoahmpStabilityInput *input,
                                 NoahmpStabilityOutput *output)
{
    const int maximum_canopy_iterations = 20;              /* NITERC */
    const int minimum_iterations_before_exit = 5;          /* ITER >= 5 */
    const double exit_temperature_change_k = 0.01;         /* ABS(DTV) <= 0.01 */
    NoahmpStabilityOutput candidate;
    NoahmpSfcdif1Inputs surface_layer_inputs;
    NoahmpSfcdif1State surface_layer_state;
    NoahmpRagrbInputs ground_canopy_inputs;
    NoahmpRagrbState ground_canopy_state;
    NoahmpAlbedoInput radiation_inputs;
    NoahmpAlbedoOutput radiation_outputs;
    NoahmpStomataInputs leaf_inputs;
    NoahmpStomataOutputs sunlit_stomata, shaded_stomata;
    NoahmpVegeFluxInputs foliage_inputs;
    NoahmpVegeFluxOutputs foliage_step;
    double vegetation_temperature_k;                   /* TV  */
    double canopy_air_temperature_k;                   /* TAH */
    double canopy_air_vapor_pressure_pa;               /* EAH */
    double ground_saturation_vapor_pressure_pa;        /* ESTG */
    double ground_saturation_slope_pa_per_k;
    double leaf_saturation_vapor_pressure_pa;          /* ESTV */
    double leaf_saturation_slope_pa_per_k;             /* DESTV */
    double canopy_to_atmosphere_sensible_heat_w_per_m2 = 0.0;  /* H  */
    double ground_to_canopy_sensible_heat_w_per_m2 = 0.0;      /* HG */
    double above_canopy_heat_resistance_s_per_m = 0.0;         /* RAHC */
    double air_density_kg_per_m3, specific_heat_air_j_per_kg_k, wind_speed_m_per_s;
    int iteration, last_iteration_flag = 0;

    if (input == NULL || output == NULL) return -1;
    memset(&candidate, 0, sizeof(candidate));
    memset(&surface_layer_state, 0, sizeof(surface_layer_state));
    memset(&ground_canopy_state, 0, sizeof(ground_canopy_state));
    surface_layer_inputs = input->surface_layer;
    ground_canopy_inputs = input->ground_canopy;
    foliage_inputs = input->canopy.foliage;
    radiation_inputs = input->canopy.radiation_stomata.radiation;
    leaf_inputs = input->canopy.radiation_stomata.canopy_air;

    air_density_kg_per_m3 = foliage_inputs.air_density_kg_per_m3;
    specific_heat_air_j_per_kg_k = foliage_inputs.cp_air_j_per_kg_k;
    wind_speed_m_per_s = surface_layer_inputs.wind_speed_m_per_s;          /* UR */
    if (!(wind_speed_m_per_s > 0.0) || !(air_density_kg_per_m3 > 0.0) ||
        !(specific_heat_air_j_per_kg_k > 0.0)) return -1;

    vegetation_temperature_k = foliage_inputs.foliage_temperature_k;
    canopy_air_temperature_k = ground_canopy_inputs.canopy_air_temperature_k;
    if (!(canopy_air_temperature_k > 0.0))
        canopy_air_temperature_k = foliage_inputs.air_temperature_k;
    canopy_air_vapor_pressure_pa = leaf_inputs.ea;

    /* Radiation: RADIATION is called by ENERGY before VEGE_FLUX, once, at the incoming TV. */
    radiation_inputs.tv_k = vegetation_temperature_k;
    if (noahmp_albedo_radiation(&radiation_inputs, &radiation_outputs) != 0) return -1;
    foliage_inputs.lai_sun = radiation_outputs.radiation_inputs.laisun;
    foliage_inputs.lai_shade = radiation_outputs.radiation_inputs.laisha;
    foliage_inputs.absorbed_shortwave_w_per_m2 = radiation_outputs.radiation.sav;

    /* ESTG at the (prescribed) ground temperature, before loop1. */
    if (noahmp_esat(foliage_inputs.ground_temperature_k, &ground_saturation_vapor_pressure_pa,
                    &ground_saturation_slope_pa_per_k) != 0) return -1;
    foliage_inputs.ground_saturation_vapor_pressure_pa = ground_saturation_vapor_pressure_pa;

    if (input->iteration_trace_stream != NULL) {
        fprintf(input->iteration_trace_stream,
                "iteration,vegetation_temperature_k,vegetation_temperature_change_k,"
                "canopy_air_temperature_k,canopy_air_vapor_pressure_pa,"
                "sunlit_stomatal_resistance_s_per_m,shaded_stomatal_resistance_s_per_m,"
                "leaf_boundary_layer_resistance_s_per_m,above_canopy_heat_resistance_s_per_m,"
                "under_canopy_heat_resistance_s_per_m,momentum_drag_coefficient,"
                "heat_exchange_coefficient,friction_velocity_m_per_s,"
                "monin_obukhov_stability_parameter,canopy_net_longwave_w_per_m2,"
                "canopy_sensible_heat_w_per_m2,transpiration_heat_w_per_m2\n");
    }

    for (iteration = 1; iteration <= maximum_canopy_iterations; ++iteration) {
        /* SFCDIF1: aerodynamic exchange between ZLVL and ZPD+Z0M, using H from the last iteration. */
        surface_layer_inputs.iteration = iteration;
        surface_layer_inputs.sensible_heat_w_per_m2 = canopy_to_atmosphere_sensible_heat_w_per_m2;
        if (noahmp_sfcdif1(&surface_layer_inputs, &surface_layer_state,
                           &candidate.surface_layer) != 0) return -1;
        above_canopy_heat_resistance_s_per_m =
            fmax(1.0, 1.0 / (candidate.surface_layer.ch * wind_speed_m_per_s));   /* RAHC */

        /* RAGRB: under-canopy resistance and leaf boundary resistance, using HG and TAH from the last iteration. */
        ground_canopy_inputs.iteration = iteration;
        ground_canopy_inputs.friction_velocity_m_per_s = candidate.surface_layer.friction_velocity_m_per_s;
        ground_canopy_inputs.ground_sensible_heat_w_per_m2 = ground_to_canopy_sensible_heat_w_per_m2;
        ground_canopy_inputs.canopy_air_temperature_k = canopy_air_temperature_k;
        if (noahmp_ragrb(&ground_canopy_inputs, &ground_canopy_state,
                         &candidate.ground_canopy) != 0) return -1;

        /* ESTV and DESTV at the current TV. */
        if (noahmp_esat(vegetation_temperature_k, &leaf_saturation_vapor_pressure_pa,
                        &leaf_saturation_slope_pa_per_k) != 0) return -1;

        /* STOMATA only on the first iteration, with the incoming TV, ESTV, EAH and this RB. */
        if (iteration == 1) {
            leaf_inputs.tv = vegetation_temperature_k;
            leaf_inputs.ei = leaf_saturation_vapor_pressure_pa;
            leaf_inputs.ea = canopy_air_vapor_pressure_pa;
            leaf_inputs.rb = candidate.ground_canopy.leaf_boundary_resistance_s_per_m;
            leaf_inputs.apar = radiation_outputs.radiation.parsun;
            if (noahmp_stomata(&input->canopy.radiation_stomata.physiology,
                                         &leaf_inputs, &sunlit_stomata) != 0) return -1;
            leaf_inputs.apar = radiation_outputs.radiation.parsha;
            if (noahmp_stomata(&input->canopy.radiation_stomata.physiology,
                                         &leaf_inputs, &shaded_stomata) != 0) return -1;
            foliage_inputs.rs_sun_s_per_m = sunlit_stomata.rs;
            foliage_inputs.rs_shade_s_per_m = shaded_stomata.rs;
        }

        /* One linearized foliage energy balance step (TAH, EAH, fluxes, DTV). */
        foliage_inputs.foliage_temperature_k = vegetation_temperature_k;
        foliage_inputs.foliage_saturation_vapor_pressure_pa = leaf_saturation_vapor_pressure_pa;
        foliage_inputs.foliage_saturation_slope_pa_per_k = leaf_saturation_slope_pa_per_k;
        foliage_inputs.rb_s_per_m = candidate.ground_canopy.leaf_boundary_resistance_s_per_m;
        foliage_inputs.rahc_s_per_m = above_canopy_heat_resistance_s_per_m;
        foliage_inputs.rawc_s_per_m = above_canopy_heat_resistance_s_per_m;        /* RAWC = RAHC */
        foliage_inputs.rahg_s_per_m = candidate.ground_canopy.ground_heat_resistance_s_per_m;
        foliage_inputs.rawg_s_per_m = candidate.ground_canopy.ground_vapor_resistance_s_per_m;
        if (noahmp_vege_flux_foliage_step(&foliage_inputs, &foliage_step) != 0) return -1;

        vegetation_temperature_k = foliage_step.updated_foliage_temperature_k;     /* TV = TV + DTV */
        canopy_air_temperature_k = foliage_step.canopy_air_temperature_k;          /* TAH (pre-update TV) */
        canopy_air_vapor_pressure_pa = foliage_step.canopy_air_vapor_pressure_pa;  /* EAH */

        /* H and HG for the next iteration's Monin-Obukhov lengths. */
        canopy_to_atmosphere_sensible_heat_w_per_m2 = air_density_kg_per_m3 * specific_heat_air_j_per_kg_k *
            (canopy_air_temperature_k - foliage_inputs.air_temperature_k) / above_canopy_heat_resistance_s_per_m;
        ground_to_canopy_sensible_heat_w_per_m2 = air_density_kg_per_m3 * specific_heat_air_j_per_kg_k *
            (foliage_inputs.ground_temperature_k - canopy_air_temperature_k) /
            candidate.ground_canopy.ground_heat_resistance_s_per_m;

        if (input->iteration_trace_stream != NULL) {
            fprintf(input->iteration_trace_stream,
                    "%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                    iteration, vegetation_temperature_k, foliage_step.foliage_temperature_increment_k,
                    canopy_air_temperature_k, canopy_air_vapor_pressure_pa,
                    foliage_inputs.rs_sun_s_per_m, foliage_inputs.rs_shade_s_per_m,
                    foliage_inputs.rb_s_per_m, above_canopy_heat_resistance_s_per_m,
                    foliage_inputs.rahg_s_per_m, candidate.surface_layer.cm, candidate.surface_layer.ch,
                    candidate.surface_layer.friction_velocity_m_per_s, surface_layer_state.moz,
                    foliage_step.net_longwave_w_per_m2, foliage_step.sensible_heat_w_per_m2,
                    foliage_step.transpiration_w_per_m2);
        }

        if (last_iteration_flag == 1) break;
        if (iteration >= minimum_iterations_before_exit &&
            fabs(foliage_step.foliage_temperature_increment_k) <= exit_temperature_change_k &&
            last_iteration_flag == 0) last_iteration_flag = 1;
    }
    if (iteration > maximum_canopy_iterations) iteration = maximum_canopy_iterations;

    /* Outputs, in the same structure as the nested solver. */
    candidate.canopy.radiation_stomata.radiation = radiation_outputs;
    candidate.canopy.radiation_stomata.sunlit = sunlit_stomata;
    candidate.canopy.radiation_stomata.shaded = shaded_stomata;
    candidate.canopy.radiation_stomata.photosynthesis_umol_m2_s =
        sunlit_stomata.psn * radiation_outputs.radiation_inputs.laisun +
        shaded_stomata.psn * radiation_outputs.radiation_inputs.laisha;
    candidate.canopy.foliage.foliage = foliage_step;
    candidate.canopy.foliage.foliage.canopy_air_temperature_k = canopy_air_temperature_k;
    candidate.canopy.foliage.foliage.canopy_air_vapor_pressure_pa = canopy_air_vapor_pressure_pa;
    candidate.canopy.foliage.foliage_temperature_k = vegetation_temperature_k;
    candidate.canopy.foliage.ground_skin_temperature_k = foliage_inputs.ground_temperature_k;
    candidate.canopy.foliage.transpiration_m_per_s =
        foliage_step.transpiration_w_per_m2 / foliage_inputs.latent_heat_j_per_kg / 1000.0;
    candidate.canopy.iterations = iteration;
    candidate.canopy.frozen_sunlit_resistance_s_per_m = sunlit_stomata.rs;
    candidate.canopy.frozen_shaded_resistance_s_per_m = shaded_stomata.rs;
    candidate.canopy.final_temperature_increment_k = foliage_step.foliage_temperature_increment_k;
    candidate.canopy.energy_balance_residual_w_per_m2 = foliage_step.energy_balance_residual_w_per_m2;
    candidate.stability_iterations = iteration;
    candidate.relative_resistance_change = 0.0;
    candidate.reference_canopy_resistance_s_per_m = above_canopy_heat_resistance_s_per_m;
    if (!isfinite(candidate.canopy.foliage.transpiration_m_per_s)) return -1;
    *output = candidate;
    return 0;
}
