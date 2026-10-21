#ifndef BALL_BERRY_HELPERS_H
#define BALL_BERRY_HELPERS_H

/* Consolidated Noah-MP canopy physics API for PET method 6. */
#include <stdio.h>
#include "noahmp_stomata.h"


/* Original: noahmp_ragrb.h */
/* Noah-MP module_sf_noahmplsm.F: VEGE_FLUX canopy wind and RAGRB. */
typedef struct {
    int iteration;                       // Noah-MP stability iteration [1..]
    double reference_wind_speed_m_per_s; // UR [m/s]
    double reference_height_m;           // ZLVL [m]
    double momentum_roughness_m;         // Z0M [m]
    double displacement_height_m;        // ZPD [m]
    double canopy_height_m;              // HCAN [m]
    double vegetation_area_index;        // VAIE [-]
    double air_density_kg_per_m3;        // RHOAIR [kg/m3]
    double ground_sensible_heat_w_per_m2;// HG [W/m2]
    double canopy_air_temperature_k;     // TAH [K]
    double ground_momentum_roughness_m;  // Z0MG [m]
    double ground_heat_roughness_m;      // Z0HG [m]
    double canopy_heat_roughness_m;      // Z0H [m]
    double friction_velocity_m_per_s;    // FV [m/s]
    double cwp;                          // CWP [-]
    double leaf_dimension_m;             // DLEAF [m]
    double mpe;                          // MPE [-]
} NoahmpRagrbInputs;
typedef struct {
    double mozg;                         // Ground stability parameter [-]
    double fhg;                          // Ground stability correction [-]
} NoahmpRagrbState;
typedef struct {
    double canopy_wind_speed_m_per_s;    // UC [m/s]
    double ground_momentum_resistance_s_per_m; // RAMG [s/m]
    double ground_heat_resistance_s_per_m; // RAHG [s/m]
    double ground_vapor_resistance_s_per_m;// RAWG [s/m]
    double leaf_boundary_resistance_s_per_m; // RB [s/m]
} NoahmpRagrbOutputs;
int noahmp_ragrb(const NoahmpRagrbInputs *inputs,
                 NoahmpRagrbState *state, NoahmpRagrbOutputs *outputs);


/* Original: noahmp_twostream.h */
/* Direct C99 translation of Noah-MP TWOSTREAM, OPT_RAD=3.
 * Array index 0=visible, 1=near-infrared; ic=0 direct, 1 diffuse.
 * Call once for each waveband and incoming-radiation type.
 */
typedef struct {
    int ib; // Spectral band: 0 or 1.
    int ic; // Incident radiation: 0 direct, 1 diffuse.
    double cosz; // Cosine of solar zenith angle.
    double vai; // Leaf plus stem area index.
    double fwet; // Wet canopy fraction.
    double temperature_k; // Canopy surface temperature (K).
    double albgrd[2]; // Direct ground albedo by band.
    double albgri[2]; // Diffuse ground albedo by band.
    double rho[2]; // Leaf and stem reflectance by band.
    double tau[2]; // Leaf and stem transmittance by band.
    double fveg; // Vegetation fraction.
    double xl; // Noah-MP leaf angle distribution parameter.
    double omegas[2]; // Snow scattering albedo by band.
    double betads; // Snow direct upscatter parameter.
    double betais; // Snow diffuse upscatter parameter.
} NoahmpTwostreamInput;
typedef struct {
    double fab[2]; // Absorbed fraction.
    double fre[2]; // Reflected fraction.
    double ftd[2]; // Downward direct transmitted fraction.
    double fti[2]; // Downward diffuse transmitted fraction.
    double frev[2]; // Vegetation-reflected fraction.
    double freg[2]; // Ground-reflected fraction.
    double gdir; // Projected leaf and stem area.
    double bgap; // Between-crown gap; OPT_RAD=3 leaves unchanged.
    double wgap; // Within-crown gap; OPT_RAD=3 leaves unchanged.
} NoahmpTwostreamOutput;
/* Return 0 on success, nonzero on invalid or singular input. */
int noahmp_twostream(const NoahmpTwostreamInput *in, NoahmpTwostreamOutput *out);


/* Original: noahmp_surrad.h */
/* Direct C99 port of SURRAD, module_sf_noahmplsm(3).F, lines 2746-2867.
 * Two wavebands: index 0 visible, index 1 near infrared.
 * All radiation fluxes are W/m2. */
typedef struct {
    double fsun, fsha, elai, vai, laisun, laisha, mpe;
    double solad[2], solai[2], fabd[2], fabi[2];
    double ftdd[2], ftid[2], ftii[2];
    double albgrd[2], albgri[2], albd[2], albi[2];
    double frevd[2], frevi[2], fregd[2], fregi[2];
} NoahmpSurradInput;
typedef struct {
    double parsun, parsha, sav, sag, fsa, fsr, fsrv, fsrg;
} NoahmpSurradOutput;
/* Return 0 on success, nonzero on invalid pointers or MPE. */
int noahmp_surrad(const NoahmpSurradInput *in, NoahmpSurradOutput *out);


/* Original: noahmp_albedo.h */
/* ALBEDO optical section + RADIATION orchestration, OPT_RAD=3.
 * Ground albedos are prescribed until SNOW_AGE, SNOWALB and GROUNDALB
 * are ported. All arrays index 0=visible, 1=near infrared. */
typedef struct {
    double cosz; // Cosine of solar zenith angle.
    double elai; // Exposed leaf area index.
    double esai; // Exposed stem area index.
    double fveg; // Vegetated fraction.
    double fwet; // Wet canopy fraction.
    double tv_k; // Foliage temperature (K).
    double rhol[2]; // Leaf reflectance.
    double rhos[2]; // Stem reflectance.
    double taul[2]; // Leaf transmittance.
    double taus[2]; // Stem transmittance.
    double albgrd[2]; // Ground direct albedo (prescribed).
    double albgri[2]; // Ground diffuse albedo (prescribed).
    double xl; // Leaf angle parameter.
    double omegas[2]; // Snow optical scattering.
    double betads; // Snow direct upscatter.
    double betais; // Snow diffuse upscatter.
    double solad[2]; // Incoming direct radiation (W/m2).
    double solai[2]; // Incoming diffuse radiation (W/m2).
} NoahmpAlbedoInput;
typedef struct {
    double fsun; // Sunlit canopy fraction.
    NoahmpSurradInput radiation_inputs; // SURRAD inputs after ALBEDO.
    NoahmpSurradOutput radiation; // Absorbed and reflected fluxes.
} NoahmpAlbedoOutput;
/* Returns 0 on success, nonzero on invalid inputs or singular optics. */
int noahmp_albedo_radiation(const NoahmpAlbedoInput *in, NoahmpAlbedoOutput *out);


/* Original: noahmp_radiation_stomata.h */
/* This is an interface increment, not a VEGE_FLUX implementation.
 * Radiation comes from ALBEDO -> TWOSTREAM -> SURRAD.
 * Canopy air vapor pressure, foliage temperature, boundary-layer
 * resistance, phenology and leaf nitrogen must be supplied by the caller.
 * Potential stomatal calculations fix BTRAN=1.
 */
typedef struct {
    NoahmpAlbedoInput radiation;
    NoahmpStomataParameters physiology;
    NoahmpStomataInputs canopy_air;
} NoahmpRadiationStomataInput;
typedef struct {
    NoahmpAlbedoOutput radiation;
    NoahmpStomataOutputs sunlit;
    NoahmpStomataOutputs shaded;
    double photosynthesis_umol_m2_s;
} NoahmpRadiationStomataOutput;
int noahmp_radiation_stomata(const NoahmpRadiationStomataInput *in,
                            NoahmpRadiationStomataOutput *out);


/* Original: noahmp_leaf_boundary_resistance.h */

/* Noah-MP RAGRB CWP branch, isolated from upstream aerodynamic calculations.
 * All inputs must be computed or supplied by the calling canopy model.
 */
typedef struct {
    double cwp;                  // Canopy wind attenuation coefficient [-]
    double vegetation_area_index; // LAI + SAI [m2/m2]
    double canopy_height_m;      // Canopy height [m]
    double fhg;                  // Ground aerodynamic stability factor [-]
    double leaf_dimension_m;     // Noah-MP DLEAF [m]
    double canopy_wind_speed_m_per_s; // Noah-MP UC [m/s]
} NoahmpLeafBoundaryInputs;

/* Returns 0 on success, -1 on invalid/nonfinite inputs. Output in s/m. */
int noahmp_leaf_boundary_resistance(const NoahmpLeafBoundaryInputs *inputs,
                                   double *resistance_s_per_m);


/* Original: noahmp_vege_flux.h */
/* VEGE_FLUX foliage energy-balance substep. All fluxes W/m2. */
typedef struct {
    double fveg;                 // Vegetated fraction [-]
    double vai;                  // Effective leaf plus stem area index [-]
    double lai_sun;              // Effective sunlit LAI [-]
    double lai_shade;            // Effective shaded LAI [-]
    double wet_fraction;         // Wet canopy fraction [-]
    double rb_s_per_m;           // Leaf boundary-layer resistance [s/m]
    double rs_sun_s_per_m;       // Sunlit stomatal resistance [s/m]
    double rs_shade_s_per_m;     // Shaded stomatal resistance [s/m]
    double rahc_s_per_m;         // Reference-to-canopy sensible resistance [s/m]
    double rahg_s_per_m;         // Ground-to-canopy sensible resistance [s/m]
    double rawc_s_per_m;         // Reference-to-canopy vapor resistance [s/m]
    double rawg_s_per_m;         // Ground-to-canopy vapor resistance [s/m]
    double rsurf_s_per_m;        // Ground surface vapor resistance [s/m]
    double air_temperature_k;    // Reference air temperature [K]
    double ground_temperature_k; // Ground temperature [K]
    double foliage_temperature_k;// Foliage temperature [K]
    double air_vapor_pressure_pa;// Reference vapor pressure [Pa]
    double ground_saturation_vapor_pressure_pa; // Ground saturation vapor pressure [Pa]
    double foliage_saturation_vapor_pressure_pa; // Foliage saturation vapor pressure [Pa]
    double foliage_saturation_slope_pa_per_k; // d(es)/dT [Pa/K]
    double air_density_kg_per_m3; // Air density [kg/m3]
    double cp_air_j_per_kg_k;   // Air heat capacity [J/kg/K]
    double psychrometric_constant_pa_per_k; // Vegetation gamma [Pa/K]
    double incoming_longwave_w_per_m2; // Incoming LW [W/m2]
    double absorbed_shortwave_w_per_m2; // SAV [W/m2]
    double vegetation_emissivity; // EMV [-]
    double ground_emissivity;     // EMG [-]
    double stefan_boltzmann_w_per_m2_k4; // SB [W/m2/K4]
    double canopy_liquid_mm;    // Intercepted liquid water [mm]
    double canopy_ice_mm;       // Intercepted ice water [mm]
    double latent_heat_j_per_kg;// LATHEAV [J/kg]
    double timestep_s;          // DT [s]
    double additional_heat_w_per_m2; // PAHV [W/m2]
} NoahmpVegeFluxInputs;
typedef struct {
    double energy_balance_residual_w_per_m2; // Exact unlinearized foliage balance [W/m2]
    double unlinearized_transpiration_w_per_m2; // Flux at input TV [W/m2]
    double canopy_air_temperature_k; // TAH [K]
    double canopy_air_vapor_pressure_pa; // EAH [Pa]
    double transpiration_w_per_m2; // TR [W/m2]
    double interception_evaporation_w_per_m2; // EVC [W/m2]
    double sensible_heat_w_per_m2; // SHC [W/m2]
    double net_longwave_w_per_m2; // IRC [W/m2]
    double foliage_temperature_increment_k; // DTV [K]
    double updated_foliage_temperature_k; // TV + DTV [K]
    double transpiration_conductance_m_per_s; // CTW [m/s]
} NoahmpVegeFluxOutputs;
/* A single foliage Newton update from VEGE_FLUX, not the full coupled solver. */
int noahmp_vege_flux_foliage_step(const NoahmpVegeFluxInputs *in,
                                 NoahmpVegeFluxOutputs *out);


/* Original: noahmp_prescribed_tg.h */
/* Prescribed ground skin temperature (TG) driver for the existing
 * Noah-MP foliage Newton substep. Aerodynamic resistances and stomatal
 * resistances are held fixed by this driver. This is not full VEGE_FLUX. */
typedef struct {
    int iterations;                    // Number of foliage Newton updates [-]
    double maximum_temperature_step_k; // Safety limit for Newton updates [K]
} NoahmpPrescribedTgOptions;
typedef struct {
    NoahmpVegeFluxOutputs foliage;    // Final foliage fluxes and canopy-air state
    double foliage_temperature_k;     // Final foliage temperature [K]
    double ground_skin_temperature_k; // Prescribed ground skin temperature [K]
    double transpiration_m_per_s; // Potential transpiration [m/s]
} NoahmpPrescribedTgOutputs;
int noahmp_esat(double temperature_k, double *esat_pa, double *slope_pa_per_k);
int noahmp_prescribed_tg_solve(const NoahmpVegeFluxInputs *inputs,
                              const NoahmpPrescribedTgOptions *options,
                              NoahmpPrescribedTgOutputs *outputs);


/* Original: noahmp_canopy_iterative.h */
/* Prescribed-TG canopy iteration. Aerodynamic resistances are supplied
 * by caller and remain fixed; this is not complete VEGE_FLUX. */
typedef enum {
    NOAHMP_STOMATA_REEVALUATE = 0, // Recalculate Ball-Berry at each trial temperature
    NOAHMP_STOMATA_REFERENCE_FROZEN = 1 // Freeze at initial foliage temperature
} NoahmpStomataUpdateMode;
typedef struct {
    NoahmpRadiationStomataInput radiation_stomata;
    NoahmpVegeFluxInputs foliage;
    NoahmpStomataUpdateMode stomata_update_mode;
    int use_prescribed_frozen_stomata; // Outer stability driver supplies reference first-iteration RSS
    double prescribed_frozen_sunlit_resistance_s_per_m; // s/m
    double prescribed_frozen_shaded_resistance_s_per_m; // s/m
    int maximum_iterations;
    double temperature_tolerance_k;
    double maximum_temperature_step_k;
} NoahmpCanopyIterativeInput;
typedef struct {
    NoahmpRadiationStomataOutput radiation_stomata;
    NoahmpPrescribedTgOutputs foliage;
    int iterations;
    double frozen_sunlit_resistance_s_per_m; // Reference mode only
    double frozen_shaded_resistance_s_per_m; // Reference mode only
    double final_temperature_increment_k;
    double energy_balance_residual_w_per_m2;
} NoahmpCanopyIterativeOutput;
/* Pure residual evaluation: no state mutation or iteration. */
int noahmp_canopy_evaluate(const NoahmpCanopyIterativeInput *in,
                           double foliage_temperature_k,
                           double *residual_w_per_m2,
                           double *transpiration_m_per_s,
                           NoahmpRadiationStomataOutput *radiation_out,
                           NoahmpVegeFluxOutputs *foliage_out);
/* Implicit derivative with respect to foliage additional heat [W/m2].
 * dTV/dQ = -1/(dF/dTV); dTR/dQ = dTR/dTV * dTV/dQ.
 * The temperature partials use centered finite differences for now.
 * This interface can later accept Tapenade-computed partial derivatives. */
int noahmp_canopy_heat_sensitivity(const NoahmpCanopyIterativeInput *in,
                                   double converged_temperature_k,
                                   double *d_temperature_d_heat_k_m2_per_w,
                                   double *d_transpiration_d_heat_m3_per_j);
int noahmp_canopy_iterative_solve(const NoahmpCanopyIterativeInput *in,
                                  NoahmpCanopyIterativeOutput *out);


/* Original: noahmp_sfcdif1.h */
/* Direct C99 translation of Noah-MP SFCDIF1, with explicit iteration state. */
typedef struct {
    double moz, fm, fh, fm2, fh2;
    double friction_velocity_m_per_s;
    int moz_sign_changes;
} NoahmpSfcdif1State;
typedef struct {
    int iteration;
    double air_temperature_k, air_density_kg_per_m3;
    double sensible_heat_w_per_m2, specific_humidity_kg_per_kg;
    double reference_height_m, displacement_height_m;
    double momentum_roughness_m, heat_roughness_m;
    double wind_speed_m_per_s, mpe;
} NoahmpSfcdif1Inputs;
typedef struct {
    double cm, ch, friction_velocity_m_per_s, ch2;
} NoahmpSfcdif1Outputs;
int noahmp_sfcdif1(const NoahmpSfcdif1Inputs *in,
                   NoahmpSfcdif1State *state, NoahmpSfcdif1Outputs *out);


/* Original: noahmp_stability_coupling.h */
/* Outer Picard coupling: SFCDIF1/RAGRB resistances updated between
 * converged scalar foliage solves. No aerodynamic terms in implicit Jacobian.
 * This is a prescribed-TG approximation, NOT full Noah-MP VEGE_FLUX.
 */
typedef struct {
    NoahmpCanopyIterativeInput canopy;
    NoahmpSfcdif1Inputs surface_layer;
    NoahmpRagrbInputs ground_canopy;
    int maximum_stability_iterations;
    double resistance_relative_tolerance;
    FILE *iteration_trace_stream;        // Optional (NULL = off): exact mode writes one CSV row per iteration
} NoahmpStabilityInput;
typedef struct {
    NoahmpCanopyIterativeOutput canopy;
    NoahmpSfcdif1Outputs surface_layer;
    NoahmpRagrbOutputs ground_canopy;
    int stability_iterations;
    double relative_resistance_change;
    double reference_canopy_resistance_s_per_m;
} NoahmpStabilityOutput;
int noahmp_stability_coupled_solve(const NoahmpStabilityInput *input,
                                  NoahmpStabilityOutput *output);

/* FLO 10/20: Noah-MP exact iteration mode. Replicates VEGE_FLUX loop1 of the
 * NWM module_sf_noahmplsm.F with the ground temperature TG prescribed:
 *   - radiation (ALBEDO/TWOSTREAM/SURRAD) evaluated once, at the initial TV;
 *   - per iteration: SFCDIF1 (with H from the previous iteration), RAGRB
 *     (with HG and TAH from the previous iteration), ESAT at TV, and ONE
 *     linearized foliage Newton step (noahmp_vege_flux_foliage_step);
 *   - STOMATA evaluated only on iteration 1, with the initial TV, ESAT(TV),
 *     the initial canopy air vapor pressure (canopy_air.ea) and iteration 1 RB;
 *   - RAHC = RAWC = MAX(1, 1/(CH*UR));
 *   - at most NITERC = 20 iterations; once ITER >= 5 and |DTV| <= 0.01 K,
 *     one more iteration is taken and the loop exits;
 *   - reported TR is the linearized TR + CTR*DESTV*DTV of the last iteration.
 * Inputs: the same NoahmpStabilityInput as noahmp_stability_coupled_solve.
 * ground_canopy.canopy_air_temperature_k is the initial TAH (if <= 0, the air
 * temperature is used). The nested-solver controls (tolerances, maximum
 * iterations, stomata_update_mode) are ignored. */
int noahmp_vege_flux_exact_solve(const NoahmpStabilityInput *input,
                                 NoahmpStabilityOutput *output);


/* Original: noahmp_canopy_aerodynamics.h */

/* Method 6 only. UC and FHG remain explicit Noah-MP upstream inputs.
 * This interface does not approximate their computation.
 */
typedef struct {
    NoahmpSfcdif1Inputs surface_layer;
    NoahmpLeafBoundaryInputs leaf_boundary;
} NoahmpCanopyAerodynamicsInputs;

typedef struct {
    NoahmpSfcdif1Outputs surface_layer;
    double leaf_boundary_resistance_s_per_m;
} NoahmpCanopyAerodynamicsOutputs;

int noahmp_canopy_aerodynamics(const NoahmpCanopyAerodynamicsInputs *inputs,
                              NoahmpSfcdif1State *state,
                              NoahmpCanopyAerodynamicsOutputs *outputs);


#endif /* BALL_BERRY_HELPERS_H */
