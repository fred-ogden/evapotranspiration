#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef PET_C
#define PET_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

//local includes
#include "../include/pet.h"
#include "../include/pet_tools.h"
#include "../include/PEtEnergyBalanceMethod.h"
#include "../include/PEtAerodynamicMethod.h"
#include "../include/PEtCombinationMethod.h"
#include "../include/PEtPriestleyTaylorMethod.h"
#include "../include/PEtPenmanMonteithMethod.h"

extern void alloc_pet_model(pet_model *model) {
    // TODO: *******************
}

extern void free_pet_model(pet_model *model) {
    // TODO: ******************
}

// ######################    RUN    ########    RUN    ########    RUN    ########    RUN    #################################
// ######################    RUN    ########    RUN    ########    RUN    ########    RUN    #################################
// ######################    RUN    ########    RUN    ########    RUN    ########    RUN    #################################
// ######################    RUN    ########    RUN    ########    RUN    ########    RUN    #################################
/* FLO 10/21: method 6 diagnostic summary, one "name,value" row per quantity, Noah-MP names.
 * Written from the model's own state, so a test driver needs only the BMI. */
static int bb6_write_diagnostic_summary(pet_model *model)
{
  FILE *summary_file = fopen(model->bb_diagnostic_summary_file, "w");
  const NoahmpStabilityOutput *o = &model->bb_last_output;
  const NoahmpSurradInput *ri = &o->canopy.radiation_stomata.radiation.radiation_inputs;
  const NoahmpSurradOutput *ro = &o->canopy.radiation_stomata.radiation.radiation;
  const BallBerryCouplingInputs *ci = &model->bb_coupling_inputs;
  const BallBerryCouplingOutputs *co = &model->bb_coupling_outputs;
  if (summary_file == NULL) {
    fprintf(stderr, "ERROR: cannot open bb_diagnostic_summary_file %s\n", model->bb_diagnostic_summary_file);
    return -1;
  }
  fprintf(summary_file, "noahmp_name,value\n");
  fprintf(summary_file, "SOLVER,%s\n", model->bb_noahmp_exact_iteration ? "noahmp_exact_iteration" : "nested_picard_newton");
  fprintf(summary_file, "ITERATIONS,%d\n", o->stability_iterations);
  fprintf(summary_file, "BTRAN,%.16e\n", ci->root_stress_factor_0_to_1);
  fprintf(summary_file, "RSURF,%.16e\n", ci->surface_resistance_s_per_m);
  fprintf(summary_file, "FSNO,%.16e\n", ci->snow_cover_fraction);
  fprintf(summary_file, "SNOWH,%.16e\n", ci->snow_depth_m);
  fprintf(summary_file, "FWET,%.16e\n", ci->canopy_wetted_fraction);
  fprintf(summary_file, "LAI,%.16e\n", model->bb_canopy.lai);
  fprintf(summary_file, "SAI,%.16e\n", model->bb_canopy.sai);
  fprintf(summary_file, "ELAI,%.16e\n", ri->elai);
  fprintf(summary_file, "ESAI,%.16e\n", ri->vai - ri->elai);
  fprintf(summary_file, "FSUN,%.16e\n", ri->fsun);
  fprintf(summary_file, "LAISUN,%.16e\n", ri->laisun);
  fprintf(summary_file, "LAISHA,%.16e\n", ri->laisha);
  fprintf(summary_file, "PARSUN,%.16e\n", ro->parsun);
  fprintf(summary_file, "PARSHA,%.16e\n", ro->parsha);
  fprintf(summary_file, "SAV,%.16e\n", ro->sav);
  fprintf(summary_file, "SAG,%.16e\n", ro->sag);
  fprintf(summary_file, "FSR,%.16e\n", ro->fsr);
  fprintf(summary_file, "RSSUN,%.16e\n", o->canopy.radiation_stomata.sunlit.rs);
  fprintf(summary_file, "RSSHA,%.16e\n", o->canopy.radiation_stomata.shaded.rs);
  fprintf(summary_file, "PSNSUN,%.16e\n", o->canopy.radiation_stomata.sunlit.psn);
  fprintf(summary_file, "PSNSHA,%.16e\n", o->canopy.radiation_stomata.shaded.psn);
  fprintf(summary_file, "TV,%.16e\n", co->vegetation_temperature_K);
  fprintf(summary_file, "TAH,%.16e\n", co->canopy_air_temperature_K);
  fprintf(summary_file, "EAH,%.16e\n", co->canopy_air_vapor_pressure_Pa);
  fprintf(summary_file, "TG,%.16e\n", o->canopy.foliage.ground_skin_temperature_k);
  fprintf(summary_file, "CM,%.16e\n", o->surface_layer.cm);
  fprintf(summary_file, "IRC,%.16e\n", o->canopy.foliage.foliage.net_longwave_w_per_m2);
  fprintf(summary_file, "SHC,%.16e\n", o->canopy.foliage.foliage.sensible_heat_w_per_m2);
  fprintf(summary_file, "EVC,%.16e\n", co->canopy_evaporation_heat_flux_W_per_m2);
  fprintf(summary_file, "TR,%.16e\n", o->canopy.foliage.foliage.transpiration_w_per_m2);
  fprintf(summary_file, "ETRAN_MM_PER_H,%.16e\n", co->canopy_transpiration_flux_m_per_s * 1000.0 * 3600.0);
  fprintf(summary_file, "RAHG,%.16e\n", co->below_canopy_aerodynamic_resistance_s_per_m);
  fprintf(summary_file, "GROUND_LW_ABSORBED,%.16e\n", co->ground_longwave_absorbed_below_canopy_W_per_m2);
  fclose(summary_file);
  return 0;
}

extern int run_pet(pet_model* model)
{
  if (model == NULL) return -1;
  if (model->bmi.verbose >2){
    printf("Running the PET model \n");
    printf("model->bmi.is_forcing_from_bmi %d \n", model->bmi.is_forcing_from_bmi);
  }
  if (model->pet_method == 6) {
    NoahmpStabilityInput input;
    NoahmpStabilityOutput output;
    /* FLO 10/21: build this step's canopy input from the BMI forcing, time and coupling inputs
     * and the parameters read at Initialize (ball_berry_method6_inputs.c). */
    if (bb6_build_canopy_input(model) != 0 || model->bb_stability_input_ready != TRUE) {
      fprintf(stderr, "ERROR: PET method 6 could not build the canopy input for this step.\n");
      return -1;
    }
    input = model->bb_stability_input;
    if (model->bb_noahmp_compatibility == TRUE) {
      input.canopy.stomata_update_mode = NOAHMP_STOMATA_REFERENCE_FROZEN;
    } else {
      input.canopy.stomata_update_mode = NOAHMP_STOMATA_REEVALUATE;
    }
    model->bb_stability_input_ready = FALSE; // Consume forcing even if the solver fails.
    // FLO 10/20: optional per-iteration trace (written by the exact solver only).
    input.iteration_trace_stream = NULL;
    if (model->bb_iteration_trace_file != NULL && model->bb_noahmp_exact_iteration == TRUE) {
      input.iteration_trace_stream = fopen(model->bb_iteration_trace_file, "w");
      if (input.iteration_trace_stream == NULL) {
        fprintf(stderr, "ERROR: cannot open bb_iteration_trace_file %s\n", model->bb_iteration_trace_file);
        return -1;
      }
    }
    int solver_status;
    if (model->bb_noahmp_exact_iteration == TRUE) {
      solver_status = noahmp_vege_flux_exact_solve(&input, &output);   // VEGE_FLUX loop1 replica
    } else {
      solver_status = noahmp_stability_coupled_solve(&input, &output); // nested Picard/Newton solver
    }
    if (input.iteration_trace_stream != NULL) fclose(input.iteration_trace_stream);
    if (solver_status != 0) return -1;
    model->bb_last_output = output;

    // FLO 10/21: BMI coupling outputs. Noah-MP names on the right; see ball_berry_state.h for
    // the area basis of each. The longwave absorbed by the under-canopy ground is the absorbed
    // part of Noah-MP's IRG (VEGE_FLUX loop2), evaluated at the final vegetation temperature:
    //   IRG = EMG*SB*TG**4 - EMG*(1-EMV)*LWDN - EMG*EMV*SB*TV**4   (+ = to atmosphere)
    {
      const NoahmpVegeFluxInputs *canopy_inputs = &input.canopy.foliage;
      const NoahmpVegeFluxOutputs *canopy_fluxes = &output.canopy.foliage.foliage;
      BallBerryCouplingOutputs *coupling_outputs = &model->bb_coupling_outputs;
      double TV   = output.canopy.foliage.foliage_temperature_k;                 // vegetation temperature
      double EMV  = canopy_inputs->vegetation_emissivity;                        // vegetation emissivity
      double EMG  = canopy_inputs->ground_emissivity;                            // ground emissivity
      double LWDN = canopy_inputs->incoming_longwave_w_per_m2;                   // downward longwave
      double stefan_boltzmann_w_per_m2_k4 = canopy_inputs->stefan_boltzmann_w_per_m2_k4; // Noah-MP SB (pet_tools.h #defines SB)
      coupling_outputs->canopy_transpiration_flux_m_per_s = output.canopy.foliage.transpiration_m_per_s;          // ETRAN
      coupling_outputs->canopy_evaporation_heat_flux_W_per_m2 = canopy_fluxes->interception_evaporation_w_per_m2; // EVC
      coupling_outputs->ground_shortwave_absorbed_below_canopy_W_per_m2 =
          output.canopy.radiation_stomata.radiation.radiation.sag;                                               // SAG
      coupling_outputs->ground_longwave_absorbed_below_canopy_W_per_m2 =
          EMG * (1.0 - EMV) * LWDN + EMG * EMV * stefan_boltzmann_w_per_m2_k4 * pow(TV, 4.0);
      coupling_outputs->canopy_air_temperature_K = canopy_fluxes->canopy_air_temperature_k;                     // TAH
      coupling_outputs->canopy_air_vapor_pressure_Pa = canopy_fluxes->canopy_air_vapor_pressure_pa;              // EAH
      coupling_outputs->below_canopy_aerodynamic_resistance_s_per_m =
          output.ground_canopy.ground_heat_resistance_s_per_m;                                                   // RAHG
      coupling_outputs->vegetation_temperature_K = TV;                                                           // TV
      coupling_outputs->canopy_latent_heat_J_per_kg = canopy_inputs->latent_heat_j_per_kg;                       // LATHEAV (FLO 10/24)
    }

    // FLO 10/21: optional diagnostic summary of this update (Noah-MP names), for comparison with
    // the Fortran NWM Noah-MP reference. Overwritten at each update: it holds the last step.
    if (model->bb_diagnostic_summary_file != NULL && bb6_write_diagnostic_summary(model) != 0) return -1;
    model->bb_transpiration_m_per_s = output.canopy.foliage.transpiration_m_per_s;
    model->pet_m_per_s = model->bb_transpiration_m_per_s;
    if (model->pet_m_per_s < 0.0) model->pet_m_per_s = 0.0;
    return 0;
  }
  if (model->pet_method < 1 || model->pet_method > 5) return -1;
  if (model->pet_params.zero_plane_displacement_height_m <= 0.0){
    fprintf(stderr, "ERROR: zero_plane_displacement_height_m must be > 0.0 m. Current value: %lf\n", model->pet_params.zero_plane_displacement_height_m);
    return -1;
  }
  // populate the evapotranspiration forcing data structure:
  //---------------------------------------------------------------------------------------------------------------
  /*
      jmframe: I think it would be better down below the setting of the model->aorc.forcings
               That way we don't have to index the larger arrays twice...
               So we would delete the first block in this "if" statement,
               And move the "else" section below the model->aorc.forcings setting block.
  */
  if (model->bmi.is_forcing_from_bmi == 0){
    model->pet_forcing.air_temperature_C = model->forcing_data_air_temperature_2m_K[model->bmi.current_step] - TK;//convert to C
    model->pet_forcing.relative_humidity_percent     = (double)-99.9; // this negative number means use specific humidity
    /* FLO 10/26: The forcing-file path previously used precipitation as
     * specific humidity. Read the SPFH field instead (kg/kg).
     */
    model->pet_forcing.specific_humidity_2m_kg_per_kg = model->forcing_data_specific_humidity_2m_kg_per_kg[model->bmi.current_step];
    model->pet_forcing.air_pressure_Pa    = model->forcing_data_surface_pressure_Pa[model->bmi.current_step];
    model->pet_forcing.wind_speed_m_per_s = hypot(model->forcing_data_u_wind_speed_10m_m_per_s[model->bmi.current_step],
                                           model->forcing_data_v_wind_speed_10m_m_per_s[model->bmi.current_step]);                 
  }
  else{
    model->pet_forcing.air_temperature_C = model->aorc.air_temperature_2m_K - TK;//convert to C
    model->pet_forcing.relative_humidity_percent     = (double)-99.9; // this negative number means use specific humidity
    model->pet_forcing.specific_humidity_2m_kg_per_kg = model->aorc.specific_humidity_2m_kg_per_kg;
    model->pet_forcing.air_pressure_Pa    = model->aorc.surface_pressure_Pa;
    model->pet_forcing.wind_speed_m_per_s = hypot(model->aorc.u_wind_speed_10m_m_per_s, model->aorc.v_wind_speed_10m_m_per_s);                 
  }

  if(model->pet_options.yes_aorc==1)
  {
    if (model->bmi.verbose >1)
        printf("YES AORC \n");
    // Validate AORC humidity measurement height
    if(model->pet_params.humidity_measurement_height_m !=2.0)
    {
        fprintf(stderr, "ERROR: humidity measurement height is not 2.0 m. Humidity adjustment not yet implemented. Current value: %lf\n", model->pet_params.humidity_measurement_height_m);
        return -1;
    }
    /* jmframe: If we are getting forcing through BMI, then we don't need this, the forcings should already be in place */
    if (model->bmi.is_forcing_from_bmi == 0){
      model->aorc.incoming_longwave_W_per_m2     =  model->forcing_data_incoming_longwave_W_per_m2[model->bmi.current_step];
      model->aorc.incoming_shortwave_W_per_m2    =  model->forcing_data_incoming_shortwave_W_per_m2[model->bmi.current_step];
      model->aorc.surface_pressure_Pa            =  model->forcing_data_surface_pressure_Pa[model->bmi.current_step];
      model->aorc.specific_humidity_2m_kg_per_kg =  model->forcing_data_specific_humidity_2m_kg_per_kg[model->bmi.current_step];
      model->aorc.air_temperature_2m_K           =  model->forcing_data_air_temperature_2m_K[model->bmi.current_step];
      model->aorc.u_wind_speed_10m_m_per_s       =  model->forcing_data_u_wind_speed_10m_m_per_s[model->bmi.current_step];
      model->aorc.v_wind_speed_10m_m_per_s       =  model->forcing_data_v_wind_speed_10m_m_per_s[model->bmi.current_step];
    }

    // jframe: not sure if this belongs here or not, but it needs to happen somewhere.
    model->pet_forcing.specific_humidity_2m_kg_per_kg =  model->aorc.specific_humidity_2m_kg_per_kg;

    model->aorc.latitude                       =  model->solar_params.latitude_degrees;
    model->aorc.longitude                      =  model->solar_params.longitude_degrees;

    // Heights validated at beginning of run_pet function
    // wind speed was measured at 10.0 m height, so we need to calculate the wind speed at 2.0m
    if (model->pet_params.wind_speed_measurement_height_m > model->pet_params.zero_plane_displacement_height_m){
        double numerator=log(2.0/model->pet_params.zero_plane_displacement_height_m);
        double denominator=log(model->pet_params.wind_speed_measurement_height_m/model->pet_params.zero_plane_displacement_height_m);
        model->pet_forcing.wind_speed_m_per_s = model->pet_forcing.wind_speed_m_per_s*numerator/denominator;  // this is the 2 m value
    }
    // Otherwise, log profile isn't valid...Could use a power law, such as
    // u(z) = u_ref * (z / z_ref)^α for some exponent α, 
    // but α varies with stability and surface roughness.
    // To avoid adding another parameter, though, we will just
    // Use observed wind speed as 2m value (already warned above)
    // wind_speed_measurement_height_m <= zero_plane_displacement_height_m. Using observed wind speed as 2m value.\n");
  
    /* FLO 10/26: The original code overwrote the configured measurement
     * height with 2 m after the first update. That prevented the next
     * forcing record from receiving its 10 m to 2 m adjustment. Retain
     * the measurement height and use 2 m only in downstream formulas.
     */
    // transfer aorc forcing data into our data structure for surface radiation calculations
    model->surf_rad_forcing.incoming_shortwave_radiation_W_per_sq_m = (double)model->aorc.incoming_shortwave_W_per_m2;
    model->surf_rad_forcing.incoming_longwave_radiation_W_per_sq_m  = (double)model->aorc.incoming_longwave_W_per_m2; 
    model->surf_rad_forcing.air_temperature_C                       = (double)model->aorc.air_temperature_2m_K-TK;

    // compute relative humidity from specific humidity...
    double saturation_vapor_pressure_Pa = calc_air_saturation_vapor_pressure_Pa(model->surf_rad_forcing.air_temperature_C);

    double actual_vapor_pressure_Pa = (double)model->aorc.specific_humidity_2m_kg_per_kg*(double)model->aorc.surface_pressure_Pa/0.622;

    model->surf_rad_forcing.relative_humidity_percent = 100.0*actual_vapor_pressure_Pa/saturation_vapor_pressure_Pa;
    // sanity check the resulting value.  Should be less than 100%.  Sometimes air can be supersaturated.
    if(100.0< model->surf_rad_forcing.relative_humidity_percent) model->surf_rad_forcing.relative_humidity_percent = 99.0;
  }

  if(model->pet_options.shortwave_radiation_provided==0)
  {
    // populate the elements of the structures needed to calculate shortwave (solar) radiation, and calculate it
    // ### OPTIONS ###
    model->solar_options.cloud_base_height_known=0;  // set to TRUE if the solar_forcing.cloud_base_height_m is known.

    calculate_solar_radiation(model);
  }
  
  // we must calculate the net radiation before calling the ET subroutine.
  if(model->pet_options.use_aerodynamic_method==0) 
  {
    if (model->bmi.verbose > 1)
      printf("calculate the net radiation before calling the PET subroutine");
    // NOTE don't call this function use_aerodynamic_method option is TRUE
    model->pet_forcing.net_radiation_W_per_sq_m=calculate_net_radiation_W_per_sq_m(model);
  }

  if(model->pet_options.use_energy_balance_method ==1)
    model->pet_m_per_s=pevapotranspiration_energy_balance_method(model);
  if(model->pet_options.use_aerodynamic_method ==1)
    model->pet_m_per_s=pevapotranspiration_aerodynamic_method(model);
  if(model->pet_options.use_combination_method ==1)
    model->pet_m_per_s=pevapotranspiration_combination_method(model);
  if(model->pet_options.use_priestley_taylor_method ==1)
    model->pet_m_per_s=pevapotranspiration_priestley_taylor_method(model);
  if(model->pet_options.use_penman_monteith_method ==1)
    model->pet_m_per_s=pevapotranspiration_penman_monteith_method(model);

  // prevent dew from forming (i.e., PET < 0)
  if(model->pet_m_per_s<0) {
    model->pet_m_per_s = 0;
  }

  if (model->bmi.verbose >=1){
    printf("\n");
    printf("_______________________________________________________________________________\n");
    if(model->pet_options.use_energy_balance_method ==1)   printf("energy balance method:\n");
    if(model->pet_options.use_aerodynamic_method ==1)      printf("aerodynamic method:\n");
    if(model->pet_options.use_combination_method ==1)      printf("combination method:\n");
    if(model->pet_options.use_priestley_taylor_method ==1) printf("Priestley-Taylor method:\n");
    if(model->pet_options.use_penman_monteith_method ==1)  printf("Penman Monteith method:\n");

    printf("calculated instantaneous potential evapotranspiration (PET) =%8.6e m/s\n",model->pet_m_per_s);
    if (model->bmi.verbose > 1)
      printf("calculated instantaneous potential evapotranspiration (PET) =%8.6lf mm/d\n",model->pet_m_per_s*86400.0*1000.0);
  
  }

  return 0;
}

//########################    SETUP    ########    SETUP    ########    SETUP    ########################################
//########################    SETUP    ########    SETUP    ########    SETUP    ########################################
//########################    SETUP    ########    SETUP    ########    SETUP    ########################################
//########################    SETUP    ########    SETUP    ########    SETUP    ########################################
void pet_setup(pet_model* model)
{

  //##########################################################
  // THE VALUE OF THESE FLAGS DETERMINE HOW THIS CODE BEHAVES.
  //##########################################################
  model->pet_options.use_energy_balance_method   = 0;
  model->pet_options.use_aerodynamic_method      = 0;
  model->pet_options.use_combination_method      = 0;
  model->pet_options.use_priestley_taylor_method = 0;
  model->pet_options.use_penman_monteith_method  = 0;
  if (model->pet_method == 1)
    model->pet_options.use_energy_balance_method   = 1;
  if (model->pet_method == 2)
    model->pet_options.use_aerodynamic_method      = 1;
  if (model->pet_method == 3)
    model->pet_options.use_combination_method      = 1;
  if (model->pet_method == 4)
    model->pet_options.use_priestley_taylor_method = 1;
  if (model->pet_method == 5)
    model->pet_options.use_penman_monteith_method  = 1;


  //###################################################################################################
  // These data now come from aorc reading/parsing function.
  
  /*
  model->aorc.incoming_longwave_W_per_m2     =  117.1;
  model->aorc.incoming_shortwave_W_per_m2    =  599.7;
  model->aorc.surface_pressure_Pa            =  101300.0;
  model->aorc.specific_humidity_2m_kg_per_kg =  0.00778;      // results in a relative humidity of 40%
  model->aorc.air_temperature_2m_K           =  25.0+TK;
  model->aorc.u_wind_speed_10m_m_per_s       =  1.54;
  model->aorc.v_wind_speed_10m_m_per_s       =  3.2;
  model->aorc.latitude                       =  37.865211;
  model->aorc.longitude                      =  -98.12345;
  model->aorc.time                           =  111111112;
  */

  // populate the evapotranspiration forcing data structure:
  // this part of code does not explicitly setting values, moved to et_wrapper_function()


  // ET forcing values that come from somewhere else...
  //---------------------------------------------------------------------------------------------------------------
  model->pet_forcing.canopy_resistance_sec_per_m   = 50.0; // TODO: from plant growth model
  model->pet_forcing.water_temperature_C           = 15.5; // TODO: from soil or lake thermal model
  model->pet_forcing.ground_heat_flux_W_per_sq_m=-10.0;    // TODO from soil thermal model.  Negative denotes downward.

  if(model->pet_options.yes_aorc!=1)
  {
    // these values are needed if we don't have incoming longwave radiation measurements.
    model->surf_rad_forcing.incoming_shortwave_radiation_W_per_sq_m     = 440.1;     // must come from somewhere
    model->surf_rad_forcing.incoming_longwave_radiation_W_per_sq_m      = -1.0e+05;  // this huge negative value tells to calc.
    model->surf_rad_forcing.air_temperature_C                           = 15.0;      // from some forcing data file
    model->surf_rad_forcing.relative_humidity_percent                   = 63.0;      // from some forcing data file
    model->surf_rad_forcing.ambient_temperature_lapse_rate_deg_C_per_km = 6.49;      // ICAO standard atmosphere lapse rate
    model->surf_rad_forcing.cloud_cover_fraction                        = 0.6;       // from some forcing data file
    model->surf_rad_forcing.cloud_base_height_m                         = 2500.0/3.281; // assumed 2500 ft.
  }

    // these values are needed if we don't have incoming longwave radiation measurements.
    model->surf_rad_forcing.incoming_shortwave_radiation_W_per_sq_m     = 440.1;     // must come from somewhere
    model->surf_rad_forcing.incoming_longwave_radiation_W_per_sq_m      = -1.0e+05;  // this huge negative value tells to calc.
    model->surf_rad_forcing.air_temperature_C                           = 15.0;      // from some forcing data file
    model->surf_rad_forcing.relative_humidity_percent                   = 63.0;      // from some forcing data file
    model->surf_rad_forcing.ambient_temperature_lapse_rate_deg_C_per_km = 6.49;      // ICAO standard atmosphere lapse rate
    model->surf_rad_forcing.cloud_cover_fraction                        = 0.6;       // from some forcing data file
    model->surf_rad_forcing.cloud_base_height_m                         = 2500.0/3.281; // assumed 2500 ft.
    
  // Surface radiation forcing parameter values that must come from other models
  //---------------------------------------------------------------------------------------------------------------
  model->surf_rad_forcing.surface_skin_temperature_C = 12.0;  // TODO from soil thermal model or vegetation model.

  if(model->pet_options.shortwave_radiation_provided==0)
  {
    // populate the elements of the structures needed to calculate shortwave (solar) radiation, and calculate it
    // ### OPTIONS ###
    model->solar_options.cloud_base_height_known=0;  // set to TRUE if the solar_forcing.cloud_base_height_m is known.

    // ### FORCING ###
    model->surf_rad_forcing.cloud_cover_fraction         =   0.5;   // THESE VALUES ARE FOR THE UNIT TEST 
    model->surf_rad_forcing.atmospheric_turbidity_factor =   2.0;   // 2.0 = clear mountain air, 5.0= smoggy air
    model->surf_rad_forcing.day_of_year                  =  208;    // THESE VALUES ARE FOR THE UNIT TEST
    model->surf_rad_forcing.zulu_time                  =  20.567; // THESE VALUES ARE FOR THE UNIT TEST
  
    // UNIT TEST RESULTS
    // CALCULATED SOLAR FLUXES
    // at time:     20.56700000 UTC
    // at site latitude: 37.250000 deg. longitude:-97.555400 deg.  elevation:303.333000 m
    // Shortwave radiation clear-sky flux calculations:
    // -above canopy/snow perpendicular to Earth-Sun line is      =964.56166277 W/m2
    // -at the top of a horizontal canopy/snow surface is:        =661.40396086 W/m2
    // Shortwave radiation clear-sky flux calculations with 0.5000 cloud cover fraction:
    // -above canopy/snow perpendicular to Earth-Sun line is      =807.82039257 W/m2
    // -at the top of a horizontal canopy/snow surface is:        =553.92581722 W/m2
    // CALCULATED ANGLES DESCRIBING VECTOR POINTING TO THE SUN
    // solar elevation angle:     43.29101185 degrees
    // solar azimuth:            225.06371958 degrees
    // local hour angle:          31.01549773 degrees
    // Number of tests passed=7 of 7.
    // UNIT TEST PASSED.
  }

  return;
}


void pet_unit_tests(pet_model* model)
{
  printf("\n #----------- BEGIN UNIT TESTS   ---------------# \n");
  printf("\n # Solar calcs assume the unit test site of configs/method_N_*.cfg: latitude 37.865211, longitude -98.12345 # \n");
  printf("\n # And that pet_doy = 208 & pet_zulu_time = 20.567\n");

  printf("\n #-----------       UNIT TEST    ---------------# \n");
  printf("solar elevation angle is %lf degrees,\n and should be: 43.168329 degrees \n",
      model->solar_results.solar_elevation_angle_degrees);
  
  printf("\n #-----------       UNIT TEST    ---------------# \n");
  printf("solar azimuth angle is %lf degrees,\n and should be: 224.010087 degrees \n",
      model->solar_results.solar_azimuth_angle_degrees);
  
  printf("\n #-----------       UNIT TEST    ---------------# \n");
  printf("solar local hour angle is %lf degrees,\n and should be: 30.447448 degrees \n",
      model->solar_results.solar_local_hour_angle_degrees);
  printf("\n #-----------       UNIT TEST    ---------------# \n");
  if (model->pet_options.use_energy_balance_method == 1)
      printf("potential ET is %8.6e m/s,\n and should be: 8.594743e-08 m/s \n", model->pet_m_per_s);
  if (model->pet_options.use_aerodynamic_method == 1)
      printf("potential ET is %8.6e m/s,\n and should be: 8.977490e-08 m/s \n", model->pet_m_per_s);
  if (model->pet_options.use_combination_method == 1)
      printf("potential ET is %8.6e m/s,\n and should be: 8.694909e-08 m/s \n", model->pet_m_per_s);
  if (model->pet_options.use_priestley_taylor_method == 1)
      printf("potential ET is %8.6e m/s,\n and should be: 8.249098e-08 m/s \n", model->pet_m_per_s);
  if (model->pet_options.use_penman_monteith_method == 1)
      printf("potential ET is %8.6e m/s,\n and should be: 7.628671e-08 m/s \n", model->pet_m_per_s);

  printf("\n #------------  END UNIT TESTS   ---------------# \n");
  printf("\n");
  printf("\n");
  return;
}

/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/* ALL THE STUFF BELOW HERE IS JUST UTILITY MEMORY AND TIME FUNCTION CODE */
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/


/*####################################################################*/
/*########################### PARSE LINE #############################*/
/*####################################################################*/
void parse_aorc_line_pet(char *theString, long *year, long *month, long *day, long *hour, long *minute, double *second,
                     struct aorc_forcing_data_pet *aorc) {
    char str[20];
    long yr, mo, da, hr, mi;
    double mm, julian, se;
    float val;
    int i, start, end, len;
    int yes_pm, wordlen;
    char theWord[150];

    len = strlen(theString);

    char *copy, *copy_to_free, *value;
    copy_to_free = copy = strdup(theString);

    // time
    value = strsep(&copy, ",");
    // TODO: handle this
    // struct tm{
    //   int tm_year;
    //   int tm_mon; 
    //   int tm_mday; 
    //   int tm_hour; 
    //   int tm_min; 
    //   int tm_sec; 
    //   int tm_isdst;
    // } t;
    struct tm t;
    time_t t_of_day;

    t.tm_year = (int)strtol(strsep(&value, "-"), NULL, 10) - 1900;
    t.tm_mon = (int)strtol(strsep(&value, "-"), NULL, 10);
    t.tm_mday = (int)strtol(strsep(&value, " "), NULL, 10);
    t.tm_hour = (int)strtol(strsep(&value, ":"), NULL, 10);
    t.tm_min = (int)strtol(strsep(&value, ":"), NULL, 10);
    t.tm_sec = (int)strtol(value, NULL, 10);
    t.tm_isdst = -1;        // Is DST on? 1 = yes, 0 = no, -1 = unknown
    aorc->time = mktime(&t);

    // APCP_surface
    value = strsep(&copy, ",");
    // Not sure what this is

    // DLWRF_surface
    value = strsep(&copy, ",");
    aorc->incoming_longwave_W_per_m2 = strtof(value, NULL);
    // DSWRF_surface
    value = strsep(&copy, ",");
    aorc->incoming_shortwave_W_per_m2 = strtof(value, NULL);
    // PRES_surface
    value = strsep(&copy, ",");
    aorc->surface_pressure_Pa = strtof(value, NULL);
    // SPFH_2maboveground
    value = strsep(&copy, ",");
    aorc->specific_humidity_2m_kg_per_kg = strtof(value, NULL);;
    // TMP_2maboveground
    value = strsep(&copy, ",");
    aorc->air_temperature_2m_K = strtof(value, NULL);
    // UGRD_10maboveground
    value = strsep(&copy, ",");
    aorc->u_wind_speed_10m_m_per_s = strtof(value, NULL);
    // VGRD_10maboveground
    value = strsep(&copy, ",");
    aorc->v_wind_speed_10m_m_per_s = strtof(value, NULL);
    // precip_rate
    value = strsep(&copy, ",");
    aorc->precip_kg_per_m2 = strtof(value, NULL);

    // Go ahead and free the duplicate copy now
    free(copy_to_free);

    return;
}

/*####################################################################*/
/*############################## GET WORD ############################*/
/*####################################################################*/
void get_word_pet(char *theString, int *start, int *end, char *theWord, int *wordlen) {
    int i, lenny, j;
    lenny = strlen(theString);

    while (theString[*start] == ' ' || theString[*start] == '\t') {
        (*start)++;
    };

    j = 0;
    for (i = (*start); i < lenny; i++) {
        if (theString[i] != ' ' && theString[i] != '\t' && theString[i] != ',' && theString[i] != ':' &&
            theString[i] != '/') {
            theWord[j] = theString[i];
            j++;
        }
        else if (theString[i]!=' ' && theString[i]!='\t' && theString[i]!=',' && theString[i]!=':' && theString[i]!='/' \
             && (theString[i]=='-' && theString[i-1]==',') )
        {
            theWord[j]=theString[i];
            j++;
        }
        else if (theString[i]!=' ' && theString[i]!='\t' && theString[i]!=',' && theString[i]!=':' && theString[i]!='/' \
             && (theString[i]=='-' && theString[i-1]=='e') )
        {
            theWord[j]=theString[i];
            j++;
        }
        else {
            break;
        }
    }
    theWord[j] = '\0';
    *start = i + 1;
    *wordlen = strlen(theWord);
    return;
}

/****************************************/
void itwo_alloc_pet(int ***array, int rows, int cols) {
    int i, frows, fcols, numgood = 0;
    int error = 0;

    if ((rows == 0) || (cols == 0)) {
        printf("Error: Attempting to allocate array of size 0\n");
        exit(1);
    }

    frows = rows + 1;  /* added one for FORTRAN numbering */
    fcols = cols + 1;  /* added one for FORTRAN numbering */

    *array = (int **) malloc(frows * sizeof(int *));
    if (*array) {
        memset((*array), 0, frows * sizeof(int *));
        for (i = 0; i < frows; i++) {
            (*array)[i] = (int *) malloc(fcols * sizeof(int));
            if ((*array)[i] == NULL) {
                error = 1;
                numgood = i;
                i = frows;
            } else
                memset((*array)[i], 0, fcols * sizeof(int));
        }
    }
    return;
}


void dtwo_alloc_pet(double ***array, int rows, int cols) {
    int i, frows, fcols, numgood = 0;
    int error = 0;

    if ((rows == 0) || (cols == 0)) {
        printf("Error: Attempting to allocate array of size 0\n");
        exit(1);
    }

    frows = rows + 1;  /* added one for FORTRAN numbering */
    fcols = cols + 1;  /* added one for FORTRAN numbering */

    *array = (double **) malloc(frows * sizeof(double *));
    if (*array) {
        memset((*array), 0, frows * sizeof(double *));
        for (i = 0; i < frows; i++) {
            (*array)[i] = (double *) malloc(fcols * sizeof(double));
            if ((*array)[i] == NULL) {
                error = 1;
                numgood = i;
                i = frows;
            } else
                memset((*array)[i], 0, fcols * sizeof(double));
        }
    }
    return;
}


void d_alloc_pet(double **var, int size) {
    size++;  /* just for safety */

    *var = (double *) malloc(size * sizeof(double));
    if (*var == NULL) {
        printf("Problem allocating memory for array in d_alloc.\n");
        return;
    } else
        memset(*var, 0, size * sizeof(double));
    return;
}

void i_alloc_pet(int **var, int size) {
    size++;  /* just for safety */

    *var = (int *) malloc(size * sizeof(int));
    if (*var == NULL) {
        printf("Problem allocating memory in i_alloc\n");
        return;
    } else
        memset(*var, 0, size * sizeof(int));
    return;
}

/*
 * convert Gregorian days to Julian date
 *
 * Modify as needed for your application.
 *
 * The Julian day starts at noon of the Gregorian day and extends
 * to noon the next Gregorian day.
 *
 */
/*
** Takes a date, and returns a Julian day. A Julian day is the number of
** days since some base date  (in the very distant past).
** Handy for getting date of x number of days after a given Julian date
** (use jdate to get that from the Gregorian date).
** Author: Robert G. Tantzen, translator: Nat Howard
** Translated from the algol original in Collected Algorithms of CACM
** (This and jdate are algorithm 199).
*/


double greg_2_jul_pet(
        long year,
        long mon,
        long day,
        long h,
        long mi,
        double se) {
    long m = mon, d = day, y = year;
    long c, ya, j;
    double seconds = h * 3600.0 + mi * 60 + se;

    if (m > 2)
        m -= 3;
    else {
        m += 9;
        --y;
    }
    c = y / 100L;
    ya = y - (100L * c);
    j = (146097L * c) / 4L + (1461L * ya) / 4L + (153L * m + 2L) / 5L + d + 1721119L;
    if (seconds < 12 * 3600.0) {
        j--;
        seconds += 12.0 * 3600.0;
    } else {
        seconds = seconds - 12.0 * 3600.0;
    }
    return (j + (seconds / 3600.0) / 24.0);
}

/* Julian date converter. Takes a julian date (the number of days since
** some distant epoch or other), and returns an int pointer to static space.
** ip[0] = month;
** ip[1] = day of month;
** ip[2] = year (actual year, like 1977, not 77 unless it was  77 a.d.);
** ip[3] = day of week (0->Sunday to 6->Saturday)
** These are Gregorian.
** Copied from Algorithm 199 in Collected algorithms of the CACM
** Author: Robert G. Tantzen, Translator: Nat Howard
**
** Modified by FLO 4/99 to account for nagging round off error 
**
*/
void calc_date_pet(double jd, long *y, long *m, long *d, long *h, long *mi,
               double *sec) {
    static int ret[4];

    long j;
    double tmp;
    double frac;

    j = (long) jd;
    frac = jd - j;

    if (frac >= 0.5) {
        frac = frac - 0.5;
        j++;
    } else {
        frac = frac + 0.5;
    }

    ret[3] = (j + 1L) % 7L;
    j -= 1721119L;
    *y = (4L * j - 1L) / 146097L;
    j = 4L * j - 1L - 146097L * *y;
    *d = j / 4L;
    j = (4L * *d + 3L) / 1461L;
    *d = 4L * *d + 3L - 1461L * j;
    *d = (*d + 4L) / 4L;
    *m = (5L * *d - 3L) / 153L;
    *d = 5L * *d - 3 - 153L * *m;
    *d = (*d + 5L) / 5L;
    *y = 100L * *y + j;
    if (*m < 10)
        *m += 3;
    else {
        *m -= 9;
        *y = *y + 1; /* Invalid use: *y++. Modified by Tony */
    }

    /* if (*m < 3) *y++; */
    /* incorrectly repeated the above if-else statement. Deleted by Tony.*/

    tmp = 3600.0 * (frac * 24.0);
    *h = (long) (tmp / 3600.0);
    tmp = tmp - *h * 3600.0;
    *mi = (long) (tmp / 60.0);
    *sec = tmp - *mi * 60.0;
}





#endif  // PET_C
