#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <ctype.h>
#include <errno.h>
#include "../include/pet.h"
#include "../include/bmi.h"
#include "../include/bmi_pet.h"

#define INPUT_VAR_NAME_COUNT 36 // FLO 10/26: 7 AORC forcings, 6 Ball-Berry state/parameters, 13 coupling inputs, 3 canopy state, 6 time inputs, top-disc soil moisture
#define OUTPUT_VAR_NAME_COUNT 16 // FLO 10/26: PET plus 15 Ball-Berry coupling outputs

// Remove leading and trailing whitespace from a string in place.
// Returns a pointer to the first non-whitespace character.
static char *trim_config_whitespace(char *text)
{
    char *end;

    while (isspace((unsigned char)*text)) {
        text++;
    }

    if (*text == '\0') {
        return text;
    }

    end = text + strlen(text) - 1;

    while (end > text && isspace((unsigned char)*end)) {
        *end = '\0';
        end--;
    }

    return text;
}

static int 
Initialize (Bmi *self, const char *cfg_file)
{
    
    pet_model *pet;
    pet = (pet_model *) self->data;
    
    int config_read_result = read_init_config_pet(pet, cfg_file);
    if (config_read_result == BMI_FAILURE)
        return BMI_FAILURE;

    /* FLO 10/21: remember the configuration file name (reported in diagnostics). */
    free(pet->config_file_path);
    pet->config_file_path = strdup(cfg_file);
    if (pet->config_file_path == NULL) return BMI_FAILURE;
    /* FLO 10/21: method 6 parses its parameters, constants and the initial values of its BMI
     * input variables from this same file, once, here (ball_berry_method6_inputs.c). Its
     * forcing, time and coupling inputs then come through BMI set_value before each update,
     * so method 6 requires forcing_file = BMI. */
    if (pet->pet_method == 6) {
        if (pet->bmi.is_forcing_from_bmi != 1) {
            fprintf(stderr, "PET method 6 requires forcing_file = BMI: the driver sets the forcing "
                            "and time with set_value (configuration %s)\n", cfg_file);
            return BMI_FAILURE;
        }
        if (bb6_initialize_from_configuration(pet, cfg_file) != 0) {
            fprintf(stderr, "PET method 6: invalid configuration in %s\n", cfg_file);
            return BMI_FAILURE;
        }
    }

    pet_setup(pet);

    if (pet->bmi.verbose >1)
        printf("BMI Initialization PET ... setup just finished \n");
    
    /*
        We might be taking forcing data from the framework via BMI.
        Or we might be reading in our own forcing data. 
        We need a logical function to determine this functionality.
        We will use the pet->bmi.is_forcing_from_bmi as the boolean option.
        That should be determined in pet_setup.
    */
    if (pet->bmi.is_forcing_from_bmi == 1)
        if (pet->bmi.verbose>1)
            printf("Using BMI to pass in forcing data, not reading in forcing from file.\n");
    if (pet->bmi.is_forcing_from_bmi == 0){
        if (pet->bmi.verbose>1)
            printf("Reading in forcing from file. %s\n", pet->forcing_file);

        // Figure out the number of lines first (also char count)
        int forcing_line_count, max_forcing_line_length;
        int count_result = read_file_line_counts_pet(pet->forcing_file, &forcing_line_count, &max_forcing_line_length);
        if (count_result == -1) {
            printf("Configured forcing file '%s' could not be opened for reading\n", pet->forcing_file);
            return BMI_FAILURE;
        }
        if (forcing_line_count <= 1 || pet->bmi.num_timesteps > forcing_line_count - 1) {
            printf("Invalid header-only forcing file '%s'\n", pet->forcing_file);
            return BMI_FAILURE;
        }
    
        // Now initialize empty arrays that depend on number of time steps
        pet->forcing_data_precip_kg_per_m2 = malloc(sizeof(double) * (pet->bmi.num_timesteps + 1));
        pet->forcing_data_surface_pressure_Pa = malloc(sizeof(double) * (pet->bmi.num_timesteps + 1));
        pet->forcing_data_time = malloc(sizeof(long) * (pet->bmi.num_timesteps + 1));
        pet->forcing_data_incoming_shortwave_W_per_m2 = malloc(sizeof(double) * (pet->bmi.num_timesteps + 1));
        pet->forcing_data_incoming_longwave_W_per_m2 = malloc(sizeof(double) * (pet->bmi.num_timesteps + 1));
        pet->forcing_data_specific_humidity_2m_kg_per_kg = malloc(sizeof(double) * (pet->bmi.num_timesteps + 1));
        pet->forcing_data_air_temperature_2m_K = malloc(sizeof(double) * (pet->bmi.num_timesteps + 1));
        pet->forcing_data_u_wind_speed_10m_m_per_s = malloc(sizeof(double) * (pet->bmi.num_timesteps + 1));
        pet->forcing_data_v_wind_speed_10m_m_per_s = malloc(sizeof(double) * (pet->bmi.num_timesteps + 1));
    
        // Now open it again to read the forcings
        if (!pet->forcing_data_precip_kg_per_m2 || !pet->forcing_data_surface_pressure_Pa ||
            !pet->forcing_data_time || !pet->forcing_data_incoming_shortwave_W_per_m2 ||
            !pet->forcing_data_incoming_longwave_W_per_m2 ||
            !pet->forcing_data_specific_humidity_2m_kg_per_kg ||
            !pet->forcing_data_air_temperature_2m_K ||
            !pet->forcing_data_u_wind_speed_10m_m_per_s ||
            !pet->forcing_data_v_wind_speed_10m_m_per_s) return BMI_FAILURE;
        FILE* ffp = fopen(pet->forcing_file, "r");
        // Ensure still exists
        if (ffp == NULL) {
            printf("Forcing file '%s' disappeared!", pet->forcing_file);
            return BMI_FAILURE;
        }
    
        // Read forcing file and parse forcings
        char line_str[max_forcing_line_length + 1];
        long year, month, day, hour, minute;
        double dsec;
        char* ret = NULL;
        // First read the header line
        ret = fgets(line_str, max_forcing_line_length + 1, ffp);
        if (ret == NULL) { fclose(ffp); return BMI_FAILURE; }

        if (pet->bmi.verbose > 2) 
            printf("the number of time steps from the forcing file is: %8.6e \n", (double)pet->bmi.num_timesteps);
    
        aorc_forcing_data_pet forcings;
        for (int i = 0; i < pet->bmi.num_timesteps; i++) {
            ret = fgets(line_str, max_forcing_line_length + 1, ffp);  // read in a line of AORC data.
            if (ret == NULL) { fclose(ffp); return BMI_FAILURE; }
            parse_aorc_line_pet(line_str, &year, &month, &day, &hour, &minute, &dsec, &forcings);
            pet->forcing_data_precip_kg_per_m2[i] = forcings.precip_kg_per_m2 * ((double)pet->bmi.time_step_size_s);
            if (pet->bmi.verbose >4)
                printf("precip %f \n", pet->forcing_data_precip_kg_per_m2[i]);
            pet->forcing_data_surface_pressure_Pa[i] = forcings.surface_pressure_Pa;
            if (pet->bmi.verbose >4)
                printf("surface pressure %f \n", pet->forcing_data_surface_pressure_Pa[i]);
            pet->forcing_data_incoming_longwave_W_per_m2[i] = forcings.incoming_longwave_W_per_m2;
            if (pet->bmi.verbose >4)
                printf("longwave %f \n", pet->forcing_data_incoming_longwave_W_per_m2[i]);
            pet->forcing_data_incoming_shortwave_W_per_m2[i] = forcings.incoming_shortwave_W_per_m2;
            if (pet->bmi.verbose >4)
                printf("shortwave %f \n", pet->forcing_data_incoming_shortwave_W_per_m2[i]);
            pet->forcing_data_specific_humidity_2m_kg_per_kg[i] = forcings.specific_humidity_2m_kg_per_kg;
            if (pet->bmi.verbose >4)
                printf("humidity %f \n", pet->forcing_data_specific_humidity_2m_kg_per_kg[i]);
            pet->forcing_data_air_temperature_2m_K[i] = forcings.air_temperature_2m_K;
            if (pet->bmi.verbose >4)
                printf("air temperature %f \n", pet->forcing_data_air_temperature_2m_K[i]);
            pet->forcing_data_u_wind_speed_10m_m_per_s[i] = forcings.u_wind_speed_10m_m_per_s;
            if (pet->bmi.verbose >4)
                printf("u wind speed %f \n", pet->forcing_data_u_wind_speed_10m_m_per_s[i]);
            pet->forcing_data_v_wind_speed_10m_m_per_s[i] = forcings.v_wind_speed_10m_m_per_s;
            if (pet->bmi.verbose >4)
                printf("v wind speed %f \n", pet->forcing_data_v_wind_speed_10m_m_per_s[i]);
    
    
            pet->forcing_data_time[i] = forcings.time;
            if (i == 0)
                pet->bmi.current_time =forcings.time;
        }
        fclose(ffp);
    }

    // Set the current time step to the first item in the forcing time series.
    // But should this be an option? Would we ever initialize to a point in the
    //     middle of a forcing file?
    pet->bmi.current_step = 0;

    return BMI_SUCCESS;
}

static int 
Update (Bmi *self)
{
    pet_model *pet = (pet_model *) self->data;
  
    if (!pet) return BMI_FAILURE;
    if (pet->bmi.verbose >1)
      printf("BMI Update PET ...\n");
  
    if (pet == NULL || pet->bmi.time_step_size_s <= 0.0 ||
        (pet->bmi.is_forcing_from_bmi == 0 &&
         pet->bmi.current_step >= pet->bmi.num_timesteps))
        return BMI_FAILURE;
    if (run_pet(pet) != 0)
        return BMI_FAILURE;

    pet->bmi.current_time_step += pet->bmi.time_step_size_s; // Seconds since start of run
    pet->bmi.current_step +=1;                            // time steps since start of run
    pet->bmi.current_time += pet->bmi.time_step_size_s;   // Seconds since 1970

    return BMI_SUCCESS;
}

// JLG: Implementing a more standard version of update_until here...

static int 
Update_until (Bmi *self, double t)
{
    // https://bmi.readthedocs.io/en/latest/#update-until
    // "the time argument can be a non-integral multiple of time steps"

    pet_model* pet = (pet_model *) self->data;
    
    double dt;
    double now;

    if(self->get_time_step (self, &dt) == BMI_FAILURE)
        return BMI_FAILURE;
    if(dt <= 0.0 || !isfinite(t))
        return BMI_FAILURE;
    if(self->get_current_time(self, &now) == BMI_FAILURE)
        return BMI_FAILURE;    

    if (t < now || (pet->bmi.is_forcing_from_bmi == 0 &&
        t > pet->bmi.num_timesteps * dt)) return BMI_FAILURE;
    /* FLO 10/26: A fractional Update previously advanced current_step by
     * one whole forcing record while advancing time by only a fraction.
     * Reject fractional requests until forcing interpolation is defined.
     */
    double steps = (t - now) / dt;
    double rounded_steps = nearbyint(steps);
    if (fabs(steps - rounded_steps) > 1.0e-9 ||
        rounded_steps < 0.0 || rounded_steps > (double)INT_MAX)
        return BMI_FAILURE;
    for (int n = 0; n < (int)rounded_steps; ++n)
        if (Update(self) != BMI_SUCCESS) return BMI_FAILURE;

  return BMI_SUCCESS;
}
/*static int
Update_until(Bmi *self, double t)
{

    pet_model* pet = ((pet_model *) self->data);
    
    // Since this model's time units are seconds, it is assumed that the param is either a valid time in seconds, a
    // relative number of time steps into the future, or invalid

    if (pet->bmi.verbose >= 1)
        printf("running update until time %lld \n", t);

    // Don't support negative parameter values
    if (t < 0.0)
        return BMI_FAILURE;

    // jmframe troubleshooting. delete this printout
    if (pet->bmi.verbose > 3)
        printf("t is greater than zero, please delete this printout\n");

    // Don't continue if current time is at or beyond end time (or we can't determine this)
    double current_time = pet->bmi.current_time_step;
    double end_time = 0;
    int current_time_result = self->get_current_time(self, &current_time);
    if (current_time_result == BMI_FAILURE)
        return BMI_FAILURE;

    // jmframe troubleshooting. delete this printout
    if (pet->bmi.verbose > 3)
        printf("passed current time result, please delete this printout\n");

    int end_time_result = self->get_end_time(self, &end_time);

    // jmframe troubleshooting. Please delete this printout
    if (pet->bmi.verbose > 3){
        printf("end time %8.6e, please delete this printout\n", &end_time);
        printf("the current time is %8.6e \n", current_time);
        printf("the end time result is %d \n", end_time_result);
    }

    if (end_time_result == BMI_FAILURE || current_time >= end_time){
        printf("BMI Failued due to end_time_result OR current time greater than end time \n");
        if (end_time_result == BMI_FAILURE)
            printf("end_time_result failed %s \n");
        if (current_time >= end_time)
            printf("current_time %8.6e is greater than end time %8.6e \n", current_time, end_time);
        return BMI_FAILURE;
    }

    // jmframe troubleshooting. Please delete this printout
    if (pet->bmi.verbose > 3)
        printf("passed end time result, please delete this printout\n");

    // Handle easy case of t == current_time by just returning success
    if (t == current_time)
        return BMI_SUCCESS;

    // jmframe troubleshooting. Please delete this printout
    if (pet->bmi.verbose > 3)
        printf("About to extract future time, please delete this printout\n");

    // First, determine if t is some future time that will be arrived at exactly after some number of future time steps
    int is_exact_future_time = (t == end_time) ? TRUE : FALSE;
    // Compare to time step endings unless obvious that t lines up (i.e., t == end_time) or doesn't (t <= current_time)
    if (is_exact_future_time == FALSE && t > current_time) {
        int future_time_step = pet->bmi.current_time_step;
        double future_time_step_time = current_time;
        while (future_time_step < pet->bmi.num_timesteps && future_time_step_time < end_time) {
            future_time_step_time += pet->bmi.time_step_size_s;
            if (future_time_step_time == t) {
                is_exact_future_time = TRUE;
                break;
            }
        }
    }

    // jmframe troubleshooting. Please delete
    if (pet->bmi.verbose > 3)
        printf("Extracted future time, please delete this printout\n");

    // If it is an exact time, advance to that time step
    if (is_exact_future_time == TRUE) {
        while (current_time < t) {
            // jmframe troubleshooting. Please delete
            printf("Updating PET with exact time\n");
            
            self->update(self);
            
            // jmframe troubleshooting. Please delete
            if (pet->bmi.verbose > 3){
                printf("current time step %f \n", pet->bmi.current_time_step);
                printf("current step %f \n", pet->bmi.current_step);
                printf("current time %f \n", pet->bmi.current_time);
            }

        }
        return BMI_SUCCESS;
    }

    // jmframe troubleshooting. Please delete
    if (pet->bmi.verbose > 3)
        printf("Checked if exact time, please delete this printout \n");

    // If t is not an exact time, it could be a number of time step forward to proceed

    // The model doesn't support partial time step value args (i.e., fractions)
    int t_int = (int) t;
    if ((t - ((double)t_int)) != 0)
        return BMI_FAILURE;

    // jmframe troubleshooting. Please delete
    if (pet->bmi.verbose > 3)
        printf("Checked setting t_int to t, please delete this line \n");

    // Keep in mind the current_time_step hasn't been processed yet (hence, using <= for this test)
    // E.g., if (unprocessed) current_time_step = 0, t = 2, num_timesteps = 2, this is valid a valid t (run 0, run 1)
    if ((pet->bmi.current_step + t_int) <= pet->bmi.num_timesteps) {
        for (int i = 0; i < t_int; i++){

            // jmframe troubleshooting. Please delete
            if (pet->bmi.verbose > 3)
                printf("Updating PET with time step\n");
            
            // Call the update function, instead of run
            // because update also advances the bmi time
            self->update(self);

            // jmframe troubleshooting. Please delete
            if (pet->bmi.verbose > 3){
                printf("current time step %f \n", pet->bmi.current_time_step);
                printf("current step %ld \n", pet->bmi.current_step);
                printf("current time %f \n", pet->bmi.current_time);
            }
        } // end time loop
        return BMI_SUCCESS;
    }

    // If we arrive here, t wasn't an exact time at end of a time step or a valid relative time step jump, so invalid.
    return BMI_FAILURE;
    
}*/


pet_model *
new_bmi_pet()
{
    pet_model *data;
    data = (pet_model*) calloc(1, sizeof(pet_model));

    return data;
}

static int
Finalize (Bmi *self)
{
    if (self == NULL || self->data == NULL)
        return BMI_FAILURE;
    pet_model *pet = (pet_model *)self->data;
    if (pet->bmi.run_unit_tests == 1)
        pet_unit_tests(pet);
    free(pet->forcing_file);
    free(pet->config_file_path);           // FLO 10/21
    free(pet->bb_iteration_trace_file);    // FLO 10/20
    free(pet->bb_diagnostic_summary_file); // FLO 10/21
    free(pet->forcing_data_precip_kg_per_m2);
    free(pet->forcing_data_surface_pressure_Pa);
    free(pet->forcing_data_time);
    free(pet->forcing_data_incoming_shortwave_W_per_m2);
    free(pet->forcing_data_incoming_longwave_W_per_m2);
    free(pet->forcing_data_specific_humidity_2m_kg_per_kg);
    free(pet->forcing_data_air_temperature_2m_K);
    free(pet->forcing_data_u_wind_speed_10m_m_per_s);
    free(pet->forcing_data_v_wind_speed_10m_m_per_s);
    free(pet);
    self->data = NULL;
    return BMI_SUCCESS;
}

//---------------------------------------------------------------------------------------------------------------------
static int 
Get_start_time (Bmi *self, double * time)
{
    // TODO: Get the real starting time, not zero. 
    // jmframe: I think this should be the time of the first line in the forcing file. In seconds since 1970
    if (!self || !self->data || !time) return BMI_FAILURE;
    *time = 0.0; /* BMI simulation-relative seconds */
    return BMI_SUCCESS;
}

//---------------------------------------------------------------------------------------------------------------------
static const char *output_var_names[OUTPUT_VAR_NAME_COUNT] = {
    "water_potential_evaporation_flux",
    "canopy_transpiration_flux_m_per_s",
    "canopy_evaporation_heat_flux_W_per_m2",
    "ground_shortwave_absorbed_below_canopy_W_per_m2",
    "ground_longwave_absorbed_below_canopy_W_per_m2",
    "canopy_air_temperature_K",
    "canopy_air_vapor_pressure_Pa",
    "below_canopy_aerodynamic_resistance_s_per_m",
    "vegetation_temperature_K",
    "canopy_latent_heat_J_per_kg",
    "canopy_freezing_rate_mm_per_s",
    "canopy_fusion_heat_flux_W_per_m2",
    "canopy_heat_storage_flux_W_per_m2",
    "surface_resistance_with_snow_s_per_m",
    "cosine_solar_zenith_angle",
    "vegetated_fraction"
};

//---------------------------------------------------------------------------------------------------------------------
static const char *output_var_types[OUTPUT_VAR_NAME_COUNT] = {
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double"
};

//---------------------------------------------------------------------------------------------------------------------
static const int output_var_item_count[OUTPUT_VAR_NAME_COUNT] = {
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1
};

//---------------------------------------------------------------------------------------------------------------------
static const char *output_var_units[OUTPUT_VAR_NAME_COUNT] = {
    "m s-1",
    "m s-1",
    "W m-2",
    "W m-2",
    "W m-2",
    "K",
    "Pa",
    "s m-1",
    "K",
    "J kg-1",
    "mm s-1",
    "W m-2",
    "W m-2",
    "s m-1",
    "1",
    "1"
};

//---------------------------------------------------------------------------------------------------------------------
static const int output_var_grids[OUTPUT_VAR_NAME_COUNT] = {
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
};

//---------------------------------------------------------------------------------------------------------------------
static const char *output_var_locations[OUTPUT_VAR_NAME_COUNT] = {
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node"
};

//---------------------------------------------------------------------------------------------------------------------
static const int input_var_item_count[INPUT_VAR_NAME_COUNT] = {
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1
};

//---------------------------------------------------------------------------------------------------------------------
static const char input_var_grids[INPUT_VAR_NAME_COUNT] = {
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
};

//---------------------------------------------------------------------------------------------------------------------
static const char *input_var_locations[INPUT_VAR_NAME_COUNT] = {
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node",
    "node"
};




//---------------------------------------------------------------------------------------------------------------------
// Don't forget to update Get_value/Get_value_at_indices 
// (and setter) implementation if these are adjusted
static const char *input_var_names[INPUT_VAR_NAME_COUNT] = {
    "land_surface_radiation~incoming~longwave__energy_flux",
    "land_surface_air__pressure",
    "atmosphere_air_water~vapor__specific_humidity",
    "land_surface_radiation~incoming~shortwave__energy_flux",
    "land_surface_air__temperature",
    "land_surface_wind__x_component_of_velocity",
    "land_surface_wind__y_component_of_velocity",
    "LAI",
    "SAI",
    "SHDFAC",
    "CWP",
    "MP",
    "VCMX25",
    "root_stress_factor_0_to_1",
    "ground_skin_temperature_K",
    "surface_resistance_s_per_m",
    "snow_cover_fraction",
    "snow_depth_m",
    "snow_albedo_direct_visible",
    "snow_albedo_direct_nir",
    "snow_albedo_diffuse_visible",
    "snow_albedo_diffuse_nir",
    "canopy_wetted_fraction",
    "canopy_liquid_water_mm",
    "canopy_ice_mm",
    "precipitation_advected_heat_to_canopy_W_per_m2",
    "vegetation_temperature_start_of_step_K",
    "canopy_air_temperature_start_of_step_K",
    "canopy_air_vapor_pressure_start_of_step_Pa",
    "time_year",
    "time_month",
    "time_day_of_month",
    "time_hour",
    "time_minute",
    "time_second",
    "soil_top_disc_moisture_content"
};

//---------------------------------------------------------------------------------------------------------------------
static const char *input_var_types[INPUT_VAR_NAME_COUNT] = {
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double",
    "double"
};

//---------------------------------------------------------------------------------------------------------------------
static const char *input_var_units[INPUT_VAR_NAME_COUNT] = {
    "W m-2",
    "Pa",
    "kg kg-1",
    "W m-2",
    "K",
    "m s-1",
    "m s-1",
    "1",
    "1",
    "1",
    "1",
    "1",
    "umol CO2 m-2 s-1",
    "1",
    "K",
    "s m-1",
    "1",
    "m",
    "1",
    "1",
    "1",
    "1",
    "1",
    "mm",
    "mm",
    "W m-2",
    "K",
    "K",
    "Pa",
    "year",
    "month",
    "d",
    "h",
    "min",
    "s",
    "m3 m-3"
};

//---------------------------------------------------------------------------------------------------------------------
/* FLO 10/21: the variables an instance exposes depend on its method. Methods 1-5 expose the
 * seven AORC forcing inputs and the one PET output, as before method 6 existed. Method 6 also
 * exposes the Ball-Berry inputs (vegetation state, coupling, time) and its eight coupling outputs.
 * The tables above list the methods 1-5 variables first, so an instance uses the first
 * active_*_variable_count entries. Before initialize, the methods 1-5 set is reported. */
#define METHODS_1_TO_5_INPUT_VARIABLE_COUNT 7
#define METHODS_1_TO_5_OUTPUT_VARIABLE_COUNT 1

static int active_input_variable_count(Bmi *self)
{
    if (self && self->data && ((pet_model *)self->data)->pet_method == 6) return INPUT_VAR_NAME_COUNT;
    return METHODS_1_TO_5_INPUT_VARIABLE_COUNT;
}

static int active_output_variable_count(Bmi *self)
{
    if (self && self->data && ((pet_model *)self->data)->pet_method == 6) return OUTPUT_VAR_NAME_COUNT;
    return METHODS_1_TO_5_OUTPUT_VARIABLE_COUNT;
}

static int is_active_variable(Bmi *self, const char *name)
{
    for (int i = 0; i < active_input_variable_count(self); ++i)
        if (strcmp(name, input_var_names[i]) == 0) return 1;
    for (int i = 0; i < active_output_variable_count(self); ++i)
        if (strcmp(name, output_var_names[i]) == 0) return 1;
    return 0;
}


//---------------------------------------------------------------------------------------------------------------------

static int Get_end_time (Bmi *self, double *time)
{
    if (!self || !self->data || !time) return BMI_FAILURE;
    pet_model *pet = (pet_model *)self->data;
    if (pet->bmi.num_timesteps < 0 || pet->bmi.time_step_size_s <= 0.0)
        return BMI_FAILURE;
    *time = pet->bmi.num_timesteps * pet->bmi.time_step_size_s;
    return BMI_SUCCESS;
}

//---------------------------------------------------------------------------------------------------------------------

static int Get_time_step (Bmi *self, double * dt)
{
    *dt = ((pet_model *) self->data)->bmi.time_step_size_s;
    return BMI_SUCCESS;
}

static int Get_time_units (Bmi *self, char * units)
{
    strncpy (units, "s", BMI_MAX_UNITS_NAME);
    return BMI_SUCCESS;
}

//---------------------------------------------------------------------------------------------------------------------

static int Get_current_time (Bmi *self, double * time)
{
    Get_start_time(self, time);
    if (((pet_model *) self->data)->bmi.verbose > 2){
        printf("Current model time step: '%8.6e'\n", ((pet_model *) self->data)->bmi.current_time_step);
    }
    *time += ((pet_model *)self->data)->bmi.current_time_step;
    return BMI_SUCCESS;
} // end Get_current_time


//---------------------------------------------------------------------------------------------------------------------
/** Count the number of values in a delimited string representing an array of values. */

static int count_delimited_values(char* string_val, char* delimiter)
{
    char *copy, *copy_to_free, *value;
    int count = 0;

    // Make duplicate to avoid changing original string
    // Then work on copy, but keep 2nd pointer to copy so that memory can be freed
    copy_to_free = copy = strdup(string_val);
    while ((value = strsep(&copy, delimiter)) != NULL)
        count++;
    free(copy_to_free);
    return count;
} //end count_delimited_values

//---------------------------------------------------------------------------------------------------------------------
int read_file_line_counts_pet(const char* file_name, int* line_count, int* max_line_length)
{

    *line_count = 0;
    *max_line_length = 0;
    int current_line_length = 0;
    FILE* fp = fopen(file_name, "r");
    // Ensure exists
    if (fp == NULL) {
        printf("File does not exist.\n Failed in function read_file_line_counts_pet\n");
        return -1;
    }
    int seen_non_whitespace = 0;
    int c; //EOF is a negative constant...and char may be either signed OR unsigned
    //depending on the compiler, system, architecture, etc.  So there are cases
    //where this loop could go infinite comparing EOF to unsigned char
    //the return of fgetc is int, and should be stored as such!
    //https://stackoverflow.com/questions/35356322/difference-between-int-and-char-in-getchar-fgetc-and-putchar-fputc
    for (c = fgetc(fp); c != EOF; c = fgetc(fp)) {
        // keep track if this line has seen any char other than space or tab
        if (c != ' ' && c != '\t' && c != '\n')
            seen_non_whitespace++;
        // Update line count, reset non-whitespace count, adjust max_line_length (if needed), and reset current line count
        if (c == '\n') {
            *line_count += 1;
            seen_non_whitespace = 0;
            if (current_line_length > *max_line_length)
                *max_line_length = current_line_length;
            current_line_length = 0;
        }
        else {
            current_line_length += 1;
        }
    }
    fclose(fp);

    // If we saw some non-whitespace char on last line, assume last line didn't have its own \n, so count needs to be
    // incremented by 1.
    if (seen_non_whitespace > 0) {
        *line_count += 1;
    }

    // Before returning, increment the max line length by 1, since the \n will be on the line also.
    *max_line_length += 1;

    return 0;
}  // end: read_file_line_counts

//---------------------------------------------------------------------------------------------------------------------
int read_init_config_pet(pet_model* model, const char* config_file)//,
{
    if (!model || !config_file) return BMI_FAILURE;
    int config_line_count, max_config_line_length;
    // Note that this determines max line length including the ending return character, if present
    int count_result = read_file_line_counts_pet(config_file, &config_line_count, &max_config_line_length);
    if (count_result == -1) {
        return BMI_FAILURE;
    }

    FILE* fp = fopen(config_file, "r");
    if (fp == NULL)
        return BMI_FAILURE;

    char* ret = NULL;

    // TODO: document config file format (<param_key>=<param_val>, where array values are comma delim strings)

    char config_line[max_config_line_length + 1];

    for (int i = 0; i < config_line_count; i++) {
        char *param_key;
        char *param_value;
        char *equals;
        char *comment;
        char *units;

        ret = fgets(config_line, max_config_line_length + 1, fp);
        if (ret == NULL) {
            fclose(fp);
            return BMI_FAILURE;
        }

        // Remove comments, including comments following a parameter.
        comment = strchr(config_line, '#');
        if (comment != NULL) {
            *comment = '\0';
        }

        // Ignore blank lines and full-line comments.
        param_key = trim_config_whitespace(config_line);
        if (*param_key == '\0') {
            continue;
        }

        // Split the parameter name and value at the first equals sign.
        equals = strchr(param_key, '=');
        if (equals == NULL) {
            fprintf(stderr, "Invalid PET configuration line: %s\n", param_key);
            fclose(fp);
            return BMI_FAILURE;
        }

        *equals = '\0';

        param_value = trim_config_whitespace(equals + 1);
        param_key = trim_config_whitespace(param_key);

        if (*param_key == '\0' || *param_value == '\0') {
            fprintf(stderr, "Missing PET configuration key or value\n");
            fclose(fp);
            return BMI_FAILURE;
        }

        // Remove optional trailing units annotations, e.g. 3600 [s].
        // Units are documentation; no unit conversion is performed.
        if (strcmp(param_key, "forcing_file") != 0 && strcmp(param_key, "bb_iteration_trace_file") != 0 && strcmp(param_key, "bb_diagnostic_summary_file") != 0 && strcmp(param_key, "forcing_csv_file") != 0) {   // FLO 10/20: string-valued keys keep their text
            units = strchr(param_value, '[');

            if (units != NULL) {
                char *units_end = strchr(units, ']');

                if (units_end == NULL ||
                    *trim_config_whitespace(units_end + 1) != '\0') {
                    fprintf(stderr, "Invalid units annotation for %s\n", param_key);
                    fclose(fp);
                    return BMI_FAILURE;
                }

                *units = '\0';
                param_value = trim_config_whitespace(param_value);

                if (*param_value == '\0') {
                    fclose(fp);
                    return BMI_FAILURE;
                }
            }
        }
        /* FLO 10/26: strtod(..., NULL) silently accepted malformed numeric
         * settings as zero. Validate every numeric configuration value before
         * assigning it, while leaving the forcing_file string untouched.
         */
        if (strcmp(param_key, "forcing_file") != 0 && strcmp(param_key, "bb_iteration_trace_file") != 0 && strcmp(param_key, "bb_diagnostic_summary_file") != 0 && strcmp(param_key, "forcing_csv_file") != 0 && strcmp(param_key, "cloud_base_height_known") != 0) {
            char *end = NULL;
            errno = 0;
            double parsed_value = strtod(param_value, &end);
            while (end && isspace((unsigned char)*end)) ++end;
            if (errno == ERANGE || !isfinite(parsed_value) || end == param_value || !end || *end != '\0') {
                fprintf(stderr, "Invalid PET configuration value for %s: %s\n",
                        param_key, param_value);
                fclose(fp);
                return BMI_FAILURE;
            }
        }
        if (strcmp(param_key, "verbose") == 0){
            model->bmi.verbose = strtod(param_value, NULL);
            if(model->bmi.verbose > 1){
                printf("printing some stuff (level > 1) for unit tests and troubleshooting \n");
            }
            if(model->bmi.verbose > 2){
                printf("printing a lot of stuff (level > 2) for unit tests and troubleshooting \n");
            }
        }
        // jmframe: this should be strtol instead of strtod
        if (strcmp(param_key, "pet_method") == 0){
            model->pet_method = strtod(param_value, NULL);
            if(model->bmi.verbose > 1){
                printf("set PET method from config file \n");
                printf("%d\n", model->pet_method);
            }
        }
        if (strcmp(param_key, "yes_aorc") == 0) {
            model->pet_options.yes_aorc = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("set aorc boolean from config file \n");
                printf("%d\n", model->pet_options.yes_aorc);
            }
            continue;
        }
        // FLO 10/20: method 6 solver selection and test-only file names.
        if (strcmp(param_key, "bb_noahmp_compatibility") == 0) {
            model->bb_noahmp_compatibility = (int)strtod(param_value, NULL);
            continue;
        }
        if (strcmp(param_key, "bb_noahmp_exact_iteration") == 0) {
            model->bb_noahmp_exact_iteration = (int)strtod(param_value, NULL);
            continue;
        }
        if (strcmp(param_key, "bb_diagnostic_summary_file") == 0) {   // FLO 10/21
            free(model->bb_diagnostic_summary_file);
            model->bb_diagnostic_summary_file = strdup(param_value);
            if (!model->bb_diagnostic_summary_file) { fclose(fp); return BMI_FAILURE; }
            continue;
        }
        if (strcmp(param_key, "bb_iteration_trace_file") == 0) {
            free(model->bb_iteration_trace_file);
            model->bb_iteration_trace_file = strdup(param_value);
            if (!model->bb_iteration_trace_file) { fclose(fp); return BMI_FAILURE; }
            continue;
        }
        if (strcmp(param_key, "forcing_file") == 0) {
            free(model->forcing_file);
            model->forcing_file = strdup(param_value);
            if (!model->forcing_file) { fclose(fp); return BMI_FAILURE; }
            if (strcmp(model->forcing_file,"BMI") == 0){
                if(model->bmi.verbose >=2)
                    printf("in pet_setup: Getting forcing values from BMI. Not reading in forcing from file. \n");
                model->bmi.is_forcing_from_bmi = 1;
            }
            if(model->bmi.verbose >=2){
                printf("set forcing file from config file \n");
                printf("%s\n", model->forcing_file);
            }
            continue;
        }
        if (strcmp(param_key, "wind_speed_measurement_height_m") == 0) {
            model->pet_params.wind_speed_measurement_height_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("set wind speed measurement height from config file \n");
                printf("%lf\n", model->pet_params.wind_speed_measurement_height_m);
            }
            continue;
        }
        if (strcmp(param_key, "humidity_measurement_height_m") == 0) {
            model->pet_params.humidity_measurement_height_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("set humidity measurement height from config file \n");
                printf("%lf\n", model->pet_params.humidity_measurement_height_m);
            }
            continue;
        }
        if (strcmp(param_key, "vegetation_height_m") == 0) {
            model->pet_params.vegetation_height_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("vegetation height from config file \n");
                printf("%lf\n", model->pet_params.vegetation_height_m);
            }
            continue;
        }
        if (strcmp(param_key, "zero_plane_displacement_height_m") == 0) {
            model->pet_params.zero_plane_displacement_height_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("zero_plane_displacement height from config file \n");
                printf("%lf\n", model->pet_params.zero_plane_displacement_height_m);
            }
            continue;
        }
        if (strcmp(param_key, "shortwave_radiation_provided") == 0) {
            model->pet_options.shortwave_radiation_provided = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("shortwave radiation provided boolean from config file \n");
                printf("%d\n", model->pet_options.shortwave_radiation_provided);
            }
            continue;
        }
        if (strcmp(param_key, "momentum_transfer_roughness_length_m") == 0) {
            model->pet_params.momentum_transfer_roughness_length_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("momentum_transfer_roughness_length_m from config file \n");
                printf("%lf\n", model->pet_params.momentum_transfer_roughness_length_m);
            }
            continue;
        }
        if (strcmp(param_key, "heat_transfer_roughness_length_m") == 0) {
            model->pet_params.heat_transfer_roughness_length_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("heat_transfer_roughness_length_m from config file \n");
                printf("%lf\n", model->pet_params.heat_transfer_roughness_length_m);
            }
            continue;
        }
        if (strcmp(param_key, "surface_longwave_emissivity") == 0) {
            model->surf_rad_params.surface_longwave_emissivity = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("surface_longwave_emissivity from config file \n");
                printf("%lf\n", model->surf_rad_params.surface_longwave_emissivity);
            }
            continue;
        }
        if (strcmp(param_key, "surface_shortwave_albedo") == 0) {
            model->surf_rad_params.surface_shortwave_albedo = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("surface_shortwave_albedo from config file \n");
                printf("%lf\n", model->surf_rad_params.surface_shortwave_albedo);
            }
            continue;
        }
        /* FLO 10/26: Removed a duplicate albedo configuration branch;
         * the first branch already handles this parameter.
         */
        if (strcmp(param_key, "latitude_degrees") == 0) {
            model->solar_params.latitude_degrees = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("latitude_degrees from config file \n");
                printf("%lf\n", model->solar_params.latitude_degrees);
            }
            continue;
        }
        if (strcmp(param_key, "longitude_degrees") == 0) {
            model->solar_params.longitude_degrees = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("longitude_degrees from config file \n");
                printf("%lf\n", model->solar_params.longitude_degrees);
            }
            continue;
        }
        if (strcmp(param_key, "site_elevation_m") == 0) {
            model->solar_params.site_elevation_m = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("site_elevation_m from config file \n");
                printf("%lf\n", model->solar_params.site_elevation_m);
            }
            continue;
        }
        if (strcmp(param_key, "time_step_size_s") == 0) {
            model->bmi.time_step_size_s = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("time_step_size_s from config file \n");
                /* FLO 10/26: time_step_size_s is an int; passing it to %g
                 * invokes undefined behavior. Use the matching integer format. */
                printf("%d\n", model->bmi.time_step_size_s);
            }
            continue;
        }
        if (strcmp(param_key, "num_timesteps") == 0) {
            model->bmi.num_timesteps = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("num_timesteps from config file \n");
                printf("%ld\n", model->bmi.num_timesteps);
            }
            continue;
        }
        if (strcmp(param_key, "run_unit_tests") == 0) {
            model->bmi.run_unit_tests = strtod(param_value, NULL);
            if(model->bmi.verbose >=2){
                printf("Running unit tests \n");
            }
            continue;
        }
        // Method 6 standalone test parameters.
        /* FLO 10/21: method 6 physics inputs (geometry, optics, physiology, LAI, SAI,
         * SHDFAC, CWP, MP, VCMX25 initial values) are now read from this same file by
         * ball_berry_reference_config.c, using the key names shared with the Fortran
         * Noah-MP reference. The former bb_* physics keys and the LAI/SAI/SHDFAC/CWP/MP/
         * VCMX25 config keys were removed; those six remain BMI input variables.
         * Only the nested-solver controls are read here. */
        if (strcmp(param_key, "bb_maximum_stability_iterations") == 0) {
            model->bb_solver.maximum_stability_iterations = strtod(param_value, NULL);
            continue;
        }
        if (strcmp(param_key, "bb_maximum_canopy_iterations") == 0) {
            model->bb_solver.maximum_canopy_iterations = strtod(param_value, NULL);
            continue;
        }
        if (strcmp(param_key, "bb_stability_relative_tolerance") == 0) {
            model->bb_solver.stability_relative_tolerance = strtod(param_value, NULL);
            continue;
        }
        if (strcmp(param_key, "bb_canopy_temperature_tolerance_k") == 0) {
            model->bb_solver.canopy_temperature_tolerance_k = strtod(param_value, NULL);
            continue;
        }
        if (strcmp(param_key, "bb_maximum_temperature_step_k") == 0) {
            model->bb_solver.maximum_temperature_step_k = strtod(param_value, NULL);
            continue;
        }
    } // end loop through config
    fclose(fp);

    /* FLO 10/26: Reject invalid timestep and forcing count at initialization
     * rather than allowing a later divide-by-zero or invalid allocation.
     */
    if (model->bmi.time_step_size_s <= 0 ||
        model->bmi.num_timesteps < 0) return BMI_FAILURE;

    // Validate parameters here if needed
    if (model->pet_params.zero_plane_displacement_height_m <= 0.0){
        // if displacement height is 0 (or less for some unknown reason...)
        // set to small value to avoid division by zero in wind profile equation
        model->pet_params.zero_plane_displacement_height_m = 0.01;
    }

    return BMI_SUCCESS;
} // end: read_init_config

static int Get_var_type (Bmi *self, const char *name, char * type)
{
    // Check to see if in output array first
    for (int i = 0; i < active_output_variable_count(self); i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            strncpy(type, output_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (int i = 0; i < active_input_variable_count(self); i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            strncpy(type, input_var_types[i], BMI_MAX_TYPE_NAME);
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    type[0] = '\0';
    return BMI_FAILURE;
}

static int Get_var_itemsize (Bmi *self, const char *name, int * size)
{
    char type[BMI_MAX_TYPE_NAME];
    int type_result = Get_var_type(self, name, type);
    if (type_result != BMI_SUCCESS) {
        return BMI_FAILURE;
    }

    if (strcmp (type, "double") == 0) {
        *size = sizeof(double);
        return BMI_SUCCESS;
    }
    else if (strcmp (type, "float") == 0) {
        *size = sizeof(float);
        return BMI_SUCCESS;
    }
    else if (strcmp (type, "int") == 0) {
        *size = sizeof(int);
        return BMI_SUCCESS;
    }
    else if (strcmp (type, "short") == 0) {
        *size = sizeof(short);
        return BMI_SUCCESS;
    }
    else if (strcmp (type, "long") == 0) {
        *size = sizeof(long);
        return BMI_SUCCESS;
    }
    else {
        *size = 0;
        return BMI_FAILURE;
    }
}



//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------
static int Get_var_location (Bmi *self, const char *name, char * location)
{
    // Check to see if in output array first
    for (int i = 0; i < active_output_variable_count(self); i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            strncpy(location, output_var_locations[i], BMI_MAX_LOCATION_NAME);
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (int i = 0; i < active_input_variable_count(self); i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            strncpy(location, input_var_locations[i], BMI_MAX_LOCATION_NAME);
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    location[0] = '\0';
    return BMI_FAILURE;
}

//------------------------------------------------------------------------------
static int Get_var_grid(Bmi *self, const char *name, int *grid)
{

    // Check to see if in output array first
    for (int i = 0; i < active_output_variable_count(self); i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            *grid = output_var_grids[i];
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (int i = 0; i < active_input_variable_count(self); i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            *grid = input_var_grids[i];
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    grid[0] = '\0';
    return BMI_FAILURE;
}

// ***********************************************************
// ********* BMI: VARIABLE GETTER & SETTER FUNCTIONS *********
// ***********************************************************
static int Get_value_ptr (Bmi *self, const char *name, void **dest)
{
    /* FLO 10/26: BMI callers may supply null arguments; avoid strcmp and
     * dereference on invalid pointers.
     */
    if (!self || !self->data || !name || !dest) return BMI_FAILURE;
    *dest = NULL;
    if (!is_active_variable(self, name)) return BMI_FAILURE;   // FLO 10/21: method 6 variables only for method 6
    void *src = NULL;
    /***********************************************************/
    /***********    OUTPUT   ***********************************/
    /***********************************************************/
    if (strcmp (name, "water_potential_evaporation_flux") == 0) {
        pet_model *pet;
        pet = (pet_model *) self->data;
        src = (void*)&pet-> pet_m_per_s;
        *dest = src;
        return BMI_SUCCESS;
    }

    /***********************************************************/
    /***********    INPUT    ***********************************/
    /***********************************************************/
    if (strcmp (name, "land_surface_radiation~incoming~longwave__energy_flux") == 0) {
        pet_model *pet;
        pet = (pet_model *) self->data;
        src = (void*)&pet->aorc.incoming_longwave_W_per_m2;
        *dest = src;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_radiation~incoming~shortwave__energy_flux") == 0) {
        pet_model *pet;
        pet = (pet_model *) self->data;
        src = (void*)&pet->aorc.incoming_shortwave_W_per_m2;
        *dest = src;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_air__pressure") == 0) {
        pet_model *pet;
        pet = (pet_model *) self->data;
        src = (void*)&pet->aorc.surface_pressure_Pa;
        *dest = src;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "atmosphere_air_water~vapor__specific_humidity") == 0) {
        pet_model *pet;
        pet = (pet_model *) self->data;
        src = (void*)&pet->aorc.specific_humidity_2m_kg_per_kg;
        *dest = src;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_air__temperature") == 0) {
        pet_model *pet;
        pet = (pet_model *) self->data;
        src = (void*)&pet->aorc.air_temperature_2m_K;
        *dest = src;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_wind__x_component_of_velocity") == 0) {
        pet_model *pet;
        pet = (pet_model *) self->data;
        src = (void*)&pet->aorc.u_wind_speed_10m_m_per_s;
        *dest = src;
        return BMI_SUCCESS;
    }
    if (strcmp (name, "land_surface_wind__y_component_of_velocity") == 0) {
        pet_model *pet;
        pet = (pet_model *) self->data;
        src = (void*)&pet->aorc.v_wind_speed_10m_m_per_s;
        *dest = src;
        return BMI_SUCCESS;
    }

    pet_model *pet = (pet_model *)self->data;
    if (strcmp(name, "LAI") == 0) { *dest = &pet->bb_canopy.lai; return BMI_SUCCESS; }
    if (strcmp(name, "SAI") == 0) { *dest = &pet->bb_canopy.sai; return BMI_SUCCESS; }
    if (strcmp(name, "SHDFAC") == 0) { *dest = &pet->bb_canopy.shdfac; return BMI_SUCCESS; }
    if (strcmp(name, "CWP") == 0) { *dest = &pet->bb_cwp; return BMI_SUCCESS; }
    if (strcmp(name, "MP") == 0) { *dest = &pet->bb_stomata.mp; return BMI_SUCCESS; }
    if (strcmp(name, "VCMX25") == 0) { *dest = &pet->bb_stomata.vcmx25; return BMI_SUCCESS; }

    /* FLO 10/21: Ball-Berry (method 6) coupling inputs [Noah-MP name] */
    BallBerryCouplingInputs *coupling = &pet->bb_coupling_inputs;
    if (strcmp(name, "root_stress_factor_0_to_1") == 0) { *dest = &coupling->root_stress_factor_0_to_1; return BMI_SUCCESS; }     /* BTRAN */
    if (strcmp(name, "ground_skin_temperature_K") == 0) { *dest = &coupling->ground_skin_temperature_K; return BMI_SUCCESS; }     /* TG */
    if (strcmp(name, "surface_resistance_s_per_m") == 0) { *dest = &coupling->surface_resistance_s_per_m; return BMI_SUCCESS; }   /* RSURF */
    if (strcmp(name, "snow_cover_fraction") == 0) { *dest = &coupling->snow_cover_fraction; return BMI_SUCCESS; }                 /* FSNO */
    if (strcmp(name, "snow_depth_m") == 0) { *dest = &coupling->snow_depth_m; return BMI_SUCCESS; }                               /* SNOWH */
    if (strcmp(name, "snow_albedo_direct_visible") == 0) { *dest = &coupling->snow_albedo_direct_visible; return BMI_SUCCESS; }   /* ALBSND(1) */
    if (strcmp(name, "snow_albedo_direct_nir") == 0) { *dest = &coupling->snow_albedo_direct_nir; return BMI_SUCCESS; }           /* ALBSND(2) */
    if (strcmp(name, "snow_albedo_diffuse_visible") == 0) { *dest = &coupling->snow_albedo_diffuse_visible; return BMI_SUCCESS; } /* ALBSNI(1) */
    if (strcmp(name, "snow_albedo_diffuse_nir") == 0) { *dest = &coupling->snow_albedo_diffuse_nir; return BMI_SUCCESS; }         /* ALBSNI(2) */
    if (strcmp(name, "canopy_wetted_fraction") == 0) { *dest = &coupling->canopy_wetted_fraction; return BMI_SUCCESS; }           /* FWET */
    if (strcmp(name, "canopy_liquid_water_mm") == 0) { *dest = &coupling->canopy_liquid_water_mm; return BMI_SUCCESS; }           /* CANLIQ */
    if (strcmp(name, "canopy_ice_mm") == 0) { *dest = &coupling->canopy_ice_mm; return BMI_SUCCESS; }                             /* CANICE */
    if (strcmp(name, "precipitation_advected_heat_to_canopy_W_per_m2") == 0) {                                                     /* PAHV, FLO 10/24 */
        *dest = &coupling->precipitation_advected_heat_to_canopy_W_per_m2; return BMI_SUCCESS; }
    /* FLO 10/24: canopy state at the start of the step [TV, TAH, EAH] */
    if (strcmp(name, "vegetation_temperature_start_of_step_K") == 0) {
        *dest = &coupling->vegetation_temperature_start_of_step_K; return BMI_SUCCESS; }
    if (strcmp(name, "canopy_air_temperature_start_of_step_K") == 0) {
        *dest = &coupling->canopy_air_temperature_start_of_step_K; return BMI_SUCCESS; }
    if (strcmp(name, "canopy_air_vapor_pressure_start_of_step_Pa") == 0) {
        *dest = &coupling->canopy_air_vapor_pressure_start_of_step_Pa; return BMI_SUCCESS; }
    if (strcmp(name, "soil_top_disc_moisture_content") == 0) {                       /* SMC(1), FLO 10/26 */
        *dest = &coupling->soil_top_disc_moisture_content; return BMI_SUCCESS; }

    /* FLO 10/21: method 6 time inputs, UTC at the start of the step */
    BallBerryTimeInputs *time_inputs = &pet->bb_time_inputs;
    if (strcmp(name, "time_year") == 0) { *dest = &time_inputs->time_year; return BMI_SUCCESS; }
    if (strcmp(name, "time_month") == 0) { *dest = &time_inputs->time_month; return BMI_SUCCESS; }
    if (strcmp(name, "time_day_of_month") == 0) { *dest = &time_inputs->time_day_of_month; return BMI_SUCCESS; }
    if (strcmp(name, "time_hour") == 0) { *dest = &time_inputs->time_hour; return BMI_SUCCESS; }
    if (strcmp(name, "time_minute") == 0) { *dest = &time_inputs->time_minute; return BMI_SUCCESS; }
    if (strcmp(name, "time_second") == 0) { *dest = &time_inputs->time_second; return BMI_SUCCESS; }

    /* FLO 10/21: Ball-Berry (method 6) coupling outputs [Noah-MP name] */
    BallBerryCouplingOutputs *coupling_outputs = &pet->bb_coupling_outputs;
    if (strcmp(name, "canopy_transpiration_flux_m_per_s") == 0) { *dest = &coupling_outputs->canopy_transpiration_flux_m_per_s; return BMI_SUCCESS; }  /* ETRAN */
    if (strcmp(name, "canopy_evaporation_heat_flux_W_per_m2") == 0) { *dest = &coupling_outputs->canopy_evaporation_heat_flux_W_per_m2; return BMI_SUCCESS; } /* EVC */
    if (strcmp(name, "ground_shortwave_absorbed_below_canopy_W_per_m2") == 0) { *dest = &coupling_outputs->ground_shortwave_absorbed_below_canopy_W_per_m2; return BMI_SUCCESS; } /* SAG */
    if (strcmp(name, "ground_longwave_absorbed_below_canopy_W_per_m2") == 0) { *dest = &coupling_outputs->ground_longwave_absorbed_below_canopy_W_per_m2; return BMI_SUCCESS; }
    if (strcmp(name, "canopy_air_temperature_K") == 0) { *dest = &coupling_outputs->canopy_air_temperature_K; return BMI_SUCCESS; }          /* TAH */
    if (strcmp(name, "canopy_air_vapor_pressure_Pa") == 0) { *dest = &coupling_outputs->canopy_air_vapor_pressure_Pa; return BMI_SUCCESS; }  /* EAH */
    if (strcmp(name, "below_canopy_aerodynamic_resistance_s_per_m") == 0) { *dest = &coupling_outputs->below_canopy_aerodynamic_resistance_s_per_m; return BMI_SUCCESS; } /* RAHG */
    if (strcmp(name, "vegetation_temperature_K") == 0) { *dest = &coupling_outputs->vegetation_temperature_K; return BMI_SUCCESS; }          /* TV */
    if (strcmp(name, "canopy_latent_heat_J_per_kg") == 0) { *dest = &coupling_outputs->canopy_latent_heat_J_per_kg; return BMI_SUCCESS; }    /* LATHEAV, FLO 10/24 */
    if (strcmp(name, "canopy_freezing_rate_mm_per_s") == 0) { *dest = &coupling_outputs->canopy_freezing_rate_mm_per_s; return BMI_SUCCESS; }
    if (strcmp(name, "canopy_fusion_heat_flux_W_per_m2") == 0) { *dest = &coupling_outputs->canopy_fusion_heat_flux_W_per_m2; return BMI_SUCCESS; }
    if (strcmp(name, "canopy_heat_storage_flux_W_per_m2") == 0) { *dest = &coupling_outputs->canopy_heat_storage_flux_W_per_m2; return BMI_SUCCESS; }
    /* FLO 10/26 */
    if (strcmp(name, "surface_resistance_with_snow_s_per_m") == 0) {
        *dest = &coupling_outputs->surface_resistance_with_snow_s_per_m; return BMI_SUCCESS; }   /* RSURF */
    if (strcmp(name, "cosine_solar_zenith_angle") == 0) {
        *dest = &coupling_outputs->cosine_solar_zenith_angle; return BMI_SUCCESS; }              /* COSZ */
    if (strcmp(name, "vegetated_fraction") == 0) {
        *dest = &coupling_outputs->vegetated_fraction; return BMI_SUCCESS; }                     /* FVEG */

    return BMI_FAILURE;
}

//------------------------------------------------------------------------------
static int Get_value_at_indices (Bmi *self, const char *name, void *dest, int * inds, int len)
{
    void *src = NULL;
    int itemsize = 0;

    /* FLO 10/26: Reject invalid BMI pointers before dereferencing them. */
    if (!self || !self->data || !name || len < 0 ||
        (len > 0 && (!inds || !dest))) return BMI_FAILURE;
    for (int i = 0; i < len; ++i) if (inds[i] != 0) return BMI_FAILURE;
    if (self->get_value_ptr(self, name, &src) == BMI_FAILURE)
        return BMI_FAILURE;

    if (self->get_var_itemsize(self, name, &itemsize) == BMI_FAILURE)
        return BMI_FAILURE;

    { /* Copy the data */
        size_t i;
        size_t offset;
        char * ptr;
        for (i=0, ptr=(char*)dest; i<len; i++, ptr+=itemsize) {
            offset = inds[i] * itemsize;
            memcpy (ptr, (char*)src + offset, itemsize);
        }
    }

    return BMI_SUCCESS;
}

//------------------------------------------------------------------------------
static int Get_value(Bmi * self, const char * name, void *dest)
{
    /* FLO 10/26: A null destination previously reached memcpy and crashed. */
    if (!self || !name || !dest) return BMI_FAILURE;
    void *src = NULL;
    int nbytes = 0;

    if (self->get_value_ptr (self, name, &src) == BMI_FAILURE)
        return BMI_FAILURE;

    if (self->get_var_nbytes (self, name, &nbytes) == BMI_FAILURE)
        return BMI_FAILURE;
    
    memcpy(dest, src, nbytes);

    return BMI_SUCCESS;
}

/* Validate Ball-Berry state and calibration inputs before changing storage.
 * The BMI pointer interface also permits direct writes, so method 6 must
 * validate the complete state again immediately before computing fluxes.
 */
static int Set_ball_berry_value(Bmi *self, const char *name, const void *value)
{
    if (!self || !self->data || !name || !value) return BMI_FAILURE;
    pet_model *pet = (pet_model *)self->data;
    double proposed_value = *(const double *)value;
    if (!isfinite(proposed_value)) return BMI_FAILURE;

    if (strcmp(name, "LAI") == 0) {
        if (ball_berry_state_set(&pet->bb_canopy, BB_STATE_LAI, proposed_value) == 0) return BMI_SUCCESS;
        return BMI_FAILURE;
    }
    if (strcmp(name, "SAI") == 0) {
        if (ball_berry_state_set(&pet->bb_canopy, BB_STATE_SAI, proposed_value) == 0) return BMI_SUCCESS;
        return BMI_FAILURE;
    }
    if (strcmp(name, "SHDFAC") == 0) {
        if (ball_berry_state_set(&pet->bb_canopy, BB_STATE_SHDFAC, proposed_value) == 0) return BMI_SUCCESS;
        return BMI_FAILURE;
    }
    /* FLO 10/21: time inputs: whole numbers in range (time_second may be fractional) */
    {
        int time_check = bb6_validate_time_input(name, proposed_value);
        if (time_check == -1) return BMI_FAILURE;
        if (time_check == 0) {
            double *time_target = NULL;
            if (self->get_value_ptr(self, name, (void **)&time_target) != BMI_SUCCESS || time_target == NULL)
                return BMI_FAILURE;
            *time_target = proposed_value;
            return BMI_SUCCESS;
        }
    }
    /* FLO 10/21: coupling inputs; validate the whole set with the proposed value, change only if valid */
    {
        BallBerryCouplingInputs proposed_coupling = pet->bb_coupling_inputs;
        double *target = NULL;
        if (strcmp(name, "root_stress_factor_0_to_1") == 0) target = &proposed_coupling.root_stress_factor_0_to_1;
        else if (strcmp(name, "ground_skin_temperature_K") == 0) target = &proposed_coupling.ground_skin_temperature_K;
        else if (strcmp(name, "surface_resistance_s_per_m") == 0) target = &proposed_coupling.surface_resistance_s_per_m;
        else if (strcmp(name, "snow_cover_fraction") == 0) target = &proposed_coupling.snow_cover_fraction;
        else if (strcmp(name, "snow_depth_m") == 0) target = &proposed_coupling.snow_depth_m;
        else if (strcmp(name, "snow_albedo_direct_visible") == 0) target = &proposed_coupling.snow_albedo_direct_visible;
        else if (strcmp(name, "snow_albedo_direct_nir") == 0) target = &proposed_coupling.snow_albedo_direct_nir;
        else if (strcmp(name, "snow_albedo_diffuse_visible") == 0) target = &proposed_coupling.snow_albedo_diffuse_visible;
        else if (strcmp(name, "snow_albedo_diffuse_nir") == 0) target = &proposed_coupling.snow_albedo_diffuse_nir;
        else if (strcmp(name, "canopy_wetted_fraction") == 0) target = &proposed_coupling.canopy_wetted_fraction;
        else if (strcmp(name, "canopy_liquid_water_mm") == 0) target = &proposed_coupling.canopy_liquid_water_mm;
        else if (strcmp(name, "canopy_ice_mm") == 0) target = &proposed_coupling.canopy_ice_mm;
        else if (strcmp(name, "precipitation_advected_heat_to_canopy_W_per_m2") == 0)
            target = &proposed_coupling.precipitation_advected_heat_to_canopy_W_per_m2;   /* FLO 10/24 */
        else if (strcmp(name, "vegetation_temperature_start_of_step_K") == 0)
            target = &proposed_coupling.vegetation_temperature_start_of_step_K;
        else if (strcmp(name, "canopy_air_temperature_start_of_step_K") == 0)
            target = &proposed_coupling.canopy_air_temperature_start_of_step_K;
        else if (strcmp(name, "canopy_air_vapor_pressure_start_of_step_Pa") == 0)
            target = &proposed_coupling.canopy_air_vapor_pressure_start_of_step_Pa;
        else if (strcmp(name, "soil_top_disc_moisture_content") == 0)
            target = &proposed_coupling.soil_top_disc_moisture_content;                     /* FLO 10/26 */
        if (target != NULL) {
            *target = proposed_value;
            if (bb6_validate_coupling_inputs(&proposed_coupling) != 0) return BMI_FAILURE;
            pet->bb_coupling_inputs = proposed_coupling;
            return BMI_SUCCESS;
        }
    }
    if (proposed_value < 0.0) return BMI_FAILURE;
    if (strcmp(name, "CWP") == 0) { pet->bb_cwp = proposed_value; return BMI_SUCCESS; }
    if (strcmp(name, "MP") == 0) { pet->bb_stomata.mp = proposed_value; return BMI_SUCCESS; }
    if (strcmp(name, "VCMX25") == 0) { pet->bb_stomata.vcmx25 = proposed_value; return BMI_SUCCESS; }
    return BMI_FAILURE;
}

static int Is_ball_berry_input(const char *name)
{
    /* Ball-Berry inputs are validated by Set_ball_berry_value: every input after the
     * seven AORC forcings in input_var_names (FLO 10/21: includes the coupling and time inputs). */
    for (int i = 7; i < INPUT_VAR_NAME_COUNT; ++i)
        if (strcmp(name, input_var_names[i]) == 0) return 1;
    return 0;
}

static int Set_value (Bmi *self, const char *name, void *array)
{
    void * dest = NULL;
    int nbytes = 0;

    /* FLO 10/26: Guard strcmp against null BMI variable names. */
    if (!self || !name || !array) return BMI_FAILURE;
    int writable = 0;
    for (int i = 0; i < active_input_variable_count(self); ++i)
        if (strcmp(name, input_var_names[i]) == 0) writable = 1;
    if (!writable || !array) return BMI_FAILURE;
    if (self->get_value_ptr(self, name, &dest) == BMI_FAILURE)
        return BMI_FAILURE;

    if (self->get_var_nbytes(self, name, &nbytes) == BMI_FAILURE)
        return BMI_FAILURE;

    if (Is_ball_berry_input(name)) return Set_ball_berry_value(self, name, array);

    memcpy (dest, array, nbytes);

    return BMI_SUCCESS;
}


//------------------------------------------------------------------------------
static int Set_value_at_indices (Bmi *self, const char *name, int * inds, int len, void *src)
{
    void * to = NULL;
    int itemsize = 0;

    /* FLO 10/26: Check the name before strcmp; NULL previously crashed. */
    if (!self || !self->data || !name || len < 0 ||
        (len > 0 && (!inds || !src))) return BMI_FAILURE;
    int writable = 0;
    for (int i = 0; i < active_input_variable_count(self); ++i)
        if (strcmp(name, input_var_names[i]) == 0) writable = 1;
    if (!writable) return BMI_FAILURE;
    for (int i = 0; i < len; ++i) if (inds[i] != 0) return BMI_FAILURE;
    if (self->get_value_ptr (self, name, &to) == BMI_FAILURE)
        return BMI_FAILURE;

    if (self->get_var_itemsize(self, name, &itemsize) == BMI_FAILURE)
        return BMI_FAILURE;

    if (Is_ball_berry_input(name)) {
        if (len == 0) return BMI_SUCCESS;
        return Set_ball_berry_value(self, name, src);
    }

    { /* Copy the data */
        size_t i;
        size_t offset;
        char * ptr;
        for (i=0, ptr=(char*)src; i<len; i++, ptr+=itemsize) {
            offset = inds[i] * itemsize;
            memcpy ((char*)to + offset, ptr, itemsize);
        }
    }
    return BMI_SUCCESS;
}

// ***********************************************************
// ************ BMI: MODEL INFORMATION FUNCTIONS *************
// ***********************************************************
static int Get_component_name (Bmi *self, char * name)
{
    strncpy (name, "Potential Evapotranspiration", BMI_MAX_COMPONENT_NAME);
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_input_item_count (Bmi *self, int * count)
{
    *count = active_input_variable_count(self);
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_input_var_names (Bmi *self, char ** names)
{
    for (int i = 0; i < active_input_variable_count(self); i++) {
        strncpy (names[i], input_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_output_item_count (Bmi *self, int * count)
{
    *count = active_output_variable_count(self);
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_output_var_names (Bmi *self, char ** names)
{
    for (int i = 0; i < active_output_variable_count(self); i++) {
        strncpy (names[i], output_var_names[i], BMI_MAX_VAR_NAME);
    }
    return BMI_SUCCESS;
}

//----------------------------------------------------------------------
static int Get_var_units (Bmi *self, const char *name, char * units)
{
    // Check to see if in output array first
    for (int i = 0; i < active_output_variable_count(self); i++) {
        if (strcmp(name, output_var_names[i]) == 0) {
            strncpy(units, output_var_units[i], BMI_MAX_UNITS_NAME);
            return BMI_SUCCESS;
        }
    }
    // Then check to see if in input array
    for (int i = 0; i < active_input_variable_count(self); i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            strncpy(units, input_var_units[i], BMI_MAX_UNITS_NAME);
            return BMI_SUCCESS;
        }
    }
    // If we get here, it means the variable name wasn't recognized
    units[0] = '\0';
    return BMI_FAILURE;
}

//----------------------------------------------------------------------
static int Get_var_nbytes (Bmi *self, const char *name, int * nbytes)
{
    /* FLO 10/26: Reject null arguments before metadata lookup. */
    if (!self || !self->data || !name || !nbytes) return BMI_FAILURE;
    int item_size;
    int item_size_result = Get_var_itemsize(self, name, &item_size);
    if (item_size_result != BMI_SUCCESS) {
        return BMI_FAILURE;
    }
    int item_count = -1;
    for (int i = 0; i < active_input_variable_count(self); i++) {
        if (strcmp(name, input_var_names[i]) == 0) {
            item_count = input_var_item_count[i];
            break;
        }
    }
    if (item_count < 1) {
        for (int i = 0; i < active_output_variable_count(self); i++) {
            if (strcmp(name, output_var_names[i]) == 0) {
                item_count = output_var_item_count[i];
                break;
            }
        }
    }
    /* FLO 10/26: Unknown item counts must fail; num_timesteps is not
     * the size of a scalar BMI variable and could overflow memcpy. */
    if (item_count < 1) return BMI_FAILURE;

    *nbytes = item_size * item_count;
    return BMI_SUCCESS;
}


// ***********************************************************
// **************** BMI: MODEL GRID FUNCTIONS ****************
// ***********************************************************
/* Grid information */
static int Get_grid_rank (Bmi *self, int grid, int * rank)
{
    if (grid == 0) {
        *rank = 0;
        return BMI_SUCCESS;
    }
    else {
        *rank = -1;
        return BMI_FAILURE;
    }
}

//----------------------------------------------------------------------
static int Get_grid_size(Bmi *self, int grid, int * size)
{
    if (grid == 0) {
        *size = 1;
        return BMI_SUCCESS;
    }
    else {
        *size = -1;
        return BMI_FAILURE;
    }
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_type (Bmi *self, int grid, char * type)
{
    int status = BMI_FAILURE;

    if (grid == 0) {
        strncpy(type, "scalar", BMI_MAX_TYPE_NAME);
        status = BMI_SUCCESS;
    }
    else {
        type[0] = '\0';
        status = BMI_FAILURE;
    }
    return status;
}

//---------------------------------------------------------------------------------------------------------------------
/* Uniform rectilinear (grid type) */
static int Get_grid_shape(Bmi *self, int grid, int *shape)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_spacing(Bmi *self, int grid, double *spacing)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_origin(Bmi *self, int grid, double *origin)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
/* Non-uniform rectilinear, curvilinear (grid type)*/
static int Get_grid_x(Bmi *self, int grid, double *x)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_y(Bmi *self, int grid, double *y)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_z(Bmi *self, int grid, double *z)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
/*Unstructured (grid type)*/
static int Get_grid_node_count(Bmi *self, int grid, int *count)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_edge_count(Bmi *self, int grid, int *count)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_face_count(Bmi *self, int grid, int *count)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_edge_nodes(Bmi *self, int grid, int *edge_nodes)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_face_edges(Bmi *self, int grid, int *face_edges)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_face_nodes(Bmi *self, int grid, int *face_nodes)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
static int Get_grid_nodes_per_face(Bmi *self, int grid, int *nodes_per_face)
{
    return BMI_FAILURE;
}

//---------------------------------------------------------------------------------------------------------------------
Bmi* 
register_bmi_pet(Bmi *model)
{
    if (model) {
        model->data = (void*)new_bmi_pet();
        model->initialize = Initialize;
        model->update = Update;
        model->update_until = Update_until;
        model->finalize = Finalize;

        model->get_component_name = Get_component_name;
        model->get_input_item_count = Get_input_item_count;
        model->get_output_item_count = Get_output_item_count;
        model->get_input_var_names = Get_input_var_names;
        model->get_output_var_names = Get_output_var_names;

        model->get_var_grid = Get_var_grid;
        model->get_var_type = Get_var_type;
        model->get_var_itemsize = Get_var_itemsize;
        model->get_var_units = Get_var_units;
        model->get_var_nbytes = Get_var_nbytes;
        model->get_var_location = Get_var_location;

        model->get_current_time = Get_current_time;
        model->get_start_time = Get_start_time;
        model->get_end_time = Get_end_time;
        model->get_time_units = Get_time_units;
        model->get_time_step = Get_time_step;

        model->get_value = Get_value;
        model->get_value_ptr = Get_value_ptr;   
        model->get_value_at_indices = Get_value_at_indices;

        model->set_value = Set_value;
        model->set_value_at_indices = Set_value_at_indices;

        model->get_grid_size = Get_grid_size;    
        model->get_grid_rank = Get_grid_rank;    
        model->get_grid_type = Get_grid_type;    
        model->get_grid_shape = Get_grid_shape;
        model->get_grid_spacing = Get_grid_spacing;
        model->get_grid_origin = Get_grid_origin;
        model->get_grid_x = Get_grid_x;
        model->get_grid_y = Get_grid_y;
        model->get_grid_z = Get_grid_z;
        model->get_grid_node_count = Get_grid_node_count;
        model->get_grid_edge_count = Get_grid_edge_count;
        model->get_grid_face_count = Get_grid_face_count;
        model->get_grid_edge_nodes = Get_grid_edge_nodes;
        model->get_grid_face_edges = Get_grid_face_edges;
        model->get_grid_face_nodes = Get_grid_face_nodes;
        model->get_grid_nodes_per_face = Get_grid_nodes_per_face;

    }

    return model;
}  // end: register_bmi_et
