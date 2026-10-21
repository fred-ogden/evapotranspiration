#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "../include/bmi.h"
#include "../include/pet.h"
#include "../include/bmi_pet.h"

#include "../extern/forcing_code/include/aorc.h"
#include "../extern/forcing_code/include/bmi_aorc.h"

/***************************************************************
    Function to pass the forcing data from AORC to PET using BMI.
    This requires a lot of getters and setters, 
    so no need to clutter up main program
***************************************************************/
static int transfer_bmi_value(Bmi *aorc_bmi_model, Bmi *pet_bmi_model,
                              const char *aorc_name, const char *pet_name)
{
    double value = 0.0;

    /* FLO 10/26: Check both BMI calls. A failed getter previously left the
     * shared buffer holding the preceding variable (surface pressure),
     * which was then incorrectly passed to PET as specific humidity.
     */
    if (aorc_bmi_model->get_value(aorc_bmi_model, aorc_name, &value) != BMI_SUCCESS) {
        fprintf(stderr, "AORC BMI get_value failed: %s\n", aorc_name);
        return BMI_FAILURE;
    }
    if (pet_bmi_model->set_value(pet_bmi_model, pet_name, &value) != BMI_SUCCESS) {
        fprintf(stderr, "PET BMI set_value failed: %s\n", pet_name);
        return BMI_FAILURE;
    }
    return BMI_SUCCESS;
}

int pass_forcing_from_aorc_to_pet(Bmi *pet_bmi_model, Bmi *aorc_bmi_model)
{
    const char *temperature_name = "land_surface_air__temperature";
    const char *pressure_name = "land_surface_air__pressure";
    const char *shortwave_name = "land_surface_radiation~incoming~shortwave__energy_flux";
    const char *longwave_name = "land_surface_radiation~incoming~longwave__energy_flux";
    const char *u_wind_name = "land_surface_wind__x_component_of_velocity";
    const char *v_wind_name = "land_surface_wind__y_component_of_velocity";
    double specific_humidity_kg_per_kg = 0.0;

    if (transfer_bmi_value(aorc_bmi_model, pet_bmi_model, temperature_name, temperature_name) != BMI_SUCCESS ||
        transfer_bmi_value(aorc_bmi_model, pet_bmi_model, pressure_name, pressure_name) != BMI_SUCCESS)
        return BMI_FAILURE;

    /* FLO 10/26: AORC exposes SPFH (kg/kg) under its legacy BMI name
     * relative_saturation. Preserve that upstream interface and supply
     * the physical specific humidity to PET under PET's correct BMI name.
     */
    if (aorc_bmi_model->get_value(aorc_bmi_model,
            "atmosphere_air_water~vapor__relative_saturation",
            &specific_humidity_kg_per_kg) != BMI_SUCCESS) {
        fprintf(stderr, "AORC BMI failed to retrieve SPFH\n");
        return BMI_FAILURE;
    }
    if (specific_humidity_kg_per_kg < 0.0 || specific_humidity_kg_per_kg > 1.0) {
        fprintf(stderr, "AORC BMI returned invalid specific humidity: %.17g\n",
                specific_humidity_kg_per_kg);
        return BMI_FAILURE;
    }
    if (pet_bmi_model->set_value(pet_bmi_model,
            "atmosphere_air_water~vapor__specific_humidity",
            &specific_humidity_kg_per_kg) != BMI_SUCCESS) {
        fprintf(stderr, "PET BMI failed to set specific humidity\n");
        return BMI_FAILURE;
    }

    if (transfer_bmi_value(aorc_bmi_model, pet_bmi_model, shortwave_name, shortwave_name) != BMI_SUCCESS ||
        transfer_bmi_value(aorc_bmi_model, pet_bmi_model, longwave_name, longwave_name) != BMI_SUCCESS ||
        transfer_bmi_value(aorc_bmi_model, pet_bmi_model, u_wind_name, u_wind_name) != BMI_SUCCESS ||
        transfer_bmi_value(aorc_bmi_model, pet_bmi_model, v_wind_name, v_wind_name) != BMI_SUCCESS)
        return BMI_FAILURE;
    return BMI_SUCCESS;
}

/***************************************************************
    FLO 10/21: calendar time for PET method 6.

    Method 6 takes the calendar time (UTC, start of the forcing
    step) through six BMI inputs, time_year ... time_second.
    The AORC BMI reports time only as seconds since the start of
    its run (get_current_time) and does not expose the calendar
    date, so the driver reads the first time stamp of the AORC
    forcing file itself and adds the AORC current time to it.
    Calendar arithmetic is done here in UTC (days from civil and
    back, proleptic Gregorian), independent of the local time
    zone. Methods 1-5 do not expose the time inputs, so the
    driver sets them only when PET lists them as inputs.
***************************************************************/
#define DRIVER_MAXIMUM_LINE_CHARACTERS 4096
static const char *time_bmi_names[6] = {
    "time_year", "time_month", "time_day_of_month", "time_hour", "time_minute", "time_second" };

// days since 1970-01-01 of a proleptic Gregorian date
static long days_from_civil(long year, long month, long day_of_month)
{
    long era, year_of_era, day_of_year, day_of_era;
    year -= (month <= 2);
    era = (year >= 0 ? year : year - 399) / 400;
    year_of_era = year - era * 400;
    day_of_year = (153 * (month + (month > 2 ? -3 : 9)) + 2) / 5 + day_of_month - 1;
    day_of_era = year_of_era * 365 + year_of_era / 4 - year_of_era / 100 + day_of_year;
    return era * 146097 + day_of_era - 719468;
}

// proleptic Gregorian date of a number of days since 1970-01-01
static void civil_from_days(long days_since_1970, long *year, long *month, long *day_of_month)
{
    long era, day_of_era, year_of_era, day_of_year, month_index;
    days_since_1970 += 719468;
    era = (days_since_1970 >= 0 ? days_since_1970 : days_since_1970 - 146096) / 146097;
    day_of_era = days_since_1970 - era * 146097;
    year_of_era = (day_of_era - day_of_era / 1460 + day_of_era / 36524 - day_of_era / 146096) / 365;
    day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    month_index = (5 * day_of_year + 2) / 153;
    *day_of_month = day_of_year - (153 * month_index + 2) / 5 + 1;
    *month = month_index + (month_index < 10 ? 3 : -9);
    *year = year_of_era + era * 400 + (*month <= 2);
}

// Value of one key=value entry of a configuration file ('#' starts a comment).
static int read_configuration_key(const char *configuration_file, const char *key, char *value_out, size_t value_size)
{
    char line_text[DRIVER_MAXIMUM_LINE_CHARACTERS];
    int found = 0;
    FILE *file_pointer = fopen(configuration_file, "r");
    if (file_pointer == NULL) return -1;
    while (fgets(line_text, sizeof(line_text), file_pointer) != NULL) {
        char *comment = strchr(line_text, '#'), *equals, *key_start, *key_end, *value_start, *value_end;
        if (comment != NULL) *comment = '\0';
        equals = strchr(line_text, '=');
        if (equals == NULL) continue;
        *equals = '\0';
        for (key_start = line_text; *key_start == ' ' || *key_start == '\t'; ++key_start) ;
        for (key_end = equals; key_end > key_start && (key_end[-1] == ' ' || key_end[-1] == '\t'); --key_end) ;
        *key_end = '\0';
        if (strcmp(key_start, key) != 0) continue;
        for (value_start = equals + 1; *value_start == ' ' || *value_start == '\t'; ++value_start) ;
        for (value_end = value_start + strlen(value_start);
             value_end > value_start && (value_end[-1] == ' ' || value_end[-1] == '\t' || value_end[-1] == '\n' || value_end[-1] == '\r');
             --value_end) ;
        *value_end = '\0';
        snprintf(value_out, value_size, "%s", value_start);
        found = 1;
    }
    fclose(file_pointer);
    return found ? 0 : -1;
}

// Seconds since 1970-01-01 00:00:00 UTC of the first time stamp in the AORC forcing file.
static int read_aorc_start_time_s(const char *aorc_configuration_file, double *start_time_since_1970_s)
{
    char forcing_file[1024], header_line[DRIVER_MAXIMUM_LINE_CHARACTERS], data_line[DRIVER_MAXIMUM_LINE_CHARACTERS];
    int year, month, day_of_month, hour, minute;
    double second;
    FILE *file_pointer;
    if (read_configuration_key(aorc_configuration_file, "forcing_file", forcing_file, sizeof(forcing_file)) != 0) {
        fprintf(stderr, "no forcing_file in %s\n", aorc_configuration_file);
        return -1;
    }
    file_pointer = fopen(forcing_file, "r");
    if (file_pointer == NULL || fgets(header_line, sizeof(header_line), file_pointer) == NULL ||
        fgets(data_line, sizeof(data_line), file_pointer) == NULL ||
        sscanf(data_line, "%d-%d-%d %d:%d:%lf", &year, &month, &day_of_month, &hour, &minute, &second) != 6) {
        fprintf(stderr, "cannot read the first time stamp (YYYY-MM-DD HH:MM:SS) of %s\n", forcing_file);
        if (file_pointer != NULL) fclose(file_pointer);
        return -1;
    }
    fclose(file_pointer);
    *start_time_since_1970_s = 86400.0 * (double)days_from_civil(year, month, day_of_month)
                             + 3600.0 * hour + 60.0 * minute + second;
    return 0;
}

// Set the six PET time inputs to a time given in seconds since 1970-01-01 00:00:00 UTC.
static int set_pet_calendar_time(Bmi *pet_bmi_model, double time_since_1970_s)
{
    double days_since_1970 = floor(time_since_1970_s / 86400.0);
    double second_of_day = time_since_1970_s - 86400.0 * days_since_1970;
    long year, month, day_of_month;
    double time_values[6];
    civil_from_days((long)days_since_1970, &year, &month, &day_of_month);
    time_values[0] = (double)year;
    time_values[1] = (double)month;
    time_values[2] = (double)day_of_month;
    time_values[3] = floor(second_of_day / 3600.0);
    time_values[4] = floor((second_of_day - 3600.0 * time_values[3]) / 60.0);
    time_values[5] = second_of_day - 3600.0 * time_values[3] - 60.0 * time_values[4];
    for (int i = 0; i < 6; ++i) {
        if (pet_bmi_model->set_value(pet_bmi_model, time_bmi_names[i], &time_values[i]) != BMI_SUCCESS) {
            fprintf(stderr, "PET BMI set_value failed: %s = %.10g\n", time_bmi_names[i], time_values[i]);
            return BMI_FAILURE;
        }
    }
    return BMI_SUCCESS;
}

// 1 if PET lists the variable among its inputs (method 6 lists the time inputs, methods 1-5 do not)
static int pet_has_input(Bmi *pet_bmi_model, const char *variable_name)
{
    char type[BMI_MAX_TYPE_NAME];
    int input_count = 0, has_input = 0;
    char **input_names;
    if (pet_bmi_model->get_var_type(pet_bmi_model, variable_name, type) != BMI_SUCCESS) return 0;
    if (pet_bmi_model->get_input_item_count(pet_bmi_model, &input_count) != BMI_SUCCESS || input_count <= 0) return 0;
    input_names = (char **)malloc(sizeof(char *) * input_count);
    if (input_names == NULL) return 0;
    for (int i = 0; i < input_count; ++i) input_names[i] = (char *)malloc(BMI_MAX_VAR_NAME);
    if (pet_bmi_model->get_input_var_names(pet_bmi_model, input_names) == BMI_SUCCESS)
        for (int i = 0; i < input_count; ++i)
            if (strcmp(input_names[i], variable_name) == 0) has_input = 1;
    for (int i = 0; i < input_count; ++i) free(input_names[i]);
    free(input_names);
    return has_input;
}

/************************************************************************
    This main program is a mock framwork.
    This is not part of BMI, but acts as the driver that calls the model.
************************************************************************/
int main(int argc, const char *argv[])
{
    Bmi *pet_bmi_model = NULL;
    Bmi *aorc_bmi_model = NULL;
    pet_model *pet = NULL;
    int result = EXIT_FAILURE;

    if (argc < 3) {
        fprintf(stderr, "Usage: %s PET_config AORC_config\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* FLO 10/26: Check allocation, registration and initialization before
     * accessing BMI model data; previously failures were ignored.
     */
    pet_bmi_model = calloc(1, sizeof(*pet_bmi_model));
    aorc_bmi_model = calloc(1, sizeof(*aorc_bmi_model));
    if (!pet_bmi_model || !aorc_bmi_model) goto cleanup;
    if (!register_bmi_pet(pet_bmi_model) || !pet_bmi_model->data) goto cleanup;
    if (!register_bmi_aorc(aorc_bmi_model) || !aorc_bmi_model->data) goto cleanup;
    if (pet_bmi_model->initialize(pet_bmi_model, argv[1]) != BMI_SUCCESS) {
        fprintf(stderr, "PET initialization failed\n");
        goto cleanup;
    }
    if (aorc_bmi_model->initialize(aorc_bmi_model, argv[2]) != BMI_SUCCESS) {
        fprintf(stderr, "AORC initialization failed\n");
        goto cleanup;
    }
    pet = (pet_model *)pet_bmi_model->data;

    /* FLO 10/21: method 6 needs the calendar time (see set_pet_calendar_time). */
    int pet_takes_calendar_time = pet_has_input(pet_bmi_model, "time_year");
    double aorc_start_time_since_1970_s = 0.0;
    if (pet_takes_calendar_time && read_aorc_start_time_s(argv[2], &aorc_start_time_since_1970_s) != 0)
        goto cleanup;

    /* FLO 10/26: Check both model updates as well as BMI transfers. A failed
     * AORC update must never lead to reusing stale forcing in PET.
     */
    for (int i = 0; i <= 100; ++i) {
        /* FLO 10/21: the AORC update about to run supplies the forcing of the step that starts at
         * the AORC current time (seconds since the start of its run). */
        double aorc_time_since_start_s = 0.0;
        if (pet_takes_calendar_time) {
            if (aorc_bmi_model->get_current_time(aorc_bmi_model, &aorc_time_since_start_s) != BMI_SUCCESS) {
                fprintf(stderr, "AORC get_current_time failed at step %d\n", i);
                goto cleanup;
            }
        }
        if (aorc_bmi_model->update(aorc_bmi_model) != BMI_SUCCESS) {
            fprintf(stderr, "AORC update failed at step %d\n", i);
            goto cleanup;
        }
        if (pass_forcing_from_aorc_to_pet(pet_bmi_model, aorc_bmi_model) != BMI_SUCCESS)
            goto cleanup;
        if (pet_takes_calendar_time &&
            set_pet_calendar_time(pet_bmi_model, aorc_start_time_since_1970_s + aorc_time_since_start_s) != BMI_SUCCESS)
            goto cleanup;
        if (pet_bmi_model->update(pet_bmi_model) != BMI_SUCCESS) {
            fprintf(stderr, "PET update failed at step %d\n", i);
            goto cleanup;
        }
        if (i == 0) continue; /* Preserve the driver's original output count. */
        if (pet_takes_calendar_time) {   /* FLO 10/21: the calendar time PET method 6 received */
            double time_values[6] = { 0.0 };
            for (int t = 0; t < 6; ++t) pet_bmi_model->get_value(pet_bmi_model, time_bmi_names[t], &time_values[t]);
            printf("calendar time (UTC, start of step) %04.0f-%02.0f-%02.0f %02.0f:%02.0f:%02.0f\n",
                   time_values[0], time_values[1], time_values[2], time_values[3], time_values[4], time_values[5]);
        }
        printf("LWDOWN after set value %lf\n", pet->aorc.incoming_longwave_W_per_m2);
        printf("SWDOWN before set value %lf\n", pet->aorc.incoming_shortwave_W_per_m2);
        printf("surface_pressure_Pa %lf \n", pet->aorc.surface_pressure_Pa);
        printf("specific_humidity_2m_kg_per_kg %lf \n", pet->aorc.specific_humidity_2m_kg_per_kg);
        printf("air_temperature_2m_K %lf \n", pet->aorc.air_temperature_2m_K);
        printf("u_wind_speed_10m_m_per_s %lf \n", pet->aorc.u_wind_speed_10m_m_per_s);
        printf("v_wind_speed_10m_m_per_s %lf \n", pet->aorc.v_wind_speed_10m_m_per_s);
        printf("potential evapotranspiration (m s-1): %8.4e\n", pet->pet_m_per_s);
        printf("--------------------- END OF TIMESTEP ---------------------------------\n");
    }
    result = EXIT_SUCCESS;

cleanup:
    /* FLO 10/26: The old driver returned early without finalizing either
     * model or freeing BMI wrappers. Clean up on success and error alike.
     */
    if (pet_bmi_model) {
        if (pet_bmi_model->data && pet_bmi_model->finalize) {
            if (pet_bmi_model->finalize(pet_bmi_model) != BMI_SUCCESS)
                result = EXIT_FAILURE;
        }
        free(pet_bmi_model);
    }
    if (aorc_bmi_model) {
        if (aorc_bmi_model->data && aorc_bmi_model->finalize) {
            if (aorc_bmi_model->finalize(aorc_bmi_model) != BMI_SUCCESS)
                result = EXIT_FAILURE;
        }
        free(aorc_bmi_model);
    }
    if (result == EXIT_SUCCESS) printf("Finalizing models\n");
    return result;
}
