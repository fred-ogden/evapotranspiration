// =========================================================================================
// main_ball_berry_bmi.c        test driver for PET method 6 (Ball-Berry), BMI only
//
// FLO 10/21: drives the PET BMI module exactly as a model engine would, using only BMI
// functions: initialize, set_value, update, get_value, finalize. The driver, not the model,
// reads the forcing file and supplies the forcing and the calendar time through set_value.
//
// Usage:
//   test_ball_berry_method6 config_file
//
// The configuration file is the shared C/Fortran configuration. The driver reads two keys
// from it itself:
//   forcing_csv_file  forcing CSV: a time stamp column "time" (YYYY-MM-DD HH:MM:SS, UTC, start of
//                     step) and AORC columns DLWRF, DSWRF, PRES, SPFH, TMP, UGRD, VGRD
//   num_timesteps     number of forcing rows to run (from the first data row)
// PET itself is configured with forcing_file = BMI.
//
// For each step it sets the seven AORC forcing variables and the six time variables, calls
// update, and prints every method 6 BMI output obtained with get_value. The coupling inputs keep
// the initial values the model read from the configuration (a stand-alone run).
// =========================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bmi.h"
#include "bmi_pet.h"

#define DRIVER_MAXIMUM_LINE_CHARACTERS 4096
#define DRIVER_MAXIMUM_COLUMNS 64

// AORC forcing columns and the PET BMI input variables they feed
static const char *forcing_column_names[7] = { "DLWRF", "DSWRF", "PRES", "SPFH", "TMP", "UGRD", "VGRD" };
static const char *forcing_bmi_names[7] = {
    "land_surface_radiation~incoming~longwave__energy_flux",
    "land_surface_radiation~incoming~shortwave__energy_flux",
    "land_surface_air__pressure",
    "atmosphere_air_water~vapor__specific_humidity",
    "land_surface_air__temperature",
    "land_surface_wind__x_component_of_velocity",
    "land_surface_wind__y_component_of_velocity" };

// method 6 outputs printed after each update
static const char *method6_output_names[13] = {
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
    "canopy_heat_storage_flux_W_per_m2" };

static void print_usage(const char *program_name)
{
    fprintf(stderr,
            "usage: %s config_file\n"
            "  config_file  shared configuration (configs/ball_berry_noahmp_comparison_config.txt);\n"
            "               PET must have forcing_file = BMI; the driver reads forcing_csv_file and num_timesteps.\n",
            program_name);
}

static char *trim_text(char *text)
{
    char *end;
    while (*text == ' ' || *text == '\t') ++text;
    end = text + strlen(text);
    while (end > text && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\n' || end[-1] == '\r')) --end;
    *end = '\0';
    return text;
}

// Value of one "key = value" entry of the configuration file ('#' starts a comment).
static int read_driver_key(const char *configuration_file, const char *key, char *value_out, size_t value_size)
{
    FILE *file_pointer = fopen(configuration_file, "r");
    char line_text[DRIVER_MAXIMUM_LINE_CHARACTERS];
    int found = 0;
    if (file_pointer == NULL) return -1;
    while (fgets(line_text, sizeof(line_text), file_pointer) != NULL) {
        char *comment = strchr(line_text, '#');
        char *equals;
        if (comment != NULL) *comment = '\0';
        equals = strchr(line_text, '=');
        if (equals == NULL) continue;
        *equals = '\0';
        if (strcmp(trim_text(line_text), key) != 0) continue;
        snprintf(value_out, value_size, "%s", trim_text(equals + 1));
        found = 1;
    }
    fclose(file_pointer);
    return found ? 0 : -1;
}

// Split a comma separated line in place.
static int split_csv_line(char *line_text, char **field_list, int maximum_fields)
{
    int field_count = 0;
    char *field_start = line_text;
    while (field_count < maximum_fields) {
        char *comma = strchr(field_start, ',');
        if (comma != NULL) *comma = '\0';
        field_list[field_count++] = trim_text(field_start);
        if (comma == NULL) break;
        field_start = comma + 1;
    }
    return field_count;
}

static int set_bmi_double(Bmi *model, const char *name, double value)
{
    if (model->set_value(model, name, &value) != BMI_SUCCESS) {
        fprintf(stderr, "ERROR: set_value failed for %s = %.10g\n", name, value);
        return -1;
    }
    return 0;
}

int main(int argc, const char *argv[])
{
    char forcing_csv_file[1024], number_text[64];
    char header_line[DRIVER_MAXIMUM_LINE_CHARACTERS], data_line[DRIVER_MAXIMUM_LINE_CHARACTERS];
    char *header_fields[DRIVER_MAXIMUM_COLUMNS], *data_fields[DRIVER_MAXIMUM_COLUMNS];
    int forcing_column_index[7], time_column_index = -1;
    int header_field_count, data_field_count, step, number_of_steps, column, variable, exit_status = 0;
    FILE *forcing_file_pointer;
    Bmi *model;

    if (argc != 2 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_usage(argv[0]);
        return 1;
    }
    if (read_driver_key(argv[1], "forcing_csv_file", forcing_csv_file, sizeof(forcing_csv_file)) != 0 ||
        read_driver_key(argv[1], "num_timesteps", number_text, sizeof(number_text)) != 0) {
        fprintf(stderr, "ERROR: forcing_csv_file and num_timesteps must be in %s\n", argv[1]);
        return 1;
    }
    number_of_steps = atoi(number_text);

    // ---- forcing file header: locate the columns by name
    forcing_file_pointer = fopen(forcing_csv_file, "r");
    if (forcing_file_pointer == NULL || fgets(header_line, sizeof(header_line), forcing_file_pointer) == NULL) {
        fprintf(stderr, "ERROR: cannot read forcing file %s\n", forcing_csv_file);
        if (forcing_file_pointer != NULL) fclose(forcing_file_pointer);
        return 1;
    }
    header_field_count = split_csv_line(header_line, header_fields, DRIVER_MAXIMUM_COLUMNS);
    for (variable = 0; variable < 7; ++variable) {
        forcing_column_index[variable] = -1;
        for (column = 0; column < header_field_count; ++column)
            if (strcmp(header_fields[column], forcing_column_names[variable]) == 0) forcing_column_index[variable] = column;
        if (forcing_column_index[variable] < 0) {
            fprintf(stderr, "ERROR: column %s not in %s\n", forcing_column_names[variable], forcing_csv_file);
            fclose(forcing_file_pointer);
            return 1;
        }
    }
    for (column = 0; column < header_field_count; ++column)
        if (strcmp(header_fields[column], "time") == 0) time_column_index = column;
    if (time_column_index < 0) {
        fprintf(stderr, "ERROR: column time not in %s\n", forcing_csv_file);
        fclose(forcing_file_pointer);
        return 1;
    }

    // ---- BMI initialize
    model = (Bmi *)malloc(sizeof(Bmi));
    if (model == NULL) { fclose(forcing_file_pointer); return 1; }
    register_bmi_pet(model);
    if (model->initialize(model, argv[1]) != BMI_SUCCESS) {
        fprintf(stderr, "ERROR: PET BMI initialize failed for %s\n", argv[1]);
        if (model->data != NULL) model->finalize(model);
        free(model);
        fclose(forcing_file_pointer);
        return 1;
    }

    // ---- time loop: set forcing and time, update, get outputs
    for (step = 0; step < number_of_steps; ++step) {
        int year, month, day_of_month, hour, minute;
        double second;
        if (fgets(data_line, sizeof(data_line), forcing_file_pointer) == NULL) {
            fprintf(stderr, "ERROR: forcing file %s has fewer than %d data rows\n", forcing_csv_file, number_of_steps);
            exit_status = 1;
            break;
        }
        data_field_count = split_csv_line(data_line, data_fields, DRIVER_MAXIMUM_COLUMNS);
        if (data_field_count != header_field_count) {
            fprintf(stderr, "ERROR: forcing row %d has %d fields, header has %d\n", step + 1, data_field_count, header_field_count);
            exit_status = 1;
            break;
        }
        // forcing, parsed in double precision
        for (variable = 0; variable < 7; ++variable)
            if (set_bmi_double(model, forcing_bmi_names[variable], strtod(data_fields[forcing_column_index[variable]], NULL)) != 0)
                { exit_status = 1; break; }
        if (exit_status != 0) break;
        // calendar time, UTC at the start of the step
        if (sscanf(data_fields[time_column_index], "%d-%d-%d %d:%d:%lf", &year, &month, &day_of_month, &hour, &minute, &second) != 6) {
            fprintf(stderr, "ERROR: cannot read the time stamp '%s'\n", data_fields[time_column_index]);
            exit_status = 1;
            break;
        }
        if (set_bmi_double(model, "time_year", (double)year) != 0 || set_bmi_double(model, "time_month", (double)month) != 0 ||
            set_bmi_double(model, "time_day_of_month", (double)day_of_month) != 0 || set_bmi_double(model, "time_hour", (double)hour) != 0 ||
            set_bmi_double(model, "time_minute", (double)minute) != 0 || set_bmi_double(model, "time_second", second) != 0)
            { exit_status = 1; break; }

        if (model->update(model) != BMI_SUCCESS) {
            fprintf(stderr, "ERROR: PET update failed at step %d\n", step + 1);
            exit_status = 1;
            break;
        }
        printf("step %d  %s\n", step + 1, data_fields[time_column_index]);
        for (variable = 0; variable < 13; ++variable) {
            double output_value;
            char units[BMI_MAX_UNITS_NAME];
            if (model->get_value(model, method6_output_names[variable], &output_value) != BMI_SUCCESS ||
                model->get_var_units(model, method6_output_names[variable], units) != BMI_SUCCESS) {
                fprintf(stderr, "ERROR: get_value failed for %s\n", method6_output_names[variable]);
                exit_status = 1;
                break;
            }
            printf("  %-50s %.12e %s\n", method6_output_names[variable], output_value, units);
        }
        if (exit_status != 0) break;
    }

    // ---- BMI finalize: the model frees everything it allocated; the driver frees the Bmi struct
    fclose(forcing_file_pointer);
    if (model->finalize(model) != BMI_SUCCESS) exit_status = 1;
    free(model);
    return exit_status;
}
