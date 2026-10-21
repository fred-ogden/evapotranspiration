// =========================================================================================
// main_unit_test_bmi.c        BMI unit test of the PET module, any of the six methods
//
// Usage:
//   run_pet_bmi_test configuration_file
//
// The configuration file is one of configs/method_N_*.cfg, with forcing_file = BMI. The test
// reads forcing_csv_file from that file itself, as a model engine would read its forcing, and
// passes the first data row to PET with set_value.
//
// The test calls every BMI function and checks what it returns, in this order:
//   1. initialize
//   2. model information      component name, input and output item counts and names
//   3. variable information   grid, itemsize, location, units, type, nbytes of every variable
//   4. grid information       rank, size, type of grid 0
//   5. time information       start, end, step, units, current time
//   6. set_value              seven forcing variables and, for method 6, six time variables;
//                             get_value must return what was set. Methods 1-5 expose 7 inputs
//                             and 1 output; method 6 exposes 31 inputs and 9 outputs
//   7. update                 must succeed and advance the current time by one time step;
//                             every output must be finite through get_value,
//                             get_value_at_indices and get_value_ptr; methods 1-5 with the unit
//                             test forcing must give the unit test PET values
//   8. set_value checks       every input accepts its own current value through set_value and
//                             set_value_at_indices; index 1 is rejected; outputs and unknown
//                             names are rejected; out-of-range method 6 inputs are rejected
//                             and leave the value unchanged
//   9. update_until           five more steps; the current time must reach the target
//  10. finalize
//
// It prints one line per failed check and a summary, and returns 0 only if every check passed.
// =========================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/pet.h"
#include "../include/bmi.h"
#include "../include/bmi_pet.h"

#define TEST_MAXIMUM_LINE_CHARACTERS 4096
#define TEST_MAXIMUM_COLUMNS 64
#define NUMBER_OF_FORCING_VARIABLES 7
#define NUMBER_OF_TIME_VARIABLES 6

static int number_of_checks = 0;
static int number_of_failed_checks = 0;

// Record one check; print it only when it fails.
static void check(int condition_is_true, const char *description, const char *variable_name)
{
    ++number_of_checks;
    if (!condition_is_true) {
        ++number_of_failed_checks;
        printf("  FAIL: %s%s%s\n", description, variable_name ? ": " : "", variable_name ? variable_name : "");
    }
}

// AORC forcing columns (a column matches NAME or NAME_suffix, e.g. TMP_2maboveground)
// and the PET BMI input variables they feed.
static const char *forcing_column_names[NUMBER_OF_FORCING_VARIABLES] = { "DLWRF", "DSWRF", "PRES", "SPFH", "TMP", "UGRD", "VGRD" };
static const char *forcing_bmi_names[NUMBER_OF_FORCING_VARIABLES] = {
    "land_surface_radiation~incoming~longwave__energy_flux",
    "land_surface_radiation~incoming~shortwave__energy_flux",
    "land_surface_air__pressure",
    "atmosphere_air_water~vapor__specific_humidity",
    "land_surface_air__temperature",
    "land_surface_wind__x_component_of_velocity",
    "land_surface_wind__y_component_of_velocity" };
static const char *time_bmi_names[NUMBER_OF_TIME_VARIABLES] = {
    "time_year", "time_month", "time_day_of_month", "time_hour", "time_minute", "time_second" };

// Expected water_potential_evaporation_flux [m s-1] for the unit test forcing
// (forcing/unit_test_forcing.csv) and site, methods 1-5: the values pet_unit_tests() prints.
static const double expected_unit_test_pet_m_per_s[6] = { 0.0, 8.594743e-08, 8.977490e-08, 8.694909e-08, 8.249098e-08, 7.628671e-08 };
#define EXPECTED_PET_RELATIVE_TOLERANCE 1.0e-6   // the expected values carry 7 significant digits

// Out-of-range values that set_value must reject (method 6 inputs only; methods 1-5 store
// the forcing as given and validate it in update).
typedef struct { const char *name; double rejected_value; } RejectedInput;
static const RejectedInput method6_rejected_inputs[] = {
    { "LAI", -1.0 }, { "SAI", -1.0 }, { "SHDFAC", 1.2 }, { "CWP", -0.1 }, { "MP", -2.0 }, { "VCMX25", -3.0 },
    { "root_stress_factor_0_to_1", 1.5 }, { "snow_cover_fraction", -0.1 }, { "canopy_wetted_fraction", 1.1 },
    { "snow_depth_m", -1.0 }, { "time_month", 13.0 }, { "time_hour", 12.5 }, { "time_day_of_month", 0.0 } };

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
static int read_configuration_key(const char *configuration_file, const char *key, char *value_out, size_t value_size)
{
    FILE *file_pointer = fopen(configuration_file, "r");
    char line_text[TEST_MAXIMUM_LINE_CHARACTERS];
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

static int column_matches(const char *column_name, const char *forcing_name)
{
    size_t length = strlen(forcing_name);
    return strncmp(column_name, forcing_name, length) == 0 && (column_name[length] == '\0' || column_name[length] == '_');
}

// First data row of the forcing file: seven forcing values and the six time values.
static int read_first_forcing_row(const char *forcing_csv_file, double forcing_values[NUMBER_OF_FORCING_VARIABLES],
                                  double time_values[NUMBER_OF_TIME_VARIABLES])
{
    char header_line[TEST_MAXIMUM_LINE_CHARACTERS], data_line[TEST_MAXIMUM_LINE_CHARACTERS];
    char *header_fields[TEST_MAXIMUM_COLUMNS], *data_fields[TEST_MAXIMUM_COLUMNS];
    int header_field_count, data_field_count, variable, column, time_column_index = -1;
    int year, month, day_of_month, hour, minute;
    double second;
    FILE *file_pointer = fopen(forcing_csv_file, "r");
    if (file_pointer == NULL || fgets(header_line, sizeof(header_line), file_pointer) == NULL ||
        fgets(data_line, sizeof(data_line), file_pointer) == NULL) {
        fprintf(stderr, "ERROR: cannot read a header and a data row from %s\n", forcing_csv_file);
        if (file_pointer != NULL) fclose(file_pointer);
        return -1;
    }
    fclose(file_pointer);
    header_field_count = split_csv_line(header_line, header_fields, TEST_MAXIMUM_COLUMNS);
    data_field_count = split_csv_line(data_line, data_fields, TEST_MAXIMUM_COLUMNS);
    if (data_field_count != header_field_count) {
        fprintf(stderr, "ERROR: the first data row of %s has %d fields, the header has %d\n",
                forcing_csv_file, data_field_count, header_field_count);
        return -1;
    }
    for (variable = 0; variable < NUMBER_OF_FORCING_VARIABLES; ++variable) {
        int column_index = -1;
        for (column = 0; column < header_field_count; ++column)
            if (column_matches(header_fields[column], forcing_column_names[variable])) column_index = column;
        if (column_index < 0) {
            fprintf(stderr, "ERROR: no column %s in %s\n", forcing_column_names[variable], forcing_csv_file);
            return -1;
        }
        forcing_values[variable] = strtod(data_fields[column_index], NULL);
    }
    for (column = 0; column < header_field_count; ++column)
        if (strcmp(header_fields[column], "time") == 0) time_column_index = column;
    if (time_column_index < 0 ||
        sscanf(data_fields[time_column_index], "%d-%d-%d %d:%d:%lf", &year, &month, &day_of_month, &hour, &minute, &second) != 6) {
        fprintf(stderr, "ERROR: no time stamp YYYY-MM-DD HH:MM:SS in column time of %s\n", forcing_csv_file);
        return -1;
    }
    time_values[0] = year; time_values[1] = month; time_values[2] = day_of_month;
    time_values[3] = hour; time_values[4] = minute; time_values[5] = second;
    return 0;
}

// Set one double input and check that get_value returns it.
static void set_and_check(Bmi *model, const char *name, double value)
{
    double value_returned = NAN;
    check(model->set_value(model, name, &value) == BMI_SUCCESS, "set_value", name);
    check(model->get_value(model, name, &value_returned) == BMI_SUCCESS && value_returned == value,
          "get_value returns the value set", name);
}

// grid, itemsize, location, units, type and nbytes of one variable
static void check_variable_information(Bmi *model, const char *name)
{
    int grid = -1, itemsize = 0, nbytes = 0;
    char type[BMI_MAX_TYPE_NAME], location[BMI_MAX_LOCATION_NAME], units[BMI_MAX_UNITS_NAME];
    check(model->get_var_grid(model, name, &grid) == BMI_SUCCESS && grid == 0, "get_var_grid is 0", name);
    check(model->get_var_type(model, name, type) == BMI_SUCCESS && strcmp(type, "double") == 0, "get_var_type is double", name);
    check(model->get_var_itemsize(model, name, &itemsize) == BMI_SUCCESS && itemsize == (int)sizeof(double), "get_var_itemsize", name);
    check(model->get_var_nbytes(model, name, &nbytes) == BMI_SUCCESS && nbytes == itemsize, "get_var_nbytes is one item", name);
    check(model->get_var_location(model, name, location) == BMI_SUCCESS && strlen(location) > 0, "get_var_location", name);
    check(model->get_var_units(model, name, units) == BMI_SUCCESS, "get_var_units", name);
}

int main(int argc, const char *argv[])
{
    const char *configuration_file;
    char forcing_csv_file[1024], method_text[16], component_name[BMI_MAX_COMPONENT_NAME];
    double forcing_values[NUMBER_OF_FORCING_VARIABLES], time_values[NUMBER_OF_TIME_VARIABLES];
    double start_time_s = NAN, end_time_s = NAN, time_step_s = NAN, current_time_s = NAN, target_time_s;
    char time_units[BMI_MAX_UNITS_NAME], grid_type[BMI_MAX_COMPONENT_NAME];
    int input_item_count = 0, output_item_count = 0, grid_rank = -1, grid_size = 0, variable, pet_method;
    char **input_names = NULL, **output_names = NULL;
    Bmi *model;

    if (argc != 2) {
        fprintf(stderr, "usage: %s configs/method_N_*.cfg   (forcing_file = BMI)\n", argv[0]);
        return 1;
    }
    configuration_file = argv[1];
    if (read_configuration_key(configuration_file, "forcing_csv_file", forcing_csv_file, sizeof(forcing_csv_file)) != 0 ||
        read_configuration_key(configuration_file, "pet_method", method_text, sizeof(method_text)) != 0) {
        fprintf(stderr, "ERROR: forcing_csv_file and pet_method must be in %s\n", configuration_file);
        return 1;
    }
    pet_method = atoi(method_text);
    if (read_first_forcing_row(forcing_csv_file, forcing_values, time_values) != 0) return 1;
    printf("BMI unit test: %s (method %d, forcing %s)\n", configuration_file, pet_method, forcing_csv_file);

    // ---- 1. initialize
    model = (Bmi *)malloc(sizeof(Bmi));
    if (model == NULL) return 1;
    register_bmi_pet(model);
    if (model->initialize(model, configuration_file) != BMI_SUCCESS) {
        printf("  FAIL: initialize\n");
        if (model->data != NULL) model->finalize(model);
        free(model);
        return 1;
    }

    // ---- 2. model information
    check(model->get_component_name(model, component_name) == BMI_SUCCESS && strlen(component_name) > 0, "get_component_name", NULL);
    check(model->get_input_item_count(model, &input_item_count) == BMI_SUCCESS && input_item_count > 0, "get_input_item_count", NULL);
    check(model->get_output_item_count(model, &output_item_count) == BMI_SUCCESS && output_item_count > 0, "get_output_item_count", NULL);
    input_names = (char **)malloc(sizeof(char *) * input_item_count);
    output_names = (char **)malloc(sizeof(char *) * output_item_count);
    for (variable = 0; variable < input_item_count; ++variable) input_names[variable] = (char *)malloc(BMI_MAX_VAR_NAME);
    for (variable = 0; variable < output_item_count; ++variable) output_names[variable] = (char *)malloc(BMI_MAX_VAR_NAME);
    check(model->get_input_var_names(model, input_names) == BMI_SUCCESS, "get_input_var_names", NULL);
    check(model->get_output_var_names(model, output_names) == BMI_SUCCESS, "get_output_var_names", NULL);
    printf("  component %s: %d inputs, %d outputs\n", component_name, input_item_count, output_item_count);

    // ---- 3. variable information
    for (variable = 0; variable < input_item_count; ++variable) check_variable_information(model, input_names[variable]);
    for (variable = 0; variable < output_item_count; ++variable) check_variable_information(model, output_names[variable]);
    {
        char type[BMI_MAX_TYPE_NAME];
        check(model->get_var_type(model, "no_such_variable", type) == BMI_FAILURE, "get_var_type rejects an unknown name", NULL);
    }

    // ---- 4. grid information
    check(model->get_grid_rank(model, 0, &grid_rank) == BMI_SUCCESS, "get_grid_rank", NULL);
    check(model->get_grid_size(model, 0, &grid_size) == BMI_SUCCESS && grid_size == 1, "get_grid_size is 1", NULL);
    check(model->get_grid_type(model, 0, grid_type) == BMI_SUCCESS, "get_grid_type", NULL);

    // ---- 5. time information
    check(model->get_start_time(model, &start_time_s) == BMI_SUCCESS, "get_start_time", NULL);
    check(model->get_end_time(model, &end_time_s) == BMI_SUCCESS && end_time_s >= start_time_s, "get_end_time >= start time", NULL);
    check(model->get_time_step(model, &time_step_s) == BMI_SUCCESS && time_step_s > 0.0, "get_time_step > 0", NULL);
    check(model->get_time_units(model, time_units) == BMI_SUCCESS && strlen(time_units) > 0, "get_time_units", NULL);
    check(model->get_current_time(model, &current_time_s) == BMI_SUCCESS && current_time_s == start_time_s,
          "get_current_time is the start time", NULL);

    // ---- 6. forcing and time through set_value
    for (variable = 0; variable < NUMBER_OF_FORCING_VARIABLES; ++variable)
        set_and_check(model, forcing_bmi_names[variable], forcing_values[variable]);
    for (variable = 0; variable < NUMBER_OF_TIME_VARIABLES; ++variable) {
        // the time inputs exist only where the method uses them (method 6)
        int name_index, is_an_input = 0;
        for (name_index = 0; name_index < input_item_count; ++name_index)
            if (strcmp(input_names[name_index], time_bmi_names[variable]) == 0) is_an_input = 1;
        if (is_an_input) set_and_check(model, time_bmi_names[variable], time_values[variable]);
        else check(model->set_value(model, time_bmi_names[variable], &time_values[variable]) == BMI_FAILURE,
                   "set_value rejects a time input this method does not expose", time_bmi_names[variable]);
    }

    // ---- 7. update and outputs
    check(model->update(model) == BMI_SUCCESS, "update", NULL);
    check(model->get_current_time(model, &current_time_s) == BMI_SUCCESS && current_time_s == start_time_s + time_step_s,
          "update advances the current time by one time step", NULL);
    for (variable = 0; variable < output_item_count; ++variable) {
        double value = NAN, value_at_index = NAN, *value_pointer = NULL;
        int index_zero = 0;
        const char *name = output_names[variable];
        check(model->get_value(model, name, &value) == BMI_SUCCESS && isfinite(value), "get_value is finite", name);
        check(model->get_value_at_indices(model, name, &value_at_index, &index_zero, 1) == BMI_SUCCESS && value_at_index == value,
              "get_value_at_indices equals get_value", name);
        check(model->get_value_ptr(model, name, (void **)&value_pointer) == BMI_SUCCESS && value_pointer != NULL && *value_pointer == value,
              "get_value_ptr points at the value", name);
        printf("  %-50s %.10e\n", name, value);
        if (pet_method >= 1 && pet_method <= 5 && strcmp(name, "water_potential_evaporation_flux") == 0 &&
            strstr(forcing_csv_file, "unit_test_forcing.csv") != NULL) {
            double expected_value = expected_unit_test_pet_m_per_s[pet_method];
            check(fabs(value - expected_value) <= EXPECTED_PET_RELATIVE_TOLERANCE * expected_value,
                  "PET equals the unit test value", name);
        }
    }

    // ---- 8. set_value checks
    for (variable = 0; variable < input_item_count; ++variable) {
        double current_value = NAN, value_returned = NAN;
        int index_zero = 0, index_one = 1;
        const char *name = input_names[variable];
        check(model->get_value(model, name, &current_value) == BMI_SUCCESS, "get_value", name);
        if (isnan(current_value)) continue;   // not yet set (e.g. ground_skin_temperature_K): nothing to write back
        check(model->set_value(model, name, &current_value) == BMI_SUCCESS, "set_value accepts the current value", name);
        check(model->set_value_at_indices(model, name, &index_zero, 1, &current_value) == BMI_SUCCESS,
              "set_value_at_indices accepts index 0", name);
        check(model->set_value_at_indices(model, name, &index_one, 1, &current_value) == BMI_FAILURE,
              "set_value_at_indices rejects index 1", name);
        check(model->get_value(model, name, &value_returned) == BMI_SUCCESS && value_returned == current_value,
              "value unchanged", name);
    }
    for (variable = 0; variable < output_item_count; ++variable) {
        double any_value = 1.0;
        check(model->set_value(model, output_names[variable], &any_value) == BMI_FAILURE, "set_value rejects an output", output_names[variable]);
    }
    {
        double any_value = 1.0;
        check(model->set_value(model, "no_such_variable", &any_value) == BMI_FAILURE, "set_value rejects an unknown name", NULL);
    }
    if (pet_method == 6) {
        for (variable = 0; variable < (int)(sizeof(method6_rejected_inputs) / sizeof(method6_rejected_inputs[0])); ++variable) {
            const char *name = method6_rejected_inputs[variable].name;
            double value_before = NAN, value_after = NAN, rejected_value = method6_rejected_inputs[variable].rejected_value;
            double not_a_number = NAN;
            model->get_value(model, name, &value_before);
            check(model->set_value(model, name, &rejected_value) == BMI_FAILURE, "set_value rejects an out-of-range value", name);
            check(model->set_value(model, name, &not_a_number) == BMI_FAILURE, "set_value rejects NaN", name);
            check(model->get_value(model, name, &value_after) == BMI_SUCCESS &&
                  (value_after == value_before || (isnan(value_after) && isnan(value_before))),
                  "a rejected value leaves the value unchanged", name);
        }
    }

    // ---- 9. update_until: five more steps with the same forcing
    target_time_s = current_time_s + 5.0 * time_step_s;
    check(model->update_until(model, target_time_s) == BMI_SUCCESS, "update_until", NULL);
    check(model->get_current_time(model, &current_time_s) == BMI_SUCCESS && current_time_s == target_time_s,
          "update_until reaches the target time", NULL);

    // ---- 10. finalize
    for (variable = 0; variable < input_item_count; ++variable) free(input_names[variable]);
    for (variable = 0; variable < output_item_count; ++variable) free(output_names[variable]);
    free(input_names);
    free(output_names);
    check(model->finalize(model) == BMI_SUCCESS, "finalize", NULL);
    free(model);

    printf("%s: %d of %d checks passed (%s)\n", number_of_failed_checks == 0 ? "PASS" : "FAIL",
           number_of_checks - number_of_failed_checks, number_of_checks, configuration_file);
    return number_of_failed_checks == 0 ? 0 : 1;
}
