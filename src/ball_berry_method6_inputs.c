// =========================================================================================
// ball_berry_method6_inputs.c
//
// PET method 6 (Noah-MP Ball-Berry canopy energy balance and transpiration): configuration
// and per-step input construction.
//
//   bb6_initialize_from_configuration  called once by BMI Initialize
//       Parses the configuration file (the same file the Fortran NWM Noah-MP reference
//       reads) into pet->bb_parameters, the method 6 physiology into pet->bb_stomata, and
//       the initial values of the BMI input variables (LAI, SAI, SHDFAC, CWP, MP, VCMX25 and
//       the coupling inputs). The parsed key table is local to this call and freed before it
//       returns, so nothing here is static or shared between model instances.
//
//   bb6_prepare_step_inputs            FLO 10/28: called by the core step at every BMI update
//       The checks, solar geometry, LAI/SAI and NaN conventions of a step, into the step inputs of
//       the physics (ball_berry_canopy_physics.c, which builds the canopy input; differentiated).
//
//   bb6_build_canopy_input             the preparation and the canopy input physics (nested solver)
//       Builds the canopy solver input (pet->bb_stability_input) from
//         - the BMI forcing variables (pet->aorc, set by the driver with set_value),
//         - the BMI time variables (pet->bb_time_inputs), used only for the monthly LAI/SAI table,
//         - the BMI coupling inputs (pet->bb_coupling_inputs),
//         - the BMI vegetation variables (pet->bb_canopy, pet->bb_cwp, pet->bb_stomata),
//         - pet->bb_parameters.
//       It reads no files.
//
// Inside the science code below the Noah-MP variable names are used for provenance; the
// mapping from the long BMI names is made at the top of bb6_build_canopy_input.
//
// Everything is set the way ENERGY and VEGE_FLUX in the NWM module set it:
//   constants  CPAIR = 1004.64, RAIR = 287.04, HVAP = 2.5104E06, HSUB = 2.8440E06,
//              SB = 5.67E-08, TFRZ = 273.16, MPE = 1.E-6
//   forcing    EAIR = QAIR*SFCPRS/(0.622+0.378*QAIR), RHOAIR = (SFCPRS-0.378*EAIR)/(RAIR*SFCTMP),
//              UR = MAX(wind speed, 1 m/s), SOLAD/SOLAI by direct and visible fractions
//   phenology  monthly LAI/SAI interpolation (DVEG = 4), burial by snow -> ELAI, ESAI
//   geometry   ZPD = 0.65*HVT (raised to SNOWH by deep snow), ZLVL = MAX(ZPD,HVT) + ZREF,
//              Z0M = Z0MVT, Z0H = Z0M, Z0MG = Z0HG = Z0*(1-FSNO) + FSNO*Z0SNO, HCAN = HVT,
//              VAIE = MIN(6, ELAI+ESAI)
//   radiation  EMV = 1 - EXP(-(ELAI+ESAI)), EMG = EG*(1-FSNO) + FSNO,
//              ground albedo = soil*(1-FSNO) + snow*FSNO (GROUNDALB)
//   physiology CO2AIR = CO2*SFCPRS, O2AIR = O2*SFCPRS, IGS = 1 if TV > TMIN else 0, BTRAN input,
//              GAMMAV = CPAIR*SFCPRS/(0.622*LATHEAV), LATHEAV from TV vs TFRZ
//   state      cold start TV = TAH = SFCTMP, EAH = EAIR, or the configured initial values
//
// Limitation (checked): only radiation_transfer_option = 3 (NWM) is ported.
// =========================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pet.h"
#include "ball_berry_canopy_physics.h"

// ---- Noah-MP constants (module_sf_noahmplsm.F, MPTABLE.TBL)

// ---------------------------------------------------------------------------------------------
// Configuration key table, local to one bb6_initialize_from_configuration call.
// ---------------------------------------------------------------------------------------------
typedef struct {
    int    entry_count;
    char (*key_list)[128];
    char (*value_list)[256];
    int    lookup_failed;          // set when a required key is missing or malformed
} ConfigurationTable;

#define CONFIGURATION_MAXIMUM_ENTRIES 512
#define CONFIGURATION_MAXIMUM_LINE_CHARACTERS 2048

static char *trim_text(char *text)
{
    char *end;
    while (*text == ' ' || *text == '\t') ++text;
    end = text + strlen(text);
    while (end > text && (end[-1] == ' ' || end[-1] == '\t' || end[-1] == '\n' || end[-1] == '\r')) --end;
    *end = '\0';
    return text;
}

// Read every "key = value" line; '#' starts a comment anywhere on a line.
static int read_configuration_table(const char *file_path, ConfigurationTable *table)
{
    FILE *file_pointer;
    char line_text[CONFIGURATION_MAXIMUM_LINE_CHARACTERS];
    table->entry_count = 0;
    table->lookup_failed = 0;
    table->key_list = malloc(sizeof(*table->key_list) * CONFIGURATION_MAXIMUM_ENTRIES);
    table->value_list = malloc(sizeof(*table->value_list) * CONFIGURATION_MAXIMUM_ENTRIES);
    if (table->key_list == NULL || table->value_list == NULL) return -1;
    file_pointer = fopen(file_path, "r");
    if (file_pointer == NULL) {
        fprintf(stderr, "ERROR: cannot open configuration file %s\n", file_path);
        return -1;
    }
    while (fgets(line_text, sizeof(line_text), file_pointer) != NULL) {
        char *comment = strchr(line_text, '#');
        char *equals, *key, *value;
        if (comment != NULL) *comment = '\0';
        equals = strchr(line_text, '=');
        if (equals == NULL) continue;
        *equals = '\0';
        key = trim_text(line_text);
        value = trim_text(equals + 1);
        if (*key == '\0' || *value == '\0') continue;
        if (table->entry_count >= CONFIGURATION_MAXIMUM_ENTRIES) { fclose(file_pointer); return -1; }
        snprintf(table->key_list[table->entry_count], sizeof(table->key_list[0]), "%s", key);
        snprintf(table->value_list[table->entry_count], sizeof(table->value_list[0]), "%s", value);
        ++table->entry_count;
    }
    fclose(file_pointer);
    return 0;
}

static void free_configuration_table(ConfigurationTable *table)
{
    free(table->key_list);
    free(table->value_list);
    table->key_list = NULL;
    table->value_list = NULL;
}

// Required numeric value; records a failure (reported per key) if missing or malformed.
static double configuration_real(ConfigurationTable *table, const char *key)
{
    int entry_index;
    for (entry_index = 0; entry_index < table->entry_count; ++entry_index) {
        if (strcmp(table->key_list[entry_index], key) == 0) {
            char *end = NULL;
            double value = strtod(table->value_list[entry_index], &end);
            if (end == table->value_list[entry_index] || !isfinite(value)) break;
            return value;
        }
    }
    fprintf(stderr, "ERROR: method 6 configuration key missing or invalid: %s\n", key);
    table->lookup_failed = 1;
    return NAN;
}

// FLO 10/26: optional numeric value: the default if the key is absent (records a failure only if
// the key is present and malformed). Used for the keys added after the NWM replica configuration.
static double configuration_real_or_default(ConfigurationTable *table, const char *key, double default_value)
{
    int entry_index;
    for (entry_index = 0; entry_index < table->entry_count; ++entry_index) {
        if (strcmp(table->key_list[entry_index], key) == 0) return configuration_real(table, key);
    }
    return default_value;
}

// ---------------------------------------------------------------------------------------------
// FLO 10/26: Noah-MP calendar and solar geometry, identical to NoahMP_snow_ablation
// (noahmp_snow_surface_physics.c, a port of NOM UtilitiesModule calc_declin_components):
// proleptic Gregorian with the 3600-year exception; slope and aspect after Corripio (2003).
// ---------------------------------------------------------------------------------------------
static int noahmp_february_days(int year)
{
    int days = 28;
    if (year % 4 == 0) {
        days = 29;
        if (year % 100 == 0) {
            days = 28;
            if (year % 400 == 0) {
                days = 29;
                if (year % 3600 == 0) days = 28;
            }
        }
    }
    return days;
}

static int noahmp_day_of_year_offset(int year, int month, int day)
{
    static const int cumulative_days[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    int offset = cumulative_days[month - 1] + day - 1;
    if (month > 2) offset = offset + (noahmp_february_days(year) - 28);
    return offset;
}

static long noahmp_year_start_day(int year)
{
    long y = (long)year - 1;
    return 365 * y + y / 4 - y / 100 + y / 400 - y / 3600;
}

static void noahmp_advance_datetime(int year, int month, int day, int hour, int minute, long minutes,
                                    int *year2, int *month2, int *day2, int *hour2, int *minute2)
{
    static const int month_days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    long total = (long)hour * 60 + minute + minutes;
    long minute_of_day = ((total % 1440) + 1440) % 1440;
    long days = noahmp_year_start_day(year) + noahmp_day_of_year_offset(year, month, day) + (total - minute_of_day) / 1440;
    long day_of_year;
    int y, m, month_length;
    y = (int)(((double)days + 0.5) / 365.2425) + 1;
    while (noahmp_year_start_day(y + 1) <= days) ++y;
    while (noahmp_year_start_day(y) > days) --y;
    day_of_year = days - noahmp_year_start_day(y);
    m = 1;
    for (;;) {
        month_length = (m == 2) ? noahmp_february_days(y) : month_days[m - 1];
        if (day_of_year < month_length) break;
        day_of_year -= month_length;
        ++m;
    }
    *year2 = y;
    *month2 = m;
    *day2 = (int)day_of_year + 1;
    *hour2 = (int)(minute_of_day / 60);
    *minute2 = (int)(minute_of_day % 60);
}

static void noahmp_calc_declin_components(int iday, int ihour, int iminute, int isecond, double latitude,
                                          double longitude, double slope, double azimuth, double *COSZ,
                                          double *COSZ_HORIZ)
{
    const double DEGRAD = 3.14159265 / 180.;
    const double DPD = 360. / 365.;
    double JULIAN, HRANG, DECLIN, TLOCTIM, SINOB, SXLONG, nvx, nvy, nvz, svx, svy, svz;
    JULIAN = (double)iday + (double)ihour / 24.;
    SINOB = sin(23.5 * DEGRAD);
    if (JULIAN >= 80.) {
        SXLONG = DPD * (JULIAN - 80.) * DEGRAD;
    } else {
        SXLONG = DPD * (JULIAN + 285.) * DEGRAD;
    }
    DECLIN = asin(SINOB * sin(SXLONG));
    TLOCTIM = (double)ihour + (double)iminute / 60.0 + (double)isecond / 3600.0 + longitude / 15.0;
    TLOCTIM = fmod(TLOCTIM + 24.0, 24.0);
    HRANG = 15. * (TLOCTIM - 12.) * DEGRAD;
    nvx = sin(azimuth * DEGRAD) * sin(slope * DEGRAD);
    nvy = -cos(azimuth * DEGRAD) * sin(slope * DEGRAD);
    nvz = cos(slope * DEGRAD);
    svx = -sin(HRANG) * cos(DECLIN);
    svy = (sin(latitude * DEGRAD) * cos(HRANG) * cos(DECLIN)) - (cos(latitude * DEGRAD) * sin(DECLIN));
    svz = (cos(latitude * DEGRAD) * cos(HRANG) * cos(DECLIN)) + (sin(latitude * DEGRAD) * sin(DECLIN));
    *COSZ = (nvx * svx) + (nvy * svy) + (nvz * svz);
    *COSZ_HORIZ = svz;
}

// COSZ and COSZ_HORIZ of this step from the BMI time inputs (start of step, UTC). Returns 0 or -1.
static int bb6_solar_geometry(const pet_model *pet, double *COSZ, double *COSZ_HORIZ)
{
    const BallBerryMethod6Parameters *parameters = &pet->bb_parameters;
    const BallBerryTimeInputs *time_inputs = &pet->bb_time_inputs;
    int year, month, day, hour, minute;
    if (isnan(time_inputs->time_year) || isnan(time_inputs->time_month) || isnan(time_inputs->time_day_of_month) ||
        isnan(time_inputs->time_hour) || isnan(time_inputs->time_minute) || isnan(time_inputs->time_second)) {
        fprintf(stderr, "ERROR: PET method 6 with solar_geometry_option = 1 needs the BMI time inputs\n");
        return -1;
    }
    year = (int)time_inputs->time_year;
    month = (int)time_inputs->time_month;
    day = (int)time_inputs->time_day_of_month;
    hour = (int)time_inputs->time_hour;
    minute = (int)time_inputs->time_minute;
    if (parameters->solar_geometry_at_end_of_step == 1) {       // NOM UtilitiesMain: itime*(dt/60) minutes
        noahmp_advance_datetime(year, month, day, hour, minute, (long)(pet->bmi.time_step_size_s / 60.0),
                                &year, &month, &day, &hour, &minute);
    }
    noahmp_calc_declin_components(noahmp_day_of_year_offset(year, month, day), hour, minute,
                                  (int)time_inputs->time_second, parameters->latitude_degrees,
                                  parameters->longitude_degrees, parameters->terrain_slope_degrees,
                                  parameters->terrain_azimuth_degrees, COSZ, COSZ_HORIZ);
    return 0;
}

// ---------------------------------------------------------------------------------------------
// Range checks used at Initialize and by BMI set_value.
// ---------------------------------------------------------------------------------------------

// Coupling inputs. TG may be NaN (not yet supplied). Returns 0 if all are valid, -1 otherwise.
int bb6_validate_coupling_inputs(const BallBerryCouplingInputs *coupling)
{
    if (coupling == NULL) return -1;
    if (!(coupling->root_stress_factor_0_to_1 >= 0.0 && coupling->root_stress_factor_0_to_1 <= 1.0)) return -1;
    if (!isnan(coupling->ground_skin_temperature_K) && !(coupling->ground_skin_temperature_K > 150.0 &&
                                                         coupling->ground_skin_temperature_K < 400.0)) return -1;
    if (!(coupling->surface_resistance_s_per_m >= 0.0) || !isfinite(coupling->surface_resistance_s_per_m)) return -1;
    if (!(coupling->snow_cover_fraction >= 0.0 && coupling->snow_cover_fraction <= 1.0)) return -1;
    if (!(coupling->snow_depth_m >= 0.0) || !isfinite(coupling->snow_depth_m)) return -1;
    if (!(coupling->snow_albedo_direct_visible >= 0.0 && coupling->snow_albedo_direct_visible < 1.0)) return -1;
    if (!(coupling->snow_albedo_direct_nir >= 0.0 && coupling->snow_albedo_direct_nir < 1.0)) return -1;
    if (!(coupling->snow_albedo_diffuse_visible >= 0.0 && coupling->snow_albedo_diffuse_visible < 1.0)) return -1;
    if (!(coupling->snow_albedo_diffuse_nir >= 0.0 && coupling->snow_albedo_diffuse_nir < 1.0)) return -1;
    if (!(coupling->canopy_wetted_fraction >= 0.0 && coupling->canopy_wetted_fraction <= 1.0)) return -1;
    if (!(coupling->canopy_liquid_water_mm >= 0.0) || !isfinite(coupling->canopy_liquid_water_mm)) return -1;
    if (!(coupling->canopy_ice_mm >= 0.0) || !isfinite(coupling->canopy_ice_mm)) return -1;
    // FLO 10/26: snow albedos may be 0 (NoahMP_snow_ablation at night, as Noah-MP ALBEDO); SMC(1) NaN = not supplied
    if (!isnan(coupling->soil_top_disc_moisture_content) &&
        !(coupling->soil_top_disc_moisture_content >= 0.0 && coupling->soil_top_disc_moisture_content <= 1.0)) return -1;
    // FLO 10/24: PAHV, any sign; NWM Noah-MP clamps it to +/-20 W m-2, the CIESM revised physics does not
    if (!isfinite(coupling->precipitation_advected_heat_to_canopy_W_per_m2) ||
        fabs(coupling->precipitation_advected_heat_to_canopy_W_per_m2) > 1.0e4) return -1;
    // FLO 10/24: start-of-step canopy state, NaN = not known (cold start)
    if (!isnan(coupling->vegetation_temperature_start_of_step_K) &&
        !(coupling->vegetation_temperature_start_of_step_K > 150.0 &&
          coupling->vegetation_temperature_start_of_step_K < 400.0)) return -1;
    if (!isnan(coupling->canopy_air_temperature_start_of_step_K) &&
        !(coupling->canopy_air_temperature_start_of_step_K > 150.0 &&
          coupling->canopy_air_temperature_start_of_step_K < 400.0)) return -1;
    if (!isnan(coupling->canopy_air_vapor_pressure_start_of_step_Pa) &&
        !(coupling->canopy_air_vapor_pressure_start_of_step_Pa >= 0.0 &&
          coupling->canopy_air_vapor_pressure_start_of_step_Pa < 2.0e4)) return -1;
    return 0;
}

// One time input: whole numbers in range (time_second may be fractional). The day of month is
// checked against the month and year at update, because the driver may set it first.
// Returns 0 if valid, -1 if invalid, 1 if the name is not a time input.
int bb6_validate_time_input(const char *time_variable_name, double value)
{
    if (time_variable_name == NULL) return -1;
    if (strcmp(time_variable_name, "time_second") == 0)
        return (value >= 0.0 && value < 60.0) ? 0 : -1;
    if (strcmp(time_variable_name, "time_year") != 0 && strcmp(time_variable_name, "time_month") != 0 &&
        strcmp(time_variable_name, "time_day_of_month") != 0 && strcmp(time_variable_name, "time_hour") != 0 &&
        strcmp(time_variable_name, "time_minute") != 0) return 1;
    if (!isfinite(value) || value != floor(value)) return -1;      // whole numbers only
    if (strcmp(time_variable_name, "time_year") == 0) return (value >= 1.0 && value <= 9999.0) ? 0 : -1;
    if (strcmp(time_variable_name, "time_month") == 0) return (value >= 1.0 && value <= 12.0) ? 0 : -1;
    if (strcmp(time_variable_name, "time_day_of_month") == 0) return (value >= 1.0 && value <= 31.0) ? 0 : -1;
    if (strcmp(time_variable_name, "time_hour") == 0) return (value >= 0.0 && value <= 23.0) ? 0 : -1;
    return (value >= 0.0 && value <= 59.0) ? 0 : -1;                  // time_minute
}

// ---------------------------------------------------------------------------------------------
// BMI Initialize: parse the configuration once.
// ---------------------------------------------------------------------------------------------
int bb6_initialize_from_configuration(pet_model *pet, const char *configuration_file_path)
{
    ConfigurationTable table;
    BallBerryMethod6Parameters *parameters;
    BallBerryCouplingInputs *coupling;
    BallBerryCouplingOutputs *coupling_outputs;
    char month_key[32];
    int month_index, status = 0;

    if (pet == NULL || configuration_file_path == NULL) return -1;
    memset(&table, 0, sizeof(table));
    if (read_configuration_table(configuration_file_path, &table) != 0) {
        free_configuration_table(&table);
        return -1;
    }
    parameters = &pet->bb_parameters;
    memset(parameters, 0, sizeof(*parameters));

    // check on an option the C method 6 does not implement
    if ((int)configuration_real(&table, "radiation_transfer_option") != 3) {
        fprintf(stderr, "ERROR: only radiation_transfer_option = 3 (NWM) is ported to C\n");
        free_configuration_table(&table);
        return -1;
    }

    // ---- parameters and constants that are not BMI variables
    parameters->use_monthly_lai_sai_table = (int)configuration_real(&table, "use_monthly_lai_sai_table");
    parameters->latitude_degrees = configuration_real(&table, "latitude_degrees");
    for (month_index = 0; month_index < 12; ++month_index) {
        snprintf(month_key, sizeof(month_key), "lai_month_%02d", month_index + 1);
        parameters->leaf_area_index_by_month[month_index] = configuration_real(&table, month_key);
        snprintf(month_key, sizeof(month_key), "sai_month_%02d", month_index + 1);
        parameters->stem_area_index_by_month[month_index] = configuration_real(&table, month_key);
    }
    parameters->canopy_top_height_m = configuration_real(&table, "canopy_top_height_m");
    parameters->canopy_bottom_height_m = configuration_real(&table, "canopy_bottom_height_m");
    parameters->momentum_roughness_length_m = configuration_real(&table, "momentum_roughness_length_m");
    parameters->forcing_height_above_canopy_top_m = configuration_real(&table, "forcing_height_above_canopy_top_m");
    parameters->characteristic_leaf_dimension_m = configuration_real(&table, "characteristic_leaf_dimension_m");
    parameters->leaf_orientation_index = configuration_real(&table, "leaf_orientation_index");
    parameters->cosine_solar_zenith_angle = configuration_real(&table, "cosine_solar_zenith_angle");
    parameters->direct_shortwave_fraction = configuration_real(&table, "direct_shortwave_fraction");
    parameters->visible_shortwave_fraction = configuration_real(&table, "visible_shortwave_fraction");
    parameters->soil_albedo_direct[0] = configuration_real(&table, "soil_albedo_direct_visible");
    parameters->soil_albedo_direct[1] = configuration_real(&table, "soil_albedo_direct_nir");
    parameters->soil_albedo_diffuse[0] = configuration_real(&table, "soil_albedo_diffuse_visible");
    parameters->soil_albedo_diffuse[1] = configuration_real(&table, "soil_albedo_diffuse_nir");
    parameters->leaf_reflectance[0] = configuration_real(&table, "leaf_reflectance_visible");
    parameters->leaf_reflectance[1] = configuration_real(&table, "leaf_reflectance_nir");
    parameters->stem_reflectance[0] = configuration_real(&table, "stem_reflectance_visible");
    parameters->stem_reflectance[1] = configuration_real(&table, "stem_reflectance_nir");
    parameters->leaf_transmittance[0] = configuration_real(&table, "leaf_transmittance_visible");
    parameters->leaf_transmittance[1] = configuration_real(&table, "leaf_transmittance_nir");
    parameters->stem_transmittance[0] = configuration_real(&table, "stem_transmittance_visible");
    parameters->stem_transmittance[1] = configuration_real(&table, "stem_transmittance_nir");
    parameters->soil_emissivity = configuration_real(&table, "soil_emissivity");
    parameters->foliage_nitrogen_percent = configuration_real(&table, "foliage_nitrogen_percent");
    parameters->minimum_temperature_for_photosynthesis_k = configuration_real(&table, "minimum_temperature_for_photosynthesis_k");
    parameters->co2_mole_fraction = configuration_real(&table, "co2_mole_fraction");
    parameters->o2_mole_fraction = configuration_real(&table, "o2_mole_fraction");
    parameters->ground_temperature_offset_from_air_k = configuration_real(&table, "ground_temperature_offset_from_air_k");
    parameters->initialize_canopy_state_from_forcing = (int)configuration_real(&table, "initialize_canopy_state_from_forcing");
    parameters->initial_vegetation_temperature_k = configuration_real(&table, "initial_vegetation_temperature_k");
    parameters->initial_canopy_air_temperature_k = configuration_real(&table, "initial_canopy_air_temperature_k");
    parameters->initial_canopy_air_vapor_pressure_pa = configuration_real(&table, "initial_canopy_air_vapor_pressure_pa");
    // FLO 10/24: canopy energy physics option and carried canopy state
    parameters->canopy_energy_physics_option = (int)configuration_real(&table, "canopy_energy_physics_option");
    parameters->carry_canopy_state_between_steps = (int)configuration_real(&table, "carry_canopy_state_between_steps");
    parameters->canopy_heat_capacity_per_unit_area_index_j_per_m2_k =
        configuration_real(&table, "canopy_heat_capacity_per_unit_area_index_J_per_m2_per_K");
    parameters->canopy_phase_change_window_k = configuration_real(&table, "canopy_phase_change_window_k");
    parameters->frozen_transpiration_ramp_k = configuration_real(&table, "frozen_transpiration_ramp_k");
    // FLO 10/26: options for the decoupled Noah-MP stack (optional keys, defaults = before 10/26)
    parameters->solar_geometry_option = (int)configuration_real_or_default(&table, "solar_geometry_option", 0.0);
    parameters->solar_geometry_at_end_of_step =
        (int)configuration_real_or_default(&table, "solar_geometry_at_end_of_step", 1.0);
    // FLO 10/27: the date of the monthly LAI/SAI table; 0 keeps the behavior before 10/27
    parameters->lai_sai_at_end_of_step = (int)configuration_real_or_default(&table, "lai_sai_at_end_of_step", 0.0);
    parameters->longitude_degrees = configuration_real_or_default(&table, "longitude_degrees", 0.0);
    parameters->terrain_slope_degrees = configuration_real_or_default(&table, "terrain_slope_degrees", 0.0);
    parameters->terrain_azimuth_degrees = configuration_real_or_default(&table, "terrain_azimuth_degrees", 0.0);
    parameters->dynamic_vegetation_option = (int)configuration_real_or_default(&table, "dynamic_vegetation_option", 0.0);
    parameters->surface_resistance_snow_blend_option =
        (int)configuration_real_or_default(&table, "surface_resistance_snow_blend_option", 0.0);
    parameters->surface_resistance_snow_s_per_m =
        configuration_real_or_default(&table, "surface_resistance_snow_s_per_m", 50.0);
    parameters->soil_albedo_option = (int)configuration_real_or_default(&table, "soil_albedo_option", 0.0);
    parameters->soil_albedo_saturated[0] = configuration_real_or_default(&table, "soil_albedo_saturated_visible", 0.09);
    parameters->soil_albedo_saturated[1] = configuration_real_or_default(&table, "soil_albedo_saturated_nir", 0.18);
    parameters->soil_albedo_dry[0] = configuration_real_or_default(&table, "soil_albedo_dry_visible", 0.18);
    parameters->soil_albedo_dry[1] = configuration_real_or_default(&table, "soil_albedo_dry_nir", 0.36);

    // ---- STOMATA physiology; MP and VCMX25 are also BMI variables (initial values here)
    pet->bb_stomata.bp = configuration_real(&table, "minimum_leaf_conductance_umol_per_m2_per_s");
    pet->bb_stomata.folnmx = configuration_real(&table, "foliage_nitrogen_at_saturation_percent");
    pet->bb_stomata.qe25 = configuration_real(&table, "quantum_efficiency_25c");
    pet->bb_stomata.kc25 = configuration_real(&table, "co2_michaelis_menten_constant_25c_pa");
    pet->bb_stomata.akc = configuration_real(&table, "q10_for_kc25");
    pet->bb_stomata.ko25 = configuration_real(&table, "o2_michaelis_menten_constant_25c_pa");
    pet->bb_stomata.ako = configuration_real(&table, "q10_for_ko25");
    pet->bb_stomata.vcmx25 = configuration_real(&table, "maximum_carboxylation_rate_25c_umol_per_m2_per_s");
    pet->bb_stomata.avcmx = configuration_real(&table, "q10_for_vcmx25");
    pet->bb_stomata.c3psn = configuration_real(&table, "c3_photosynthesis_flag");
    pet->bb_stomata.mp = configuration_real(&table, "conductance_photosynthesis_slope");

    // ---- initial values of the BMI vegetation variables LAI, SAI, SHDFAC, CWP
    pet->bb_canopy.lai = configuration_real(&table, "leaf_area_index");
    pet->bb_canopy.sai = configuration_real(&table, "stem_area_index");
    pet->bb_canopy.shdfac = configuration_real(&table, "green_vegetation_fraction");
    pet->bb_cwp = configuration_real(&table, "canopy_wind_parameter");

    // ---- initial values of the BMI coupling inputs (same key names as the BMI variables)
    coupling = &pet->bb_coupling_inputs;
    coupling->root_stress_factor_0_to_1 = configuration_real(&table, "root_stress_factor_0_to_1");
    coupling->ground_skin_temperature_K = NAN;   // until set through BMI: air temperature + offset
    coupling->surface_resistance_s_per_m = configuration_real(&table, "surface_resistance_s_per_m");
    coupling->snow_cover_fraction = configuration_real(&table, "snow_cover_fraction");
    coupling->snow_depth_m = configuration_real(&table, "snow_depth_m");
    coupling->snow_albedo_direct_visible = configuration_real(&table, "snow_albedo_direct_visible");
    coupling->snow_albedo_direct_nir = configuration_real(&table, "snow_albedo_direct_nir");
    coupling->snow_albedo_diffuse_visible = configuration_real(&table, "snow_albedo_diffuse_visible");
    coupling->snow_albedo_diffuse_nir = configuration_real(&table, "snow_albedo_diffuse_nir");
    coupling->canopy_wetted_fraction = configuration_real(&table, "canopy_wetted_fraction");
    coupling->canopy_liquid_water_mm = configuration_real(&table, "canopy_liquid_water_mm");
    coupling->canopy_ice_mm = configuration_real(&table, "canopy_ice_mm");
    coupling->precipitation_advected_heat_to_canopy_W_per_m2 =
        configuration_real(&table, "precipitation_advected_heat_to_canopy_W_per_m2");   // FLO 10/24
    coupling->vegetation_temperature_start_of_step_K = NAN;      // FLO 10/24: not known until the first update
    coupling->canopy_air_temperature_start_of_step_K = NAN;
    coupling->canopy_air_vapor_pressure_start_of_step_Pa = NAN;
    coupling->soil_top_disc_moisture_content =                    // FLO 10/26: optional, see soil_albedo_option
        configuration_real_or_default(&table, "soil_top_disc_moisture_content", NAN);

    // ---- time inputs are set by the driver before each update; undefined until then
    pet->bb_time_inputs.time_year = NAN;
    pet->bb_time_inputs.time_month = NAN;
    pet->bb_time_inputs.time_day_of_month = NAN;
    pet->bb_time_inputs.time_hour = NAN;
    pet->bb_time_inputs.time_minute = NAN;
    pet->bb_time_inputs.time_second = NAN;

    // ---- outputs are undefined until the first update
    coupling_outputs = &pet->bb_coupling_outputs;
    coupling_outputs->canopy_transpiration_flux_m_per_s = NAN;
    coupling_outputs->canopy_evaporation_heat_flux_W_per_m2 = NAN;
    coupling_outputs->ground_shortwave_absorbed_below_canopy_W_per_m2 = NAN;
    coupling_outputs->ground_longwave_absorbed_below_canopy_W_per_m2 = NAN;
    coupling_outputs->canopy_air_temperature_K = NAN;
    coupling_outputs->canopy_air_vapor_pressure_Pa = NAN;
    coupling_outputs->below_canopy_aerodynamic_resistance_s_per_m = NAN;
    coupling_outputs->vegetation_temperature_K = NAN;
    coupling_outputs->canopy_latent_heat_J_per_kg = NAN;   // FLO 10/24
    coupling_outputs->canopy_freezing_rate_mm_per_s = NAN;
    coupling_outputs->canopy_fusion_heat_flux_W_per_m2 = NAN;
    coupling_outputs->canopy_heat_storage_flux_W_per_m2 = NAN;
    coupling_outputs->surface_resistance_with_snow_s_per_m = NAN;   // FLO 10/26
    coupling_outputs->cosine_solar_zenith_angle = NAN;
    coupling_outputs->vegetated_fraction = NAN;

    if (table.lookup_failed) status = -1;
    free_configuration_table(&table);           // the parsed text is no longer needed
    if (status != 0) return -1;
    // FLO 10/24: the revised canopy energy physics is built on the nested (fully converged) solver;
    // the exact solver stays a line-by-line replica of NWM VEGE_FLUX
    if (parameters->canopy_energy_physics_option != 1 && parameters->canopy_energy_physics_option != 2) {
        fprintf(stderr, "ERROR: canopy_energy_physics_option must be 1 (Noah-MP) or 2 (revised)\n");
        return -1;
    }
    if (parameters->carry_canopy_state_between_steps != 0 && parameters->carry_canopy_state_between_steps != 1) {
        fprintf(stderr, "ERROR: carry_canopy_state_between_steps must be 0 or 1\n");
        return -1;
    }
    if (parameters->canopy_energy_physics_option == 2) {
        if (pet->bb_noahmp_exact_iteration == 1) {
            fprintf(stderr, "ERROR: canopy_energy_physics_option = 2 needs bb_noahmp_exact_iteration = 0 "
                            "(the exact solver is the NWM replica)\n");
            return -1;
        }
        if (parameters->carry_canopy_state_between_steps != 1) {
            fprintf(stderr, "ERROR: canopy_energy_physics_option = 2 needs carry_canopy_state_between_steps = 1 "
                            "(the heat storage term needs the TV of the previous step)\n");
            return -1;
        }
        if (!(parameters->canopy_heat_capacity_per_unit_area_index_j_per_m2_k >= 0.0) ||
            !(parameters->canopy_phase_change_window_k > 0.0) || !(parameters->frozen_transpiration_ramp_k > 0.0)) {
            fprintf(stderr, "ERROR: canopy_energy_physics_option = 2 needs canopy heat capacity >= 0, "
                            "canopy_phase_change_window_k > 0 and frozen_transpiration_ramp_k > 0\n");
            return -1;
        }
    }
    if ((parameters->solar_geometry_option != 0 && parameters->solar_geometry_option != 1) ||
        (parameters->solar_geometry_at_end_of_step != 0 && parameters->solar_geometry_at_end_of_step != 1) ||
        (parameters->lai_sai_at_end_of_step != 0 && parameters->lai_sai_at_end_of_step != 1) ||
        (parameters->initialize_canopy_state_from_forcing < 0 || parameters->initialize_canopy_state_from_forcing > 2) ||
        (parameters->dynamic_vegetation_option != 0 && parameters->dynamic_vegetation_option != 1) ||
        (parameters->surface_resistance_snow_blend_option != 0 && parameters->surface_resistance_snow_blend_option != 1) ||
        (parameters->soil_albedo_option != 0 && parameters->soil_albedo_option != 1) ||
        !(parameters->surface_resistance_snow_s_per_m > 0.0) || !(fabs(parameters->longitude_degrees) <= 360.0) ||
        !(parameters->terrain_slope_degrees >= 0.0 && parameters->terrain_slope_degrees <= 90.0)) {
        fprintf(stderr, "ERROR: solar_geometry_option, solar_geometry_at_end_of_step, lai_sai_at_end_of_step, "
                        "dynamic_vegetation_option, "
                        "surface_resistance_snow_blend_option and soil_albedo_option must be 0 or 1; "
                        "initialize_canopy_state_from_forcing 0, 1 or 2; surface_resistance_snow_s_per_m > 0; "
                        "longitude and terrain slope in range\n");
        return -1;
    }
    if (ball_berry_state_validate(&pet->bb_canopy) != 0) {
        fprintf(stderr, "ERROR: invalid leaf_area_index, stem_area_index or green_vegetation_fraction\n");
        return -1;
    }
    if (bb6_validate_coupling_inputs(coupling) != 0) {
        fprintf(stderr, "ERROR: invalid initial value of a coupling input in the configuration file\n");
        return -1;
    }
    return 0;
}

// ---------------------------------------------------------------------------------------------
// Monthly LAI/SAI interpolation, PHENOLOGY DVEG = 1, 3, 4 branch, from the BMI time inputs.
// JULIAN = day of year - 1 + fraction of day (UTC), T = 12*DAY/YEARLEN.
// ---------------------------------------------------------------------------------------------
static int interpolate_monthly_lai_sai(const BallBerryMethod6Parameters *parameters,
                                       const BallBerryTimeInputs *time_inputs,
                                       double *leaf_area_index, double *stem_area_index)
{
    int days_in_month[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
    int year, month, day_of_month, month_index, day_of_year, year_length_days;
    int month_index_1, month_index_2;
    double julian_day_since_jan1, day_value, month_position, weight_month_1, weight_month_2;

    if (isnan(time_inputs->time_year) || isnan(time_inputs->time_month) || isnan(time_inputs->time_day_of_month) ||
        isnan(time_inputs->time_hour) || isnan(time_inputs->time_minute) || isnan(time_inputs->time_second)) {
        fprintf(stderr, "ERROR: PET method 6 with use_monthly_lai_sai_table = 1 needs the BMI time inputs "
                        "time_year, time_month, time_day_of_month, time_hour, time_minute, time_second\n");
        return -1;
    }
    year = (int)time_inputs->time_year;
    month = (int)time_inputs->time_month;
    day_of_month = (int)time_inputs->time_day_of_month;
    if ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) days_in_month[1] = 29;
    if (day_of_month > days_in_month[month - 1]) {
        fprintf(stderr, "ERROR: PET method 6 time input: day %d does not exist in month %d of %d\n",
                day_of_month, month, year);
        return -1;
    }
    year_length_days = 0;
    for (month_index = 0; month_index < 12; ++month_index) year_length_days += days_in_month[month_index];
    day_of_year = day_of_month;
    for (month_index = 0; month_index < month - 1; ++month_index) day_of_year += days_in_month[month_index];
    julian_day_since_jan1 = (double)(day_of_year - 1) + time_inputs->time_hour / 24.0 +
                            time_inputs->time_minute / 1440.0 + time_inputs->time_second / 86400.0;

    day_value = julian_day_since_jan1;
    if (parameters->latitude_degrees < 0.0)          // southern hemisphere: shifted by half a year
        day_value = fmod(julian_day_since_jan1 + 0.5 * year_length_days, (double)year_length_days);
    month_position = 12.0 * day_value / (double)year_length_days;     // T
    month_index_1 = (int)(month_position + 0.5);                       // IT1, Fortran INTEGER truncation
    month_index_2 = month_index_1 + 1;                                 // IT2
    weight_month_1 = (month_index_1 + 0.5) - month_position;           // WT1
    weight_month_2 = 1.0 - weight_month_1;                             // WT2
    if (month_index_1 < 1) month_index_1 = 12;
    if (month_index_2 > 12) month_index_2 = 1;
    *leaf_area_index = weight_month_1 * parameters->leaf_area_index_by_month[month_index_1 - 1] +
                       weight_month_2 * parameters->leaf_area_index_by_month[month_index_2 - 1];
    *stem_area_index = weight_month_1 * parameters->stem_area_index_by_month[month_index_1 - 1] +
                       weight_month_2 * parameters->stem_area_index_by_month[month_index_2 - 1];
    return 0;
}

// ---------------------------------------------------------------------------------------------
// BMI update: build the canopy solver input for this step. Returns 0 on success, -1 on error.
// ---------------------------------------------------------------------------------------------
// FLO 10/28: the preparation of one step: the checks, the solar geometry, the leaf and stem area index
// (monthly table or BMI values), and the NaN conventions of the BMI inputs resolved into flags, in a
// BallBerryStepInputs for the physics (ball_berry_canopy_physics.c, differentiated by Tapenade). The checks
// and messages are those of bb6_build_canopy_input() before 10/28, in the same order. Returns 0 or -1.
int bb6_prepare_step_inputs(pet_model *pet, BallBerryStepInputs *step)
{
    const BallBerryMethod6Parameters *parameters;
    const BallBerryCouplingInputs *coupling;
    double COSZ, COSZ_HORIZ;
    double leaf_area_index, stem_area_index;

    if (pet == NULL || step == NULL) return -1;
    pet->bb_stability_input_ready = 0;
    parameters = &pet->bb_parameters;
    coupling = &pet->bb_coupling_inputs;
    if (bb6_validate_coupling_inputs(coupling) != 0) {
        fprintf(stderr, "ERROR: PET method 6 coupling input out of range\n");
        return -1;
    }

    // ---- BMI forcing (set by the driver)
    if (!(pet->aorc.air_temperature_2m_K > 150.0 && pet->aorc.air_temperature_2m_K < 400.0) ||
        !(pet->aorc.surface_pressure_Pa > 0.0) ||
        !(pet->aorc.specific_humidity_2m_kg_per_kg >= 0.0 && pet->aorc.specific_humidity_2m_kg_per_kg < 1.0) ||
        !(pet->aorc.incoming_shortwave_W_per_m2 >= 0.0) || !(pet->aorc.incoming_longwave_W_per_m2 >= 0.0) ||
        !isfinite(pet->aorc.u_wind_speed_10m_m_per_s) || !isfinite(pet->aorc.v_wind_speed_10m_m_per_s)) {
        fprintf(stderr, "ERROR: PET method 6 forcing missing or out of range (set it through BMI before update)\n");
        return -1;
    }

    // ---- FLO 10/26: GROUNDALB from the top-disc soil moisture needs it supplied
    if (parameters->soil_albedo_option == 1 && isnan(coupling->soil_top_disc_moisture_content)) {
        fprintf(stderr, "ERROR: PET method 6 with soil_albedo_option = 1 needs soil_top_disc_moisture_content\n");
        return -1;
    }

    // ---- FLO 10/26: solar geometry. Option 1: COSZ and COSZ_HORIZ from the time inputs (the physics
    //      corrects the shortwave to the sloped terrain); option 0: the configured COSZ.
    COSZ_HORIZ = 0.0;
    if (parameters->solar_geometry_option == 1) {
        if (bb6_solar_geometry(pet, &COSZ, &COSZ_HORIZ) != 0) return -1;
    } else {
        COSZ = parameters->cosine_solar_zenith_angle;
    }
    pet->bb_coupling_outputs.cosine_solar_zenith_angle = COSZ;

    // ---- leaf and stem area index: monthly table (written back to BMI) or current BMI values
    if (parameters->use_monthly_lai_sai_table == 1) {
        // FLO 10/27: dated at the BMI time inputs, or at start + DT (lai_sai_at_end_of_step = 1, NOM)
        BallBerryTimeInputs lai_sai_date = pet->bb_time_inputs;
        if (parameters->lai_sai_at_end_of_step == 1 && !isnan(lai_sai_date.time_year) &&
            !isnan(lai_sai_date.time_month) && !isnan(lai_sai_date.time_day_of_month) &&
            !isnan(lai_sai_date.time_hour) && !isnan(lai_sai_date.time_minute)) {
            int year, month, day, hour, minute;
            noahmp_advance_datetime((int)lai_sai_date.time_year, (int)lai_sai_date.time_month,
                                    (int)lai_sai_date.time_day_of_month, (int)lai_sai_date.time_hour,
                                    (int)lai_sai_date.time_minute, (long)(pet->bmi.time_step_size_s / 60.0),
                                    &year, &month, &day, &hour, &minute);
            lai_sai_date.time_year = (double)year;
            lai_sai_date.time_month = (double)month;
            lai_sai_date.time_day_of_month = (double)day;
            lai_sai_date.time_hour = (double)hour;
            lai_sai_date.time_minute = (double)minute;
        }
        if (interpolate_monthly_lai_sai(parameters, &lai_sai_date, &leaf_area_index, &stem_area_index) != 0)
            return -1;
        pet->bb_canopy.lai = leaf_area_index;
        pet->bb_canopy.sai = stem_area_index;
    } else {
        leaf_area_index = pet->bb_canopy.lai;
        stem_area_index = pet->bb_canopy.sai;
    }
    if (ball_berry_state_validate(&pet->bb_canopy) != 0) return -1;

    // ---- the step inputs
    step->air_temperature_K = pet->aorc.air_temperature_2m_K;                    // SFCTMP
    step->surface_pressure_Pa = pet->aorc.surface_pressure_Pa;                    // SFCPRS
    step->specific_humidity_kg_per_kg = pet->aorc.specific_humidity_2m_kg_per_kg;   // QAIR
    step->shortwave_down_W_per_m2 = pet->aorc.incoming_shortwave_W_per_m2;        // SOLDN
    step->longwave_down_W_per_m2 = pet->aorc.incoming_longwave_W_per_m2;          // LWDN
    step->wind_x_m_per_s = pet->aorc.u_wind_speed_10m_m_per_s;
    step->wind_y_m_per_s = pet->aorc.v_wind_speed_10m_m_per_s;
    step->cosine_solar_zenith_angle = COSZ;
    step->cosine_solar_zenith_angle_horizontal = COSZ_HORIZ;
    step->leaf_area_index = leaf_area_index;
    step->stem_area_index = stem_area_index;
    step->green_vegetation_fraction = pet->bb_canopy.shdfac;                      // FVEG = SHDFAC
    step->canopy_wind_parameter = pet->bb_cwp;                                    // CWP
    step->stomata = pet->bb_stomata;                                              // includes MP, VCMX25
    step->coupling = *coupling;
    step->ground_temperature_known = isnan(coupling->ground_skin_temperature_K) ? 0 : 1;
    step->vegetation_temperature_known = isnan(coupling->vegetation_temperature_start_of_step_K) ? 0 : 1;
    step->canopy_air_temperature_known = isnan(coupling->canopy_air_temperature_start_of_step_K) ? 0 : 1;
    step->canopy_air_vapor_pressure_known = isnan(coupling->canopy_air_vapor_pressure_start_of_step_Pa) ? 0 : 1;
    step->time_step_s = pet->bmi.time_step_size_s;                               // DT
    return 0;
}

// FLO 10/28: the canopy input results kept in the model: the coupling outputs the canopy input sets and
// the no-canopy flag and values (as bb6_build_canopy_input() set them before 10/28).
void bb6_store_canopy_input_results(pet_model *pet, const BallBerryStepResults *results)
{
    pet->bb_coupling_outputs.cosine_solar_zenith_angle = results->coupling_outputs.cosine_solar_zenith_angle;
    pet->bb_coupling_outputs.vegetated_fraction = results->coupling_outputs.vegetated_fraction;
    pet->bb_coupling_outputs.surface_resistance_with_snow_s_per_m =
        results->coupling_outputs.surface_resistance_with_snow_s_per_m;
    pet->bb_no_canopy_this_step = results->no_canopy_this_step;
    if (results->no_canopy_this_step == 1) {
        pet->bb_no_canopy_ground_shortwave_W_per_m2 = results->no_canopy_ground_shortwave_W_per_m2;
        pet->bb_no_canopy_ground_emissivity = results->no_canopy_ground_emissivity;
        pet->bb_no_canopy_air_temperature_K = results->no_canopy_air_temperature_K;
        pet->bb_no_canopy_air_vapor_pressure_Pa = results->no_canopy_air_vapor_pressure_Pa;
        pet->bb_no_canopy_vegetation_temperature_K = results->no_canopy_vegetation_temperature_K;
    }
}

// FLO 10/28: the controls of the nested solver in the canopy input (ignored by the exact solver);
// defaults if not configured.
void bb6_set_nested_solver_controls(pet_model *pet)
{
    const BallBerrySolverControls *solver = &pet->bb_solver;
    NoahmpStabilityInput *input = &pet->bb_stability_input;
    input->canopy.maximum_iterations = (solver->maximum_canopy_iterations >= 1.0) ? (int)solver->maximum_canopy_iterations : 60;
    input->canopy.temperature_tolerance_k = (solver->canopy_temperature_tolerance_k > 0.0) ? solver->canopy_temperature_tolerance_k : 1.0e-5;
    input->canopy.maximum_temperature_step_k = (solver->maximum_temperature_step_k > 0.0) ? solver->maximum_temperature_step_k : 100.0;
    input->maximum_stability_iterations = (solver->maximum_stability_iterations >= 2.0) ? (int)solver->maximum_stability_iterations : 40;
    input->resistance_relative_tolerance = (solver->stability_relative_tolerance > 0.0) ? solver->stability_relative_tolerance : 1.0e-5;
}

// The canopy input of this step (pet->bb_stability_input), for the nested solver: the preparation and
// the canopy input physics (FLO 10/28; the same result as before). Returns 0 or -1.
int bb6_build_canopy_input(pet_model *pet)
{
    BallBerryStepInputs step;
    BallBerryStepResults results;
    if (bb6_prepare_step_inputs(pet, &step) != 0) return -1;
    bb6_canopy_input_physics(&pet->bb_parameters, &step, &pet->bb_stability_input, &results);
    bb6_store_canopy_input_results(pet, &results);
    if (results.no_canopy_this_step == 1) return 0;
    bb6_set_nested_solver_controls(pet);
    pet->bb_stability_input_ready = 1;
    return 0;
}
