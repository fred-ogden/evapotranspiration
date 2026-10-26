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
//   bb6_build_canopy_input             called by run_pet at every BMI update
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

// ---- Noah-MP constants (module_sf_noahmplsm.F, MPTABLE.TBL)
static const double noahmp_specific_heat_air_j_per_kg_k = 1004.64;       // CPAIR
static const double noahmp_gas_constant_dry_air_j_per_kg_k = 287.04;     // RAIR
static const double noahmp_latent_heat_vaporization_j_per_kg = 2.5104e6; // HVAP
static const double noahmp_latent_heat_sublimation_j_per_kg = 2.8440e6;  // HSUB
static const double noahmp_stefan_boltzmann_w_per_m2_k4 = 5.67e-8;       // SB
static const double noahmp_freezing_point_k = 273.16;                    // TFRZ
static const double noahmp_minimum_positive_value = 1.0e-6;              // MPE
static const double noahmp_bare_soil_roughness_m = 0.002;                // Z0 in ENERGY
static const double noahmp_snow_surface_roughness_m = 0.002;             // Z0SNO
static const double noahmp_minimum_wind_speed_m_per_s = 1.0;             // UR >= 1 in ENERGY
static const double noahmp_maximum_effective_area_index = 6.0;           // VAIE = MIN(6, VAI)

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
    if (!(coupling->snow_albedo_direct_visible > 0.0 && coupling->snow_albedo_direct_visible < 1.0)) return -1;
    if (!(coupling->snow_albedo_direct_nir > 0.0 && coupling->snow_albedo_direct_nir < 1.0)) return -1;
    if (!(coupling->snow_albedo_diffuse_visible > 0.0 && coupling->snow_albedo_diffuse_visible < 1.0)) return -1;
    if (!(coupling->snow_albedo_diffuse_nir > 0.0 && coupling->snow_albedo_diffuse_nir < 1.0)) return -1;
    if (!(coupling->canopy_wetted_fraction >= 0.0 && coupling->canopy_wetted_fraction <= 1.0)) return -1;
    if (!(coupling->canopy_liquid_water_mm >= 0.0) || !isfinite(coupling->canopy_liquid_water_mm)) return -1;
    if (!(coupling->canopy_ice_mm >= 0.0) || !isfinite(coupling->canopy_ice_mm)) return -1;
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
int bb6_build_canopy_input(pet_model *pet)
{
    NoahmpStabilityInput *input;
    const BallBerryMethod6Parameters *parameters;
    const BallBerrySolverControls *solver;
    const BallBerryCouplingInputs *coupling;
    double air_temperature_k, surface_pressure_pa, specific_humidity_kg_per_kg, wind_speed_m_per_s;
    double downward_shortwave_w_per_m2, downward_longwave_w_per_m2;
    double air_vapor_pressure_pa, air_density_kg_per_m3;
    double leaf_area_index, stem_area_index, vegetation_area_index, green_vegetation_fraction;
    double zero_plane_displacement_m, reference_height_m;
    double initial_vegetation_temperature_k, initial_canopy_air_temperature_k;
    double initial_canopy_air_vapor_pressure_pa, latent_heat_vegetation_j_per_kg;
    // Noah-MP names used in the science code below
    double BTRAN, TG, RSURF, FSNO, SNOWH, FWET, CANLIQ, CANICE;
    double ALBSND[2], ALBSNI[2];
    double HVT, HVB, DB, FB, SNOWHC, ELAI, ESAI, ZPDG, Z0MG, EG;
    int band;

    if (pet == NULL) return -1;
    pet->bb_stability_input_ready = 0;
    parameters = &pet->bb_parameters;
    solver = &pet->bb_solver;
    coupling = &pet->bb_coupling_inputs;
    if (bb6_validate_coupling_inputs(coupling) != 0) {
        fprintf(stderr, "ERROR: PET method 6 coupling input out of range\n");
        return -1;
    }
    input = &pet->bb_stability_input;
    memset(input, 0, sizeof(*input));

    // ---- map the BMI coupling inputs (long names) to the Noah-MP names
    BTRAN  = coupling->root_stress_factor_0_to_1;   // soil moisture factor for stomatal resistance
    RSURF  = coupling->surface_resistance_s_per_m;  // ground surface resistance to evaporation (s/m)
    FSNO   = coupling->snow_cover_fraction;         // fraction of the ground covered by snow
    SNOWH  = coupling->snow_depth_m;                // snow depth (m)
    ALBSND[0] = coupling->snow_albedo_direct_visible;   // snow albedo, direct, visible
    ALBSND[1] = coupling->snow_albedo_direct_nir;       // snow albedo, direct, near infrared
    ALBSNI[0] = coupling->snow_albedo_diffuse_visible;  // snow albedo, diffuse, visible
    ALBSNI[1] = coupling->snow_albedo_diffuse_nir;      // snow albedo, diffuse, near infrared
    FWET   = coupling->canopy_wetted_fraction;      // wetted or snowed fraction of the canopy
    CANLIQ = coupling->canopy_liquid_water_mm;      // intercepted liquid water (mm)
    CANICE = coupling->canopy_ice_mm;               // intercepted ice (mm)

    // ---- BMI forcing (set by the driver), prepared as Noah-MP ATM and ENERGY do
    air_temperature_k = pet->aorc.air_temperature_2m_K;                       // SFCTMP
    surface_pressure_pa = pet->aorc.surface_pressure_Pa;                       // SFCPRS
    specific_humidity_kg_per_kg = pet->aorc.specific_humidity_2m_kg_per_kg;   // QAIR
    downward_shortwave_w_per_m2 = pet->aorc.incoming_shortwave_W_per_m2;       // SOLDN
    downward_longwave_w_per_m2 = pet->aorc.incoming_longwave_W_per_m2;         // LWDN
    if (!(air_temperature_k > 150.0 && air_temperature_k < 400.0) || !(surface_pressure_pa > 0.0) ||
        !(specific_humidity_kg_per_kg >= 0.0 && specific_humidity_kg_per_kg < 1.0) ||
        !(downward_shortwave_w_per_m2 >= 0.0) || !(downward_longwave_w_per_m2 >= 0.0) ||
        !isfinite(pet->aorc.u_wind_speed_10m_m_per_s) || !isfinite(pet->aorc.v_wind_speed_10m_m_per_s)) {
        fprintf(stderr, "ERROR: PET method 6 forcing missing or out of range (set it through BMI before update)\n");
        return -1;
    }
    wind_speed_m_per_s = fmax(hypot(pet->aorc.u_wind_speed_10m_m_per_s, pet->aorc.v_wind_speed_10m_m_per_s),
                              noahmp_minimum_wind_speed_m_per_s);              // UR
    air_vapor_pressure_pa = specific_humidity_kg_per_kg * surface_pressure_pa /
                            (0.622 + 0.378 * specific_humidity_kg_per_kg);     // EAIR
    air_density_kg_per_m3 = (surface_pressure_pa - 0.378 * air_vapor_pressure_pa) /
                            (noahmp_gas_constant_dry_air_j_per_kg_k * air_temperature_k);   // RHOAIR

    // ---- leaf and stem area index: monthly table (written back to BMI) or current BMI values
    if (parameters->use_monthly_lai_sai_table == 1) {
        if (interpolate_monthly_lai_sai(parameters, &pet->bb_time_inputs, &leaf_area_index, &stem_area_index) != 0)
            return -1;
        pet->bb_canopy.lai = leaf_area_index;
        pet->bb_canopy.sai = stem_area_index;
    } else {
        leaf_area_index = pet->bb_canopy.lai;
        stem_area_index = pet->bb_canopy.sai;
    }
    if (ball_berry_state_validate(&pet->bb_canopy) != 0) return -1;
    // PHENOLOGY small-value checks, applied to the working copies
    if (stem_area_index < 0.05) stem_area_index = 0.0;
    if (leaf_area_index < 0.05 || stem_area_index == 0.0) leaf_area_index = 0.0;

    // ---- burial of the canopy by snow, as in PHENOLOGY (ELAI, ESAI = exposed indices)
    HVT = parameters->canopy_top_height_m;
    HVB = parameters->canopy_bottom_height_m;
    DB = fmin(fmax(SNOWH - HVB, 0.0), HVT - HVB);
    FB = DB / fmax(1.0e-06, HVT - HVB);
    if (HVT > 0.0 && HVT <= 1.0) {                  // short vegetation (MB: 1.0 and 0.2 in Noah-MP)
        SNOWHC = HVT * exp(-SNOWH / 0.2);
        if (SNOWHC <= SNOWH) FB = 1.0;
        else FB = SNOWH / SNOWHC;
    }
    ELAI = leaf_area_index * (1.0 - FB);
    ESAI = stem_area_index * (1.0 - FB);
    if (ESAI < 0.05) ESAI = 0.0;
    if (ELAI < 0.05 || ESAI == 0.0) ELAI = 0.0;
    leaf_area_index = ELAI;                         // from here on, the exposed indices
    stem_area_index = ESAI;
    vegetation_area_index = leaf_area_index + stem_area_index;
    if (!(vegetation_area_index > 0.0)) {
        fprintf(stderr, "ERROR: PET method 6 needs exposed LAI + SAI > 0 (canopy buried or bare)\n");
        return -1;
    }
    green_vegetation_fraction = pet->bb_canopy.shdfac;   // FVEG

    // ---- ENERGY geometry: ZPD = 0.65*HVT, raised to SNOWH by deep snow; ZLVL = MAX(ZPD,HVT) + ZREF,
    //      raised if the snow surface (ZPDG = SNOWH) reaches it; Z0MG blends bare soil and snow.
    zero_plane_displacement_m = 0.65 * HVT;
    if (SNOWH > zero_plane_displacement_m) zero_plane_displacement_m = SNOWH;
    reference_height_m = fmax(zero_plane_displacement_m, HVT) + parameters->forcing_height_above_canopy_top_m;
    ZPDG = SNOWH;
    if (ZPDG >= reference_height_m) reference_height_m = ZPDG + parameters->forcing_height_above_canopy_top_m;
    Z0MG = noahmp_bare_soil_roughness_m * (1.0 - FSNO) + FSNO * noahmp_snow_surface_roughness_m;

    // ---- initial canopy state: cold start (forcing or configured values) ...
    if (parameters->initialize_canopy_state_from_forcing == 1) {
        initial_vegetation_temperature_k = air_temperature_k;
        initial_canopy_air_temperature_k = air_temperature_k;
        initial_canopy_air_vapor_pressure_pa = air_vapor_pressure_pa;
    } else {
        initial_vegetation_temperature_k = parameters->initial_vegetation_temperature_k;
        initial_canopy_air_temperature_k = parameters->initial_canopy_air_temperature_k;
        initial_canopy_air_vapor_pressure_pa = parameters->initial_canopy_air_vapor_pressure_pa;
    }
    // FLO 10/24: ... unless the start-of-step state is known (carried from the previous step or set)
    if (!isnan(coupling->vegetation_temperature_start_of_step_K))
        initial_vegetation_temperature_k = coupling->vegetation_temperature_start_of_step_K;
    if (!isnan(coupling->canopy_air_temperature_start_of_step_K))
        initial_canopy_air_temperature_k = coupling->canopy_air_temperature_start_of_step_K;
    if (!isnan(coupling->canopy_air_vapor_pressure_start_of_step_Pa))
        initial_canopy_air_vapor_pressure_pa = coupling->canopy_air_vapor_pressure_start_of_step_Pa;
    latent_heat_vegetation_j_per_kg = (initial_vegetation_temperature_k > noahmp_freezing_point_k) ?
        noahmp_latent_heat_vaporization_j_per_kg : noahmp_latent_heat_sublimation_j_per_kg;   // LATHEAV

    // ---- radiation (ALBEDO, TWOSTREAM, SURRAD)
    NoahmpAlbedoInput *radiation = &input->canopy.radiation_stomata.radiation;
    radiation->cosz = parameters->cosine_solar_zenith_angle;
    radiation->elai = leaf_area_index;
    radiation->esai = stem_area_index;
    radiation->fveg = green_vegetation_fraction;
    radiation->fwet = FWET;
    radiation->tv_k = initial_vegetation_temperature_k;
    radiation->xl = parameters->leaf_orientation_index;
    radiation->omegas[0] = 0.8;     // OMEGAS, MPTABLE rad_parameters (weighted by FWET below TFRZ)
    radiation->omegas[1] = 0.4;
    radiation->betads = 0.5;        // BETADS
    radiation->betais = 0.5;        // BETAIS
    for (band = 0; band < 2; ++band) {       // 0 = visible, 1 = near infrared
        radiation->rhol[band] = parameters->leaf_reflectance[band];
        radiation->rhos[band] = parameters->stem_reflectance[band];
        radiation->taul[band] = parameters->leaf_transmittance[band];
        radiation->taus[band] = parameters->stem_transmittance[band];
        // GROUNDALB: ALBGRD = ALBSOD*(1-FSNO) + ALBSND*FSNO, ALBGRI = ALBSOI*(1-FSNO) + ALBSNI*FSNO
        radiation->albgrd[band] = parameters->soil_albedo_direct[band] * (1.0 - FSNO) + ALBSND[band] * FSNO;
        radiation->albgri[band] = parameters->soil_albedo_diffuse[band] * (1.0 - FSNO) + ALBSNI[band] * FSNO;
    }
    radiation->solad[0] = downward_shortwave_w_per_m2 * parameters->direct_shortwave_fraction *
                          parameters->visible_shortwave_fraction;
    radiation->solad[1] = downward_shortwave_w_per_m2 * parameters->direct_shortwave_fraction *
                          (1.0 - parameters->visible_shortwave_fraction);
    radiation->solai[0] = downward_shortwave_w_per_m2 * (1.0 - parameters->direct_shortwave_fraction) *
                          parameters->visible_shortwave_fraction;
    radiation->solai[1] = downward_shortwave_w_per_m2 * (1.0 - parameters->direct_shortwave_fraction) *
                          (1.0 - parameters->visible_shortwave_fraction);

    // ---- physiology (STOMATA): current BMI values, including MP and VCMX25
    input->canopy.radiation_stomata.physiology = pet->bb_stomata;

    NoahmpStomataInputs *canopy_air = &input->canopy.radiation_stomata.canopy_air;
    canopy_air->mpe = noahmp_minimum_positive_value;
    canopy_air->foln = parameters->foliage_nitrogen_percent;
    canopy_air->igs = (initial_vegetation_temperature_k > parameters->minimum_temperature_for_photosynthesis_k) ? 1.0 : 0.0;
    canopy_air->ea = initial_canopy_air_vapor_pressure_pa;   // EAH seen by STOMATA on iteration 1
    canopy_air->sfctmp = air_temperature_k;
    canopy_air->sfcpres = surface_pressure_pa;
    canopy_air->co2 = parameters->co2_mole_fraction * surface_pressure_pa;   // CO2AIR
    canopy_air->o2 = parameters->o2_mole_fraction * surface_pressure_pa;     // O2AIR
    canopy_air->btran = BTRAN;                                               // root stress factor

    // ---- foliage energy balance
    NoahmpVegeFluxInputs *foliage = &input->canopy.foliage;
    foliage->fveg = green_vegetation_fraction;
    foliage->vai = vegetation_area_index;
    foliage->wet_fraction = FWET;
    foliage->rsurf_s_per_m = RSURF;
    foliage->air_temperature_k = air_temperature_k;
    // TG from soil freeze-thaw once supplied through BMI; before that, air temperature + offset
    if (isnan(coupling->ground_skin_temperature_K))
        TG = air_temperature_k + parameters->ground_temperature_offset_from_air_k;
    else
        TG = coupling->ground_skin_temperature_K;
    foliage->ground_temperature_k = TG;
    foliage->foliage_temperature_k = initial_vegetation_temperature_k;
    foliage->air_vapor_pressure_pa = air_vapor_pressure_pa;
    foliage->air_density_kg_per_m3 = air_density_kg_per_m3;
    foliage->cp_air_j_per_kg_k = noahmp_specific_heat_air_j_per_kg_k;
    foliage->psychrometric_constant_pa_per_k = noahmp_specific_heat_air_j_per_kg_k * surface_pressure_pa /
                                               (0.622 * latent_heat_vegetation_j_per_kg);   // GAMMAV
    foliage->incoming_longwave_w_per_m2 = downward_longwave_w_per_m2;
    foliage->vegetation_emissivity = 1.0 - exp(-(leaf_area_index + stem_area_index));        // EMV
    EG = parameters->soil_emissivity;                                                         // EG(1)
    foliage->ground_emissivity = EG * (1.0 - FSNO) + 1.0 * FSNO;                              // EMG
    foliage->stefan_boltzmann_w_per_m2_k4 = noahmp_stefan_boltzmann_w_per_m2_k4;
    foliage->canopy_liquid_mm = CANLIQ;
    foliage->canopy_ice_mm = CANICE;
    foliage->latent_heat_j_per_kg = latent_heat_vegetation_j_per_kg;
    foliage->timestep_s = pet->bmi.time_step_size_s;                                          // DT
    foliage->additional_heat_w_per_m2 = coupling->precipitation_advected_heat_to_canopy_W_per_m2;   // PAHV (FLO 10/24)

    // ---- FLO 10/24: revised canopy energy physics (option 2). With option 1 these stay 0 (memset above),
    //      which is exactly the Noah-MP calculation.
    if (parameters->canopy_energy_physics_option == 2) {
        // C_v = FVEG * c * (ELAI+ESAI) + heat capacity of the intercepted water (CANLIQ, CANICE per ground area)
        foliage->canopy_heat_capacity_j_per_m2_k =
            green_vegetation_fraction * parameters->canopy_heat_capacity_per_unit_area_index_j_per_m2_k *
            vegetation_area_index + (4.188e6 * CANLIQ + 2.094e6 * CANICE) / 1000.0;   // CWAT, CICE (J m-3 K-1)
        foliage->previous_foliage_temperature_k = initial_vegetation_temperature_k;
        foliage->phase_change_window_k = parameters->canopy_phase_change_window_k;
        foliage->canopy_water_total_mm = CANLIQ + CANICE;
        foliage->canopy_liquid_start_mm = CANLIQ;
        foliage->frozen_transpiration_ramp_k = parameters->frozen_transpiration_ramp_k;
        foliage->evaporation_phase_from_latent_heat = 1;
        input->canopy.hold_radiation_at_start_temperature = 1;
        input->noahmp_outer_iteration_rule = 1;     // never fails in stable, calm air (see ball_berry_helpers.h)
    }

    // ---- nested solver controls (ignored by the exact solver); defaults if not configured
    input->canopy.maximum_iterations = (solver->maximum_canopy_iterations >= 1.0) ? (int)solver->maximum_canopy_iterations : 60;
    input->canopy.temperature_tolerance_k = (solver->canopy_temperature_tolerance_k > 0.0) ? solver->canopy_temperature_tolerance_k : 1.0e-5;
    input->canopy.maximum_temperature_step_k = (solver->maximum_temperature_step_k > 0.0) ? solver->maximum_temperature_step_k : 100.0;
    input->maximum_stability_iterations = (solver->maximum_stability_iterations >= 2.0) ? (int)solver->maximum_stability_iterations : 40;
    input->resistance_relative_tolerance = (solver->stability_relative_tolerance > 0.0) ? solver->stability_relative_tolerance : 1.0e-5;

    // ---- SFCDIF1 surface layer
    NoahmpSfcdif1Inputs *surface_layer = &input->surface_layer;
    surface_layer->air_temperature_k = air_temperature_k;
    surface_layer->air_density_kg_per_m3 = air_density_kg_per_m3;
    surface_layer->specific_humidity_kg_per_kg = specific_humidity_kg_per_kg;
    surface_layer->reference_height_m = reference_height_m;                    // ZLVL
    surface_layer->displacement_height_m = zero_plane_displacement_m;          // ZPD
    surface_layer->momentum_roughness_m = parameters->momentum_roughness_length_m;   // Z0M = Z0MVT
    surface_layer->heat_roughness_m = parameters->momentum_roughness_length_m;       // Z0H = Z0M
    surface_layer->wind_speed_m_per_s = wind_speed_m_per_s;
    surface_layer->mpe = noahmp_minimum_positive_value;

    // ---- RAGRB under canopy and leaf boundary layer
    NoahmpRagrbInputs *ground_canopy = &input->ground_canopy;
    ground_canopy->reference_wind_speed_m_per_s = wind_speed_m_per_s;
    ground_canopy->reference_height_m = reference_height_m;
    ground_canopy->momentum_roughness_m = parameters->momentum_roughness_length_m;
    ground_canopy->displacement_height_m = zero_plane_displacement_m;
    ground_canopy->canopy_height_m = HVT;                                         // HCAN = HVT
    ground_canopy->vegetation_area_index = fmin(noahmp_maximum_effective_area_index, vegetation_area_index);  // VAIE
    ground_canopy->air_density_kg_per_m3 = air_density_kg_per_m3;
    ground_canopy->canopy_air_temperature_k = initial_canopy_air_temperature_k;
    ground_canopy->ground_momentum_roughness_m = Z0MG;
    ground_canopy->ground_heat_roughness_m = Z0MG;                                // Z0HG = Z0MG
    ground_canopy->canopy_heat_roughness_m = parameters->momentum_roughness_length_m;   // Z0H = Z0M
    ground_canopy->leaf_dimension_m = parameters->characteristic_leaf_dimension_m;     // DLEAF
    ground_canopy->cwp = pet->bb_cwp;                                             // BMI variable CWP
    ground_canopy->mpe = noahmp_minimum_positive_value;

    pet->bb_stability_input_ready = 1;
    return 0;
}
