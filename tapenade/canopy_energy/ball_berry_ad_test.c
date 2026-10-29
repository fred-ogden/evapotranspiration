/*
 * tapenade/canopy_energy/ball_berry_ad_test.c - tests of the Tapenade differentiation of one PET method 6
 * (Ball-Berry canopy energy balance) step with the VEGE_FLUX loop1 replica (ball_berry_head.h).
 *
 * Author: Fred L. Ogden, Ph.D., P.E.
 *         NOAA/National Weather Service
 *
 * This software was developed by a U.S. Government employee as part of
 * official duties and is not subject to U.S. copyright protection under
 * 17 U.S.C. Section 105.
 *
 * Licensed under the Apache License, Version 2.0.
 * SPDX-License-Identifier: Apache-2.0
 */

/*
 * The inputs are the ones dcfe-noahmp-stack gave Ball-Berry for one divide (build/stack_driver --core-path
 * --dump-slot-inputs, the canopy_energy slot: cat-18560 of the 01123000 basin over the whole realization
 * period; data/).  Ball-Berry is initialized through its BMI from the divide's configuration, as the stack
 * does.  Each hour the dumped inputs are set into the model as the stack's adapter does
 * (tapenade/stack/stack_slot_ball_berry.c of dcfe-noahmp-stack); the preparation (bb6_prepare_step_inputs)
 * is run on a copy of the model, the head's independents are packed from the step inputs it gives, and
 * the model is stepped with ball_berry_method6_core_step(), as the stack's adapter steps it.
 *
 * Experiments (letters, any order):
 *   p  the head against the production step at every hour, bitwise, every dependent
 *   t  the tangent against central differences at the test states, along each independent; step sizes
 *      scanned; kinks (forward and backward differences disagree at every step size in some dependent,
 *      relative to its own size) reported apart, with whether the tangent equals one side's difference
 *      there; the tangent and central differences are compared over the dependents the step resolves
 *      (each changes by more than 1e-9 of its value); directions below finite-difference resolution (no
 *      dependent resolved) tested only for agreement (the tangent must not predict a resolvable change
 *      either); directions where every step would cross a bound of the domain (zero snow or wetted
 *      fraction, BTRAN = 1, ...) not tested; the tangent code's primal checked bitwise against the head; a
 *      non-finite derivative is a failure
 *   a  the reverse Jacobian against the tangent Jacobian, and the dot-product test
 * Built with -DBB_AD_PRIMAL_ONLY (no generated code), t reports the central differences'
 * self-consistency across step sizes, and a is skipped.
 *
 * Test states: the first hour (after a 48-hour spin-up) of each regime (below), and every 5000th hour.
 *
 * Usage: ball_berry_ad_test CONFIG CANOPY_ENERGY_CSV RESULTS_DIRECTORY EXPERIMENTS
 *
 * FLO 10/28: first version (from CFE3.1's tapenade/cfe_hour/cfe_hour_ad_test.c).  ASCII only.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bmi.h"
#include "bmi_pet.h"
#include "pet.h"
#include "ball_berry_head.h"
#ifndef BB_AD_PRIMAL_ONLY
#include "ball_berry_adapters.h"
#endif

#define MAXIMUM_LINE 16384
#define MAXIMUM_COLUMNS 64
#define MAXIMUM_TEST_POINTS 64
#define SPIN_UP_HOURS 48
#define REGULAR_TEST_STRIDE 5000
#define STEP_SIZE_COUNT 5
#define KINK_THRESHOLD 1.0e-3
#define KINK_SIGNIFICANCE 1.0e-3      /* a dependent's one-sided differences count for the kink test when
                                         they exceed this fraction of the column's largest */
#define ONE_SIDED_THRESHOLD 1.0e-4    /* at a kink, the tangent should equal one side's difference */
#define RESOLUTION_THRESHOLD 1.0e-9   /* a dependent is resolved by a central difference when it changes by more
                                         than this fraction of its value */
#define TANGENT_PASS_THRESHOLD 1.0e-4
#define REVERSE_PASS_THRESHOLD 1.0e-10
#define REVERSE_ENTRIES_REPORTED 12   /* at a point that fails a, the Jacobian entries printed */
#define NO_LOWER_BOUND (-1.0e300)
#define NO_UPPER_BOUND (1.0e300)
#define LATENT_HEAT_SUBLIMATION_J_PER_KG 2.8440e6
#define NIN BB_HEAD_INDEPENDENT_COUNT
#define NOUT BB_HEAD_DEPENDENT_COUNT

static const double relative_step_size_list[STEP_SIZE_COUNT] = { 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7 };

#define COMPONENT_COUNT 1
static const char *const component_names[COMPONENT_COUNT] = { "canopy" };

static const char *const independent_names[NIN] = {
    "TV_start", "TAH_start", "EAH_start", "air_temperature", "surface_pressure", "specific_humidity",
    "shortwave_down", "longwave_down", "wind_x", "wind_y", "BTRAN", "TG", "RSURF", "FSNO", "SNOWH",
    "ALBSND_visible", "ALBSND_nir", "ALBSNI_visible", "ALBSNI_nir", "FWET", "CANLIQ", "CANICE", "PAHV", "SMC_1",
    "MP", "VCMX25", "BP", "QE25", "CWP", "SHDFAC", "Z0M", "DLEAF" };
static const double typical_scale[NIN] = {
    1.0, 1.0, 10.0, 1.0, 100.0, 1.0e-3,
    10.0, 10.0, 1.0, 1.0, 0.1, 1.0, 10.0, 0.1, 0.1,
    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 0.1,
    1.0, 10.0, 100.0, 0.01, 0.1, 0.1, 0.1, 0.01 };
/* the bounds of each independent's domain: a central difference must not cross them */
static const double lower_bound[NIN] = {
    NO_LOWER_BOUND, NO_LOWER_BOUND, 0.0, NO_LOWER_BOUND, 0.0, 0.0,
    0.0, 0.0, NO_LOWER_BOUND, NO_LOWER_BOUND, 0.0, NO_LOWER_BOUND, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NO_LOWER_BOUND, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
static const double upper_bound[NIN] = {
    NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND,
    NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND, 1.0, NO_UPPER_BOUND, NO_UPPER_BOUND, 1.0, NO_UPPER_BOUND,
    1.0, 1.0, 1.0, 1.0, 1.0, NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND, 1.0,
    NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND, NO_UPPER_BOUND, 1.0, NO_UPPER_BOUND, NO_UPPER_BOUND };
static const char *const dependent_names[NOUT] = {
    "ETRAN", "EVC", "SAG", "ground_longwave_absorbed", "TAH", "EAH", "RAHG", "TV", "LATHEAV",
    "fusion_heat", "heat_storage", "RSURF_with_snow", "FVEG", "PET" };

typedef struct {
    int component;
    char time_text[32];
    char regime[64];
    BallBerryMethod6Parameters parameters;
    BallBerryStepInputs step;
    double independent[NIN];
} TestPoint;

/* ======================================================================================
 * CSV
 * ====================================================================================== */
typedef struct {
    FILE *file;
    int column_count;
    char names[MAXIMUM_COLUMNS][96];
    char time_text[32];
    double values[MAXIMUM_COLUMNS];
} CsvReader;

static void csv_open(CsvReader *reader, const char *path)
{
    char line[MAXIMUM_LINE], *token;
    reader->file = fopen(path, "r");
    if (reader->file == NULL || fgets(line, sizeof(line), reader->file) == NULL) {
        fprintf(stderr, "ERROR: cannot read %s\n", path);
        exit(1);
    }
    reader->column_count = 0;
    for (token = strtok(line, ",\r\n"); token != NULL && reader->column_count < MAXIMUM_COLUMNS; token = strtok(NULL, ",\r\n"))
        snprintf(reader->names[reader->column_count++], sizeof(reader->names[0]), "%s", token);
}

static int csv_next(CsvReader *reader)
{
    char line[MAXIMUM_LINE], *token;
    int index = 0;
    if (fgets(line, sizeof(line), reader->file) == NULL) return 0;
    for (token = strtok(line, ",\r\n"); token != NULL && index < reader->column_count; token = strtok(NULL, ",\r\n")) {
        if (index == 0) snprintf(reader->time_text, sizeof(reader->time_text), "%s", token);
        else reader->values[index] = strtod(token, NULL);
        ++index;
    }
    return 1;
}

static double csv_value(const CsvReader *reader, const char *name)
{
    int index;
    for (index = 1; index < reader->column_count; ++index) {
        if (strcmp(reader->names[index], name) == 0) return reader->values[index];
    }
    fprintf(stderr, "ERROR: column %s not in the dump\n", name);
    exit(1);
}

/* ======================================================================================
 * the model: initialize (as the stack); the hourly inputs (as its adapter)
 * ====================================================================================== */
static Bmi *initialize_model(const char *configuration)
{
    Bmi *model = (Bmi *)calloc(1, sizeof(Bmi));
    if (model == NULL) exit(1);
    register_bmi_pet(model);
    if (model->initialize(model, configuration) != BMI_SUCCESS) {
        fprintf(stderr, "ERROR: Ball-Berry initialize failed for %s\n", configuration);
        exit(1);
    }
    return model;
}

/* the slot inputs into the model, as stack_slot_canopy_energy_run() sets them */
static void set_hour_inputs(pet_model *pet, const CsvReader *row)
{
    BallBerryCouplingInputs *coupling = &pet->bb_coupling_inputs;
    BallBerryTimeInputs *time_inputs = &pet->bb_time_inputs;
    time_inputs->time_year = csv_value(row, "time_year");
    time_inputs->time_month = csv_value(row, "time_month");
    time_inputs->time_day_of_month = csv_value(row, "time_day_of_month");
    time_inputs->time_hour = csv_value(row, "time_hour");
    time_inputs->time_minute = csv_value(row, "time_minute");
    time_inputs->time_second = csv_value(row, "time_second");
    pet->aorc.incoming_longwave_W_per_m2 = csv_value(row, "longwave_down_W_per_m2");
    pet->aorc.surface_pressure_Pa = csv_value(row, "air_pressure_Pa");
    pet->aorc.specific_humidity_2m_kg_per_kg = csv_value(row, "specific_humidity_kg_per_kg");
    pet->aorc.incoming_shortwave_W_per_m2 = csv_value(row, "shortwave_down_W_per_m2");
    pet->aorc.air_temperature_2m_K = csv_value(row, "air_temperature_K");
    pet->aorc.u_wind_speed_10m_m_per_s = csv_value(row, "wind_x_m_per_s");
    pet->aorc.v_wind_speed_10m_m_per_s = csv_value(row, "wind_y_m_per_s");
    coupling->root_stress_factor_0_to_1 = csv_value(row, "root_stress_factor_0_to_1");
    coupling->ground_skin_temperature_K = csv_value(row, "ground_skin_temperature_K");
    coupling->surface_resistance_s_per_m = csv_value(row, "surface_resistance_s_per_m");
    coupling->snow_cover_fraction = csv_value(row, "snow_cover_fraction");
    coupling->snow_depth_m = csv_value(row, "snow_depth_m");
    coupling->snow_albedo_direct_visible = csv_value(row, "snow_albedo_direct_visible");
    coupling->snow_albedo_direct_nir = csv_value(row, "snow_albedo_direct_nir");
    coupling->snow_albedo_diffuse_visible = csv_value(row, "snow_albedo_diffuse_visible");
    coupling->snow_albedo_diffuse_nir = csv_value(row, "snow_albedo_diffuse_nir");
    coupling->canopy_wetted_fraction = csv_value(row, "canopy_wetted_fraction");
    coupling->canopy_liquid_water_mm = csv_value(row, "canopy_liquid_water_mm");
    coupling->canopy_ice_mm = csv_value(row, "canopy_ice_mm");
    coupling->precipitation_advected_heat_to_canopy_W_per_m2 =
        csv_value(row, "precipitation_advected_heat_to_canopy_W_per_m2");
    coupling->soil_top_disc_moisture_content = csv_value(row, "soil_top_disc_moisture_content");
    coupling->vegetation_temperature_start_of_step_K = csv_value(row, "vegetation_temperature_start_of_step_K");
}

static void pack_independents(const BallBerryMethod6Parameters *p, const BallBerryStepInputs *s, double x[NIN])
{
    x[BB_HEAD_IN_VEGETATION_TEMPERATURE_START] = s->coupling.vegetation_temperature_start_of_step_K;
    x[BB_HEAD_IN_CANOPY_AIR_TEMPERATURE_START] = s->coupling.canopy_air_temperature_start_of_step_K;
    x[BB_HEAD_IN_CANOPY_AIR_VAPOR_PRESSURE_START] = s->coupling.canopy_air_vapor_pressure_start_of_step_Pa;
    x[BB_HEAD_IN_AIR_TEMPERATURE] = s->air_temperature_K;
    x[BB_HEAD_IN_SURFACE_PRESSURE] = s->surface_pressure_Pa;
    x[BB_HEAD_IN_SPECIFIC_HUMIDITY] = s->specific_humidity_kg_per_kg;
    x[BB_HEAD_IN_SHORTWAVE_DOWN] = s->shortwave_down_W_per_m2;
    x[BB_HEAD_IN_LONGWAVE_DOWN] = s->longwave_down_W_per_m2;
    x[BB_HEAD_IN_WIND_X] = s->wind_x_m_per_s;
    x[BB_HEAD_IN_WIND_Y] = s->wind_y_m_per_s;
    x[BB_HEAD_IN_ROOT_STRESS_FACTOR] = s->coupling.root_stress_factor_0_to_1;
    x[BB_HEAD_IN_GROUND_TEMPERATURE] = s->coupling.ground_skin_temperature_K;
    x[BB_HEAD_IN_SURFACE_RESISTANCE] = s->coupling.surface_resistance_s_per_m;
    x[BB_HEAD_IN_SNOW_COVER_FRACTION] = s->coupling.snow_cover_fraction;
    x[BB_HEAD_IN_SNOW_DEPTH] = s->coupling.snow_depth_m;
    x[BB_HEAD_IN_SNOW_ALBEDO_DIRECT_VISIBLE] = s->coupling.snow_albedo_direct_visible;
    x[BB_HEAD_IN_SNOW_ALBEDO_DIRECT_NIR] = s->coupling.snow_albedo_direct_nir;
    x[BB_HEAD_IN_SNOW_ALBEDO_DIFFUSE_VISIBLE] = s->coupling.snow_albedo_diffuse_visible;
    x[BB_HEAD_IN_SNOW_ALBEDO_DIFFUSE_NIR] = s->coupling.snow_albedo_diffuse_nir;
    x[BB_HEAD_IN_CANOPY_WETTED_FRACTION] = s->coupling.canopy_wetted_fraction;
    x[BB_HEAD_IN_CANOPY_LIQUID_WATER] = s->coupling.canopy_liquid_water_mm;
    x[BB_HEAD_IN_CANOPY_ICE] = s->coupling.canopy_ice_mm;
    x[BB_HEAD_IN_PRECIPITATION_ADVECTED_HEAT] = s->coupling.precipitation_advected_heat_to_canopy_W_per_m2;
    x[BB_HEAD_IN_SOIL_TOP_DISC_MOISTURE] = s->coupling.soil_top_disc_moisture_content;
    x[BB_HEAD_IN_STOMATAL_SLOPE] = s->stomata.mp;
    x[BB_HEAD_IN_MAXIMUM_CARBOXYLATION_RATE] = s->stomata.vcmx25;
    x[BB_HEAD_IN_MINIMUM_LEAF_CONDUCTANCE] = s->stomata.bp;
    x[BB_HEAD_IN_QUANTUM_EFFICIENCY] = s->stomata.qe25;
    x[BB_HEAD_IN_CANOPY_WIND_PARAMETER] = s->canopy_wind_parameter;
    x[BB_HEAD_IN_GREEN_VEGETATION_FRACTION] = s->green_vegetation_fraction;
    x[BB_HEAD_IN_MOMENTUM_ROUGHNESS_LENGTH] = p->momentum_roughness_length_m;
    x[BB_HEAD_IN_LEAF_DIMENSION] = p->characteristic_leaf_dimension_m;
}

/* production's outputs after the step, in the head's dependent order */
static void pack_production_dependents(const pet_model *pet, double y[NOUT])
{
    const BallBerryCouplingOutputs *o = &pet->bb_coupling_outputs;
    y[BB_HEAD_OUT_TRANSPIRATION_FLUX] = o->canopy_transpiration_flux_m_per_s;
    y[BB_HEAD_OUT_CANOPY_EVAPORATION_HEAT_FLUX] = o->canopy_evaporation_heat_flux_W_per_m2;
    y[BB_HEAD_OUT_GROUND_SHORTWAVE_ABSORBED] = o->ground_shortwave_absorbed_below_canopy_W_per_m2;
    y[BB_HEAD_OUT_GROUND_LONGWAVE_ABSORBED] = o->ground_longwave_absorbed_below_canopy_W_per_m2;
    y[BB_HEAD_OUT_CANOPY_AIR_TEMPERATURE] = o->canopy_air_temperature_K;
    y[BB_HEAD_OUT_CANOPY_AIR_VAPOR_PRESSURE] = o->canopy_air_vapor_pressure_Pa;
    y[BB_HEAD_OUT_BELOW_CANOPY_RESISTANCE] = o->below_canopy_aerodynamic_resistance_s_per_m;
    y[BB_HEAD_OUT_VEGETATION_TEMPERATURE] = o->vegetation_temperature_K;
    y[BB_HEAD_OUT_CANOPY_LATENT_HEAT] = o->canopy_latent_heat_J_per_kg;
    y[BB_HEAD_OUT_CANOPY_FUSION_HEAT_FLUX] = o->canopy_fusion_heat_flux_W_per_m2;
    y[BB_HEAD_OUT_CANOPY_HEAT_STORAGE_FLUX] = o->canopy_heat_storage_flux_W_per_m2;
    y[BB_HEAD_OUT_SURFACE_RESISTANCE_WITH_SNOW] = o->surface_resistance_with_snow_s_per_m;
    y[BB_HEAD_OUT_VEGETATED_FRACTION] = o->vegetated_fraction;
    y[BB_HEAD_OUT_PET] = pet->pet_m_per_s;
}

/* ======================================================================================
 * the head and its derivatives through one interface
 * ====================================================================================== */
static void run_head(const TestPoint *point, const double independent[], double dependent[])
{
    double independent_copy[NIN];
    int k;
    for (k = 0; k < NIN; ++k) independent_copy[k] = independent[k];
    ball_berry_head(independent_copy, dependent, &point->parameters, &point->step);
}

#ifndef BB_AD_PRIMAL_ONLY
static void run_tangent(const TestPoint *point, const double independent[], const double seed[], double dependent[],
                        double tangent[])
{
    ball_berry_tangent(independent, seed, dependent, tangent, &point->parameters, &point->step);
}

static void run_adjoint(const TestPoint *point, const double independent[], const double dependent_adjoint[],
                        double independent_adjoint[])
{
    ball_berry_adjoint(independent, dependent_adjoint, independent_adjoint, &point->parameters, &point->step);
}
#endif

#ifndef BB_AD_PRIMAL_ONLY
static double normalized_difference(const double a[], const double b[], int count)
{
    double largest_difference = 0.0, scale = 0.0;
    int i;
    for (i = 0; i < count; ++i) {
        if (fabs(a[i] - b[i]) > largest_difference) largest_difference = fabs(a[i] - b[i]);
        if (fabs(a[i]) > scale) scale = fabs(a[i]);
        if (fabs(b[i]) > scale) scale = fabs(b[i]);
    }
    if (scale == 0.0) return 0.0;
    return largest_difference / scale;
}
#endif

/* normalized_difference() over the dependents with mask[i] set (1e300, no comparison, when none is) */
static double masked_normalized_difference(const double a[], const double b[], const int mask[], int count)
{
    double largest_difference = 0.0, scale = 0.0;
    int i, compared = 0;
    for (i = 0; i < count; ++i) {
        if (!mask[i]) continue;
        ++compared;
        if (fabs(a[i] - b[i]) > largest_difference) largest_difference = fabs(a[i] - b[i]);
        if (fabs(a[i]) > scale) scale = fabs(a[i]);
        if (fabs(b[i]) > scale) scale = fabs(b[i]);
    }
    if (compared == 0) return 1.0e300;
    if (scale == 0.0) return 0.0;
    return largest_difference / scale;
}

/* the kink measure of one step size: the largest relative disagreement of the forward and backward
 * differences over the dependents, each relative to its own size (floored at KINK_SIGNIFICANCE of the
 * column's largest), so a kink in a small dependent is not hidden by a large smooth one */
static double kink_mismatch(const double forward[], const double backward[], const int mask[], int count)
{
    double column_scale = 0.0, largest = 0.0, size, mismatch;
    int i, compared = 0;
    for (i = 0; i < count; ++i) {
        if (!mask[i]) continue;
        ++compared;
        if (fabs(forward[i]) > column_scale) column_scale = fabs(forward[i]);
        if (fabs(backward[i]) > column_scale) column_scale = fabs(backward[i]);
    }
    if (compared == 0) return 1.0e300;      /* this step resolves nothing: it does not count */
    if (column_scale == 0.0) return 0.0;
    for (i = 0; i < count; ++i) {
        if (!mask[i]) continue;
        size = fmax(fmax(fabs(forward[i]), fabs(backward[i])), KINK_SIGNIFICANCE * column_scale);
        mismatch = fabs(forward[i] - backward[i]) / size;
        if (mismatch > largest) largest = mismatch;
    }
    return largest;
}

/* ---- regimes: the labels that hold this hour; a test state is the first hour of a label */
#define MAXIMUM_LABELS 16
typedef struct {
    const char *label[MAXIMUM_LABELS];
    int count;
} LabelList;

static void add_label(LabelList *labels, int condition, const char *label)
{
    if (condition && labels->count < MAXIMUM_LABELS) labels->label[labels->count++] = label;
}

/* the labels of an hour, from its step inputs, the production outputs, and the solution and trace of the
 * exact solver (the trace's Monin-Obukhov parameter shows an SFCDIF1 neutral reset: 0 after nonzero) */
static void hour_labels(const BallBerryStepInputs *step, const pet_model *pet, const NoahmpVegeFluxExactTrace *trace,
                        LabelList *labels)
{
    const BallBerryCouplingOutputs *o = &pet->bb_coupling_outputs;
    int row, neutral_reset = 0, nonzero_seen = 0;
    for (row = 0; row < trace->count; ++row) {
        if (trace->record[row].monin_obukhov_stability_parameter != 0.0) nonzero_seen = 1;
        else if (nonzero_seen && trace->record[row].iteration > 1) neutral_reset = 1;
    }
    labels->count = 0;
    add_label(labels, step->cosine_solar_zenith_angle > 0.0 && step->shortwave_down_W_per_m2 > 0.0, "day");
    add_label(labels, !(step->cosine_solar_zenith_angle > 0.0), "night");
    add_label(labels, step->coupling.root_stress_factor_0_to_1 < 1.0, "stressed (BTRAN < 1)");
    add_label(labels, step->coupling.canopy_wetted_fraction > 0.0, "wet canopy");
    add_label(labels, step->coupling.canopy_wetted_fraction > 0.0 &&
              step->coupling.canopy_liquid_water_mm + step->coupling.canopy_ice_mm < 0.05,
              "wet canopy under 0.05 mm of water");
    add_label(labels, o->canopy_latent_heat_J_per_kg == LATENT_HEAT_SUBLIMATION_J_PER_KG, "frozen canopy (HSUB)");
    add_label(labels, step->coupling.snow_cover_fraction > 0.0, "snow on the ground");
    add_label(labels, pet->bb_no_canopy_this_step == 0 && pet->bb_last_output.stability_iterations == 20,
              "loop1 at 20 iterations");
    add_label(labels, neutral_reset, "SFCDIF1 neutral reset");
    add_label(labels, step->coupling.surface_resistance_s_per_m >= 1.0e6, "dry bare soil resistance");
    add_label(labels, pet->bb_no_canopy_this_step == 1, "no exposed canopy");
}

#define REGIME_SLOTS 24
typedef struct {
    const char *found[REGIME_SLOTS];
    int found_count;
} RegimeList;

static const char *first_new_label(RegimeList *list, const LabelList *labels)
{
    int k, l, seen;
    for (l = 0; l < labels->count; ++l) {
        seen = 0;
        for (k = 0; k < list->found_count; ++k) {
            if (strcmp(list->found[k], labels->label[l]) == 0) seen = 1;
        }
        if (!seen) {
            if (list->found_count < REGIME_SLOTS) list->found[list->found_count++] = labels->label[l];
            return labels->label[l];
        }
    }
    return NULL;
}

/* ======================================================================================
 * experiment t
 * ====================================================================================== */
static void experiment_tangent(const TestPoint *points, int point_count, FILE *csv, int *failed_total)
{
    int point_index, j, s, i;
    double worst_overall = 0.0;
    int kinks_overall = 0, unresolved_overall = 0, failed_overall = 0, pairs_overall = 0, primal_mismatch_overall = 0;
    int at_bound_overall = 0, non_finite_overall = 0;
#ifndef BB_AD_PRIMAL_ONLY
    int kinks_one_sided_overall = 0;    /* kinks where the tangent equals one side's difference */
#endif
    for (point_index = 0; point_index < point_count; ++point_index) {
        const TestPoint *point = &points[point_index];
        double base[NOUT], worst = 0.0;
        int kinks = 0, unresolved_count = 0, failed = 0, primal_mismatches = 0, at_bound_count = 0;
        const char *worst_name = "-";
        run_head(point, point->independent, base);
        for (j = 0; j < NIN; ++j) {
            double central[STEP_SIZE_COUNT][NOUT], perturbed[NIN];
            int resolved[STEP_SIZE_COUNT][NOUT];    /* the dependents each step resolves */
            double plus[NOUT], minus[NOUT], forward[NOUT], backward[NOUT];
            double value = point->independent[j], step, best_error = 1.0e300, best_step = 0.0, kink_measure = 1.0e300;
            double largest_relative_change = 0.0;
            int best_s = 0, unresolved, non_finite = 0, first_valid_s = -1, previous_valid_s = -1;
            double first_valid_step = 0.0;
#ifndef BB_AD_PRIMAL_ONLY
            double seed[NIN], tangent[NOUT], primal_from_tangent[NOUT];
            double one_sided_error = 1.0e300;   /* the tangent against the nearer of forward and backward */
            for (i = 0; i < NIN; ++i) seed[i] = 0.0;
            seed[j] = 1.0;
            run_tangent(point, point->independent, seed, primal_from_tangent, tangent);
            for (i = 0; i < NOUT; ++i) {
                if (!isfinite(tangent[i])) non_finite = 1;
            }
            for (i = 0; i < NOUT; ++i) {
                if (memcmp(&primal_from_tangent[i], &base[i], sizeof(double)) != 0) {
                    if (primal_mismatches == 0)
                        printf("     point %d: the tangent code's primal differs from the head, %s: %.17g vs %.17g\n",
                               point_index, dependent_names[i], primal_from_tangent[i], base[i]);
                    ++primal_mismatches;
                    break;
                }
            }
#endif
            for (s = 0; s < STEP_SIZE_COUNT; ++s) {
                double mismatch;
                step = relative_step_size_list[s] * fmax(fabs(value), typical_scale[j]);
                if (value - step < lower_bound[j] || value + step > upper_bound[j]) continue;   /* would cross a bound */
                if (first_valid_s < 0) {
                    first_valid_s = s;
                    first_valid_step = step;
                }
                for (i = 0; i < NIN; ++i) perturbed[i] = point->independent[i];
                perturbed[j] = value + step;
                run_head(point, perturbed, plus);
                perturbed[j] = value - step;
                run_head(point, perturbed, minus);
                for (i = 0; i < NOUT; ++i) {
                    double relative_change = fabs(plus[i] - minus[i]) / fmax(fabs(base[i]), 1.0e-300);
                    central[s][i] = (plus[i] - minus[i]) / (2.0 * step);
                    forward[i] = (plus[i] - base[i]) / step;
                    backward[i] = (base[i] - minus[i]) / step;
                    resolved[s][i] = (relative_change > RESOLUTION_THRESHOLD);
                    if (s == first_valid_s && relative_change > largest_relative_change)
                        largest_relative_change = relative_change;
                }
                mismatch = kink_mismatch(forward, backward, resolved[s], NOUT);
                if (mismatch < kink_measure) kink_measure = mismatch;
#ifdef BB_AD_PRIMAL_ONLY
                if (previous_valid_s >= 0) {
                    int both[NOUT];
                    double error;
                    for (i = 0; i < NOUT; ++i) both[i] = resolved[s][i] && resolved[previous_valid_s][i];
                    error = masked_normalized_difference(central[s], central[previous_valid_s], both, NOUT);
                    if (error < best_error) {
                        best_error = error;
                        best_step = relative_step_size_list[s];
                        best_s = s;
                    }
                }
#else
                one_sided_error = fmin(one_sided_error,
                                       fmin(masked_normalized_difference(tangent, forward, resolved[s], NOUT),
                                            masked_normalized_difference(tangent, backward, resolved[s], NOUT)));
                {
                    double error = masked_normalized_difference(tangent, central[s], resolved[s], NOUT);
                    if (error < best_error) {
                        best_error = error;
                        best_step = relative_step_size_list[s];
                        best_s = s;
                    }
                }
#endif
                previous_valid_s = s;
            }
            (void)first_valid_step;     /* used by the tangent's resolution check only */
            (void)best_s;               /* used by the tangent's failure report only */
            (void)previous_valid_s;     /* used by the primal-only self-consistency only */
            unresolved = (largest_relative_change <= RESOLUTION_THRESHOLD);
            ++pairs_overall;
            fprintf(csv, "%s,%d,%s,\"%s\",%s,%.17g,%.0e,%.6e,%.6e,%.6e,%d\n", component_names[point->component],
                    point_index, point->time_text, point->regime, independent_names[j], value, best_step, best_error,
                    kink_measure, largest_relative_change, non_finite);
#ifdef BB_AD_PRIMAL_ONLY
            /* fewer than two step sizes resolve a dependent inside the domain: nothing to compare */
            if (!unresolved && best_error >= 1.0e300) first_valid_s = -1;
#endif
            if (first_valid_s < 0) {
                ++at_bound_count;               /* every step crosses the bound: one-sided, not tested */
            } else if (non_finite) {
                ++failed;
                ++non_finite_overall;
                printf("     point %d, %s: non-finite tangent\n", point_index, independent_names[j]);
            } else if (unresolved) {
                /* (before the kink test: a change of a unit in the last place flips forward and backward
                 * differences at every step size, which is not a kink) */
                ++unresolved_count;
#ifndef BB_AD_PRIMAL_ONLY
                {   /* no resolvable change: the tangent must not predict one either */
                    double largest_predicted_change = 0.0;
                    for (i = 0; i < NOUT; ++i) {
                        double predicted = fabs(tangent[i]) * first_valid_step / fmax(fabs(base[i]), 1.0e-300);
                        if (predicted > largest_predicted_change) largest_predicted_change = predicted;
                    }
                    if (largest_predicted_change > 10.0 * RESOLUTION_THRESHOLD) {
                        ++failed;
                        printf("     point %d, %s: no resolvable change by finite differences, but the tangent predicts a\n"
                               "       relative change of %.3e over the step\n", point_index, independent_names[j],
                               largest_predicted_change);
                    }
                }
#endif
            } else if (kink_measure > KINK_THRESHOLD) {
                ++kinks;
#ifndef BB_AD_PRIMAL_ONLY
                if (one_sided_error <= ONE_SIDED_THRESHOLD) {
                    ++kinks_one_sided_overall;
                } else {
                    printf("     point %d, %s: kink; the tangent equals neither one-sided difference (best %.3e)\n",
                           point_index, independent_names[j], one_sided_error);
                }
#endif
            } else {
                if (best_error > worst) {
                    worst = best_error;
                    worst_name = independent_names[j];
                }
                if (best_error > TANGENT_PASS_THRESHOLD) {
                    ++failed;
#ifndef BB_AD_PRIMAL_ONLY
                    printf("     point %d, %s (value %.17g), step %.0e:\n", point_index, independent_names[j], value,
                           best_step);
                    for (i = 0; i < NOUT; ++i) {
                        if (resolved[best_s][i] && fabs(tangent[i] - central[best_s][i]) > 0.0)
                            printf("       %-26s tangent %.10e  central %.10e\n", dependent_names[i], tangent[i],
                                   central[best_s][i]);
                    }
#else
                    printf("     point %d, %s (value %.17g): central differences change by %.3e between step sizes\n",
                           point_index, independent_names[j], value, best_error);
#endif
                }
            }
        }
        printf("  %2d %-10s %s %-33s worst smooth error %.2e (%s), %d kinks, %d unresolved, %d at a bound, "
               "%d above threshold\n", point_index, component_names[point->component], point->time_text, point->regime,
               worst, worst_name, kinks, unresolved_count, at_bound_count, failed);
        at_bound_overall += at_bound_count;
        if (worst > worst_overall) worst_overall = worst;
        kinks_overall += kinks;
        unresolved_overall += unresolved_count;
        failed_overall += failed;
        primal_mismatch_overall += primal_mismatches;
    }
    printf("  t: %d point-independent pairs, %d kinks, %d below finite-difference resolution and %d at a bound of\n"
           "     the domain (zero snow, wetted fraction, BTRAN = 1, ...: one-sided; not tested), worst error away from kinks %.2e;\n"
           "     %d above %.0e; %d non-finite tangents\n", pairs_overall, kinks_overall, unresolved_overall,
           at_bound_overall, worst_overall, failed_overall, TANGENT_PASS_THRESHOLD, non_finite_overall);
#ifndef BB_AD_PRIMAL_ONLY
    printf("  t: at %d of the %d kinks the tangent equals one side's difference (within %.0e, the step-size scan)\n",
           kinks_one_sided_overall, kinks_overall, ONE_SIDED_THRESHOLD);
    printf("  t: the tangent code's primal differs from the head in %d point-independent pairs\n", primal_mismatch_overall);
    printf("%s\n", failed_overall == 0 && primal_mismatch_overall == 0 ? "PASS t" : "FAIL t");
#endif
    *failed_total += failed_overall;
}

#ifndef BB_AD_PRIMAL_ONLY
static double next_random(unsigned long long *state)
{
    *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
    return ((double)(*state >> 11) / 9007199254740992.0) * 2.0 - 1.0;
}

static void experiment_reverse(const TestPoint *points, int point_count)
{
    static double jacobian_tangent[NOUT][NIN], jacobian_reverse[NOUT][NIN];
    int point_index, i, j, failed = 0;
    unsigned long long random_state = 20261028ULL;
    double worst_jacobian = 0.0, worst_dot = 0.0;
    for (point_index = 0; point_index < point_count; ++point_index) {
        const TestPoint *point = &points[point_index];
        double seed[NIN], tangent[NOUT], primal[NOUT], ybar[NOUT], xbar[NIN], column_t[NOUT], column_r[NOUT];
        double jacobian_error = 0.0, left = 0.0, right = 0.0, dot_error;
        for (j = 0; j < NIN; ++j) {
            for (i = 0; i < NIN; ++i) seed[i] = 0.0;
            seed[j] = 1.0;
            run_tangent(point, point->independent, seed, primal, tangent);
            for (i = 0; i < NOUT; ++i) jacobian_tangent[i][j] = tangent[i];
        }
        for (i = 0; i < NOUT; ++i) {
            int k;
            for (k = 0; k < NOUT; ++k) ybar[k] = 0.0;
            ybar[i] = 1.0;
            run_adjoint(point, point->independent, ybar, xbar);
            for (j = 0; j < NIN; ++j) jacobian_reverse[i][j] = xbar[j];
        }
        for (j = 0; j < NIN; ++j) {
            double error;
            for (i = 0; i < NOUT; ++i) {
                column_t[i] = jacobian_tangent[i][j];
                column_r[i] = jacobian_reverse[i][j];
            }
            error = normalized_difference(column_t, column_r, NOUT);
            if (!(error <= 1.0)) error = 1.0;      /* NaN */
            if (error > jacobian_error) jacobian_error = error;
        }
        for (j = 0; j < NIN; ++j) seed[j] = next_random(&random_state);
        for (i = 0; i < NOUT; ++i) ybar[i] = next_random(&random_state);
        run_tangent(point, point->independent, seed, primal, tangent);
        run_adjoint(point, point->independent, ybar, xbar);
        for (i = 0; i < NOUT; ++i) left += ybar[i] * tangent[i];
        for (j = 0; j < NIN; ++j) right += xbar[j] * seed[j];
        dot_error = fabs(left - right) / fmax(fmax(fabs(left), fabs(right)), 1.0e-300);
        if (!(dot_error <= 1.0)) dot_error = 1.0;
        printf("  %2d %-10s %s %-33s Jacobian reverse vs tangent %.2e, dot product %.2e\n", point_index,
               component_names[point->component], point->time_text, point->regime, jacobian_error, dot_error);
        if (jacobian_error > REVERSE_PASS_THRESHOLD) {
            /* the Jacobian entries that differ most, relative to the largest of their column, so the failing
             * part of the reverse code can be found */
            int reported_count;
            static int entry_reported[NOUT][NIN];
            for (i = 0; i < NOUT; ++i) {
                for (j = 0; j < NIN; ++j) entry_reported[i][j] = 0;
            }
            for (reported_count = 0; reported_count < REVERSE_ENTRIES_REPORTED; ++reported_count) {
                int worst_i = -1, worst_j = -1;
                double worst_entry = REVERSE_PASS_THRESHOLD;
                for (j = 0; j < NIN; ++j) {
                    double column_largest = 0.0;
                    for (i = 0; i < NOUT; ++i) {
                        column_largest = fmax(column_largest, fabs(jacobian_tangent[i][j]));
                        column_largest = fmax(column_largest, fabs(jacobian_reverse[i][j]));
                    }
                    for (i = 0; i < NOUT; ++i) {
                        double entry_error = fabs(jacobian_tangent[i][j] - jacobian_reverse[i][j]) /
                                             fmax(column_largest, 1.0e-300);
                        if (!entry_reported[i][j] && entry_error > worst_entry) {
                            worst_entry = entry_error;
                            worst_i = i;
                            worst_j = j;
                        }
                    }
                }
                if (worst_i >= 0) {
                    entry_reported[worst_i][worst_j] = 1;
                    printf("       d %-26s / d %-24s tangent %+.10e  reverse %+.10e\n", dependent_names[worst_i],
                           independent_names[worst_j], jacobian_tangent[worst_i][worst_j],
                           jacobian_reverse[worst_i][worst_j]);
                }
            }
        }
        if (jacobian_error > worst_jacobian) worst_jacobian = jacobian_error;
        if (dot_error > worst_dot) worst_dot = dot_error;
        if (jacobian_error > REVERSE_PASS_THRESHOLD || dot_error > REVERSE_PASS_THRESHOLD) ++failed;
    }
    printf("  a: worst Jacobian difference %.2e, worst dot-product difference %.2e (pass below %.0e)\n%s\n", worst_jacobian,
           worst_dot, REVERSE_PASS_THRESHOLD, failed == 0 ? "PASS a" : "FAIL a");
}
#endif

/* ======================================================================================
 * main
 * ====================================================================================== */
int main(int argc, char **argv)
{
    static TestPoint points[MAXIMUM_TEST_POINTS];
    static TestPoint point;
    static pet_model prepared;              /* a copy of the model, prepared for the step */
    static NoahmpStabilityInput label_input;
    static NoahmpStabilityOutput label_solution;
    static NoahmpVegeFluxExactTrace label_trace;
    BallBerryStepResults label_results;
    Bmi *model;
    pet_model *pet;
    CsvReader rows;
    RegimeList regimes;
    LabelList labels;
    const char *experiments;
    char path[2048];
    FILE *csv;
    int hour = 0, point_count = 0, mismatch_hours = 0, reported = 0, failed_total = 0;

    if (argc != 5) {
        fprintf(stderr, "usage: %s CONFIG CANOPY_ENERGY_CSV RESULTS_DIRECTORY EXPERIMENTS\n", argv[0]);
        return 2;
    }
    experiments = argv[4];
    memset(&regimes, 0, sizeof(regimes));
    model = initialize_model(argv[1]);
    pet = (pet_model *)model->data;
    if (pet->pet_method != 6 || pet->bb_noahmp_exact_iteration != 1) {
        fprintf(stderr, "ERROR: the configuration is not PET method 6 with the exact solver\n");
        return 1;
    }
    csv_open(&rows, argv[2]);

    while (csv_next(&rows)) {
        double head_dependents[NOUT], production_dependents[NOUT];
        const char *regime;
        set_hour_inputs(pet, &rows);
        /* the preparation, on a copy of the model */
        prepared = *pet;
        memset(&point, 0, sizeof(point));
        if (bb6_prepare_step_inputs(&prepared, &point.step) != 0) {
            fprintf(stderr, "ERROR: the preparation failed at %s\n", rows.time_text);
            return 1;
        }
        point.component = 0;
        snprintf(point.time_text, sizeof(point.time_text), "%s", rows.time_text);
        point.parameters = pet->bb_parameters;
        pack_independents(&point.parameters, &point.step, point.independent);
        /* the exact solver's solution and trace, for the regime labels */
        ball_berry_method6_physics_step(&point.parameters, &point.step, &label_input, &label_solution, &label_trace,
                                        &label_results);

        if (ball_berry_method6_core_step(pet) != 0) {
            fprintf(stderr, "ERROR: ball_berry_method6_core_step failed at %s\n", rows.time_text);
            return 1;
        }
        pet->bmi.current_time_step += pet->bmi.time_step_size_s;   /* as the BMI update() and the stack's adapter */
        pet->bmi.current_step += 1;
        pet->bmi.current_time += pet->bmi.time_step_size_s;
        pack_production_dependents(pet, production_dependents);

        if (strchr(experiments, 'p') != NULL) {
            run_head(&point, point.independent, head_dependents);
            if (memcmp(head_dependents, production_dependents, sizeof(head_dependents)) != 0) {
                if (!reported) {
                    int i;
                    for (i = 0; i < NOUT; ++i) {
                        if (memcmp(&head_dependents[i], &production_dependents[i], sizeof(double)) != 0) {
                            printf("  %s, %s: head %.17g, production %.17g\n", rows.time_text, dependent_names[i],
                                   head_dependents[i], production_dependents[i]);
                            break;
                        }
                    }
                    reported = 1;
                }
                ++mismatch_hours;
            }
        }
        hour_labels(&point.step, pet, &label_trace, &labels);
        if (hour >= SPIN_UP_HOURS) {
            int take = 0;
            regime = first_new_label(&regimes, &labels);
            if (regime != NULL) {
                snprintf(point.regime, sizeof(point.regime), "%s", regime);
                take = 1;
            } else if (hour % REGULAR_TEST_STRIDE == 0) {
                snprintf(point.regime, sizeof(point.regime), "every %d hours", REGULAR_TEST_STRIDE);
                take = 1;
            }
            if (take && point_count < MAXIMUM_TEST_POINTS) points[point_count++] = point;
        }
        ++hour;
    }
    printf("ball_berry_ad_test: %d hours of %s, %d test states\n", hour, argv[2], point_count);
    {
        int k;
        printf("  regimes found:");
        for (k = 0; k < regimes.found_count; ++k) printf("%s %s", k == 0 ? "" : ";", regimes.found[k]);
        printf("\n");
    }

    if (strchr(experiments, 'p') != NULL) {
        printf("\n==== p: the head against the production step, every hour, every dependent\n");
        printf("  %d hours with differences\n%s\n", mismatch_hours, mismatch_hours == 0 ? "PASS p" : "FAIL p");
    }
    if (strchr(experiments, 't') != NULL) {
        snprintf(path, sizeof(path), "%s/tangent_errors.csv", argv[3]);
        csv = fopen(path, "w");
        if (csv == NULL) {
            fprintf(stderr, "ERROR: cannot open %s\n", path);
            return 1;
        }
        fprintf(csv, "component,point,time,regime,independent,value,best_relative_step,best_error,kink_forward_backward,"
                     "largest_relative_change_at_largest_step,non_finite_tangent\n");
#ifdef BB_AD_PRIMAL_ONLY
        printf("\n==== t: central differences only (built without the generated code): the error is the change of the\n"
               "     central difference between successive step sizes\n");
#else
        printf("\n==== t: tangent against central differences (column-normalized; pass below %.0e)\n", TANGENT_PASS_THRESHOLD);
#endif
        experiment_tangent(points, point_count, csv, &failed_total);
        fclose(csv);
        printf("  details: %s\n", path);
    }
#ifndef BB_AD_PRIMAL_ONLY
    if (strchr(experiments, 'a') != NULL) {
        printf("\n==== a: reverse against tangent: full Jacobian, and the dot-product test\n");
        experiment_reverse(points, point_count);
    }
#else
    if (strchr(experiments, 'a') != NULL) printf("\n==== a: skipped (built without the generated code)\n");
#endif
    model->finalize(model);
    return 0;
}
