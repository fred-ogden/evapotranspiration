// =========================================================================================
// test_ball_berry_noahmp_options.c
//
// Tests of the PET method 6 options of 10/26 for the decoupled Noah-MP stack, driven through
// the BMI with configs/ball_berry_noahmp_cat-18560.cfg. ctest ball_berry_noahmp_options.
//
//   1. solar geometry from the time inputs: COSZ near noon in July at 42 N, the canopy
//      transpires and the ground absorbs shortwave; at night COSZ < 0, the update succeeds and
//      nothing is absorbed (Noah-MP ALBEDO: IF(COSZ <= 0) GOTO 100)
//   2. surface resistance with snow (OPT_RSF = 4): 1/(FSNO/RSURF_SNOW + (1-FSNO)/RSURF)
//   3. vegetated fraction: DVEG = 1 floor of 0.05
//   4. canopy buried by snow (ELAI + ESAI = 0): no transpiration, FVEG = 0, TAH = air temperature,
//      SAG = the shortwave absorbed by the snow-free and snow ground without a canopy
//   5. soil albedo from the top-disc soil moisture (GROUNDALB): wetter soil absorbs more
//
// Usage: test_ball_berry_noahmp_options configs/ball_berry_noahmp_cat-18560.cfg
//
// FLO 10/26: first version.
// =========================================================================================
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/pet.h"
#include "../include/bmi.h"
#include "../include/bmi_pet.h"

static int number_of_checks = 0;
static int number_of_failed_checks = 0;

static void check(int condition_is_true, const char *description)
{
    ++number_of_checks;
    if (!condition_is_true) {
        ++number_of_failed_checks;
        printf("  FAIL: %s\n", description);
    }
}

static void set_double(Bmi *model, const char *name, double value)
{
    if (model->set_value(model, name, &value) != BMI_SUCCESS) {
        fprintf(stderr, "ERROR: set_value %s = %g failed\n", name, value);
        exit(1);
    }
}

static double get_double(Bmi *model, const char *name)
{
    double value = NAN;
    if (model->get_value(model, name, &value) != BMI_SUCCESS) {
        fprintf(stderr, "ERROR: get_value %s failed\n", name);
        exit(1);
    }
    return value;
}

// One hour of 2015-07-01 (UTC hour as given): clear summer day or night forcing.
static void set_forcing(Bmi *model, int hour_utc, double shortwave_w_per_m2)
{
    set_double(model, "land_surface_radiation~incoming~longwave__energy_flux", 350.0);
    set_double(model, "land_surface_air__pressure", 100000.0);
    set_double(model, "atmosphere_air_water~vapor__specific_humidity", 0.010);
    set_double(model, "land_surface_radiation~incoming~shortwave__energy_flux", shortwave_w_per_m2);
    set_double(model, "land_surface_air__temperature", 296.0);
    set_double(model, "land_surface_wind__x_component_of_velocity", 2.5);
    set_double(model, "land_surface_wind__y_component_of_velocity", 0.5);
    set_double(model, "time_year", 2015.0);
    set_double(model, "time_month", 7.0);
    set_double(model, "time_day_of_month", 1.0);
    set_double(model, "time_hour", (double)hour_utc);
    set_double(model, "time_minute", 0.0);
    set_double(model, "time_second", 0.0);
    set_double(model, "ground_skin_temperature_K", 297.0);
    set_double(model, "soil_top_disc_moisture_content", 0.25);
}

int main(int argument_count, char **argument_list)
{
    Bmi model;
    double cosz, ground_shortwave, transpiration, wet_ground_shortwave, dry_ground_shortwave;

    if (argument_count != 2) {
        fprintf(stderr, "usage: %s configs/ball_berry_noahmp_cat-18560.cfg\n", argument_list[0]);
        return 1;
    }
    memset(&model, 0, sizeof(model));
    register_bmi_pet(&model);
    if (model.initialize(&model, argument_list[1]) != BMI_SUCCESS) {
        fprintf(stderr, "ERROR: initialize failed\n");
        return 1;
    }

    // ---- 1. solar geometry: 16 UTC is about 11:10 local solar time at 72 W (end of step 17 UTC)
    printf("1. solar geometry from the time inputs\n");
    set_forcing(&model, 16, 850.0);
    check(model.update(&model) == BMI_SUCCESS, "daytime update");
    cosz = get_double(&model, "cosine_solar_zenith_angle");
    ground_shortwave = get_double(&model, "ground_shortwave_absorbed_below_canopy_W_per_m2");
    transpiration = get_double(&model, "canopy_transpiration_flux_m_per_s");
    printf("   noon: COSZ %.4f, SAG %.2f W m-2, ETRAN %.3e m s-1\n", cosz, ground_shortwave, transpiration);
    check(cosz > 0.85 && cosz < 1.0, "COSZ near noon in July at 42 N (sun about 70 degrees high)");
    check(ground_shortwave > 0.0 && transpiration > 0.0, "daytime absorption and transpiration");
    set_forcing(&model, 6, 0.0);
    check(model.update(&model) == BMI_SUCCESS, "night update (COSZ <= 0) succeeds");
    printf("   night: COSZ %.4f, SAG %.2f W m-2\n", get_double(&model, "cosine_solar_zenith_angle"),
           get_double(&model, "ground_shortwave_absorbed_below_canopy_W_per_m2"));
    check(get_double(&model, "cosine_solar_zenith_angle") < 0.0, "night COSZ < 0");
    check(get_double(&model, "ground_shortwave_absorbed_below_canopy_W_per_m2") == 0.0, "nothing absorbed at night");

    // ---- 2. surface resistance with snow
    printf("2. surface resistance with snow (OPT_RSF = 4)\n");
    set_forcing(&model, 16, 850.0);
    set_double(&model, "surface_resistance_s_per_m", 200.0);
    set_double(&model, "snow_cover_fraction", 0.5);
    model.update(&model);
    check(fabs(get_double(&model, "surface_resistance_with_snow_s_per_m") - 80.0) < 1e-9,
          "RSURF = 1/(0.5/50 + 0.5/200) = 80 s/m");
    set_double(&model, "snow_cover_fraction", 0.0);

    // ---- 3. vegetated fraction
    printf("3. vegetated fraction (DVEG = 1)\n");
    check(fabs(get_double(&model, "vegetated_fraction") - 0.8) < 1e-15, "FVEG = SHDFAC = 0.8");
    set_double(&model, "SHDFAC", 0.01);
    model.update(&model);
    check(fabs(get_double(&model, "vegetated_fraction") - 0.05) < 1e-15, "FVEG floor 0.05");
    set_double(&model, "SHDFAC", 0.8);

    // ---- 4. canopy buried by snow (20 m of snow over a 16 m canopy, a code test)
    printf("4. canopy buried by snow\n");
    set_double(&model, "snow_depth_m", 20.0);
    set_double(&model, "snow_cover_fraction", 1.0);
    set_double(&model, "snow_albedo_direct_visible", 0.9);
    set_double(&model, "snow_albedo_direct_nir", 0.6);
    set_double(&model, "snow_albedo_diffuse_visible", 0.9);
    set_double(&model, "snow_albedo_diffuse_nir", 0.6);
    check(model.update(&model) == BMI_SUCCESS, "update with the canopy buried");
    check(get_double(&model, "canopy_transpiration_flux_m_per_s") == 0.0, "no transpiration");
    check(get_double(&model, "vegetated_fraction") == 0.0, "FVEG = 0");
    check(get_double(&model, "canopy_air_temperature_K") == 296.0, "TAH = air temperature");
    ground_shortwave = get_double(&model, "ground_shortwave_absorbed_below_canopy_W_per_m2");
    printf("   SAG on snow, no canopy: %.2f W m-2\n", ground_shortwave);
    check(ground_shortwave > 0.0 && ground_shortwave < 0.4 * 850.0, "snow absorbs a small part of the shortwave");

    // ---- 5. soil albedo from soil moisture (snow-free, no canopy)
    printf("5. soil albedo from the top-disc soil moisture\n");
    set_double(&model, "snow_cover_fraction", 0.0);
    set_double(&model, "soil_top_disc_moisture_content", 0.05);
    model.update(&model);
    dry_ground_shortwave = get_double(&model, "ground_shortwave_absorbed_below_canopy_W_per_m2");
    set_double(&model, "soil_top_disc_moisture_content", 0.35);
    model.update(&model);
    wet_ground_shortwave = get_double(&model, "ground_shortwave_absorbed_below_canopy_W_per_m2");
    printf("   SAG dry %.2f, wet %.2f W m-2\n", dry_ground_shortwave, wet_ground_shortwave);
    check(wet_ground_shortwave > dry_ground_shortwave, "wetter soil absorbs more");

    model.finalize(&model);
    printf("%d of %d checks passed\n", number_of_checks - number_of_failed_checks, number_of_checks);
    return (number_of_failed_checks == 0) ? 0 : 1;
}
