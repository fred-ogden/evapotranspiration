// =========================================================================================
// test_ball_berry_canopy_energy.c
//
// Tests of the REVISED canopy energy physics of PET method 6 (canopy_energy_physics_option = 2,
// configs/ball_berry_revised_canopy_energy.cfg), driven through the BMI. ctest
// ball_berry_revised_canopy_energy. Prints one line per failed check and a summary; returns 0
// only if every check passed.
//
//   1. energy closure: at the solution, SAV - IRC - SHC - EVC - TR + PAHV + fusion - storage = 0
//   2. carried state: the next step starts from this step's TV, TAH, EAH (BMI *_start_of_step_*)
//   3. freezing night: liquid on a canopy that cools below TFRZ freezes; the fusion heat holds TV
//      near TFRZ (the same night with a dry canopy is colder); water is conserved
//   4. melting: ice on a canopy warming above TFRZ melts; the canopy is held near TFRZ
//   5. frozen canopy: transpiration ramps to almost zero below TFRZ
//   6. smoothness: TV and the freezing rate change smoothly with the advected heat PAHV through
//      the freeze (no jumps), as gradient-based parameter estimation needs
//
// Usage: test_ball_berry_canopy_energy configs/ball_berry_revised_canopy_energy.cfg
//
// FLO 10/24: first version.
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
static const double freezing_point_k = 273.16;              // TFRZ
static const double latent_heat_fusion_j_per_kg = 0.3336e6; // HFUS

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

// A winter or summer hour, 2016-01-14 or 2016-07-01 at the given hour (monthly LAI table).
typedef struct {
    int    month, hour;
    double shortwave_w_per_m2, longwave_w_per_m2, air_temperature_k, relative_humidity, wind_m_per_s;
} HourForcing;

static void set_forcing(Bmi *model, const HourForcing *f)
{
    // specific humidity from relative humidity over water (Noah-MP ESAT form is not needed here)
    double saturation_vapor_pressure_pa = 611.2 * exp(17.67 * (f->air_temperature_k - 273.15) /
                                                      (f->air_temperature_k - 29.65));
    double vapor_pressure_pa = f->relative_humidity * saturation_vapor_pressure_pa;
    double surface_pressure_pa = 100000.0;
    double specific_humidity = 0.622 * vapor_pressure_pa / (surface_pressure_pa - 0.378 * vapor_pressure_pa);
    set_double(model, "land_surface_radiation~incoming~longwave__energy_flux", f->longwave_w_per_m2);
    set_double(model, "land_surface_air__pressure", surface_pressure_pa);
    set_double(model, "atmosphere_air_water~vapor__specific_humidity", specific_humidity);
    set_double(model, "land_surface_radiation~incoming~shortwave__energy_flux", f->shortwave_w_per_m2);
    set_double(model, "land_surface_air__temperature", f->air_temperature_k);
    set_double(model, "land_surface_wind__x_component_of_velocity", f->wind_m_per_s);
    set_double(model, "land_surface_wind__y_component_of_velocity", 0.0);
    set_double(model, "time_year", 2016.0);
    set_double(model, "time_month", (double)f->month);
    if (f->month == 1) {
        set_double(model, "time_day_of_month", 14.0);
    } else {
        set_double(model, "time_day_of_month", 1.0);
    }
    set_double(model, "time_hour", (double)f->hour);
    set_double(model, "time_minute", 0.0);
    set_double(model, "time_second", 0.0);
}

static void set_canopy_water(Bmi *model, double liquid_mm, double ice_mm, double wetted_fraction)
{
    set_double(model, "canopy_liquid_water_mm", liquid_mm);
    set_double(model, "canopy_ice_mm", ice_mm);
    set_double(model, "canopy_wetted_fraction", wetted_fraction);
}

static void set_start_state(Bmi *model, double vegetation_temperature_k, double air_temperature_k,
                            double vapor_pressure_pa)
{
    set_double(model, "vegetation_temperature_start_of_step_K", vegetation_temperature_k);
    set_double(model, "canopy_air_temperature_start_of_step_K", air_temperature_k);
    set_double(model, "canopy_air_vapor_pressure_start_of_step_Pa", vapor_pressure_pa);
}

// net heating of the canopy at the solution, from the unlinearized fluxes at the final TV
static double energy_residual(Bmi *model)
{
    pet_model *pet = (pet_model *)model->data;
    return pet->bb_last_output.canopy.energy_balance_residual_w_per_m2;
}

int main(int argument_count, char **argument_list)
{
    Bmi model;
    char description[256];
    // marginal night: the canopy's net cooling is comparable to the fusion heat of the intercepted water
    HourForcing winter_night = { 1, 6, 0.0, 285.0, 272.3, 0.95, 1.5 };
    HourForcing winter_day = { 1, 14, 450.0, 280.0, 280.0, 0.60, 3.0 };
    HourForcing summer_day = { 7, 12, 600.0, 361.6, 288.0, 0.70, 2.4 };
    double tv_wet, tv_dry, freeze_rate, fusion_heat, melt_rate, tv_melt;
    int index;

    if (argument_count != 2) {
        fprintf(stderr, "usage: %s configs/ball_berry_revised_canopy_energy.cfg\n", argument_list[0]);
        return 1;
    }
    memset(&model, 0, sizeof(model));
    register_bmi_pet(&model);          /* allocates model.data */
    if (model.initialize(&model, argument_list[1]) != BMI_SUCCESS) {
        fprintf(stderr, "ERROR: initialize failed\n");
        return 1;
    }

    // ---- 1 and 2: closure and carried state, two summer steps
    printf("1-2. energy closure and carried canopy state\n");
    set_forcing(&model, &summer_day);
    set_canopy_water(&model, 0.0, 0.0, 0.0);
    check(model.update(&model) == BMI_SUCCESS, "summer update 1");
    snprintf(description, sizeof(description), "energy closure at the solution (%.2e W m-2)", energy_residual(&model));
    check(fabs(energy_residual(&model)) < 1.0e-4, description);
    check(get_double(&model, "vegetation_temperature_start_of_step_K") == get_double(&model, "vegetation_temperature_K") &&
          get_double(&model, "canopy_air_temperature_start_of_step_K") == get_double(&model, "canopy_air_temperature_K") &&
          get_double(&model, "canopy_air_vapor_pressure_start_of_step_Pa") == get_double(&model, "canopy_air_vapor_pressure_Pa"),
          "end-of-step TV, TAH, EAH written to the start-of-step inputs");
    {
        double tv_first = get_double(&model, "vegetation_temperature_K");
        double storage_second;
        check(model.update(&model) == BMI_SUCCESS, "summer update 2 (same forcing)");
        storage_second = get_double(&model, "canopy_heat_storage_flux_W_per_m2");
        printf("  TV %.4f K then %.4f K; heat storage in step 2 %.4f W m-2\n", tv_first,
               get_double(&model, "vegetation_temperature_K"), storage_second);
        check(fabs(storage_second - 1.0e-3 * 0.0) < 1.0 &&
              fabs(get_double(&model, "vegetation_temperature_K") - tv_first) < 0.05,
              "same forcing again: TV nearly steady, small heat storage");
    }

    // ---- 3: freezing night, wet against dry canopy, from the same start state
    printf("3. freezing night\n");
    set_forcing(&model, &winter_night);
    set_start_state(&model, 273.6, 272.0, 400.0);
    set_canopy_water(&model, 0.0, 0.0, 0.0);
    check(model.update(&model) == BMI_SUCCESS, "winter night, dry canopy");
    tv_dry = get_double(&model, "vegetation_temperature_K");
    set_start_state(&model, 273.6, 272.0, 400.0);
    set_canopy_water(&model, 0.40, 0.0, 0.60);
    check(model.update(&model) == BMI_SUCCESS, "winter night, 0.4 mm of liquid on the canopy");
    tv_wet = get_double(&model, "vegetation_temperature_K");
    freeze_rate = get_double(&model, "canopy_freezing_rate_mm_per_s");
    fusion_heat = get_double(&model, "canopy_fusion_heat_flux_W_per_m2");
    printf("  TV dry %.3f K, wet %.3f K; frozen %.4f mm in the hour, fusion heat %.2f W m-2, closure %.1e\n",
           tv_dry, tv_wet, freeze_rate * 3600.0, fusion_heat, energy_residual(&model));
    check(freeze_rate > 0.0 && freeze_rate * 3600.0 <= 0.40 + 1.0e-12, "liquid freezes, at most what is there");
    check(fabs(fusion_heat - latent_heat_fusion_j_per_kg * freeze_rate) < 1.0e-9, "fusion heat = HFUS * rate");
    check(tv_wet > tv_dry + 0.1, "the fusion heat keeps the wet canopy warmer than the dry one");
    check(fabs(tv_wet - freezing_point_k) < 1.0, "the freezing canopy stays within 1 K of TFRZ");
    check(fabs(energy_residual(&model)) < 1.0e-4, "energy closure with freezing");

    // ---- 4: melting day
    printf("4. melting day\n");
    set_forcing(&model, &winter_day);
    set_start_state(&model, 272.8, 273.5, 500.0);
    set_canopy_water(&model, 0.0, 1.0, 0.30);
    check(model.update(&model) == BMI_SUCCESS, "winter day, 1 mm of ice on the canopy");
    melt_rate = -get_double(&model, "canopy_freezing_rate_mm_per_s");
    tv_melt = get_double(&model, "vegetation_temperature_K");
    printf("  melted %.4f mm in the hour, TV %.3f K, air 280 K, closure %.1e\n", melt_rate * 3600.0, tv_melt,
           energy_residual(&model));
    check(melt_rate > 0.0 && melt_rate * 3600.0 <= 1.0 + 1.0e-12, "ice melts, at most what is there");
    check(fabs(energy_residual(&model)) < 1.0e-4, "energy closure with melting");

    // ---- 5: transpiration of a frozen canopy, against the same hour 8 K warmer
    printf("5. frozen-canopy transpiration ramp\n");
    {
        HourForcing cold_sunny = { 1, 13, 400.0, 240.0, 266.0, 0.60, 2.0 };
        HourForcing mild_sunny = { 1, 13, 400.0, 240.0, 284.0, 0.60, 2.0 };
        double transpiration_cold, transpiration_mild;
        set_forcing(&model, &cold_sunny);
        set_start_state(&model, 266.5, 266.0, 300.0);
        set_canopy_water(&model, 0.0, 0.0, 0.0);
        check(model.update(&model) == BMI_SUCCESS, "cold sunny hour");
        transpiration_cold = get_double(&model, "canopy_transpiration_flux_m_per_s");
        printf("  cold: TV %.2f K, transpiration %.3e m s-1\n", get_double(&model, "vegetation_temperature_K"),
               transpiration_cold);
        set_forcing(&model, &mild_sunny);
        set_start_state(&model, 284.5, 284.0, 800.0);
        check(model.update(&model) == BMI_SUCCESS, "mild sunny hour");
        transpiration_mild = get_double(&model, "canopy_transpiration_flux_m_per_s");
        printf("  mild: TV %.2f K, transpiration %.3e m s-1\n", get_double(&model, "vegetation_temperature_K"),
               transpiration_mild);
        check(get_double(&model, "vegetation_temperature_K") > freezing_point_k, "mild hour is above TFRZ");
        check(fabs(transpiration_cold) < 0.01 * transpiration_mild, "frozen canopy: transpiration below 1 percent");
    }

    // ---- 6: smoothness through the freeze: TV and the freezing rate against PAHV
    printf("6. smooth response to the advected heat through the freeze\n");
    {
        const int number_of_points = 161;
        double pahv, tv[161], rate[161], largest_tv_step = 0.0, largest_rate_step = 0.0;
        double largest_tv_curvature = 0.0;
        int monotone = 1, all_converged = 1;
        for (index = 0; index < number_of_points; ++index) {
            pahv = -40.0 + 0.5 * index;
            set_forcing(&model, &winter_night);
            set_start_state(&model, 273.6, 272.0, 400.0);
            set_canopy_water(&model, 0.40, 0.0, 0.60);
            set_double(&model, "precipitation_advected_heat_to_canopy_W_per_m2", pahv);
            if (model.update(&model) != BMI_SUCCESS) {
                all_converged = 0;
                break;
            }
            tv[index] = get_double(&model, "vegetation_temperature_K");
            rate[index] = get_double(&model, "canopy_freezing_rate_mm_per_s");
        }
        set_double(&model, "precipitation_advected_heat_to_canopy_W_per_m2", 0.0);
        check(all_converged, "solver converges for every PAHV from -40 to +40 W m-2");
        if (all_converged) {
            for (index = 1; index < number_of_points; ++index) {
                if (!(tv[index] > tv[index - 1])) monotone = 0;
                if (fabs(tv[index] - tv[index - 1]) > largest_tv_step) largest_tv_step = fabs(tv[index] - tv[index - 1]);
                if (fabs(rate[index] - rate[index - 1]) * 3600.0 > largest_rate_step)
                    largest_rate_step = fabs(rate[index] - rate[index - 1]) * 3600.0;
                if (index < number_of_points - 1) {
                    double curvature = fabs(tv[index + 1] - 2.0 * tv[index] + tv[index - 1]);
                    if (curvature > largest_tv_curvature) largest_tv_curvature = curvature;
                }
            }
            printf("  PAHV -40..+40 W m-2: TV %.3f..%.3f K, frozen %.4f..%.4f mm/h\n", tv[0], tv[number_of_points - 1],
                   rate[0] * 3600.0, rate[number_of_points - 1] * 3600.0);
            printf("  per 0.5 W m-2: largest TV change %.2e K, largest freezing change %.2e mm/h, "
                   "largest TV second difference %.2e K\n", largest_tv_step, largest_rate_step, largest_tv_curvature);
            check(monotone, "TV increases monotonically with PAHV");
            check(largest_tv_step < 0.05 && largest_rate_step < 0.01 && largest_tv_curvature < 0.01,
                  "no jumps in TV or the freezing rate");
        }
    }

    model.finalize(&model);
    printf("%d checks, %d failed\n", number_of_checks, number_of_failed_checks);
    if (number_of_failed_checks == 0) {
        return 0;
    } else {
        return 1;
    }
}
