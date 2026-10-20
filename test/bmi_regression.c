/* FLO 10/26: Exercise invalid BMI arguments and timestep behavior so
 * interface regressions fail independently of PET forcing files. */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../include/bmi.h"
#include "../include/bmi_pet.h"
#include "../include/pet.h"

#define CHECK(expression) do { if (!(expression)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #expression); \
    failures++; } } while (0)

int main(int argc, char **argv)
{
    Bmi model = {0};
    int failures = 0;
    int nbytes = 0;
    int index = 0;
    int invalid_index = 1;
    double value = 0.0;
    double timestep = 0.0;
    double now = 0.0;
    const char *humidity = "atmosphere_air_water~vapor__specific_humidity";
    const char *output = "water_potential_evaporation_flux";

    if (argc != 2 || !register_bmi_pet(&model)) {
        fprintf(stderr, "Usage: bmi_regression CONFIG\n");
        return 1;
    }
    if (model.initialize(&model, argv[1]) != BMI_SUCCESS) {
        fprintf(stderr, "FAIL: initialize\n");
        model.finalize(&model);
        return 1;
    }
    CHECK(model.get_var_nbytes(&model, humidity, &nbytes) == BMI_SUCCESS);
    CHECK(nbytes == (int)sizeof(double));
    CHECK(model.get_var_nbytes(&model, output, &nbytes) == BMI_SUCCESS);
    CHECK(nbytes == (int)sizeof(double));
    CHECK(model.get_var_nbytes(&model, "unknown", &nbytes) == BMI_FAILURE);
    CHECK(model.get_var_nbytes(&model, humidity, NULL) == BMI_FAILURE);
    CHECK(model.get_value(&model, humidity, NULL) == BMI_FAILURE);
    CHECK(model.get_value_at_indices(&model, humidity, &value,
                                      &invalid_index, 1) == BMI_FAILURE);
    CHECK(model.get_value_at_indices(&model, NULL, &value,
                                      &index, 1) == BMI_FAILURE);
    CHECK(model.set_value_at_indices(&model, NULL, &index,
                                      1, &value) == BMI_FAILURE);
    CHECK(model.set_value_at_indices(&model, humidity, &invalid_index,
                                      1, &value) == BMI_FAILURE);
    CHECK(model.set_value(&model, output, &value) == BMI_FAILURE);
    value = 0.0088;
    CHECK(model.set_value(&model, humidity, &value) == BMI_SUCCESS);
    value = 0.0;
    CHECK(model.get_value_at_indices(&model, humidity, &value,
                                      &index, 1) == BMI_SUCCESS);
    CHECK(fabs(value - 0.0088) < 1e-12);
    CHECK(model.get_time_step(&model, &timestep) == BMI_SUCCESS);
    CHECK(model.get_current_time(&model, &now) == BMI_SUCCESS);
    CHECK(model.update_until(&model, now + timestep / 2.0) == BMI_FAILURE);
    {
        double after = -1.0;
        CHECK(model.get_current_time(&model, &after) == BMI_SUCCESS);
        CHECK(after == now);
    }
    /* FLO 10/26: The previous run_pet code changed the configured wind
     * measurement height after its first update. Verify that all five PET
     * methods give the same answer for physically equivalent 10 m and 2 m
     * winds, and that repeated 10 m updates do not lose normalization.
     */
    {
        pet_model *pet = (pet_model *)model.data;
        const double original_height = pet->pet_params.wind_speed_measurement_height_m;
        const double displacement = pet->pet_params.zero_plane_displacement_height_m;
        const double speed_at_10m = 4.0;
        const double speed_at_2m = speed_at_10m *
            log(2.0 / displacement) / log(10.0 / displacement);
        pet->aorc.air_temperature_2m_K = 295.0;
        pet->aorc.surface_pressure_Pa = 100000.0;
        pet->aorc.specific_humidity_2m_kg_per_kg = 0.008;
        pet->aorc.incoming_shortwave_W_per_m2 = 600.0;
        pet->aorc.incoming_longwave_W_per_m2 = 350.0;
        pet->aorc.v_wind_speed_10m_m_per_s = 0.0;
        for (int method = 1; method <= 5; ++method) {
            double reference_pet, repeated_pet;
            pet->pet_method = method;
            pet->pet_params.wind_speed_measurement_height_m = 10.0;
            pet->aorc.u_wind_speed_10m_m_per_s = speed_at_10m;
            CHECK(run_pet(pet) == 0);
            reference_pet = pet->pet_m_per_s;
            CHECK(isfinite(reference_pet));
            CHECK(pet->pet_params.wind_speed_measurement_height_m == 10.0);
            CHECK(fabs(pet->pet_forcing.wind_speed_m_per_s - speed_at_2m) < 1e-12);
            CHECK(run_pet(pet) == 0);
            repeated_pet = pet->pet_m_per_s;
            CHECK(fabs(repeated_pet - reference_pet) < 1e-12);
            pet->pet_params.wind_speed_measurement_height_m = 2.0;
            pet->aorc.u_wind_speed_10m_m_per_s = speed_at_2m;
            CHECK(run_pet(pet) == 0);
            CHECK(fabs(pet->pet_m_per_s - reference_pet) < 1e-12);
            CHECK(pet->pet_params.wind_speed_measurement_height_m == 2.0);
        }
        pet->pet_params.wind_speed_measurement_height_m = original_height;
    }
    CHECK(model.finalize(&model) == BMI_SUCCESS);
    if (failures) {
        fprintf(stderr, "BMI regression: %d failures\n", failures);
        return 1;
    }
    puts("BMI regression: PASS");
    return 0;
}
