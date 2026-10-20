#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
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

    /* FLO 10/26: Check both model updates as well as BMI transfers. A failed
     * AORC update must never lead to reusing stale forcing in PET.
     */
    for (int i = 0; i <= 100; ++i) {
        if (aorc_bmi_model->update(aorc_bmi_model) != BMI_SUCCESS) {
            fprintf(stderr, "AORC update failed at step %d\n", i);
            goto cleanup;
        }
        if (pass_forcing_from_aorc_to_pet(pet_bmi_model, aorc_bmi_model) != BMI_SUCCESS)
            goto cleanup;
        if (pet_bmi_model->update(pet_bmi_model) != BMI_SUCCESS) {
            fprintf(stderr, "PET update failed at step %d\n", i);
            goto cleanup;
        }
        if (i == 0) continue; /* Preserve the driver's original output count. */
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
