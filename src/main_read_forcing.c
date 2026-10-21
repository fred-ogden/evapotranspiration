 #include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../include/bmi.h"
#include "../include/pet.h"
#include "../include/bmi_pet.h"

// Standalone test-only forcing adapter; not part of the BMI interface.
extern int bb6_prepare_test_forcing(pet_model *pet);

/************************************************************************
    This main program is a mock framwork.
    This is not part of BMI, but acts as the driver that calls the model.
************************************************************************/
int
 main(int argc, const char *argv[])
{

  /************************************************************************
      A configuration file is required for running this model through BMI
  ************************************************************************/
  if(argc<=1){
//    printf("make sure to include a path to the PET config file\n");
    exit(1);
  }

  /************************************************************************
      allocating memory to store the entire BMI structure for PET
  ************************************************************************/
//  printf("allocating memory to store entire BMI PET structure\n");
  Bmi *pet_bmi_model = (Bmi *) malloc(sizeof(Bmi));

  /************************************************************************
      Registering the BMI model for PET
  ************************************************************************/
//  printf("Registering BMI PET model\n");
  register_bmi_pet(pet_bmi_model);

  /************************************************************************
      Initializing the BMI model for PET
  ************************************************************************/
//  printf("Initializeing BMI PET model\n");
  const char *cfg_file_pet = argv[1];
  if (pet_bmi_model->initialize(pet_bmi_model, cfg_file_pet) != BMI_SUCCESS) {
    fprintf(stderr, "PET initialization failed: %s\n", cfg_file_pet);
    return 1;
  }

  /************************************************************************
      Defining the PET data as a new variable, so we can use its contents
  ************************************************************************/
//  printf("Get the information from the configuration here in Main\n");
  pet_model *pet;
  pet = (pet_model *) pet_bmi_model->data;

  /************************************************************************
      Defining the PET data as a new variable, so we can use its contents
  ************************************************************************/
//  printf("Updating PET model\n");
  if (pet->pet_method == 6) {
    if (bb6_prepare_test_forcing(pet) != 0 ||
        pet_bmi_model->update(pet_bmi_model) != BMI_SUCCESS) {
      fprintf(stderr, "Method 6 PET update failed\n");
      return 1;
    }
    printf("potential transpiration (m s-1): %.12e\n", pet->pet_m_per_s);
  } else if (pet_bmi_model->update(pet_bmi_model) != BMI_SUCCESS) {
    fprintf(stderr, "PET update failed\n");
    return 1;
  }

//  printf("looping through and calling update\n");
  int i=0;
  for (i = 0; i < 1; i++){
    // The unit test only runs one time step.
    if (pet->bmi.run_unit_tests == 0 && pet->pet_method != 6){
      if (pet_bmi_model->update(pet_bmi_model) != BMI_SUCCESS) {
        fprintf(stderr, "PET update failed\n");
        return 1;
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
  }

//  printf("Finalizing models\n");
  pet_bmi_model->finalize(pet_bmi_model);

  return 0;
}  /************    END MAIN PROGRAM    ***********************************************/
