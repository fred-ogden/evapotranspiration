/*
 * tapenade/canopy_energy/ball_berry_head.c - the head Tapenade differentiates: one PET method 6 step (see
 * ball_berry_head.h).
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
 * The passive parameters and step inputs are copied into local structures field by field (the generated
 * copies of src/ball_berry_struct_copies.c: Tapenade 3.16 does not handle whole-structure assignment),
 * the independents are written over them, one step of the physics runs, and the outputs are packed into
 * the dependents.
 *
 * FLO 10/28: first version.  ASCII only.
 */
#include "ball_berry_head.h"
#include "ball_berry_struct_copies.h"

void ball_berry_head(double independent[BB_HEAD_INDEPENDENT_COUNT], double dependent[BB_HEAD_DEPENDENT_COUNT],
                     const BallBerryMethod6Parameters *passive_parameters, const BallBerryStepInputs *passive_step)
{
    BallBerryMethod6Parameters parameters;
    BallBerryStepInputs step;
    NoahmpStabilityInput canopy_input;
    NoahmpStabilityOutput solution;
    NoahmpVegeFluxExactTrace trace;
    BallBerryStepResults results;
    int index, status;

    bb_copy_BallBerryMethod6Parameters(passive_parameters, &parameters);
    bb_copy_BallBerryStepInputs(passive_step, &step);

    /* ---- the independents over the passive values */
    step.coupling.vegetation_temperature_start_of_step_K = independent[BB_HEAD_IN_VEGETATION_TEMPERATURE_START];
    step.coupling.canopy_air_temperature_start_of_step_K = independent[BB_HEAD_IN_CANOPY_AIR_TEMPERATURE_START];
    step.coupling.canopy_air_vapor_pressure_start_of_step_Pa = independent[BB_HEAD_IN_CANOPY_AIR_VAPOR_PRESSURE_START];
    step.air_temperature_K = independent[BB_HEAD_IN_AIR_TEMPERATURE];
    step.surface_pressure_Pa = independent[BB_HEAD_IN_SURFACE_PRESSURE];
    step.specific_humidity_kg_per_kg = independent[BB_HEAD_IN_SPECIFIC_HUMIDITY];
    step.shortwave_down_W_per_m2 = independent[BB_HEAD_IN_SHORTWAVE_DOWN];
    step.longwave_down_W_per_m2 = independent[BB_HEAD_IN_LONGWAVE_DOWN];
    step.wind_x_m_per_s = independent[BB_HEAD_IN_WIND_X];
    step.wind_y_m_per_s = independent[BB_HEAD_IN_WIND_Y];
    step.coupling.root_stress_factor_0_to_1 = independent[BB_HEAD_IN_ROOT_STRESS_FACTOR];
    step.coupling.ground_skin_temperature_K = independent[BB_HEAD_IN_GROUND_TEMPERATURE];
    step.coupling.surface_resistance_s_per_m = independent[BB_HEAD_IN_SURFACE_RESISTANCE];
    step.coupling.snow_cover_fraction = independent[BB_HEAD_IN_SNOW_COVER_FRACTION];
    step.coupling.snow_depth_m = independent[BB_HEAD_IN_SNOW_DEPTH];
    step.coupling.snow_albedo_direct_visible = independent[BB_HEAD_IN_SNOW_ALBEDO_DIRECT_VISIBLE];
    step.coupling.snow_albedo_direct_nir = independent[BB_HEAD_IN_SNOW_ALBEDO_DIRECT_NIR];
    step.coupling.snow_albedo_diffuse_visible = independent[BB_HEAD_IN_SNOW_ALBEDO_DIFFUSE_VISIBLE];
    step.coupling.snow_albedo_diffuse_nir = independent[BB_HEAD_IN_SNOW_ALBEDO_DIFFUSE_NIR];
    step.coupling.canopy_wetted_fraction = independent[BB_HEAD_IN_CANOPY_WETTED_FRACTION];
    step.coupling.canopy_liquid_water_mm = independent[BB_HEAD_IN_CANOPY_LIQUID_WATER];
    step.coupling.canopy_ice_mm = independent[BB_HEAD_IN_CANOPY_ICE];
    step.coupling.precipitation_advected_heat_to_canopy_W_per_m2 = independent[BB_HEAD_IN_PRECIPITATION_ADVECTED_HEAT];
    step.coupling.soil_top_disc_moisture_content = independent[BB_HEAD_IN_SOIL_TOP_DISC_MOISTURE];
    step.stomata.mp = independent[BB_HEAD_IN_STOMATAL_SLOPE];
    step.stomata.vcmx25 = independent[BB_HEAD_IN_MAXIMUM_CARBOXYLATION_RATE];
    step.stomata.bp = independent[BB_HEAD_IN_MINIMUM_LEAF_CONDUCTANCE];
    step.stomata.qe25 = independent[BB_HEAD_IN_QUANTUM_EFFICIENCY];
    step.canopy_wind_parameter = independent[BB_HEAD_IN_CANOPY_WIND_PARAMETER];
    step.green_vegetation_fraction = independent[BB_HEAD_IN_GREEN_VEGETATION_FRACTION];
    parameters.momentum_roughness_length_m = independent[BB_HEAD_IN_MOMENTUM_ROUGHNESS_LENGTH];
    parameters.characteristic_leaf_dimension_m = independent[BB_HEAD_IN_LEAF_DIMENSION];
    parameters.surface_resistance_snow_s_per_m = independent[BB_HEAD_IN_SNOW_SURFACE_RESISTANCE];

    /* ---- one step */
    status = ball_berry_method6_physics_step(&parameters, &step, &canopy_input, &solution, &trace, &results);

    /* ---- the dependents (0 if the step failed; the tests use states where it does not) */
    for (index = 0; index < BB_HEAD_DEPENDENT_COUNT; ++index) dependent[index] = 0.0;
    if (status == 0) {
        dependent[BB_HEAD_OUT_TRANSPIRATION_FLUX] = results.coupling_outputs.canopy_transpiration_flux_m_per_s;
        dependent[BB_HEAD_OUT_CANOPY_EVAPORATION_HEAT_FLUX] = results.coupling_outputs.canopy_evaporation_heat_flux_W_per_m2;
        dependent[BB_HEAD_OUT_GROUND_SHORTWAVE_ABSORBED] =
            results.coupling_outputs.ground_shortwave_absorbed_below_canopy_W_per_m2;
        dependent[BB_HEAD_OUT_GROUND_LONGWAVE_ABSORBED] =
            results.coupling_outputs.ground_longwave_absorbed_below_canopy_W_per_m2;
        dependent[BB_HEAD_OUT_CANOPY_AIR_TEMPERATURE] = results.coupling_outputs.canopy_air_temperature_K;
        dependent[BB_HEAD_OUT_CANOPY_AIR_VAPOR_PRESSURE] = results.coupling_outputs.canopy_air_vapor_pressure_Pa;
        dependent[BB_HEAD_OUT_BELOW_CANOPY_RESISTANCE] = results.coupling_outputs.below_canopy_aerodynamic_resistance_s_per_m;
        dependent[BB_HEAD_OUT_VEGETATION_TEMPERATURE] = results.coupling_outputs.vegetation_temperature_K;
        dependent[BB_HEAD_OUT_CANOPY_LATENT_HEAT] = results.coupling_outputs.canopy_latent_heat_J_per_kg;
        dependent[BB_HEAD_OUT_CANOPY_FUSION_HEAT_FLUX] = results.coupling_outputs.canopy_fusion_heat_flux_W_per_m2;
        dependent[BB_HEAD_OUT_CANOPY_HEAT_STORAGE_FLUX] = results.coupling_outputs.canopy_heat_storage_flux_W_per_m2;
        dependent[BB_HEAD_OUT_SURFACE_RESISTANCE_WITH_SNOW] =
            results.coupling_outputs.surface_resistance_with_snow_s_per_m;
        dependent[BB_HEAD_OUT_VEGETATED_FRACTION] = results.coupling_outputs.vegetated_fraction;
        dependent[BB_HEAD_OUT_PET] = results.pet_m_per_s;
    }
}
