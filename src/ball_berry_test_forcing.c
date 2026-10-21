// TEST ONLY: populate method 6 from configuration and the current AORC record.
#include <math.h>
#include <string.h>
#include "pet.h"

int bb6_prepare_test_forcing(pet_model *pet)
{
    NoahmpStabilityInput *input;
    double temperature_k, pressure_pa, humidity, wind, density, vapor_pa, sw, lw;
    int step, band;
    if (pet == NULL || pet->bmi.is_forcing_from_bmi != 0) return -1;
    pet->bb_stability_input_ready = 0;
    step = pet->bmi.current_step;
    if (step < 0 || step >= pet->bmi.num_timesteps) return -1;
    input = &pet->bb_stability_input;
    memset(input, 0, sizeof(*input));
    BallBerryTestParameters *cfg = &pet->bb_test;
    // Reject missing or inconsistent test geometry and canopy parameters.
    if (ball_berry_state_validate(&pet->bb_canopy) != 0 ||
        !isfinite(pet->bb_cwp) || pet->bb_cwp <= 0.0 ||
        !isfinite(pet->bb_stomata.mp) || pet->bb_stomata.mp <= 0.0 ||
        !isfinite(pet->bb_stomata.vcmx25) || pet->bb_stomata.vcmx25 <= 0.0 ||
        cfg->canopy_height_m <= 0.0 ||
        cfg->wind_reference_height_m <= cfg->displacement_height_m + cfg->momentum_roughness_m ||
        cfg->momentum_roughness_m <= 0.0 || cfg->heat_roughness_m <= 0.0 ||
        cfg->cosz <= 0.0 || cfg->cosz > 1.0 ||
        cfg->direct_shortwave_fraction < 0.0 || cfg->direct_shortwave_fraction > 1.0 ||
        cfg->visible_shortwave_fraction < 0.0 || cfg->visible_shortwave_fraction > 1.0 ||
        cfg->maximum_stability_iterations < 1.0 ||
        cfg->maximum_canopy_iterations < 1.0) return -1;
    input->maximum_stability_iterations = (int)cfg->maximum_stability_iterations;
    input->resistance_relative_tolerance = cfg->stability_relative_tolerance;
    input->canopy.maximum_iterations = (int)cfg->maximum_canopy_iterations;
    input->canopy.temperature_tolerance_k = cfg->canopy_temperature_tolerance_k;
    input->canopy.maximum_temperature_step_k = cfg->maximum_temperature_step_k;
    input->canopy.radiation_stomata.radiation.xl = cfg->xl;
    input->canopy.radiation_stomata.radiation.cosz = cfg->cosz;
    for (band = 0; band < 2; ++band) {
        input->canopy.radiation_stomata.radiation.rhol[band] = cfg->leaf_reflectance;
        input->canopy.radiation_stomata.radiation.rhos[band] = cfg->stem_reflectance;
        input->canopy.radiation_stomata.radiation.taul[band] = cfg->leaf_transmittance;
        input->canopy.radiation_stomata.radiation.taus[band] = cfg->stem_transmittance;
        if (band == 0) {
            input->canopy.radiation_stomata.radiation.albgrd[band] = cfg->ground_albedo_visible;
            input->canopy.radiation_stomata.radiation.albgri[band] = cfg->ground_albedo_visible;
        } else {
            input->canopy.radiation_stomata.radiation.albgrd[band] = cfg->ground_albedo_nir;
            input->canopy.radiation_stomata.radiation.albgri[band] = cfg->ground_albedo_nir;
        }
    }
    input->canopy.radiation_stomata.physiology.bp = cfg->bp;
    input->canopy.radiation_stomata.physiology.folnmx = cfg->folnmx;
    input->canopy.radiation_stomata.physiology.qe25 = cfg->qe25;
    input->canopy.radiation_stomata.physiology.kc25 = cfg->kc25;
    input->canopy.radiation_stomata.physiology.akc = cfg->akc;
    input->canopy.radiation_stomata.physiology.ko25 = cfg->ko25;
    input->canopy.radiation_stomata.physiology.ako = cfg->ako;
    input->canopy.radiation_stomata.physiology.avcmx = cfg->avcmx;
    input->canopy.radiation_stomata.physiology.c3psn = cfg->c3psn;
    input->canopy.radiation_stomata.canopy_air.mpe = 1.0e-6;
    input->canopy.radiation_stomata.canopy_air.foln = 1.0;
    input->canopy.radiation_stomata.canopy_air.igs = 1.0;
    input->canopy.foliage.rsurf_s_per_m = cfg->surface_resistance_s_per_m;
    input->canopy.foliage.cp_air_j_per_kg_k = 1004.0;
    input->canopy.foliage.psychrometric_constant_pa_per_k = 66.0;
    input->canopy.foliage.vegetation_emissivity = cfg->vegetation_emissivity;
    input->canopy.foliage.ground_emissivity = cfg->ground_emissivity;
    input->canopy.foliage.stefan_boltzmann_w_per_m2_k4 = 5.670374419e-8;
    input->canopy.foliage.latent_heat_j_per_kg = 2450000.0;
    input->ground_canopy.canopy_height_m = cfg->canopy_height_m;
    input->ground_canopy.reference_height_m = cfg->wind_reference_height_m;
    input->surface_layer.reference_height_m = cfg->wind_reference_height_m;
    input->ground_canopy.displacement_height_m = cfg->displacement_height_m;
    input->surface_layer.displacement_height_m = cfg->displacement_height_m;
    input->ground_canopy.momentum_roughness_m = cfg->momentum_roughness_m;
    input->surface_layer.momentum_roughness_m = cfg->momentum_roughness_m;
    input->surface_layer.heat_roughness_m = cfg->heat_roughness_m;
    input->ground_canopy.canopy_heat_roughness_m = cfg->heat_roughness_m;
    input->ground_canopy.ground_momentum_roughness_m = cfg->ground_momentum_roughness_m;
    input->ground_canopy.ground_heat_roughness_m = cfg->ground_heat_roughness_m;
    input->ground_canopy.leaf_dimension_m = cfg->leaf_dimension_m;
    input->ground_canopy.mpe = 1.0e-8;
    input->surface_layer.mpe = 1.0e-8;
    /* AORC data for this forcing record. */
    temperature_k = pet->forcing_data_air_temperature_2m_K[step];
    pressure_pa = pet->forcing_data_surface_pressure_Pa[step];
    humidity = pet->forcing_data_specific_humidity_2m_kg_per_kg[step];
    wind = hypot(pet->forcing_data_u_wind_speed_10m_m_per_s[step],
                 pet->forcing_data_v_wind_speed_10m_m_per_s[step]);
    sw = pet->forcing_data_incoming_shortwave_W_per_m2[step];
    lw = pet->forcing_data_incoming_longwave_W_per_m2[step];
    if (!isfinite(temperature_k) || !isfinite(pressure_pa) || !isfinite(humidity) ||
        !isfinite(wind) || !isfinite(sw) || !isfinite(lw) ||
        temperature_k < 200.0 || pressure_pa <= 0.0 || humidity < 0.0 || humidity >= 1.0 ||
        sw < 0.0 || lw < 0.0) return -1;
    if (wind < cfg->minimum_wind_speed_m_per_s) wind = cfg->minimum_wind_speed_m_per_s; // Test-only wind floor.
    density = pressure_pa / (287.05 * temperature_k * (1.0 + 0.61 * humidity));
    vapor_pa = humidity * pressure_pa / (0.622 + 0.378 * humidity);
    input->canopy.radiation_stomata.radiation.elai = pet->bb_canopy.lai;
    input->canopy.radiation_stomata.radiation.esai = pet->bb_canopy.sai;
    input->canopy.radiation_stomata.radiation.fveg = pet->bb_canopy.shdfac;
    input->canopy.radiation_stomata.physiology.mp = pet->bb_stomata.mp;
    input->canopy.radiation_stomata.physiology.vcmx25 = pet->bb_stomata.vcmx25;
    input->canopy.radiation_stomata.canopy_air.ea = vapor_pa;
    input->canopy.radiation_stomata.canopy_air.sfctmp = temperature_k;
    input->canopy.radiation_stomata.canopy_air.sfcpres = pressure_pa;
    input->canopy.radiation_stomata.canopy_air.co2 = cfg->co2_mole_fraction * pressure_pa;
    input->canopy.radiation_stomata.canopy_air.o2 = cfg->o2_mole_fraction * pressure_pa;
    input->canopy.foliage.fveg = pet->bb_canopy.shdfac;
    input->canopy.foliage.vai = pet->bb_canopy.lai + pet->bb_canopy.sai;
    input->canopy.foliage.air_temperature_k = temperature_k;
    input->canopy.foliage.ground_temperature_k = temperature_k + cfg->ground_temperature_offset_k;
    input->canopy.foliage.foliage_temperature_k = temperature_k + cfg->ground_temperature_offset_k;
    input->canopy.foliage.air_vapor_pressure_pa = vapor_pa;
    input->canopy.foliage.air_density_kg_per_m3 = density;
    input->canopy.foliage.incoming_longwave_w_per_m2 = lw;
    input->surface_layer.air_temperature_k = temperature_k;
    input->surface_layer.air_density_kg_per_m3 = density;
    input->surface_layer.specific_humidity_kg_per_kg = humidity;
    input->surface_layer.wind_speed_m_per_s = wind;
    input->ground_canopy.reference_wind_speed_m_per_s = wind;
    input->ground_canopy.air_density_kg_per_m3 = density;
    input->ground_canopy.vegetation_area_index = pet->bb_canopy.lai + pet->bb_canopy.sai;
    input->ground_canopy.cwp = pet->bb_cwp;
    input->canopy.foliage.timestep_s = pet->bmi.time_step_size_s;
    for (int band = 0; band < 2; ++band) {
        input->canopy.radiation_stomata.radiation.solad[band] = 0.35 * sw;
        input->canopy.radiation_stomata.radiation.solai[band] = cfg->visible_shortwave_fraction * (1.0 - cfg->direct_shortwave_fraction) * sw;
     }
    // TEST ONLY: fixed solar geometry and shortwave partition from config.
    input->canopy.radiation_stomata.radiation.cosz = cfg->cosz;
    pet->bb_stability_input_ready = 1;
    return 0;
}
