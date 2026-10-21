! ==========================================================================================
! noahmp_canopy_reference_driver.f90
!
! Purpose:
!   Stage 2 reference: one time step of the NWM Noah-MP vegetated canopy solution
!   with Ball-Berry stomata (OPT_CRS = 1), driven exactly as ENERGY drives it:
!     ATM-style forcing preparation -> RADIATION (ALBEDO, TWOSTREAM, SURRAD) ->
!     emissivities, psychrometric constants -> VEGE_FLUX (SFCDIF1, RAGRB, STOMATA).
!   Differences from a full NOAHMP_SFLX call (all deliberate, all listed in the output):
!     - ground albedos imposed (HARNESS_FIX_GROUND_ALBEDO), not from soil color/moisture;
!     - ground skin temperature TG prescribed (HARNESS_FIX_GROUND_TEMPERATURE), so the
!       VEGE_FLUX ground loop does not update TG and GH is the implied residual;
!     - BTRAN, RSURF, RHSUR imposed (no soil state);
!     - cold-start canopy state TV, TAH, EAH (prognostic in Noah-MP).
!
! Usage:
!   noahmp_canopy_reference_driver_r4 config_file output_directory
!
! Outputs (precision tag r4 or r8 from the build):
!   noahmp_canopy_reference_summary_<tag>.csv   name, value, units, description
!   noahmp_canopy_iteration_trace_<tag>.csv     one row per canopy (loop1) iteration
! ==========================================================================================
program noahmp_canopy_reference_driver

  use module_sf_noahmplsm
  use harness_text_utilities
  use noahmp_reference_setup
  implicit none

  ! ---- Noah-MP dimensions as used by NWM
  integer, parameter :: number_of_snow_discs_maximum = 3   ! [NSNOW]
  integer, parameter :: number_of_soil_discs = 4           ! [NSOIL]
  integer, parameter :: actual_number_of_snow_discs = 0    ! [ISNOW] no snow
  integer, parameter :: surface_type_soil = 1              ! [IST]
  integer, parameter :: land_ice_flag_off = 0              ! [ICE]
  integer, parameter :: grid_index_i = 1, grid_index_j = 1 ! [ILOC, JLOC]

  real,    parameter :: bare_soil_roughness_length_m = 0.002 ! [Z0] hard-coded in ENERGY

  character(len=512) :: configuration_file_path, output_directory_path, forcing_file_path
  character(len=512) :: summary_file_path, trace_file_path
  character(len=128) :: forcing_time_string
  character(len=2)   :: precision_tag
  integer :: summary_file_unit, trace_file_unit
  integer :: usgs_vegetation_type, radiation_transfer_option, year_length_days
  integer :: initialize_canopy_state_from_forcing

  type(noahmp_parameters) :: noahmp_parameter_set

  ! ---- forcing and derived atmospheric state
  real :: time_step_size_s                         ! [DT]
  real :: downward_longwave_w_per_m2               ! [LWDN]
  real :: downward_shortwave_w_per_m2              ! [SOLDN]
  real :: surface_pressure_pa                      ! [SFCPRS], also [PSFC]
  real :: specific_humidity_kg_per_kg              ! [QAIR]
  real :: air_temperature_k                        ! [SFCTMP]
  real :: eastward_wind_m_per_s                    ! [UU]
  real :: northward_wind_m_per_s                   ! [VV]
  real :: wind_speed_reference_height_m_per_s      ! [UR] >= 1 m/s
  real :: potential_temperature_k                  ! [THAIR]
  real :: air_vapor_pressure_pa                    ! [EAIR]
  real :: air_density_kg_per_m3                    ! [RHOAIR]
  real :: direct_shortwave_fraction, visible_shortwave_fraction
  real :: incoming_direct_shortwave_w_per_m2(2)    ! [SOLAD] 1 = vis, 2 = nir
  real :: incoming_diffuse_shortwave_w_per_m2(2)   ! [SOLAI]
  real :: cosine_solar_zenith_angle                ! [COSZ]
  real :: co2_partial_pressure_pa                  ! [CO2AIR]
  real :: o2_partial_pressure_pa                   ! [O2AIR]
  real :: julian_day_since_jan1, latitude_degrees

  ! ---- vegetation and surface state
  real :: leaf_area_index, stem_area_index         ! [LAI, SAI]
  real :: exposed_leaf_area_index                  ! [ELAI] after burial by snow
  real :: exposed_stem_area_index                  ! [ESAI] after burial by snow
  real :: snow_free_soil_albedo_direct(2), snow_free_soil_albedo_diffuse(2)   ! imposed [ALBSOD, ALBSOI]
  real :: imposed_snow_albedo_direct(2), imposed_snow_albedo_diffuse(2)       ! imposed [ALBSND, ALBSNI]
  real :: ground_longwave_absorbed_below_canopy_w_per_m2   ! EMG*((1-EMV)*LWDN + EMV*SB*TV**4)
  real :: below_canopy_aerodynamic_resistance_s_per_m      ! [RAHG] = 1/CHUC
  real :: prescribed_ground_temperature_k                  ! [TG] as passed in; VEGE_FLUX may reset TG under snow
  real :: vegetation_area_index                    ! [VAI]
  real :: green_vegetation_fraction                ! [FVEG]
  real :: growing_season_index                     ! [IGS]
  real :: zero_plane_displacement_m                ! [ZPD]
  real :: momentum_roughness_length_m              ! [Z0M]
  real :: ground_momentum_roughness_length_m       ! [Z0MG]
  real :: forcing_height_above_canopy_top_m        ! [ZREF]
  real :: reference_height_m                       ! [ZLVL]
  real :: canopy_wind_parameter                    ! [CWP]
  real :: vegetation_emissivity                    ! [EMV]
  real :: ground_emissivity                        ! [EMG]
  real :: snow_cover_fraction                      ! [FSNO]
  real :: snow_depth_m                             ! [SNOWH]
  real :: wetted_canopy_fraction                   ! [FWET]
  real :: canopy_liquid_water_mm, canopy_ice_mm    ! [CANLIQ, CANICE]
  real :: root_stress_factor_0_to_1                ! [BTRAN]
  real :: foliage_nitrogen_percent                 ! [FOLN]
  real :: surface_resistance_s_per_m               ! [RSURF]
  real :: ground_surface_relative_humidity         ! [RHSUR]
  real :: latent_heat_vegetation_j_per_kg          ! [LATHEAV]
  real :: latent_heat_ground_j_per_kg              ! [LATHEAG]
  real :: psychrometric_constant_vegetation_pa_per_k ! [GAMMAV]
  real :: psychrometric_constant_ground_pa_per_k   ! [GAMMAG]
  real :: soil_moisture_volumetric(number_of_soil_discs) ! [SMC] only feeds GROUNDALB, which is overridden
  real :: snow_soil_temperature_k(-number_of_snow_discs_maximum+1:number_of_soil_discs)     ! [STC] unused, TG fixed
  real :: snow_soil_thermal_conductivity_w_per_m_per_k(-number_of_snow_discs_maximum+1:number_of_soil_discs) ! [DF] unused
  real :: snow_soil_disc_thickness_m(-number_of_snow_discs_maximum+1:number_of_soil_discs)  ! [DZSNSO] unused

  ! ---- prognostic canopy state (inout to VEGE_FLUX)
  real :: vegetation_temperature_k                 ! [TV]
  real :: canopy_air_temperature_k                 ! [TAH]
  real :: canopy_air_vapor_pressure_pa             ! [EAH]
  real :: ground_temperature_k                     ! [TG]
  real :: momentum_drag_coefficient                ! [CM]
  real :: heat_exchange_coefficient                ! [CH]
  real :: initial_vegetation_temperature_k, initial_canopy_air_temperature_k
  real :: initial_canopy_air_vapor_pressure_pa, initial_leaf_saturation_vapor_pressure_pa

  ! ---- radiation outputs
  real :: sunlit_canopy_fraction                   ! [FSUN]
  real :: sunlit_leaf_area_index                   ! [LAISUN]
  real :: shaded_leaf_area_index                   ! [LAISHA]
  real :: absorbed_par_per_sunlit_lai_w_per_m2     ! [PARSUN]
  real :: absorbed_par_per_shaded_lai_w_per_m2     ! [PARSHA]
  real :: shortwave_absorbed_by_vegetation_w_per_m2 ! [SAV]
  real :: shortwave_absorbed_by_ground_w_per_m2    ! [SAG]
  real :: total_absorbed_shortwave_w_per_m2        ! [FSA]
  real :: total_reflected_shortwave_w_per_m2       ! [FSR]
  real :: vegetation_reflected_shortwave_w_per_m2  ! [FSRV]
  real :: ground_reflected_shortwave_w_per_m2      ! [FSRG]
  real :: between_canopy_gap_fraction, within_canopy_gap_fraction ! [BGAP, WGAP]
  real :: snow_albedo_direct(2), snow_albedo_diffuse(2)           ! [ALBSND, ALBSNI]
  real :: snow_albedo_last_step, snow_age_nondimensional          ! [ALBOLD, TAUSS]

  ! ---- VEGE_FLUX outputs
  real :: wind_stress_eastward_n_per_m2, wind_stress_northward_n_per_m2 ! [TAUXV, TAUYV]
  real :: ground_net_longwave_w_per_m2             ! [IRG]  + to atmosphere
  real :: canopy_net_longwave_w_per_m2             ! [IRC]
  real :: ground_sensible_heat_w_per_m2            ! [SHG]
  real :: canopy_sensible_heat_w_per_m2            ! [SHC]
  real :: ground_evaporation_heat_w_per_m2         ! [EVG]
  real :: canopy_evaporation_heat_w_per_m2         ! [EVC]
  real :: transpiration_heat_w_per_m2              ! [TR]
  real :: ground_heat_flux_w_per_m2                ! [GH]  + into soil
  real :: two_meter_air_temperature_k              ! [T2MV]
  real :: sunlit_photosynthesis_umol_per_m2_per_s  ! [PSNSUN]
  real :: shaded_photosynthesis_umol_per_m2_per_s  ! [PSNSHA]
  real :: sunlit_stomatal_resistance_s_per_m       ! [RSSUN]
  real :: shaded_stomatal_resistance_s_per_m       ! [RSSHA]
  real :: two_meter_specific_humidity_kg_per_kg    ! [Q2V]
  real :: two_meter_heat_conductance_m_per_s       ! [CAH2]
  real :: leaf_exchange_coefficient                ! [CHLEAF]
  real :: under_canopy_exchange_coefficient        ! [CHUC]
  real :: surface_mixing_ratio_kg_per_kg           ! [QSFC] inout
  real :: canopy_energy_balance_residual_w_per_m2
  real :: transpiration_rate_mm_per_h

  ! ---- dummies required by the VEGE_FLUX interface
  real :: cloud_water_mixing_ratio, grid_spacing_m, lowest_layer_thickness_m
  real :: precipitation_heat_canopy_w_per_m2, precipitation_heat_ground_w_per_m2
  real :: water_vapor_mixing_ratio_kg_per_kg
  logical :: vegetated_surface_flag

  real :: temperature_celsius, saturation_vapor_pressure_water_pa, saturation_vapor_pressure_ice_pa
  real :: d_saturation_vapor_pressure_water_pa_per_k, d_saturation_vapor_pressure_ice_pa_per_k
  real :: ground_temperature_offset_from_air_k

  ! ------------------------------------------------------------------------------------------
  ! command line and configuration
  ! ------------------------------------------------------------------------------------------
  if (command_argument_count() /= 2) then
    write(*,'(a)') 'usage: noahmp_canopy_reference_driver_<r4|r8> config_file output_directory'
    stop 1
  end if
  call get_command_argument(1, configuration_file_path)
  call get_command_argument(2, output_directory_path)

  if (kind(1.0) == 8) then
    precision_tag = 'r8'
  else
    precision_tag = 'r4'
  end if

  call read_key_value_configuration_file(configuration_file_path)
  forcing_file_path = 'unset'
  call read_forcing_file_path_from_configuration(configuration_file_path, forcing_file_path)
  call read_forcing_csv_first_data_row(forcing_file_path)
  forcing_time_string = get_forcing_time_string()

  usgs_vegetation_type      = get_required_integer_value('usgs_vegetation_type')
  radiation_transfer_option = get_required_integer_value('radiation_transfer_option')
  time_step_size_s          = get_required_real_value('time_step_size_s')
  latitude_degrees              = get_required_real_value('latitude_degrees')

  call fill_noahmp_parameters_from_configuration(noahmp_parameter_set)
  call set_nwm_noahmp_physics_options(radiation_transfer_option)

  ! ------------------------------------------------------------------------------------------
  ! forcing, prepared as in ATM (THAIR = SFCTMP because PAIR = SFCPRS there)
  ! ------------------------------------------------------------------------------------------
  downward_longwave_w_per_m2  = get_forcing_column_value('DLWRF')
  downward_shortwave_w_per_m2 = get_forcing_column_value('DSWRF')
  surface_pressure_pa         = get_forcing_column_value('PRES')
  specific_humidity_kg_per_kg = get_forcing_column_value('SPFH')
  air_temperature_k           = get_forcing_column_value('TMP')
  eastward_wind_m_per_s       = get_forcing_column_value('UGRD')
  northward_wind_m_per_s      = get_forcing_column_value('VGRD')

  potential_temperature_k = air_temperature_k
  air_vapor_pressure_pa   = specific_humidity_kg_per_kg*surface_pressure_pa / &
                            (0.622 + 0.378*specific_humidity_kg_per_kg)
  air_density_kg_per_m3   = (surface_pressure_pa - 0.378*air_vapor_pressure_pa) / &
                            (287.04*air_temperature_k)                ! RAIR = 287.04
  water_vapor_mixing_ratio_kg_per_kg = specific_humidity_kg_per_kg   ! [Q2] passed to VEGE_FLUX, not used by it

  cosine_solar_zenith_angle  = get_required_real_value('cosine_solar_zenith_angle')
  direct_shortwave_fraction  = get_required_real_value('direct_shortwave_fraction')
  visible_shortwave_fraction = get_required_real_value('visible_shortwave_fraction')
  incoming_direct_shortwave_w_per_m2(1)  = downward_shortwave_w_per_m2*direct_shortwave_fraction* &
                                           visible_shortwave_fraction
  incoming_direct_shortwave_w_per_m2(2)  = downward_shortwave_w_per_m2*direct_shortwave_fraction* &
                                           (1.0 - visible_shortwave_fraction)
  incoming_diffuse_shortwave_w_per_m2(1) = downward_shortwave_w_per_m2*(1.0 - direct_shortwave_fraction)* &
                                           visible_shortwave_fraction
  incoming_diffuse_shortwave_w_per_m2(2) = downward_shortwave_w_per_m2*(1.0 - direct_shortwave_fraction)* &
                                           (1.0 - visible_shortwave_fraction)

  ! wind speed at reference height, UR >= 1 m/s (ENERGY)
  wind_speed_reference_height_m_per_s = max(sqrt(eastward_wind_m_per_s**2 + northward_wind_m_per_s**2), 1.0)

  ! CO2 and O2 partial pressures, as the NWM driver forms them
  co2_partial_pressure_pa = noahmp_parameter_set%CO2 * surface_pressure_pa
  o2_partial_pressure_pa  = noahmp_parameter_set%O2  * surface_pressure_pa

  ! ------------------------------------------------------------------------------------------
  ! ground temperature, initial canopy state
  ! ------------------------------------------------------------------------------------------
  ground_temperature_offset_from_air_k = get_required_real_value('ground_temperature_offset_from_air_k')
  ground_temperature_k = air_temperature_k + ground_temperature_offset_from_air_k

  initialize_canopy_state_from_forcing = get_required_integer_value('initialize_canopy_state_from_forcing')
  if (initialize_canopy_state_from_forcing == 1) then
    vegetation_temperature_k     = air_temperature_k
    canopy_air_temperature_k     = air_temperature_k
    canopy_air_vapor_pressure_pa = air_vapor_pressure_pa
  else
    vegetation_temperature_k     = get_required_real_value('initial_vegetation_temperature_k')
    canopy_air_temperature_k     = get_required_real_value('initial_canopy_air_temperature_k')
    canopy_air_vapor_pressure_pa = get_required_real_value('initial_canopy_air_vapor_pressure_pa')
  end if
  initial_vegetation_temperature_k     = vegetation_temperature_k
  initial_canopy_air_temperature_k     = canopy_air_temperature_k
  initial_canopy_air_vapor_pressure_pa = canopy_air_vapor_pressure_pa

  ! leaf internal vapor pressure STOMATA sees on iteration 1 (VEGE_FLUX: TDC then ESAT)
  temperature_celsius = min(50.0, max(-50.0, vegetation_temperature_k - 273.16))
  call ESAT(temperature_celsius, saturation_vapor_pressure_water_pa, saturation_vapor_pressure_ice_pa, &
            d_saturation_vapor_pressure_water_pa_per_k, d_saturation_vapor_pressure_ice_pa_per_k)
  if (temperature_celsius > 0.0) then
    initial_leaf_saturation_vapor_pressure_pa = saturation_vapor_pressure_water_pa
  else
    initial_leaf_saturation_vapor_pressure_pa = saturation_vapor_pressure_ice_pa
  end if

  momentum_drag_coefficient = 0.0
  heat_exchange_coefficient = 0.0

  ! ------------------------------------------------------------------------------------------
  ! phenology: LAI, SAI from the monthly table, growing season index from TV > TMIN
  ! ------------------------------------------------------------------------------------------
  call convert_time_string_to_noahmp_julian_day(forcing_time_string, julian_day_since_jan1, &
                                                year_length_days)
  if (get_required_integer_value('use_monthly_lai_sai_table') == 1) then
    call interpolate_monthly_lai_sai_as_noahmp_phenology(noahmp_parameter_set, latitude_degrees, &
           julian_day_since_jan1, year_length_days, leaf_area_index, stem_area_index)
  else
    leaf_area_index = get_required_real_value('leaf_area_index')
    stem_area_index = get_required_real_value('stem_area_index')
    ! same small-value checks PHENOLOGY applies
    if (stem_area_index < 0.05) stem_area_index = 0.0
    if (leaf_area_index < 0.05 .or. stem_area_index == 0.0) leaf_area_index = 0.0
  end if
  ! FLO 10/21: burial of the canopy by snow (PHENOLOGY); ELAI and ESAI are used from here on
  snow_depth_m = get_required_real_value('snow_depth_m')
  call apply_snow_burial_as_noahmp_phenology(noahmp_parameter_set, snow_depth_m, &
         leaf_area_index, stem_area_index, exposed_leaf_area_index, exposed_stem_area_index)
  if (vegetation_temperature_k > noahmp_parameter_set%TMIN) then
    growing_season_index = 1.0
  else
    growing_season_index = 0.0
  end if

  ! ------------------------------------------------------------------------------------------
  ! ENERGY preliminaries (soil surface type); FLO 10/21: snow and canopy water from the configuration
  ! ------------------------------------------------------------------------------------------
  green_vegetation_fraction = get_required_real_value('green_vegetation_fraction')
  vegetation_area_index     = exposed_leaf_area_index + exposed_stem_area_index   ! VAI = ELAI + ESAI
  vegetated_surface_flag    = (vegetation_area_index > 0.0)
  snow_cover_fraction       = get_required_real_value('snow_cover_fraction')      ! [FSNO]
  wetted_canopy_fraction    = get_required_real_value('canopy_wetted_fraction')   ! [FWET]
  canopy_liquid_water_mm    = get_required_real_value('canopy_liquid_water_mm')   ! [CANLIQ]
  canopy_ice_mm             = get_required_real_value('canopy_ice_mm')            ! [CANICE]

  ! ground roughness: IST = 1, Z0MG = Z0*(1-FSNO) + FSNO*Z0SNO
  ground_momentum_roughness_length_m = bare_soil_roughness_length_m*(1.0 - snow_cover_fraction) + &
                                       snow_cover_fraction*noahmp_parameter_set%Z0SNO
  momentum_roughness_length_m = noahmp_parameter_set%Z0MVT
  zero_plane_displacement_m   = 0.65*noahmp_parameter_set%HVT
  if (snow_depth_m > zero_plane_displacement_m) zero_plane_displacement_m = snow_depth_m   ! ENERGY: IF(SNOWH.GT.ZPD) ZPD = SNOWH
  forcing_height_above_canopy_top_m = get_required_real_value('forcing_height_above_canopy_top_m')
  reference_height_m = max(zero_plane_displacement_m, noahmp_parameter_set%HVT) + &
                       forcing_height_above_canopy_top_m
  ! ENERGY: ZPDG = SNOWH; IF(ZPDG >= ZLVL) ZLVL = ZPDG + ZREF
  if (snow_depth_m >= reference_height_m) reference_height_m = snow_depth_m + forcing_height_above_canopy_top_m
  canopy_wind_parameter = noahmp_parameter_set%CWPVT

  ! RADIATION inputs not otherwise needed
  soil_moisture_volumetric = 0.3
  snow_albedo_last_step    = 0.65
  snow_age_nondimensional  = 0.0

  ! imposed snow-free soil albedos and snow albedos; GROUNDALB blends them by FSNO (harness patch)
  snow_free_soil_albedo_direct(1)  = get_required_real_value('soil_albedo_direct_visible')
  snow_free_soil_albedo_direct(2)  = get_required_real_value('soil_albedo_direct_nir')
  snow_free_soil_albedo_diffuse(1) = get_required_real_value('soil_albedo_diffuse_visible')
  snow_free_soil_albedo_diffuse(2) = get_required_real_value('soil_albedo_diffuse_nir')
  imposed_snow_albedo_direct(1)    = get_required_real_value('snow_albedo_direct_visible')
  imposed_snow_albedo_direct(2)    = get_required_real_value('snow_albedo_direct_nir')
  imposed_snow_albedo_diffuse(1)   = get_required_real_value('snow_albedo_diffuse_visible')
  imposed_snow_albedo_diffuse(2)   = get_required_real_value('snow_albedo_diffuse_nir')
  HARNESS_FIX_GROUND_ALBEDO = .true.
  HARNESS_SOIL_ALBEDO_DIRECT(1:2)  = snow_free_soil_albedo_direct(1:2)
  HARNESS_SOIL_ALBEDO_DIFFUSE(1:2) = snow_free_soil_albedo_diffuse(1:2)
  HARNESS_SNOW_ALBEDO_DIRECT(1:2)  = imposed_snow_albedo_direct(1:2)
  HARNESS_SNOW_ALBEDO_DIFFUSE(1:2) = imposed_snow_albedo_diffuse(1:2)

  call RADIATION(noahmp_parameter_set, usgs_vegetation_type, surface_type_soil, land_ice_flag_off, &
                 number_of_soil_discs, &
                 0.0, 0.0, time_step_size_s, cosine_solar_zenith_angle, snow_depth_m, &
                 ground_temperature_k, vegetation_temperature_k, snow_cover_fraction, 0.0, &
                 wetted_canopy_fraction, &
                 exposed_leaf_area_index, exposed_stem_area_index, soil_moisture_volumetric, &
                 incoming_direct_shortwave_w_per_m2, incoming_diffuse_shortwave_w_per_m2, &
                 green_vegetation_fraction, grid_index_i, grid_index_j, &
                 snow_albedo_last_step, snow_age_nondimensional, &
                 sunlit_canopy_fraction, sunlit_leaf_area_index, shaded_leaf_area_index, &
                 absorbed_par_per_sunlit_lai_w_per_m2, absorbed_par_per_shaded_lai_w_per_m2, &
                 shortwave_absorbed_by_vegetation_w_per_m2, shortwave_absorbed_by_ground_w_per_m2, &
                 total_reflected_shortwave_w_per_m2, total_absorbed_shortwave_w_per_m2, &
                 vegetation_reflected_shortwave_w_per_m2, ground_reflected_shortwave_w_per_m2, &
                 snow_albedo_direct, snow_albedo_diffuse, &
                 between_canopy_gap_fraction, within_canopy_gap_fraction)

  ! emissivities (ENERGY)
  vegetation_emissivity = 1.0 - exp(-(exposed_leaf_area_index + exposed_stem_area_index)/1.0)
  ground_emissivity     = noahmp_parameter_set%EG(surface_type_soil)*(1.0 - snow_cover_fraction) + &
                          1.0*snow_cover_fraction

  root_stress_factor_0_to_1         = get_required_real_value('root_stress_factor_0_to_1')
  foliage_nitrogen_percent          = get_required_real_value('foliage_nitrogen_percent')
  surface_resistance_s_per_m        = get_required_real_value('surface_resistance_s_per_m')
  ground_surface_relative_humidity  = get_required_real_value('ground_surface_relative_humidity')

  ! psychrometric constants (ENERGY): HVAP = 2.5104E06, HSUB = 2.8440E06, CPAIR = 1004.64
  if (vegetation_temperature_k > 273.16) then
    latent_heat_vegetation_j_per_kg = 2.5104e06
  else
    latent_heat_vegetation_j_per_kg = 2.8440e06
  end if
  psychrometric_constant_vegetation_pa_per_k = 1004.64*surface_pressure_pa/(0.622*latent_heat_vegetation_j_per_kg)
  if (ground_temperature_k > 273.16) then
    latent_heat_ground_j_per_kg = 2.5104e06
  else
    latent_heat_ground_j_per_kg = 2.8440e06
  end if
  psychrometric_constant_ground_pa_per_k = 1004.64*surface_pressure_pa/(0.622*latent_heat_ground_j_per_kg)

  ! VEGE_FLUX interface dummies (soil terms unused with TG prescribed)
  snow_soil_temperature_k = ground_temperature_k
  snow_soil_thermal_conductivity_w_per_m_per_k = 1.0
  snow_soil_disc_thickness_m = 0.0
  snow_soil_disc_thickness_m(1:number_of_soil_discs) = (/ 0.1, 0.3, 0.6, 1.0 /)
  cloud_water_mixing_ratio = 0.0
  grid_spacing_m           = 1000.0
  lowest_layer_thickness_m = 20.0
  precipitation_heat_canopy_w_per_m2 = 0.0
  precipitation_heat_ground_w_per_m2 = 0.0
  surface_mixing_ratio_kg_per_kg = specific_humidity_kg_per_kg

  ! ------------------------------------------------------------------------------------------
  ! per-iteration trace
  ! ------------------------------------------------------------------------------------------
  trace_file_path = trim(output_directory_path)//'/noahmp_canopy_iteration_trace_'//precision_tag//'.csv'
  open(newunit=trace_file_unit, file=trim(trace_file_path), status='replace', action='write')
  write(trace_file_unit,'(a)') 'iteration,vegetation_temperature_k,vegetation_temperature_change_k,'// &
       'canopy_air_temperature_k,canopy_air_vapor_pressure_pa,sunlit_stomatal_resistance_s_per_m,'// &
       'shaded_stomatal_resistance_s_per_m,sunlit_photosynthesis_umol_per_m2_per_s,'// &
       'shaded_photosynthesis_umol_per_m2_per_s,leaf_boundary_layer_resistance_s_per_m,'// &
       'above_canopy_heat_resistance_s_per_m,above_canopy_vapor_resistance_s_per_m,'// &
       'under_canopy_heat_resistance_s_per_m,under_canopy_vapor_resistance_s_per_m,'// &
       'momentum_drag_coefficient,heat_exchange_coefficient,friction_velocity_m_per_s,'// &
       'monin_obukhov_stability_parameter,momentum_stability_correction,heat_stability_correction,'// &
       'under_canopy_heat_stability_correction,canopy_net_longwave_w_per_m2,'// &
       'canopy_sensible_heat_w_per_m2,canopy_evaporation_heat_w_per_m2,transpiration_heat_w_per_m2,'// &
       'leaf_to_canopy_air_heat_conductance_m_per_s,transpiration_conductance_m_per_s,'// &
       'canopy_air_to_reference_heat_conductance_m_per_s,ground_to_canopy_air_heat_conductance_m_per_s,'// &
       'leaf_saturation_vapor_pressure_pa,leaf_saturation_vapor_pressure_slope_pa_per_k'
  prescribed_ground_temperature_k = ground_temperature_k
  HARNESS_TRACE_UNIT    = trace_file_unit
  HARNESS_TRACE_ENABLED = .true.
  HARNESS_FIX_GROUND_TEMPERATURE = .true.

  ! ------------------------------------------------------------------------------------------
  ! vegetated canopy energy balance (ENERGY passes TGV = TG, CMV = CM, CHV = CH)
  ! ------------------------------------------------------------------------------------------
  call VEGE_FLUX(noahmp_parameter_set, number_of_snow_discs_maximum, number_of_soil_discs, &
                 actual_number_of_snow_discs, usgs_vegetation_type, vegetated_surface_flag, &
                 time_step_size_s, shortwave_absorbed_by_vegetation_w_per_m2, &
                 shortwave_absorbed_by_ground_w_per_m2, downward_longwave_w_per_m2, &
                 wind_speed_reference_height_m_per_s, &
                 eastward_wind_m_per_s, northward_wind_m_per_s, air_temperature_k, &
                 potential_temperature_k, specific_humidity_kg_per_kg, &
                 air_vapor_pressure_pa, air_density_kg_per_m3, snow_depth_m, vegetation_area_index, &
                 psychrometric_constant_vegetation_pa_per_k, psychrometric_constant_ground_pa_per_k, &
                 wetted_canopy_fraction, sunlit_leaf_area_index, shaded_leaf_area_index, &
                 canopy_wind_parameter, snow_soil_disc_thickness_m, &
                 reference_height_m, zero_plane_displacement_m, momentum_roughness_length_m, &
                 green_vegetation_fraction, &
                 ground_momentum_roughness_length_m, vegetation_emissivity, ground_emissivity, &
                 canopy_liquid_water_mm, snow_cover_fraction, &
                 canopy_ice_mm, snow_soil_temperature_k, snow_soil_thermal_conductivity_w_per_m_per_k, &
                 sunlit_stomatal_resistance_s_per_m, shaded_stomatal_resistance_s_per_m, &
                 surface_resistance_s_per_m, latent_heat_vegetation_j_per_kg, &
                 latent_heat_ground_j_per_kg, absorbed_par_per_sunlit_lai_w_per_m2, &
                 absorbed_par_per_shaded_lai_w_per_m2, growing_season_index, &
                 foliage_nitrogen_percent, co2_partial_pressure_pa, o2_partial_pressure_pa, &
                 root_stress_factor_0_to_1, surface_pressure_pa, &
                 ground_surface_relative_humidity, grid_index_i, grid_index_j, &
                 water_vapor_mixing_ratio_kg_per_kg, precipitation_heat_canopy_w_per_m2, &
                 precipitation_heat_ground_w_per_m2, &
                 canopy_air_vapor_pressure_pa, canopy_air_temperature_k, vegetation_temperature_k, &
                 ground_temperature_k, momentum_drag_coefficient, &
                 heat_exchange_coefficient, grid_spacing_m, lowest_layer_thickness_m, &
                 wind_stress_eastward_n_per_m2, wind_stress_northward_n_per_m2, &
                 ground_net_longwave_w_per_m2, canopy_net_longwave_w_per_m2, ground_sensible_heat_w_per_m2, &
                 canopy_sensible_heat_w_per_m2, ground_evaporation_heat_w_per_m2, &
                 canopy_evaporation_heat_w_per_m2, transpiration_heat_w_per_m2, ground_heat_flux_w_per_m2, &
                 two_meter_air_temperature_k, sunlit_photosynthesis_umol_per_m2_per_s, &
                 shaded_photosynthesis_umol_per_m2_per_s, &
                 cloud_water_mixing_ratio, surface_mixing_ratio_kg_per_kg, surface_pressure_pa, &
                 two_meter_specific_humidity_kg_per_kg, two_meter_heat_conductance_m_per_s, &
                 leaf_exchange_coefficient, under_canopy_exchange_coefficient)

  close(trace_file_unit)

  canopy_energy_balance_residual_w_per_m2 = shortwave_absorbed_by_vegetation_w_per_m2 - &
       (canopy_net_longwave_w_per_m2 + canopy_sensible_heat_w_per_m2 + &
        canopy_evaporation_heat_w_per_m2 + transpiration_heat_w_per_m2)
  ! FLO 10/21: coupling outputs; RAHG from CHUC = 1/RAHG, longwave absorbed by the under-canopy
  ! ground from the IRG expression in VEGE_FLUX loop2 at the final TV (SB = 5.67E-08)
  below_canopy_aerodynamic_resistance_s_per_m = 1.0/under_canopy_exchange_coefficient
  ground_longwave_absorbed_below_canopy_w_per_m2 = ground_emissivity*(1.0 - vegetation_emissivity)* &
       downward_longwave_w_per_m2 + ground_emissivity*vegetation_emissivity*5.67e-08*vegetation_temperature_k**4
  ! Noah-MP: ETRAN = FCTR / LATHEAV (kg/m2/s = mm/s)
  transpiration_rate_mm_per_h = transpiration_heat_w_per_m2/latent_heat_vegetation_j_per_kg*3600.0

  ! ------------------------------------------------------------------------------------------
  ! summary output
  ! ------------------------------------------------------------------------------------------
  summary_file_path = trim(output_directory_path)//'/noahmp_canopy_reference_summary_'//precision_tag//'.csv'
  open(newunit=summary_file_unit, file=trim(summary_file_path), status='replace', action='write')
  write(summary_file_unit,'(a)') 'noahmp_name,value,units,description'
  call write_text_row('PRECISION', precision_tag, 'build: r4 = 32-bit REAL as NWM; r8 = -fdefault-real-8')
  call write_text_row('TIME', trim(forcing_time_string), 'forcing time stamp')
  call write_integer_row('VEGTYP', usgs_vegetation_type, '-', 'USGS vegetation type')
  call write_integer_row('OPT_RAD', radiation_transfer_option, '-', 'radiation transfer option')
  call write_integer_row('OPT_CRS', OPT_CRS, '-', 'stomatal resistance option (1 = Ball-Berry)')
  call write_integer_row('OPT_SFC', OPT_SFC, '-', 'surface layer drag option (1 = SFCDIF1)')
  call write_real_row('JULIAN', julian_day_since_jan1, 'day', 'days since Jan 1 00Z')
  ! forcing and atmosphere
  call write_real_row('SFCTMP', air_temperature_k, 'K', 'air temperature at reference height')
  call write_real_row('SFCPRS', surface_pressure_pa, 'Pa', 'surface pressure')
  call write_real_row('QAIR', specific_humidity_kg_per_kg, 'kg/kg', 'specific humidity')
  call write_real_row('EAIR', air_vapor_pressure_pa, 'Pa', 'air vapor pressure')
  call write_real_row('RHOAIR', air_density_kg_per_m3, 'kg/m3', 'air density')
  call write_real_row('LWDN', downward_longwave_w_per_m2, 'W/m2', 'downward longwave')
  call write_real_row('SOLDN', downward_shortwave_w_per_m2, 'W/m2', 'downward shortwave')
  call write_real_row('SOLAD_VIS', incoming_direct_shortwave_w_per_m2(1), 'W/m2', 'direct visible')
  call write_real_row('SOLAD_NIR', incoming_direct_shortwave_w_per_m2(2), 'W/m2', 'direct near infrared')
  call write_real_row('SOLAI_VIS', incoming_diffuse_shortwave_w_per_m2(1), 'W/m2', 'diffuse visible')
  call write_real_row('SOLAI_NIR', incoming_diffuse_shortwave_w_per_m2(2), 'W/m2', 'diffuse near infrared')
  call write_real_row('UR', wind_speed_reference_height_m_per_s, 'm/s', 'wind speed at ZLVL (>= 1)')
  call write_real_row('COSZ', cosine_solar_zenith_angle, '-', 'cosine solar zenith angle')
  call write_real_row('CO2AIR', co2_partial_pressure_pa, 'Pa', 'CO2 partial pressure')
  call write_real_row('O2AIR', o2_partial_pressure_pa, 'Pa', 'O2 partial pressure')
  ! surface and vegetation
  call write_real_row('HVT', noahmp_parameter_set%HVT, 'm', 'canopy top height')
  call write_real_row('ZPD', zero_plane_displacement_m, 'm', 'zero plane displacement = 0.65*HVT')
  call write_real_row('Z0M', momentum_roughness_length_m, 'm', 'momentum roughness (Z0H = Z0M in VEGE_FLUX)')
  call write_real_row('Z0MG', ground_momentum_roughness_length_m, 'm', 'ground roughness (Z0HG = Z0MG)')
  call write_real_row('ZLVL', reference_height_m, 'm', 'reference height above ground')
  call write_real_row('CWP', canopy_wind_parameter, '-', 'canopy wind parameter')
  call write_real_row('LAI', leaf_area_index, 'm2/m2', 'leaf area index')
  call write_real_row('SAI', stem_area_index, 'm2/m2', 'stem area index')
  call write_real_row('ELAI', exposed_leaf_area_index, 'm2/m2', 'leaf area index above the snow')
  call write_real_row('ESAI', exposed_stem_area_index, 'm2/m2', 'stem area index above the snow')
  call write_real_row('FSNO', snow_cover_fraction, '-', 'snow cover fraction (imposed)')
  call write_real_row('SNOWH', snow_depth_m, 'm', 'snow depth (imposed)')
  call write_real_row('FWET', wetted_canopy_fraction, '-', 'wetted fraction of canopy (imposed)')
  call write_real_row('CANLIQ', canopy_liquid_water_mm, 'mm', 'intercepted liquid water (imposed)')
  call write_real_row('VAI', vegetation_area_index, 'm2/m2', 'LAI + SAI')
  call write_real_row('FVEG', green_vegetation_fraction, '-', 'green vegetation fraction')
  call write_real_row('IGS', growing_season_index, '-', 'growing season index')
  call write_real_row('EMV', vegetation_emissivity, '-', 'vegetation emissivity 1-exp(-VAI)')
  call write_real_row('EMG', ground_emissivity, '-', 'ground emissivity')
  call write_real_row('BTRAN', root_stress_factor_0_to_1, '-', 'root stress factor (imposed)')
  call write_real_row('FOLN', foliage_nitrogen_percent, '%', 'foliage nitrogen')
  call write_real_row('RSURF', surface_resistance_s_per_m, 's/m', 'ground surface resistance (imposed)')
  call write_real_row('RHSUR', ground_surface_relative_humidity, '-', 'ground pore air relative humidity (imposed)')
  call write_real_row('GAMMAV', psychrometric_constant_vegetation_pa_per_k, 'Pa/K', 'psychrometric constant for vegetation')
  call write_real_row('GAMMAG', psychrometric_constant_ground_pa_per_k, 'Pa/K', 'psychrometric constant for ground')
  call write_real_row('LATHEAV', latent_heat_vegetation_j_per_kg, 'J/kg', 'latent heat for vegetation')
  ! initial canopy state and the inputs STOMATA sees on iteration 1
  call write_real_row('TV_INITIAL', initial_vegetation_temperature_k, 'K', 'initial vegetation temperature')
  call write_real_row('TAH_INITIAL', initial_canopy_air_temperature_k, 'K', 'initial canopy air temperature')
  call write_real_row('EAH_INITIAL', initial_canopy_air_vapor_pressure_pa, 'Pa', 'initial canopy air vapor pressure')
  call write_real_row('ESTV_INITIAL', initial_leaf_saturation_vapor_pressure_pa, 'Pa', &
                      'EI passed to STOMATA (esat at initial TV); RB on iteration 1 is in the trace file')
  ! radiation
  call write_real_row('FSUN', sunlit_canopy_fraction, '-', 'sunlit fraction of canopy')
  call write_real_row('LAISUN', sunlit_leaf_area_index, 'm2/m2', 'sunlit leaf area index')
  call write_real_row('LAISHA', shaded_leaf_area_index, 'm2/m2', 'shaded leaf area index')
  call write_real_row('PARSUN', absorbed_par_per_sunlit_lai_w_per_m2, 'W/m2', 'absorbed PAR per unit sunlit LAI')
  call write_real_row('PARSHA', absorbed_par_per_shaded_lai_w_per_m2, 'W/m2', 'absorbed PAR per unit shaded LAI')
  call write_real_row('SAV', shortwave_absorbed_by_vegetation_w_per_m2, 'W/m2', 'shortwave absorbed by vegetation')
  call write_real_row('SAG', shortwave_absorbed_by_ground_w_per_m2, 'W/m2', 'shortwave absorbed by ground')
  call write_real_row('FSA', total_absorbed_shortwave_w_per_m2, 'W/m2', 'total absorbed shortwave')
  call write_real_row('FSR', total_reflected_shortwave_w_per_m2, 'W/m2', 'total reflected shortwave')
  call write_real_row('BGAP', between_canopy_gap_fraction, '-', 'between canopy gap fraction (OPT_RAD = 1 only)')
  call write_real_row('WGAP', within_canopy_gap_fraction, '-', 'within canopy gap fraction (OPT_RAD = 1 only)')
  ! canopy solution
  call write_real_row('TV', vegetation_temperature_k, 'K', 'vegetation temperature')
  call write_real_row('TAH', canopy_air_temperature_k, 'K', 'canopy air temperature')
  call write_real_row('EAH', canopy_air_vapor_pressure_pa, 'Pa', 'canopy air vapor pressure')
  call write_real_row('TG', prescribed_ground_temperature_k, 'K', 'ground skin temperature (prescribed)')
  call write_real_row('TG_AFTER_VEGE_FLUX', ground_temperature_k, 'K', &
                      'TG returned by VEGE_FLUX: OPT_STC = 3 resets it when SNOWH > 0.05 m and TG > TFRZ')
  call write_real_row('RSSUN', sunlit_stomatal_resistance_s_per_m, 's/m', 'sunlit stomatal resistance')
  call write_real_row('RSSHA', shaded_stomatal_resistance_s_per_m, 's/m', 'shaded stomatal resistance')
  call write_real_row('PSNSUN', sunlit_photosynthesis_umol_per_m2_per_s, 'umol/m2/s', 'sunlit leaf photosynthesis')
  call write_real_row('PSNSHA', shaded_photosynthesis_umol_per_m2_per_s, 'umol/m2/s', 'shaded leaf photosynthesis')
  call write_real_row('CM', momentum_drag_coefficient, '-', 'momentum drag coefficient')
  call write_real_row('CH', heat_exchange_coefficient, 'm/s', 'returned CH = CAH = canopy air to ZLVL conductance')
  call write_real_row('CHLEAF', leaf_exchange_coefficient, 'm/s', 'leaf exchange coefficient (CVH)')
  call write_real_row('CHUC', under_canopy_exchange_coefficient, 'm/s', 'under canopy exchange coefficient')
  call write_real_row('IRC', canopy_net_longwave_w_per_m2, 'W/m2', 'canopy net longwave (+ to atm)')
  call write_real_row('SHC', canopy_sensible_heat_w_per_m2, 'W/m2', 'canopy sensible heat (+ to atm)')
  call write_real_row('EVC', canopy_evaporation_heat_w_per_m2, 'W/m2', 'canopy evaporation heat (+ to atm)')
  call write_real_row('TR', transpiration_heat_w_per_m2, 'W/m2', 'transpiration heat (+ to atm)')
  call write_real_row('IRG', ground_net_longwave_w_per_m2, 'W/m2', 'ground net longwave under canopy')
  call write_real_row('SHG', ground_sensible_heat_w_per_m2, 'W/m2', 'ground sensible heat under canopy')
  call write_real_row('EVG', ground_evaporation_heat_w_per_m2, 'W/m2', 'ground evaporation heat under canopy')
  call write_real_row('GH', ground_heat_flux_w_per_m2, 'W/m2', 'ground heat flux implied by prescribed TG (+ into soil)')
  call write_real_row('T2MV', two_meter_air_temperature_k, 'K', '2 m air temperature over vegetation')
  call write_real_row('Q2V', two_meter_specific_humidity_kg_per_kg, 'kg/kg', '2 m specific humidity over vegetation')
  call write_real_row('TAUXV', wind_stress_eastward_n_per_m2, 'N/m2', 'eastward wind stress')
  call write_real_row('TAUYV', wind_stress_northward_n_per_m2, 'N/m2', 'northward wind stress')
  call write_real_row('CANOPY_RESIDUAL', canopy_energy_balance_residual_w_per_m2, 'W/m2', &
                      'SAV - (IRC + SHC + EVC + TR)')
  call write_real_row('ETRAN_MM_PER_H', transpiration_rate_mm_per_h, 'mm/h', 'transpiration = TR/LATHEAV')
  call write_real_row('RAHG', below_canopy_aerodynamic_resistance_s_per_m, 's/m', 'under-canopy resistance = 1/CHUC')
  call write_real_row('GROUND_LW_ABSORBED', ground_longwave_absorbed_below_canopy_w_per_m2, 'W/m2', &
                      'absorbed part of IRG: EMG*((1-EMV)*LWDN + EMV*SB*TV**4)')
  close(summary_file_unit)

  write(*,'(a,a)') 'wrote ', trim(summary_file_path)
  write(*,'(a,a)') 'wrote ', trim(trace_file_path)

contains

  subroutine write_real_row(noahmp_name, value_real, units_text, description_text)
    character(len=*), intent(in) :: noahmp_name, units_text, description_text
    real,             intent(in) :: value_real
    write(summary_file_unit,'(a,",",es24.16e3,",",a,",",a)') trim(noahmp_name), value_real, &
         trim(units_text), trim(description_text)
  end subroutine write_real_row

  subroutine write_integer_row(noahmp_name, value_integer, units_text, description_text)
    character(len=*), intent(in) :: noahmp_name, units_text, description_text
    integer,          intent(in) :: value_integer
    write(summary_file_unit,'(a,",",i0,",",a,",",a)') trim(noahmp_name), value_integer, &
         trim(units_text), trim(description_text)
  end subroutine write_integer_row

  subroutine write_text_row(noahmp_name, value_text, description_text)
    character(len=*), intent(in) :: noahmp_name, value_text, description_text
    write(summary_file_unit,'(a,",",a,",-,",a)') trim(noahmp_name), trim(value_text), &
         trim(description_text)
  end subroutine write_text_row

  ! forcing_file is a string key; read it from the shared configuration file
  subroutine read_forcing_file_path_from_configuration(config_path, forcing_path)
    character(len=*), intent(in)  :: config_path
    character(len=*), intent(out) :: forcing_path
    character(len=2048) :: line_text
    integer :: file_unit, io_status, equals_position, comment_position
    forcing_path = ' '
    open(newunit=file_unit, file=trim(config_path), status='old', action='read')
    do
      read(file_unit,'(a)', iostat=io_status) line_text
      if (io_status /= 0) exit
      line_text = adjustl(line_text)
      if (line_text(1:17) /= 'forcing_csv_file ' .and. line_text(1:17) /= 'forcing_csv_file=') cycle   ! FLO 10/21
      equals_position = index(line_text, '=')
      forcing_path = adjustl(line_text(equals_position+1:))
      comment_position = index(forcing_path, '#')
      if (comment_position > 0) forcing_path(comment_position:) = ' '
      forcing_path = adjustl(forcing_path)
      exit
    end do
    close(file_unit)
    if (len_trim(forcing_path) == 0) stop 'ERROR: forcing_csv_file not found in configuration'
    ! FLO 10/21: a relative forcing_file path is relative to the working directory, as in
    ! the C PET reader, so both codes read the same file when run from the project root.
  end subroutine read_forcing_file_path_from_configuration

end program noahmp_canopy_reference_driver

! ---- stubs for the WRF utilities referenced by the Noah-MP module
subroutine wrf_message(message_text)
  character(len=*), intent(in) :: message_text
  write(*,'(a)') trim(message_text)
end subroutine wrf_message

subroutine wrf_error_fatal(message_text)
  character(len=*), intent(in) :: message_text
  write(*,'(a,a)') 'FATAL: ', trim(message_text)
  stop 2
end subroutine wrf_error_fatal
