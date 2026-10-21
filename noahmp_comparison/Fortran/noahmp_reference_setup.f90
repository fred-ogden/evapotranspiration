! ==========================================================================================
! noahmp_reference_setup.f90
!
! Purpose:
!   Shared setup for the Noah-MP Ball-Berry reference drivers:
!     - fill the noahmp_parameters derived type from the key = value configuration,
!       setting every field touched by RADIATION, VEGE_FLUX, RAGRB, SFCDIF1 and STOMATA;
!     - set the Noah-MP physics options to the NWM configuration;
!     - convert the forcing time stamp to the Noah-MP JULIAN day and year length.
! ==========================================================================================
module noahmp_reference_setup

  use module_sf_noahmplsm
  use harness_text_utilities
  implicit none
  private

  public :: fill_noahmp_parameters_from_configuration
  public :: set_nwm_noahmp_physics_options
  public :: convert_time_string_to_noahmp_julian_day
  public :: interpolate_monthly_lai_sai_as_noahmp_phenology
  public :: apply_snow_burial_as_noahmp_phenology

contains

  ! ----------------------------------------------------------------------------------------
  ! Fill the Noah-MP parameter structure.  Fields not used by the routines driven here
  ! are set to zero (or to the table value where cheap), so nothing is left undefined.
  ! ----------------------------------------------------------------------------------------
  subroutine fill_noahmp_parameters_from_configuration(noahmp_parameter_set)
    type(noahmp_parameters), intent(out) :: noahmp_parameter_set

    ! ---- surface type flags (USGS table)
    noahmp_parameter_set%URBAN_FLAG = .false.
    noahmp_parameter_set%ISWATER    = 16
    noahmp_parameter_set%ISBARREN   = 19
    noahmp_parameter_set%ISICE      = 24
    noahmp_parameter_set%ISCROP     = 2
    noahmp_parameter_set%EBLFOREST  = 13

    ! ---- canopy structure
    noahmp_parameter_set%CH2OP  = 0.1                ! max intercepted water per unit LAI+SAI (mm), unused with no canopy water
    noahmp_parameter_set%DLEAF  = get_required_real_value('characteristic_leaf_dimension_m')
    noahmp_parameter_set%Z0MVT  = get_required_real_value('momentum_roughness_length_m')
    noahmp_parameter_set%HVT    = get_required_real_value('canopy_top_height_m')
    noahmp_parameter_set%HVB    = get_required_real_value('canopy_bottom_height_m')
    noahmp_parameter_set%DEN    = 0.10               ! not used in this Noah-MP version
    noahmp_parameter_set%RC     = get_required_real_value('tree_crown_radius_m')
    noahmp_parameter_set%MFSNO  = 2.5                ! snow only
    noahmp_parameter_set%XL     = get_required_real_value('leaf_orientation_index')
    noahmp_parameter_set%CWPVT  = get_required_real_value('canopy_wind_parameter')

    noahmp_parameter_set%LAIM(1)  = get_required_real_value('lai_month_01')
    noahmp_parameter_set%LAIM(2)  = get_required_real_value('lai_month_02')
    noahmp_parameter_set%LAIM(3)  = get_required_real_value('lai_month_03')
    noahmp_parameter_set%LAIM(4)  = get_required_real_value('lai_month_04')
    noahmp_parameter_set%LAIM(5)  = get_required_real_value('lai_month_05')
    noahmp_parameter_set%LAIM(6)  = get_required_real_value('lai_month_06')
    noahmp_parameter_set%LAIM(7)  = get_required_real_value('lai_month_07')
    noahmp_parameter_set%LAIM(8)  = get_required_real_value('lai_month_08')
    noahmp_parameter_set%LAIM(9)  = get_required_real_value('lai_month_09')
    noahmp_parameter_set%LAIM(10) = get_required_real_value('lai_month_10')
    noahmp_parameter_set%LAIM(11) = get_required_real_value('lai_month_11')
    noahmp_parameter_set%LAIM(12) = get_required_real_value('lai_month_12')
    noahmp_parameter_set%SAIM(1)  = get_required_real_value('sai_month_01')
    noahmp_parameter_set%SAIM(2)  = get_required_real_value('sai_month_02')
    noahmp_parameter_set%SAIM(3)  = get_required_real_value('sai_month_03')
    noahmp_parameter_set%SAIM(4)  = get_required_real_value('sai_month_04')
    noahmp_parameter_set%SAIM(5)  = get_required_real_value('sai_month_05')
    noahmp_parameter_set%SAIM(6)  = get_required_real_value('sai_month_06')
    noahmp_parameter_set%SAIM(7)  = get_required_real_value('sai_month_07')
    noahmp_parameter_set%SAIM(8)  = get_required_real_value('sai_month_08')
    noahmp_parameter_set%SAIM(9)  = get_required_real_value('sai_month_09')
    noahmp_parameter_set%SAIM(10) = get_required_real_value('sai_month_10')
    noahmp_parameter_set%SAIM(11) = get_required_real_value('sai_month_11')
    noahmp_parameter_set%SAIM(12) = get_required_real_value('sai_month_12')

    ! ---- optical properties: index 1 = visible, 2 = near infrared
    noahmp_parameter_set%RHOL(1) = get_required_real_value('leaf_reflectance_visible')
    noahmp_parameter_set%RHOL(2) = get_required_real_value('leaf_reflectance_nir')
    noahmp_parameter_set%RHOS(1) = get_required_real_value('stem_reflectance_visible')
    noahmp_parameter_set%RHOS(2) = get_required_real_value('stem_reflectance_nir')
    noahmp_parameter_set%TAUL(1) = get_required_real_value('leaf_transmittance_visible')
    noahmp_parameter_set%TAUL(2) = get_required_real_value('leaf_transmittance_nir')
    noahmp_parameter_set%TAUS(1) = get_required_real_value('stem_transmittance_visible')
    noahmp_parameter_set%TAUS(2) = get_required_real_value('stem_transmittance_nir')

    ! ---- photosynthesis and stomata
    noahmp_parameter_set%C3PSN  = get_required_real_value('c3_photosynthesis_flag')
    noahmp_parameter_set%KC25   = get_required_real_value('co2_michaelis_menten_constant_25c_pa')
    noahmp_parameter_set%AKC    = get_required_real_value('q10_for_kc25')
    noahmp_parameter_set%KO25   = get_required_real_value('o2_michaelis_menten_constant_25c_pa')
    noahmp_parameter_set%AKO    = get_required_real_value('q10_for_ko25')
    noahmp_parameter_set%VCMX25 = get_required_real_value('maximum_carboxylation_rate_25c_umol_per_m2_per_s')
    noahmp_parameter_set%AVCMX  = get_required_real_value('q10_for_vcmx25')
    noahmp_parameter_set%BP     = get_required_real_value('minimum_leaf_conductance_umol_per_m2_per_s')
    noahmp_parameter_set%MP     = get_required_real_value('conductance_photosynthesis_slope')
    noahmp_parameter_set%QE25   = get_required_real_value('quantum_efficiency_25c')
    noahmp_parameter_set%AQE    = 1.0
    noahmp_parameter_set%FOLNMX = get_required_real_value('foliage_nitrogen_at_saturation_percent')
    noahmp_parameter_set%TMIN   = get_required_real_value('minimum_temperature_for_photosynthesis_k')

    ! ---- carbon model and Jarvis fields: not used by the Ball-Berry path, table values for type 15
    noahmp_parameter_set%SLA    = 80.0
    noahmp_parameter_set%DILEFC = 0.80
    noahmp_parameter_set%DILEFW = 0.20
    noahmp_parameter_set%FRAGR  = 0.10
    noahmp_parameter_set%LTOVRC = 0.5
    noahmp_parameter_set%RMF25  = 3.00
    noahmp_parameter_set%RMS25  = 0.80
    noahmp_parameter_set%RMR25  = 0.03
    noahmp_parameter_set%ARM    = 2.0
    noahmp_parameter_set%MRP    = 0.30
    noahmp_parameter_set%WRRAT  = 30.0
    noahmp_parameter_set%WDPOOL = 1.0
    noahmp_parameter_set%TDLEF  = 268.0
    noahmp_parameter_set%NROOT  = 4
    noahmp_parameter_set%RGL    = 30.0
    noahmp_parameter_set%RSMIN  = 125.0
    noahmp_parameter_set%HS     = 51.93
    noahmp_parameter_set%TOPT   = 298.0
    noahmp_parameter_set%RSMAX  = 5000.0
    noahmp_parameter_set%SLAREA = 0.0223
    noahmp_parameter_set%EPS    = 0.0

    ! ---- rad section (soil albedos are overridden by the imposed ground albedos)
    noahmp_parameter_set%ALBSAT = (/ 0.15, 0.30 /)
    noahmp_parameter_set%ALBDRY = (/ 0.27, 0.54 /)
    noahmp_parameter_set%ALBICE = (/ 0.80, 0.55 /)
    noahmp_parameter_set%ALBLAK = (/ 0.60, 0.40 /)
    noahmp_parameter_set%OMEGAS = (/ 0.8, 0.4 /)   ! used in TWOSTREAM, weighted by FWET (= 0 here)
    noahmp_parameter_set%BETADS = 0.5
    noahmp_parameter_set%BETAIS = 0.5
    noahmp_parameter_set%EG(1)  = get_required_real_value('soil_emissivity')
    noahmp_parameter_set%EG(2)  = 0.98

    ! ---- global section
    noahmp_parameter_set%CO2          = get_required_real_value('co2_mole_fraction')
    noahmp_parameter_set%O2           = get_required_real_value('o2_mole_fraction')
    noahmp_parameter_set%TIMEAN       = 10.5
    noahmp_parameter_set%FSATMX       = 0.38
    noahmp_parameter_set%Z0SNO        = 0.002
    noahmp_parameter_set%SSI          = 0.03
    noahmp_parameter_set%SNOW_RET_FAC = 5.0e-5
    noahmp_parameter_set%SWEMX        = 1.0
    noahmp_parameter_set%TAU0         = 1.0e6
    noahmp_parameter_set%GRAIN_GROWTH = 5000.0
    noahmp_parameter_set%EXTRA_GROWTH = 10.0
    noahmp_parameter_set%DIRT_SOOT    = 0.3
    noahmp_parameter_set%BATS_COSZ    = 2.0
    noahmp_parameter_set%BATS_VIS_NEW = 0.95
    noahmp_parameter_set%BATS_NIR_NEW = 0.65
    noahmp_parameter_set%BATS_VIS_AGE = 0.2
    noahmp_parameter_set%BATS_NIR_AGE = 0.5
    noahmp_parameter_set%BATS_VIS_DIR = 0.4
    noahmp_parameter_set%BATS_NIR_DIR = 0.4
    noahmp_parameter_set%RSURF_SNOW   = 50.0
    noahmp_parameter_set%RSURF_EXP    = 5.0
    noahmp_parameter_set%IMPERV       = 0.0
    noahmp_parameter_set%SCAMAX       = 1.0
    noahmp_parameter_set%SWE_LIMIT    = 2000.0

    ! ---- soil and general parameters: not used with a prescribed ground temperature
    noahmp_parameter_set%BEXP   = 4.0
    noahmp_parameter_set%SMCDRY = 0.05
    noahmp_parameter_set%SMCWLT = 0.05
    noahmp_parameter_set%SMCREF = 0.3
    noahmp_parameter_set%SMCMAX = 0.45
    noahmp_parameter_set%PSISAT = 0.2
    noahmp_parameter_set%DKSAT  = 1.0e-5
    noahmp_parameter_set%DWSAT  = 1.0e-5
    noahmp_parameter_set%QUARTZ = 0.5
    noahmp_parameter_set%F1     = 0.0
    noahmp_parameter_set%SLOPE  = 0.1
    noahmp_parameter_set%CSOIL  = 2.0e6
    noahmp_parameter_set%ZBOT   = -8.0
    noahmp_parameter_set%CZIL   = 0.1
    noahmp_parameter_set%REFDK  = 2.0e-6
    noahmp_parameter_set%REFKDT = 3.0
    noahmp_parameter_set%KDT    = 0.0
    noahmp_parameter_set%FRZX   = 0.0
  end subroutine fill_noahmp_parameters_from_configuration

  ! ----------------------------------------------------------------------------------------
  ! NWM Noah-MP physics options relevant to the vegetated canopy path.
  ! ----------------------------------------------------------------------------------------
  subroutine set_nwm_noahmp_physics_options(radiation_transfer_option)
    integer, intent(in) :: radiation_transfer_option
    DVEG     = 4      ! table LAI, maximum green vegetation fraction
    OPT_CRS  = 1      ! Ball-Berry
    OPT_BTR  = 1      ! Noah soil moisture factor (BTRAN prescribed here anyway)
    OPT_SFC  = 1      ! Monin-Obukhov (SFCDIF1)
    OPT_RAD  = radiation_transfer_option
    OPT_ALB  = 1      ! BATS snow albedo (no snow here)
    OPT_STC  = 3      ! NWM temperature time scheme (only acts with snow)
    OPT_RSF  = 4      ! not used: RSURF imposed
    OPT_FRZ  = 1
    OPT_INF  = 1
    OPT_RUN  = 3
    OPT_SNF  = 1
    OPT_TBOT = 2
    OPT_SOIL = 1
    OPT_PEDO = 1
    OPT_CROP = 0
    OPT_IMPERV = 9
  end subroutine set_nwm_noahmp_physics_options

  ! ----------------------------------------------------------------------------------------
  ! Time stamp "YYYY-MM-DD HH:MM:SS" -> Noah-MP JULIAN (days since Jan 1 00Z, fractional)
  ! and YEARLEN (365 or 366), as the HRLDAS/NWM driver supplies them.
  ! ----------------------------------------------------------------------------------------
  subroutine convert_time_string_to_noahmp_julian_day(time_string, julian_day_since_jan1, &
                                                      year_length_days)
    character(len=*), intent(in)  :: time_string
    real,             intent(out) :: julian_day_since_jan1
    integer,          intent(out) :: year_length_days
    integer :: year_number, month_number, day_of_month, hour_of_day, minute_of_hour
    integer :: second_of_minute, month_index, day_of_year
    integer :: days_in_month(12)
    logical :: is_leap_year

    read(time_string(1:4),*) year_number
    read(time_string(6:7),*) month_number
    read(time_string(9:10),*) day_of_month
    read(time_string(12:13),*) hour_of_day
    read(time_string(15:16),*) minute_of_hour
    read(time_string(18:19),*) second_of_minute

    is_leap_year = (mod(year_number,4) == 0 .and. mod(year_number,100) /= 0) .or. &
                   (mod(year_number,400) == 0)
    days_in_month = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
    if (is_leap_year) days_in_month(2) = 29
    year_length_days = sum(days_in_month)

    day_of_year = day_of_month
    do month_index = 1, month_number - 1
      day_of_year = day_of_year + days_in_month(month_index)
    end do
    julian_day_since_jan1 = real(day_of_year - 1) + real(hour_of_day)/24.0 &
                          + real(minute_of_hour)/1440.0 + real(second_of_minute)/86400.0
  end subroutine convert_time_string_to_noahmp_julian_day

  ! ----------------------------------------------------------------------------------------
  ! Monthly LAI/SAI interpolation copied from PHENOLOGY (DVEG = 1, 3, 4 branch), followed
  ! by the same small-value checks.  No snow burial (SNOWH = 0), so ELAI = LAI, ESAI = SAI.
  ! ----------------------------------------------------------------------------------------
  subroutine interpolate_monthly_lai_sai_as_noahmp_phenology(noahmp_parameter_set, latitude_degrees, &
                   julian_day_since_jan1, year_length_days, leaf_area_index, stem_area_index)
    type(noahmp_parameters), intent(in)  :: noahmp_parameter_set
    real,                    intent(in)  :: latitude_degrees
    real,                    intent(in)  :: julian_day_since_jan1
    integer,                 intent(in)  :: year_length_days
    real,                    intent(out) :: leaf_area_index
    real,                    intent(out) :: stem_area_index
    real    :: day_value, month_position, weight_month_1, weight_month_2
    integer :: month_index_1, month_index_2

    if (latitude_degrees >= 0.0) then
      day_value = julian_day_since_jan1
    else
      day_value = mod(julian_day_since_jan1 + (0.5*year_length_days), real(year_length_days))
    end if
    month_position = 12.0 * day_value / real(year_length_days)
    month_index_1  = int(month_position + 0.5)     ! Fortran integer assignment truncates, as in Noah-MP
    month_index_2  = month_index_1 + 1
    weight_month_1 = (month_index_1 + 0.5) - month_position
    weight_month_2 = 1.0 - weight_month_1
    if (month_index_1 < 1)  month_index_1 = 12
    if (month_index_2 > 12) month_index_2 = 1

    leaf_area_index = weight_month_1*noahmp_parameter_set%LAIM(month_index_1) + &
                      weight_month_2*noahmp_parameter_set%LAIM(month_index_2)
    stem_area_index = weight_month_1*noahmp_parameter_set%SAIM(month_index_1) + &
                      weight_month_2*noahmp_parameter_set%SAIM(month_index_2)

    if (stem_area_index < 0.05) stem_area_index = 0.0
    if (leaf_area_index < 0.05 .or. stem_area_index == 0.0) leaf_area_index = 0.0
  end subroutine interpolate_monthly_lai_sai_as_noahmp_phenology

  ! ----------------------------------------------------------------------------------------
  ! FLO 10/21: burial of the canopy by snow, copied from PHENOLOGY ("buried by snow" block and
  ! the ESAI/ELAI checks that follow it, CROPTYPE = 0). Returns the exposed indices ELAI, ESAI.
  ! ----------------------------------------------------------------------------------------
  subroutine apply_snow_burial_as_noahmp_phenology(noahmp_parameter_set, snow_depth_m, &
                   leaf_area_index, stem_area_index, exposed_leaf_area_index, exposed_stem_area_index)
    type(noahmp_parameters), intent(in)  :: noahmp_parameter_set
    real,                    intent(in)  :: snow_depth_m              ! [SNOWH]
    real,                    intent(in)  :: leaf_area_index           ! [LAI]
    real,                    intent(in)  :: stem_area_index           ! [SAI]
    real,                    intent(out) :: exposed_leaf_area_index   ! [ELAI]
    real,                    intent(out) :: exposed_stem_area_index   ! [ESAI]
    real :: buried_depth_m                                            ! [DB]
    real :: buried_fraction                                           ! [FB]
    real :: snow_height_for_short_vegetation_m                        ! [SNOWHC]

    buried_depth_m  = min( max(snow_depth_m - noahmp_parameter_set%HVB, 0.), &
                           noahmp_parameter_set%HVT - noahmp_parameter_set%HVB )
    buried_fraction = buried_depth_m / max(1.e-06, noahmp_parameter_set%HVT - noahmp_parameter_set%HVB)
    if (noahmp_parameter_set%HVT > 0. .and. noahmp_parameter_set%HVT <= 1.0) then
      snow_height_for_short_vegetation_m = noahmp_parameter_set%HVT*exp(-snow_depth_m/0.2)
      if (snow_height_for_short_vegetation_m <= snow_depth_m) then
        buried_fraction = 1.
      else
        buried_fraction = snow_depth_m/snow_height_for_short_vegetation_m
      end if
    end if
    exposed_leaf_area_index = leaf_area_index*(1.-buried_fraction)
    exposed_stem_area_index = stem_area_index*(1.-buried_fraction)
    if (exposed_stem_area_index < 0.05) exposed_stem_area_index = 0.0
    if (exposed_leaf_area_index < 0.05 .or. exposed_stem_area_index == 0.0) exposed_leaf_area_index = 0.0
  end subroutine apply_snow_burial_as_noahmp_phenology

end module noahmp_reference_setup
