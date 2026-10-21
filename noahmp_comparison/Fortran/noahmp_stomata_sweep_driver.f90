! ==========================================================================================
! noahmp_stomata_sweep_driver.f90
!
! Purpose:
!   Stage 1 reference: call the NWM Noah-MP STOMATA routine (Ball-Berry, Collatz C3/C4
!   photosynthesis) directly over a grid of leaf-level inputs, isolating the stomata
!   kernel from the canopy energy balance.  Every input STOMATA receives is written on
!   each row, so the C code can be driven with exactly the same values.
!
!   Swept:  absorbed PAR per unit LAI   APAR  =   0, 25, 50, ..., 500 W/m2  (21 values)
!           foliage temperature         TV    = 278.16 .. 308.16 K by 5 K    (7 values)
!           canopy air vapor pressure   EA    = 0.2, 0.5, 0.8, 1.0 * EI      (4 values)
!           leaf boundary resistance    RB    = 20, 50, 100 s/m              (3 values)
!   Fixed from configuration and forcing: SFCTMP, SFCPRS, CO2, O2, BTRAN, FOLN, and all
!   Ball-Berry parameters.  IGS = 1 (growing season).  EI = esat(TV) from Noah-MP ESAT,
!   evaluated as VEGE_FLUX does (TDC limits, water branch above TFRZ = 273.16 K).
!
! Usage:
!   noahmp_stomata_sweep_driver_<r4|r8> config_file output_directory
!
! Output:
!   noahmp_stomata_sweep_<tag>.csv
! ==========================================================================================
program noahmp_stomata_sweep_driver

  use module_sf_noahmplsm
  use harness_text_utilities
  use noahmp_reference_setup
  implicit none

  integer, parameter :: number_of_apar_values = 21
  integer, parameter :: number_of_temperature_values = 7
  integer, parameter :: number_of_vapor_pressure_fractions = 4
  integer, parameter :: number_of_boundary_resistance_values = 3
  integer, parameter :: usgs_vegetation_type = 15
  integer, parameter :: grid_index_i = 1, grid_index_j = 1
  real,    parameter :: minimum_positive_value = 1.0e-6     ! [MPE] as passed by VEGE_FLUX

  real :: vapor_pressure_fraction_list(number_of_vapor_pressure_fractions)
  real :: boundary_resistance_list_s_per_m(number_of_boundary_resistance_values)

  character(len=512) :: configuration_file_path, output_directory_path, output_file_path
  character(len=2)   :: precision_tag
  integer :: output_file_unit
  integer :: apar_index, temperature_index, vapor_pressure_index, boundary_resistance_index
  type(noahmp_parameters) :: noahmp_parameter_set

  real :: air_temperature_k, surface_pressure_pa
  real :: co2_partial_pressure_pa, o2_partial_pressure_pa
  real :: root_stress_factor_0_to_1, foliage_nitrogen_percent, growing_season_index
  real :: absorbed_par_per_unit_lai_w_per_m2, foliage_temperature_k
  real :: leaf_internal_vapor_pressure_pa, canopy_air_vapor_pressure_pa
  real :: leaf_boundary_layer_resistance_s_per_m
  real :: stomatal_resistance_s_per_m, photosynthesis_umol_per_m2_per_s
  real :: temperature_celsius, saturation_vapor_pressure_water_pa, saturation_vapor_pressure_ice_pa
  real :: d_saturation_vapor_pressure_water_pa_per_k, d_saturation_vapor_pressure_ice_pa_per_k
  character(len=512) :: forcing_file_path

  vapor_pressure_fraction_list     = (/ 0.2, 0.5, 0.8, 1.0 /)
  boundary_resistance_list_s_per_m = (/ 20.0, 50.0, 100.0 /)

  if (command_argument_count() /= 2) then
    write(*,'(a)') 'usage: noahmp_stomata_sweep_driver_<r4|r8> config_file output_directory'
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
  call resolve_forcing_file_path(configuration_file_path, forcing_file_path)
  call read_forcing_csv_first_data_row(forcing_file_path)
  call fill_noahmp_parameters_from_configuration(noahmp_parameter_set)
  call set_nwm_noahmp_physics_options(get_required_integer_value('radiation_transfer_option'))

  air_temperature_k   = get_forcing_column_value('TMP')
  surface_pressure_pa = get_forcing_column_value('PRES')
  co2_partial_pressure_pa = noahmp_parameter_set%CO2 * surface_pressure_pa
  o2_partial_pressure_pa  = noahmp_parameter_set%O2  * surface_pressure_pa
  root_stress_factor_0_to_1 = get_required_real_value('root_stress_factor_0_to_1')
  foliage_nitrogen_percent    = get_required_real_value('foliage_nitrogen_percent')
  growing_season_index        = 1.0

  output_file_path = trim(output_directory_path)//'/noahmp_stomata_sweep_'//precision_tag//'.csv'
  open(newunit=output_file_unit, file=trim(output_file_path), status='replace', action='write')
  write(output_file_unit,'(a)') 'absorbed_par_per_unit_lai_w_per_m2,foliage_temperature_k,'// &
       'leaf_internal_vapor_pressure_pa,canopy_air_vapor_pressure_pa,'// &
       'leaf_boundary_layer_resistance_s_per_m,reference_air_temperature_k,surface_pressure_pa,'// &
       'co2_partial_pressure_pa,o2_partial_pressure_pa,root_stress_factor_0_to_1,'// &
       'foliage_nitrogen_percent,growing_season_index,'// &
       'stomatal_resistance_s_per_m,photosynthesis_umol_co2_per_m2_per_s'

  do boundary_resistance_index = 1, number_of_boundary_resistance_values
    leaf_boundary_layer_resistance_s_per_m = boundary_resistance_list_s_per_m(boundary_resistance_index)
    do temperature_index = 1, number_of_temperature_values
      foliage_temperature_k = 278.16 + 5.0*real(temperature_index - 1)
      temperature_celsius = min(50.0, max(-50.0, foliage_temperature_k - 273.16))
      call ESAT(temperature_celsius, saturation_vapor_pressure_water_pa, saturation_vapor_pressure_ice_pa, &
                d_saturation_vapor_pressure_water_pa_per_k, d_saturation_vapor_pressure_ice_pa_per_k)
      if (temperature_celsius > 0.0) then
        leaf_internal_vapor_pressure_pa = saturation_vapor_pressure_water_pa
      else
        leaf_internal_vapor_pressure_pa = saturation_vapor_pressure_ice_pa
      end if
      do vapor_pressure_index = 1, number_of_vapor_pressure_fractions
        canopy_air_vapor_pressure_pa = vapor_pressure_fraction_list(vapor_pressure_index) * &
                                       leaf_internal_vapor_pressure_pa
        do apar_index = 1, number_of_apar_values
          absorbed_par_per_unit_lai_w_per_m2 = 25.0*real(apar_index - 1)

          call STOMATA(noahmp_parameter_set, usgs_vegetation_type, minimum_positive_value, &
                       absorbed_par_per_unit_lai_w_per_m2, foliage_nitrogen_percent, &
                       grid_index_i, grid_index_j, &
                       foliage_temperature_k, leaf_internal_vapor_pressure_pa, &
                       canopy_air_vapor_pressure_pa, air_temperature_k, surface_pressure_pa, &
                       o2_partial_pressure_pa, co2_partial_pressure_pa, growing_season_index, &
                       root_stress_factor_0_to_1, leaf_boundary_layer_resistance_s_per_m, &
                       stomatal_resistance_s_per_m, photosynthesis_umol_per_m2_per_s)

          write(output_file_unit,'(13(es24.16e3,","),es24.16e3)') &
               absorbed_par_per_unit_lai_w_per_m2, foliage_temperature_k, &
               leaf_internal_vapor_pressure_pa, canopy_air_vapor_pressure_pa, &
               leaf_boundary_layer_resistance_s_per_m, air_temperature_k, surface_pressure_pa, &
               co2_partial_pressure_pa, o2_partial_pressure_pa, root_stress_factor_0_to_1, &
               foliage_nitrogen_percent, growing_season_index, &
               stomatal_resistance_s_per_m, photosynthesis_umol_per_m2_per_s
        end do
      end do
    end do
  end do
  close(output_file_unit)
  write(*,'(a,a)') 'wrote ', trim(output_file_path)

contains

  subroutine resolve_forcing_file_path(config_path, forcing_path)
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
  end subroutine resolve_forcing_file_path

end program noahmp_stomata_sweep_driver

subroutine wrf_message(message_text)
  character(len=*), intent(in) :: message_text
  write(*,'(a)') trim(message_text)
end subroutine wrf_message

subroutine wrf_error_fatal(message_text)
  character(len=*), intent(in) :: message_text
  write(*,'(a,a)') 'FATAL: ', trim(message_text)
  stop 2
end subroutine wrf_error_fatal
