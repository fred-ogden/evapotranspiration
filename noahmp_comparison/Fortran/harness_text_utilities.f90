! ==========================================================================================
! harness_text_utilities.f90
!
! Purpose:
!   Small text utilities for the Noah-MP Ball-Berry reference drivers:
!     - read a "key = value" configuration file (lines starting with # are comments,
!       text after a # on a value line is ignored), and look up real or integer values
!       by key with an error stop if a required key is missing;
!     - read the header row and one data row of a comma separated forcing file and
!       look up a column value by its header name.
!
! All REAL variables are default REAL, so they follow -fdefault-real-8 when the
! reference is built in double precision, exactly as the Noah-MP module does.
! ==========================================================================================
module harness_text_utilities

  implicit none
  private

  public :: read_key_value_configuration_file
  public :: get_required_real_value
  public :: get_required_integer_value
  public :: read_forcing_csv_first_data_row
  public :: get_forcing_column_value
  public :: get_forcing_time_string

  integer, parameter :: maximum_line_length_characters = 2048
  integer, parameter :: maximum_number_of_configuration_entries = 500
  integer, parameter :: maximum_number_of_forcing_columns = 100

  character(len=128) :: configuration_key_list(maximum_number_of_configuration_entries)
  character(len=256) :: configuration_value_list(maximum_number_of_configuration_entries)
  integer            :: number_of_configuration_entries = 0

  character(len=64)  :: forcing_column_name_list(maximum_number_of_forcing_columns)
  character(len=128) :: forcing_column_value_list(maximum_number_of_forcing_columns)
  integer            :: number_of_forcing_columns = 0

contains

  ! ----------------------------------------------------------------------------------------
  ! Read every "key = value" line of a configuration file into the module lists.
  ! ----------------------------------------------------------------------------------------
  subroutine read_key_value_configuration_file(configuration_file_path)
    character(len=*), intent(in) :: configuration_file_path
    character(len=maximum_line_length_characters) :: line_text
    integer :: file_unit, io_status, equals_position, comment_position

    open(newunit=file_unit, file=trim(configuration_file_path), status='old', action='read', &
         iostat=io_status)
    if (io_status /= 0) then
      write(*,'(a,a)') 'ERROR: cannot open configuration file ', trim(configuration_file_path)
      stop 1
    end if

    number_of_configuration_entries = 0
    do
      read(file_unit,'(a)', iostat=io_status) line_text
      if (io_status /= 0) exit
      line_text = adjustl(line_text)
      if (len_trim(line_text) == 0) cycle
      if (line_text(1:1) == '#') cycle
      equals_position = index(line_text, '=')
      if (equals_position == 0) cycle
      number_of_configuration_entries = number_of_configuration_entries + 1
      if (number_of_configuration_entries > maximum_number_of_configuration_entries) then
        write(*,'(a)') 'ERROR: too many configuration entries'
        stop 1
      end if
      configuration_key_list(number_of_configuration_entries) = &
           adjustl(line_text(1:equals_position-1))
      configuration_value_list(number_of_configuration_entries) = &
           adjustl(line_text(equals_position+1:))
      ! strip trailing comment on the value
      comment_position = index(configuration_value_list(number_of_configuration_entries), '#')
      if (comment_position > 0) then
        configuration_value_list(number_of_configuration_entries)(comment_position:) = ' '
      end if
    end do
    close(file_unit)
  end subroutine read_key_value_configuration_file

  ! ----------------------------------------------------------------------------------------
  ! Return the index of a configuration key, or stop if it is missing.
  ! ----------------------------------------------------------------------------------------
  integer function find_configuration_entry_index(key_name)
    character(len=*), intent(in) :: key_name
    integer :: entry_index
    find_configuration_entry_index = 0
    do entry_index = 1, number_of_configuration_entries
      if (trim(configuration_key_list(entry_index)) == trim(key_name)) then
        find_configuration_entry_index = entry_index
        return
      end if
    end do
    write(*,'(a,a)') 'ERROR: required configuration key not found: ', trim(key_name)
    stop 1
  end function find_configuration_entry_index

  real function get_required_real_value(key_name)
    character(len=*), intent(in) :: key_name
    integer :: entry_index, io_status
    entry_index = find_configuration_entry_index(key_name)
    read(configuration_value_list(entry_index),*, iostat=io_status) get_required_real_value
    if (io_status /= 0) then
      write(*,'(a,a)') 'ERROR: cannot read a real value for key ', trim(key_name)
      stop 1
    end if
  end function get_required_real_value

  integer function get_required_integer_value(key_name)
    character(len=*), intent(in) :: key_name
    integer :: entry_index, io_status
    entry_index = find_configuration_entry_index(key_name)
    read(configuration_value_list(entry_index),*, iostat=io_status) get_required_integer_value
    if (io_status /= 0) then
      write(*,'(a,a)') 'ERROR: cannot read an integer value for key ', trim(key_name)
      stop 1
    end if
  end function get_required_integer_value

  ! ----------------------------------------------------------------------------------------
  ! Split a comma separated line into fields (no quoting support, none needed here).
  ! ----------------------------------------------------------------------------------------
  subroutine split_comma_separated_line(line_text, field_list, number_of_fields)
    character(len=*), intent(in)  :: line_text
    character(len=*), intent(out) :: field_list(:)
    integer,          intent(out) :: number_of_fields
    integer :: field_start_position, comma_position, line_length_characters

    line_length_characters = len_trim(line_text)
    number_of_fields = 0
    field_start_position = 1
    do
      comma_position = index(line_text(field_start_position:line_length_characters), ',')
      number_of_fields = number_of_fields + 1
      if (comma_position == 0) then
        field_list(number_of_fields) = adjustl(line_text(field_start_position:line_length_characters))
        exit
      end if
      field_list(number_of_fields) = &
           adjustl(line_text(field_start_position:field_start_position+comma_position-2))
      field_start_position = field_start_position + comma_position
    end do
  end subroutine split_comma_separated_line

  ! ----------------------------------------------------------------------------------------
  ! Read the header row and the first data row of the forcing CSV file.
  ! ----------------------------------------------------------------------------------------
  subroutine read_forcing_csv_first_data_row(forcing_file_path)
    character(len=*), intent(in) :: forcing_file_path
    character(len=maximum_line_length_characters) :: header_line_text, data_line_text
    integer :: file_unit, io_status, number_of_data_fields

    open(newunit=file_unit, file=trim(forcing_file_path), status='old', action='read', &
         iostat=io_status)
    if (io_status /= 0) then
      write(*,'(a,a)') 'ERROR: cannot open forcing file ', trim(forcing_file_path)
      stop 1
    end if
    read(file_unit,'(a)', iostat=io_status) header_line_text
    if (io_status /= 0) stop 'ERROR: cannot read forcing header'
    read(file_unit,'(a)', iostat=io_status) data_line_text
    if (io_status /= 0) stop 'ERROR: cannot read forcing data row'
    close(file_unit)

    call split_comma_separated_line(header_line_text, forcing_column_name_list, &
                                    number_of_forcing_columns)
    call split_comma_separated_line(data_line_text, forcing_column_value_list, &
                                    number_of_data_fields)
    if (number_of_data_fields /= number_of_forcing_columns) then
      stop 'ERROR: forcing data row and header have different numbers of columns'
    end if
  end subroutine read_forcing_csv_first_data_row

  real function get_forcing_column_value(column_name)
    character(len=*), intent(in) :: column_name
    integer :: column_index, io_status
    do column_index = 1, number_of_forcing_columns
      if (trim(forcing_column_name_list(column_index)) == trim(column_name)) then
        read(forcing_column_value_list(column_index),*, iostat=io_status) get_forcing_column_value
        if (io_status /= 0) then
          write(*,'(a,a)') 'ERROR: cannot read forcing value for column ', trim(column_name)
          stop 1
        end if
        return
      end if
    end do
    write(*,'(a,a)') 'ERROR: forcing column not found: ', trim(column_name)
    stop 1
  end function get_forcing_column_value

  ! first column is the time stamp, e.g. 2016-07-01 12:00:00
  function get_forcing_time_string() result(time_string)
    character(len=128) :: time_string
    time_string = forcing_column_value_list(1)
  end function get_forcing_time_string

end module harness_text_utilities
