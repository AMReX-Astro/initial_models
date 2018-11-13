!! Read a MESA file and return an array containing the required data.
module mesa_reader

contains

  !> @brief reads a mesa file and puts data in an array
  !!
  !! Given the MESA input file `filename`, this shall read in the file, store all
  !! variables listed in `var_names_model` in the array `base_state` (along with all species
  !! required for the network) and store the radius of each point in `base_r`.
  !! All variables are stored in cgs units.
  !!
  !! @param[in]     filename        name of MESA file
  !! @param[inout]  base_state      allocatable real array where the data shall be stored
  !! @param[inout]  base_r          allocatable real array where radius shall be stored
  !! @param[in]     var_names_model list of variables we want to read from file and store in
  !!                                base_state, NOT including species
  !! @param[in]     nvars_model     Number of variables in var_names_model
  !! @param[out]    npts_file       Number of data points in file that shall be stored in
  !!                                base_state
  subroutine read_mesa(filename, base_state, base_r, var_names_model, nvars_model, npts_file)

    use bl_types
    use amrex_constants_module
    use amrex_error_module
    use network
    use fundamental_constants_module

    implicit none

    integer, parameter :: MAX_VARNAME_LENGTH=80

    character (len=100), intent(in) :: filename
    real(kind=dp_t), allocatable, intent(inout) :: base_state(:,:)
    real(kind=dp_t), allocatable, intent(inout) :: base_r(:)
    character (len=MAX_VARNAME_LENGTH), intent(in) :: var_names_model(:)
    integer, intent(in) :: nvars_model
    integer, intent(out) :: npts_file

    integer :: nparams_file, nvars_file, status, ipos
    character (len=5000) :: header_line
    real(kind=dp_t), allocatable :: vars_stored(:), params_stored(:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: paramnames_stored(:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: var_names_file(:)
    integer, allocatable :: var_indices_file(:)
    logical :: found
    integer :: i, j, k, n

1000 format (1x, 30(g26.16, 1x))

    open(99,file=filename)

    ! going to first do a pass through the file to count the line numbers
    npts_file = 0
    do
       read (99,*,iostat=status)
       if (status > 0) then
          write(*,*) "Something went wrong :("
       else if (status < 0) then
          exit
       else
          npts_file = npts_file + 1
       endif
    enddo

    ! subtract first two header rows and go back to start of file
    npts_file = npts_file - 6

    print *, npts_file, '    points found in the initial model file'

    rewind(99)

    ! the first line gives the column numbers of the parameters. From this we can get the number of parameters
    read (99, '(a2000)') header_line
    ! find last space in line
    ipos = index(trim(header_line), ' ', back=.true.)
    header_line = trim(adjustl(header_line(ipos:)))
    read (header_line, *) nparams_file

    print *, nparams_file, ' parameters found in the initial model file'

    allocate (params_stored(nparams_file))
    allocate (paramnames_stored(nparams_file))

    ! now read in the names of the variables
    read (99, '(a2000)') header_line
    ipos = 1
    do i = 1, nparams_file
       header_line = trim(adjustl(header_line(ipos:)))
       ipos = index(header_line, ' ') + 1
       paramnames_stored(i) = trim(adjustl(header_line(:ipos)))
    enddo

    ! and finally read and store the parameters
    read(99,*) (params_stored(j), j = 1, nparams_file)

    ! space then read in variable numbers
    read(99, *)

    ! get number of variables
    read (99, '(a5000)') header_line
    ! find last space in line
    ipos = index(trim(header_line), ' ', back=.true.)
    header_line = trim(adjustl(header_line(ipos:)))
    read (header_line, *) nvars_file

    print *, nvars_file, ' variables found in the initial model file'

    allocate (vars_stored(nvars_file))
    allocate (var_names_file(nvars_file))
    allocate (var_indices_file(nvars_file))

    var_indices_file(:) = -1

    ! now read in the names of the variables
    read (99, '(a5000)') header_line
    ipos = 1
    do i = 1, nvars_file
       header_line = trim(adjustl(header_line(ipos:)))
       ipos = index(header_line, ' ') + 1
       var_names_file(i) = trim(adjustl(header_line(:ipos)))

       ! create map of file indices to model indices
       do j = 1, nvars_model
          if (var_names_file(i) == var_names_model(j)) then
             var_indices_file(i) = j
             exit
          endif
       enddo

       ! map species as well
       if (var_indices_file(i) == -1) then
          k = network_species_index(var_names_file(i))
          if (k > 0) then
             var_indices_file(i) = nvars_model + k
          endif
       endif
    enddo

    ! check that all var_names_model have been found
    do i = 1, nvars_model
       found = .false.
       do j = 1, nvars_file
          if (var_names_model(i) == trim(var_names_file(j))) then
             found = .true.
             exit
          endif
       enddo
       if (.not. found) then
          call amrex_error("file does not contain variable: ", var_names_model(i))
       endif
    enddo

    ! allocate storage for the model data
    allocate (base_state(npts_file, nvars_model+nspec))
    allocate (base_r(npts_file))

    do i = 1, npts_file
       read(99, *) (vars_stored(j), j = 1, nvars_file)

       ! need to reverse the inputs file here

       n = npts_file - i + 1

       base_state(n,:) = ZERO

       do j = 1, nvars_file
          if (var_names_file(j) == "logR") then

             base_r(n) = R_solar*10**vars_stored(j)

          else if (var_indices_file(j) .ge. 0) then
             base_state(n, var_indices_file(j)) = vars_stored(j)
          endif

       enddo

    enddo

    print *, base_state(1,:)

    open (unit=50, file="model.orig", status="unknown")

    write (50,*) "# initial model as read in"

    do i = 1, npts_file
       write (50,1000) base_r(i), (base_state(i,j), j = 1, nvars_model)
    enddo

    close (50)


  end subroutine read_mesa

end module mesa_reader
