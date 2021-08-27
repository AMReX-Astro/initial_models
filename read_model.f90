module model_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use amrex_error_module
  use network
  use fundamental_constants_module

  implicit none

  integer, parameter :: MAX_VARNAME_LENGTH=80

  ! define convenient indices for the scalars
  integer, parameter :: nvar = 5 + nspec
  integer, parameter :: idens = 1, &
                        itemp = 2, &
                        ipres = 3, &
                        ientr = 4, &
                        iyef  = 5, &    ! this is ye -- we add "f" for file to not clash with the network
                        ispec = 6

contains

  subroutine read_file(filename, base_state, base_r, npts_model)

    implicit none


    character (len=100), intent(in) :: filename
    real(kind=rt), allocatable, intent(inout) :: base_state(:,:)
    real(kind=rt), allocatable, intent(inout) :: base_r(:)
    integer, intent(out) :: npts_model

    real(kind=rt), allocatable :: vars_stored(:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: varnames_stored(:)

    integer :: nvars_model_file, npts_model_file
    integer :: status, ipos
    character (len=5000) :: header_line
    logical :: found
    integer :: i, j, k, n

1000 format (1x, 30(g26.16, 1x))

    open(99,file=filename)

    ! the first line has the number of points in the model
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) npts_model_file

    print *, npts_model_file, '    points found in the initial model file'

    npts_model = npts_model_file

    ! now read in the number of variables
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) nvars_model_file

    print *, nvars_model_file, ' variables found in the initial model file'

    allocate (vars_stored(nvars_model_file))
    allocate (varnames_stored(nvars_model_file))

    ! now read in the names of the variables
    do i = 1, nvars_model_file
       read (99, '(a256)') header_line
       ipos = index(header_line, '#') + 1
       varnames_stored(i) = trim(adjustl(header_line(ipos:)))
    enddo

    ! allocate storage for the model data
    allocate (base_state(npts_model_file, nvar))
    allocate (base_r(npts_model_file))

    do i = 1, npts_model_file
       read(99,*) base_r(i), (vars_stored(j), j = 1, nvars_model_file)

       base_state(i,:) = ZERO

       do j = 1, nvars_model_file

          found = .false.

          select case (trim(varnames_stored(j)))

          case ("density")
             base_state(i,idens) = vars_stored(j)
             found = .true.

          case ("temperature")
             base_state(i,itemp) = vars_stored(j)
             found = .true.

          case ("pressure")
             base_state(i,ipres) = vars_stored(j)
             found = .true.

          case ("ye", "Ye")
             base_state(i,iyef) = vars_stored(j)
             found = .true.

          case default

             ! check if they are species
             n = network_species_index(trim(varnames_stored(j)))
             if (n > 0) then
                base_state(i,ispec-1+n) = vars_stored(j)
                found = .true.
             endif

          end select

          if (.NOT. found .and. i == 1) then
             print *, 'ERROR: variable not found: ', varnames_stored(j)
          endif

       enddo

    enddo

  end subroutine read_file

end module model_module
