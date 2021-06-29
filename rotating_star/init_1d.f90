!!  Take an initial model from a Lagrangian code and put it onto
!!  a uniform grid and make sure that it is happy with the EOS in
!!  our code.  The output is a .hse file that can be read directly
!!  by Maestro.
!!
!!  The model is placed into HSE by the following differencing:
!!
!!   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
!!
!!
!!  We take the temperature structure directly from the original
!!  initial model.  We adjust the density and pressure according to
!!  HSE using the EOS.
!!
!!***

module model_params
  !! Set up the model parameters

  use amrex_fort_module, only: rt => amrex_real

  use network

  integer, parameter :: nx = 6400

  ! define convenient indices for the scalars
  integer, parameter :: idens = 1, &
       itemp = 2, &
       ipres = 3, &
       ientr = 4, &
       ienuc = 5, & 
       ispec = 6
    !    isndspd = 5, &
    !    imass = 6, &
    !    igradr = 7, &
    !    igrav = 8, &
    !    ibrunt = 9, &
    !    iconv_vel = 10, &
    !    ispec = 11

  real (kind=rt), parameter :: TOL = 1.d-10

  real (kind=rt), allocatable :: xzn_hse(:), xznl(:), xznr(:)
  real (kind=rt), allocatable :: M_enclosed(:)
  real (kind=rt), allocatable :: model_mesa_hse(:,:)
  real (kind=rt), allocatable :: model_isentropic_hse(:,:)
  real (kind=rt), allocatable :: model_hybrid_hse(:,:)
  real (kind=rt), allocatable :: model_conservative(:,:)
  real (kind=rt), allocatable :: model_convective(:,:)
  real (kind=rt), allocatable :: model_ad_excess(:,:)
  real (kind=rt), allocatable :: model_ledoux(:,:)
  real (kind=rt), allocatable :: entropy_want(:)
  real (kind=rt), allocatable :: model_ener(:)

  integer, parameter :: MAX_ITER = 250, AD_ITER = 2500

  integer, parameter :: MAX_VARNAME_LENGTH=80

  real(kind=rt), parameter :: anelastic_cutoff = 9.d4  ! this is for diagnostics only -- not used in the HSEing
  real (kind=rt), parameter :: smallx = 1.d-10

!   character (len=100), parameter :: model_file = "18m_500_s_rot_b_eq_1.dat"
!   logical, parameter :: mesa = .true.

  character (len=100), parameter :: model_file = "15m_500_sec.txt"
  logical, parameter :: mesa = .false.

  real (kind=rt), parameter :: xmin = 0.d0, xmax = 1.75d10 !1.732050808d10
  real (kind=rt), parameter :: delx = (xmax - xmin) / dble(nx)

  real (kind=rt), save :: low_density_cutoff =1.d-7

  ! temp_fluff_cutoff is the density below which we hold the temperature
  ! constant for the MESA model

  ! MAESTRO
  ! temp_fluff_cutoff = 1.d-4

  ! CASTRO
  real (kind=rt), save :: temp_fluff_cutoff = 2.d-7
  real (kind=rt), save :: temp_fluff = 1.d5

  character(len=MAX_VARNAME_LENGTH), allocatable :: varnames_stored(:)

  integer, parameter :: nvar = ispec - 1 + nspec

end module model_params


module initial_model_module


contains


  function interpolate(r, npts, model_r, model_var)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    ! given the array of model coordinates (model_r), and variable (model_var),
    ! find the value of model_var at point r using linear interpolation.
    ! Eventually, we can do something fancier here.

    real(kind=rt) :: interpolate
    real(kind=rt), intent(in) :: r
    integer :: npts
    real(kind=rt), dimension(npts) :: model_r, model_var

    real(kind=rt) :: slope
    real(kind=rt) :: minvar, maxvar

    integer :: i, id

    ! find the location in the coordinate array where we want to interpolate
    do i = 1, npts
       if (model_r(i) >= r) exit
    enddo

    id = i

    if (id == 1) then

       slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
       interpolate = slope*(r - model_r(id)) + model_var(id)

       ! safety check to make sure interpolate lies within the bounding points
       !minvar = min(model_var(id+1),model_var(id))
       !maxvar = max(model_var(id+1),model_var(id))
       !interpolate = max(interpolate,minvar)
       !interpolate = min(interpolate,maxvar)

    else

       slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
       interpolate = slope*(r - model_r(id)) + model_var(id)

       ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_var(id),model_var(id-1))
       maxvar = max(model_var(id),model_var(id-1))
       interpolate = max(interpolate,minvar)
       interpolate = min(interpolate,maxvar)

    endif

    return

  end function interpolate

  ! function conservative_interpolate(r, npts, model_r, model_var)
  subroutine conservative_interpolate(interpolated, r, npts_model, model_r, model_var, dx, status, iloc)

    !     given the array of model coordinates (model_r), and variable (model_var),
    !     find the value of model_var at point r (var_r) using linear interpolation.
    !     Eventually, we can do something fancier here.

    use amrex_fort_module, only: rt => amrex_real
    use amrex_error_module

    implicit none

    real(kind=rt)        , intent(in   ) :: r, dx
    integer         , intent(in   ) :: npts_model
    real(kind=rt)        , intent(in   ) :: model_r(npts_model), model_var(npts_model)
    integer, intent(in), optional   :: iloc
    real(kind=rt)        , intent(out  ) :: interpolated
    integer         , intent(out  ) :: status

    ! Local variables
    integer                         :: max_iter = 5
    integer                         :: i, n, n_boxes
    real(kind=rt)                        :: rel_error = 1.d0
    real(kind=rt)                        :: delta = 1.d-4
    real(kind=rt)                        :: summ, rm
    ! real(kind=rt)               :: centered_interpolate


    interpolated = centered_interpolate(r, npts_model, model_r, model_var, iloc)

    status = 0

    do n = 1, max_iter
       if (rel_error <= delta) exit

       summ = 0.0d0
       n_boxes = 2**n

       do i = 1, n_boxes
          rm = r - 0.5 * dx + dx * (float(i-1) + 0.5d0) / float(n_boxes)
          summ = summ + 1.0d0 / float(n_boxes) * centered_interpolate(rm, npts_model, model_r, model_var)
       enddo

       rel_error = abs(summ - interpolated) / abs(interpolated)

       interpolated = summ

    enddo

    if (rel_error >  delta) status = 1

  end subroutine conservative_interpolate


  function centered_interpolate(r, npts_model, model_r, model_var, iloc) result(interpolated)

    !     given the array of model coordinates (model_r), and variable (model_var),
    !     find the value of model_var at point r (var_r) using linear interpolation.
    !     Eventually, we can do something fancier here.

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(kind=rt)                        :: interpolated
    real(kind=rt)        , intent(in   ) :: r
    integer         , intent(in   ) :: npts_model
    real(kind=rt)        , intent(in   ) :: model_r(npts_model), model_var(npts_model)
    integer, intent(in), optional   :: iloc

    ! Local variables
    integer                         :: id
    real(kind=rt)                        :: slope,minvar,maxvar

    !     find the location in the coordinate array where we want to interpolate
    if (present(iloc)) then
       id = iloc
    else
       call locate_sub(r, npts_model, model_r, id)
    end if

    if (id .eq. 1) then

       slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
       interpolated = slope*(r - model_r(id)) + model_var(id)

       ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_var(id+1),model_var(id))
       maxvar = max(model_var(id+1),model_var(id))
       interpolated = max(interpolated,minvar)
       interpolated = min(interpolated,maxvar)

    else if (id .eq. npts_model) then

       slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
       interpolated = slope*(r - model_r(id)) + model_var(id)

       ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_var(id),model_var(id-1))
       maxvar = max(model_var(id),model_var(id-1))
       interpolated = max(interpolated,minvar)
       interpolated = min(interpolated,maxvar)

    else

       slope = 0.5d0 *( (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id)) + &
            (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1)) )
       interpolated = slope*(r - model_r(id)) + model_var(id)

       ! ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_var(id+1),model_var(id), model_var(id-1))
       maxvar = max(model_var(id+1),model_var(id), model_var(id-1))
       interpolated = max(interpolated,minvar)
       interpolated = min(interpolated,maxvar)


    endif

  end function centered_interpolate

  subroutine locate_sub(x, n, xs, loc)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: n
    real(rt), intent(in   ) :: x, xs(n)
    integer,  intent(  out) :: loc

    integer :: ilo, ihi, imid

    !$gpu

    if (x .le. xs(1)) then
       loc = 1
    else if (x .gt. xs(n-1)) then
       loc = n
    else

       ilo = 1
       ihi = n-1

       do while (ilo+1 .ne. ihi)
          imid = (ilo+ihi)/2
          if (x .le. xs(imid)) then
             ihi = imid
          else
             ilo = imid
          end if
       end do

       loc = ihi

    end if

  end subroutine locate_sub

  subroutine write_model(model_name, model_state)

    !! Write data stored in `model_state` array to file

    use model_params

    implicit none

    character (len=100), intent(in) :: model_name
    ! integer, intent(in) :: nx, nvar
    real(kind=rt), intent(in) :: model_state(nx, nvar)

    character(len=100) :: outfile
    integer :: ipos, i, n
    character (len=8) num

1000 format (1x, 30(g26.16, 1x))
1001 format(a, i5)
1002 format(a)
1003 format(a,a)

    ipos = index(model_file, '.dat')
    if (ipos .eq. 0) then 
        ipos = index(model_file, '.txt')
    endif
    outfile = model_file(1:ipos-1) // '.' // trim(adjustl(model_name))

    write(num,'(i8)') nx

    outfile = trim(outfile) // '.' // trim(adjustl(num))

    print *, 'writing ', trim(adjustl(model_name)), ' model to ', trim(outfile)

    open (unit=50, file=outfile, status="unknown")

    write (50,1001) "# npts = ", nx
    write (50,1001) "# num of variables = ", 3 + nspec
    write (50,1002) "# density"
    write (50,1002) "# temperature"
    write (50,1002) "# pressure"
    ! write (50,1002) "# conv_vel"

    do n = 1, nspec
       write (50,1003) "# ", spec_names(n)
    enddo

    do i = 1, nx

        if (mesa) then

            ! write (50,1000) xzn_hse(i), model_state(i,idens), model_state(i,itemp), &
            !         model_state(i,ipres), model_state(i,iconv_vel), &
            !         (model_state(i,ispec-1+n), n=1,nspec)

        else

            write (50,1000) xzn_hse(i), model_state(i,idens), model_state(i,itemp), &
                 model_state(i,ipres), &
                 (model_state(i,ispec-1+n), n=1,nspec)
     

        endif

    enddo

    close (unit=50)


  end subroutine write_model

  subroutine read_file(filename, base_state, base_r, var_names_model, nvars_model, npts_file) 

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module
    use amrex_error_module
    use network
    use fundamental_constants_module

    implicit none

    integer, parameter :: MAX_VARNAME_LENGTH=80

    character (len=100), intent(in) :: filename
    real(kind=rt), allocatable, intent(inout) :: base_state(:,:)
    real(kind=rt), allocatable, intent(inout) :: base_r(:)
    character (len=MAX_VARNAME_LENGTH), intent(in) :: var_names_model(:)
    integer, intent(in) :: nvars_model
    integer, intent(out) :: npts_file

    integer :: nparams_file, nvars_file, status, ipos
    character (len=5000) :: header_line
    real(kind=rt), allocatable :: vars_stored(:), params_stored(:)
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

    rewind(99)

     ! the second line gives the number of variables
    read(99, *)
    read(99, '(a2000)') header_line
    ! find last space in line
    ipos = index(trim(header_line), ' ', back=.true.)
    header_line = trim(adjustl(header_line(ipos:)))
    read (header_line, *) nvars_file

    print *, nvars_file, ' variables found in the initial model file'

    allocate (vars_stored(0:nvars_file))
    allocate (var_names_file(nvars_file))
    allocate (var_indices_file(nvars_file))

    ! subtract header rows 
    npts_file = npts_file - nvars_file - 2

    print *, npts_file, '    points found in the initial model file'

    var_indices_file(:) = -1

    ! now read in the names of the variables
    ipos = 1
    do i = 1, nvars_file
       read (99, '(a5000)') header_line
    !    header_line = trim(adjustl(header_line(ipos:)))
    !    ipos = index(header_line, ' ') + 1
       var_names_file(i) = trim(adjustl(header_line))

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

    ! allocate storage for the model data
    allocate (base_state(npts_file, nvars_model+nspec))
    allocate (base_r(npts_file))

    do i = 1, npts_file
       read(99, *) (vars_stored(j), j = 0, nvars_file)

       base_state(i,:) = ZERO

       base_r(i) = vars_stored(0)

       do j = 1, nvars_file
          if (var_indices_file(j) .ge. 0) then
             base_state(i, var_indices_file(j)) = vars_stored(j)
          endif

       enddo

    enddo

  end subroutine read_file

  
  subroutine init_1d() bind(C, name="init_1d")

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module
    use amrex_error_module
    use eos_module, only: eos, eos_init
    use eos_type_module, only : eos_t, eos_input_rt, eos_input_re, eos_input_rp
    use network
    use fundamental_constants_module
    use extern_probin_module, only: use_eos_coulomb
    use mesa_reader
    use model_params

    implicit none

    integer :: i, k, n

    real (kind=rt) :: dens_zone, temp_zone, pres_zone, entropy, temp_min, temp_max

    real (kind=rt) :: A, B, dAdT, dAdrho, dBdT, dBdrho
    real (kind=rt) :: dpd, dpt, dsd, dst

    real (kind=rt) :: central_density

    real (kind=rt) :: p_want, drho, dtemp, chiT, chirho, p_old, eos_p_old,  ledoux

    real (kind=rt) :: g_zone, ad_error, ad_tol, adiabatic_excess, err_min, err_max

    real (kind=rt) :: max_hse_error, dpdr, rhog

    integer :: iter, iter_dens, status

    integer :: igood

    logical :: converged_hse, converged_central_density, fluff, isentropic, print_n

    real (kind=rt) :: max_temp

    integer :: index_hse_fluff = 1

    real (kind=rt), dimension(nspec) :: xn
    integer :: npts_model

    real(kind=rt), allocatable :: base_state(:,:), base_r(:), base_ener(:)
    real(kind=rt), allocatable :: base_ad(:), model_ad(:)
    real(kind=rt), allocatable :: base_led(:), model_led(:)

    integer :: ipos
    character (len=500) :: outfile, model_name

    real(kind=rt) :: summ

    integer :: ibegin
    integer :: i_isentropic
    real(kind=rt) :: M_enclosed_anel
    real(kind=rt) :: grav_ener, M_center
    real(kind=rt) :: eint_hybrid

    type (eos_t) :: eos_state

    !===========================================================================
    ! Create a 1-d uniform grid that is identical to the mesh that we are
    ! mapping onto, and then we want to force it into HSE on that mesh.
    !===========================================================================

    ! allocate storage
    allocate(xzn_hse(nx))
    allocate(xznl(nx))
    allocate(xznr(nx))
    allocate(model_mesa_hse(nx,nvar))
    allocate(model_isentropic_hse(nx,nvar))
    allocate(model_hybrid_hse(nx,nvar))
    allocate(model_convective(nx,nvar))
    allocate(M_enclosed(nx))
    allocate(entropy_want(nx))

    ! compute the coordinates of the new gridded function

    do i = 1, nx
       xznl(i) = xmin + (dble(i) - 1.0d0)*delx
       xznr(i) = xmin + (dble(i))*delx
       xzn_hse(i) = 0.5d0*(xznl(i) + xznr(i))
    enddo

    !===========================================================================
    ! read in the MESA model
    !===========================================================================

    if (mesa) then 

       allocate(varnames_stored(nvar-nspec))
       varnames_stored(idens) = "rho"
       varnames_stored(itemp) = "temperature"
       varnames_stored(ipres) = "pressure"
       varnames_stored(ientr) = "entropy"
       ! varnames_stored(isndspd) = "csound"
       ! varnames_stored(imass) = "mass"
       ! varnames_stored(igradr) = "gradT"
       ! varnames_stored(igrav) = "grav"
       ! varnames_stored(ibrunt) = "brunt_N2"
       ! varnames_stored(iconv_vel) = "conv_vel"

       call read_mesa(model_file, base_state, base_r, varnames_stored, nvar-nspec, npts_model)

    else 

       allocate(varnames_stored(nvar-nspec))
       varnames_stored(idens) = "dens"
       varnames_stored(itemp) = "temp"
       varnames_stored(ipres) = "pres"
       varnames_stored(ientr) = "entr"
       varnames_stored(ienuc) = "enuc"

       call read_file(model_file, base_state, base_r, varnames_stored, nvar-nspec, npts_model)

    endif

    !===========================================================================
    ! put the model onto our new uniform grid
    !===========================================================================

    igood = -1

    do i = 1, nx

       do n = 1, nvar

          if (xzn_hse(i) < base_r(npts_model)) then

             model_mesa_hse(i,n) = interpolate(xzn_hse(i), npts_model, &
                  base_r, base_state(:,n))

             igood = i
          else
             model_mesa_hse(i,n) = base_state(npts_model,n)
          endif

       enddo

       ! make sure that the species (mass fractions) summ to 1
       summ = 0.0d0
       do n = ispec,ispec-1+nspec
          model_mesa_hse(i,n) = max(model_mesa_hse(i,n),smallx)
          summ = summ + model_mesa_hse(i,n)
       enddo

       do n = ispec,ispec-1+nspec
          model_mesa_hse(i,n) = model_mesa_hse(i,n)/summ
       enddo

    enddo

    print *, "model_hse(1,idens) =", model_mesa_hse(1,idens)

    open (unit=30, file="model.uniform", status="unknown")

1000 format (1x, 30(g26.16, 1x))

    write (30,*) "# initial model just after putting onto a uniform grid"

    do i = 1, nx

       if (mesa) then 

          !  write (30,1000) xzn_hse(i), model_mesa_hse(i,idens), model_mesa_hse(i,itemp), &
          !       model_mesa_hse(i,ipres), model_mesa_hse(i,iconv_vel), (model_mesa_hse(i,ispec-1+n), n=1,nspec)

       else 

          write (30,1000) xzn_hse(i), model_mesa_hse(i,idens), model_mesa_hse(i,itemp), &
               model_mesa_hse(i,ipres), (model_mesa_hse(i,ispec-1+n), n=1,nspec)

       endif

    enddo

    close (unit=30)

    !===========================================================================
    ! iterate to find the central density
    !===========================================================================

    ! because the MESA model likely begins at a larger radius than our first
    ! HSE model zone, simple interpolation will not do a good job.  We want to
    ! integrate in from the zone that best matches the first MESA model zone,
    ! assumming HSE and constant entropy.

    ! find the zone in the uniformly gridded model that corresponds to the
    ! first zone of the original model
    ibegin = -1

    do i = 1, nx
       if (xzn_hse(i) >= base_r(1)) then
          ibegin = i
          exit
       endif
    enddo

    ! store the central density.  We will iterate until the central density
    ! converges
    central_density = model_mesa_hse(1,idens)
    print *, 'interpolated central density = ', central_density

    do iter_dens = 1, max_iter

       ! compute the enclosed mass
       M_enclosed(1) = FOUR3RD*M_PI*delx**3*model_mesa_hse(1,idens)

       do i = 2, ibegin
          M_enclosed(i) = M_enclosed(i-1) + &
               FOUR3RD*M_PI*(xznr(i) - xznl(i)) * &
               (xznr(i)**2 +xznl(i)*xznr(i) + xznl(i)**2)*model_mesa_hse(i,idens)
       enddo

       ! now start at ibegin and integrate inward
       eos_state%T     = model_mesa_hse(ibegin,itemp)
       eos_state%rho   = model_mesa_hse(ibegin,idens)
       eos_state%xn(:) = model_mesa_hse(ibegin,ispec:nvar)

       call eos(eos_input_rt, eos_state)

       model_mesa_hse(ibegin,ipres) = eos_state%p
       !  model_mesa_hse(ibegin,ientr) = eos_state%s
       !  model_mesa_hse(ibegin,isndspd) = &
       !   sqrt(eos_state%gam1*eos_state%p/eos_state%rho)

       entropy_want(:) = eos_state%s

       do i = ibegin-1, 1, -1

          ! as the initial guess for the temperature and density, use
          ! the previous zone
          dens_zone = model_mesa_hse(i+1,idens)
          temp_zone = model_mesa_hse(i+1,itemp)
          xn(:) = model_mesa_hse(i,ispec:nvar)

          ! compute the gravitational acceleration on the interface between zones
          ! i and i+1
          g_zone = -Gconst*M_enclosed(i)/(xznr(i)*xznr(i))

          !-----------------------------------------------------------------------
          ! iteration loop
          !-----------------------------------------------------------------------

          ! start off the Newton loop by saying that the zone has not converged
          converged_hse = .FALSE.

          do iter = 1, MAX_ITER

             p_want = model_mesa_hse(i+1,ipres) - &
                  delx*0.5d0*(dens_zone + model_mesa_hse(i+1,idens))*g_zone

             ! now we have two functions to zero:
             !   A = p_want - p(rho,T)
             !   B = entropy_want - s(rho,T)
             ! We use a two dimensional Taylor expansion and find the deltas
             ! for both density and temperature

             ! (t, rho) -> (p, s)

             eos_state%T     = temp_zone
             eos_state%rho   = dens_zone
             eos_state%xn(:) = xn(:)

             call eos(eos_input_rt, eos_state)

             entropy = eos_state%s
             pres_zone = eos_state%p

             dpt = eos_state%dpdt
             dpd = eos_state%dpdr
             dst = eos_state%dsdt
             dsd = eos_state%dsdr

             A = p_want - pres_zone
             B = entropy_want(i) - entropy

             dAdT = -dpt
             dAdrho = -0.5d0*delx*g_zone - dpd
             dBdT = -dst
             dBdrho = -dsd

             dtemp = (B - (dBdrho/dAdrho)*A)/ &
                  ((dBdrho/dAdrho)*dAdT - dBdT)

             drho = -(A + dAdT*dtemp)/dAdrho

             dens_zone = max(0.9d0*dens_zone, &
                  min(dens_zone + drho, 1.1d0*dens_zone))

             temp_zone = max(0.9d0*temp_zone, &
                  min(temp_zone + dtemp, 1.1d0*temp_zone))

             ! if (A < TOL .and. B < ETOL) then
             if (abs(drho) < TOL*dens_zone .and. abs(dtemp) < TOL*temp_zone) then
                converged_hse = .TRUE.
                exit
             endif

          enddo

          if (.NOT. converged_hse) then

             print *, 'Error zone', i, ' did not converge in init_1d'
             print *, 'integrate up'
             print *, 'dens_zone, temp_zone = ', dens_zone, temp_zone
             print *, "p_want = ", p_want
             print *, "drho = ", drho
             call amrex_error('Error: HSE non-convergence')

          endif

          ! call the EOS one more time for this zone and then go on to the next
          ! (t, rho) -> (p, s)

          eos_state%T     = temp_zone
          eos_state%rho   = dens_zone
          eos_state%xn(:) = xn(:)

          call eos(eos_input_rt, eos_state)

          pres_zone = eos_state%p

          dpd = eos_state%dpdr

          ! update the thermodynamics in this zone
          model_mesa_hse(i,idens) = dens_zone
          model_mesa_hse(i,itemp) = temp_zone
          model_mesa_hse(i,ipres) = pres_zone
          model_mesa_hse(i,ientr) = eos_state%s
          ! model_mesa_hse(i,isndspd) = &
          !      sqrt(eos_state%gam1*eos_state%p/eos_state%rho)

       enddo

       if (abs(model_mesa_hse(1,idens) - central_density) < TOL*central_density) then
          converged_central_density = .true.
          exit
       endif

       central_density = model_mesa_hse(1,idens)

    enddo

    if (.NOT. converged_central_density) then
       print *, 'ERROR: central density iterations did not converge'
       call amrex_error('ERROR: non-convergence')
    endif

    print *, 'converged central density = ', model_mesa_hse(1,idens)
    print *, ' '

    !===========================================================================
    ! compute the full HSE model using our new central density and temperature,
    ! and the temperature structure as dictated by the MESA model.
    !===========================================================================

    print *, 'putting MESA model into HSE on our grid...'

    ! compute the enclosed mass
    M_enclosed(1) = FOUR3RD*M_PI*delx**3*model_mesa_hse(1,idens)

    fluff = .FALSE.

    do i = 2, nx

       ! use previous zone as initial guess for T, rho
       dens_zone = model_mesa_hse(i-1,idens)
       temp_zone = model_mesa_hse(i-1,itemp)

       xn(:) = model_mesa_hse(i,ispec:nvar)

       ! compute the gravitational acceleration on the interface between zones
       ! i-1 and i
       g_zone = -Gconst*M_enclosed(i-1)/(xznr(i-1)*xznr(i-1))

       !-----------------------------------------------------------------------
       ! iteration loop
       !-----------------------------------------------------------------------

       converged_hse = .FALSE.

       if (.not. fluff) then

          do iter = 1, MAX_ITER

             ! HSE differencing
             p_want = model_mesa_hse(i-1,ipres) + &
                  delx*0.5d0*(dens_zone + model_mesa_hse(i-1,idens))*g_zone

             temp_zone = model_mesa_hse(i,itemp)

             if (model_mesa_hse(i-1,idens) .lt. temp_fluff_cutoff) then
                temp_zone = temp_fluff
             end if

             ! (t, rho) -> (p)
             eos_state%T     = temp_zone
             eos_state%rho   = dens_zone
             eos_state%xn(:) = xn(:)

             call eos(eos_input_rt, eos_state)

             pres_zone = eos_state%p

             dpd = eos_state%dpdr
             drho = (p_want - pres_zone)/(dpd - 0.5d0*delx*g_zone)

             dens_zone = max(0.9d0*dens_zone, &
                  min(dens_zone + drho, 1.1d0*dens_zone))

             if (abs(drho) < TOL*dens_zone) then
                converged_hse = .TRUE.
                exit
             endif

             if (dens_zone < low_density_cutoff) then
                dens_zone = low_density_cutoff
                temp_zone = temp_fluff
                converged_hse = .TRUE.
                fluff = .TRUE.
                index_hse_fluff = i
                exit

             endif

          enddo

          if (.NOT. converged_hse) then

             print *, 'Error zone', i, ' did not converge in init_1d'
             print *, 'integrate up'
             print *, 'dens_zone, temp_zone = ', dens_zone, temp_zone
             print *, "p_want = ", p_want
             print *, "drho = ", drho
             call amrex_error('Error: HSE non-convergence')

          endif

          if (temp_zone < temp_fluff) then
             temp_zone = temp_fluff
          endif

       else
          dens_zone = low_density_cutoff
          temp_zone = temp_fluff
       endif

       ! call the EOS one more time for this zone and then go on to the next
       ! (t, rho) -> (p)

       eos_state%T     = temp_zone
       eos_state%rho   = dens_zone
       eos_state%xn(:) = xn(:)

       call eos(eos_input_rt, eos_state)

       pres_zone = eos_state%p

       ! update the thermodynamics in this zone
       model_mesa_hse(i,idens) = dens_zone
       model_mesa_hse(i,itemp) = temp_zone
       model_mesa_hse(i,ipres) = pres_zone
       model_mesa_hse(i,ientr) = eos_state%s
       !  model_mesa_hse(i,isndspd) = &
       !       sqrt(eos_state%gam1*eos_state%p/eos_state%rho)

       M_enclosed(i) = M_enclosed(i-1) + &
            FOUR3RD*M_PI*(xznr(i) - xznl(i)) * &
            (xznr(i)**2 +xznl(i)*xznr(i) + xznl(i)**2)*model_mesa_hse(i,idens)

    enddo



    !===========================================================================
    ! output
    !===========================================================================

    !---------------------------------------------------------------------------
    ! MESA model
    !---------------------------------------------------------------------------

    model_name = 'hse'

    call write_model(model_name, model_mesa_hse)

    ! compute the enclosed mass
    M_enclosed(1) = FOUR3RD*M_PI*delx**3*model_mesa_hse(1,idens)

    do i = 2, nx
       M_enclosed(i) = M_enclosed(i-1) + &
            FOUR3RD*M_PI*(xznr(i) - xznl(i)) * &
            (xznr(i)**2 +xznl(i)*xznr(i) + xznl(i)**2)*model_mesa_hse(i,idens)
    enddo

    print *, 'summ mass = ', real(M_enclosed(nx)), ' g = ', real(M_enclosed(nx)/M_solar), ' solar masses'

    ! compute the maximum HSE error
    max_hse_error = -1.d30

    do i = 2, nx-1
       g_zone = -Gconst*M_enclosed(i-1)/xznr(i-1)**2
       dpdr = (model_mesa_hse(i,ipres) - model_mesa_hse(i-1,ipres))/delx
       rhog = HALF*(model_mesa_hse(i,idens) + model_mesa_hse(i-1,idens))*g_zone

       if (dpdr /= ZERO .and. model_mesa_hse(i+1,idens) > low_density_cutoff) then
          max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))
       endif
    enddo

    print *, 'maximum HSE error = ', max_hse_error
    print *, ' '

  end subroutine init_1d

end module initial_model_module
