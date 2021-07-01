!!  Take an initial model from a Lagrangian code and put it onto
!!  a uniform grid and make sure that it is happy with the EOS in
!!  our code.  The output is a .hse file that can be read directly
!!  by Maestro.
!!
!!  The model is placed into HSE by the following differencing:
!!
!!   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
!!
!!  We take the temperature structure directly from the original
!!  initial model.  For the composition, we interpolate Ye and X
!!  from the initial model.  If we are in an NSE region, then we
!!  take Ye and call the NSE table to get X.  If we are not in NSE,
!!  then we use X to compute Ye.  Since we build with USE_NSE = TRUE,
!!  the EOS always requires Ye and abar through the aux data.
!!  We adjust the density and pressure according to HSE using the EOS.
!!
!!***

module model_params
  !! Set up the model parameters

  use amrex_fort_module, only: rt => amrex_real
  use extern_probin_module
  use network

  ! define convenient indices for the scalars
  integer, parameter :: nvar = 5 + nspec
  integer, parameter :: idens = 1, &
                        itemp = 2, &
                        ipres = 3, &
                        ientr = 4, &
                        iyef  = 5, &    ! this is ye -- we add "f" for file to not clash with the network
                        ispec = 6

  real (kind=rt), parameter :: TOL = 1.d-10

  real (kind=rt), allocatable :: xzn_hse(:), xznl(:), xznr(:)
  real (kind=rt), allocatable :: M_enclosed(:)
  real (kind=rt), allocatable :: model_mesa_hse(:,:)
  real (kind=rt), allocatable :: model_hybrid_hse(:,:)
  real (kind=rt), allocatable :: entropy_want(:)

  integer, parameter :: MAX_ITER = 250, AD_ITER = 2500

  integer, parameter :: MAX_VARNAME_LENGTH=80

  real (kind=rt), parameter :: smallx = 1.d-10

  real (kind=rt) :: delx

  real (kind=rt), save :: low_density_cutoff =1.d-7

  ! temp_fluff_cutoff is the density below which we hold the temperature
  ! constant for the MESA model

  ! MAESTRO
  ! temp_fluff_cutoff = 1.d-4

  ! CASTRO
  real (kind=rt), save :: temp_fluff_cutoff = 2.d-7
  real (kind=rt), save :: temp_fluff = 1.d5

end module model_params


subroutine set_aux(eos_state)

  use eos_type_module
  use eos_module
  use network
  use nse_module
  use nse_check_module
  use burn_type_module

  type(eos_t), intent(inout) :: eos_state

  type(burn_t) :: burn_state

  real(rt) :: abar_pass, dq_pass, dyedt_pass
  integer :: nse_check

  call eos_to_burn(eos_state, burn_state)

  call in_nse(burn_state, nse_check)

  if (nse_check == 1) then

     ! we are in NSE, so leave ye alone, but get abar

     call nse_interp(eos_state % T, eos_state % rho, eos_state % aux(iye), &
                     abar_pass, dq_pass, dyedt_pass, eos_state % xn(:))

     eos_state % aux(iabar) = abar_pass
  else
     eos_state % aux(iye) = sum(eos_state % xn * zion * aion_inv)
     eos_state % aux(iabar) = 1.0d0 / sum(eos_state % xn * aion_inv)
  end if

end subroutine set_aux

module initial_model_module


contains

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
    write (50,1001) "# num of variables = ", 4 + nspec
    write (50,1002) "# density"
    write (50,1002) "# temperature"
    write (50,1002) "# pressure"
    write (50,1002) "# ye"

    ! write (50,1002) "# conv_vel"

    do n = 1, nspec
       write (50,1003) "# ", spec_names(n)
    enddo

    do i = 1, nx

       write (50,1000) xzn_hse(i), model_state(i,idens), model_state(i,itemp), &
            model_state(i,ipres), model_state(i,iyef), &
            (model_state(i,ispec-1+n), n=1,nspec)
    enddo

    close (unit=50)


  end subroutine write_model

  subroutine read_file(filename, base_state, base_r, npts_model)

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module
    use amrex_error_module
    use network
    use fundamental_constants_module
    use model_params
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

    open(99,file=model_file)

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

          case ("ye")
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


    open (unit=50, file="model.orig", status="unknown")

    write (50,*) "# initial model as read in"

    do i = 1, npts_model_file
       write (50,1000) base_r(i), (base_state(i,j), j = 1, nvar)
    enddo

    close (50)

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
    use model_params
    use interpolate_module
    use nse_module
    use nse_check_module
    use burn_type_module

    implicit none

    integer :: i, k, n

    real (kind=rt) :: dens_zone, temp_zone, pres_zone, entropy, temp_min, temp_max

    real (kind=rt) :: A, B, dAdT, dAdrho, dBdT, dBdrho
    real (kind=rt) :: dpd, dpt, dsd, dst

    real (kind=rt) :: central_density

    real (kind=rt) :: p_want, drho, dtemp, p_old, eos_p_old,  ledoux

    real (kind=rt) :: g_zone, err_min, err_max

    real (kind=rt) :: max_hse_error, dpdr, rhog

    real (kind=rt) :: abar_pass, dq_pass, dyedt_pass
    integer :: nse_check

    integer :: iter, iter_dens, status

    integer :: igood

    logical :: converged_hse, converged_central_density, fluff, isentropic

    integer :: index_hse_fluff = 1

    real (kind=rt), dimension(nspec) :: xn
    real (kind=rt) :: ye

    integer :: npts_model

    real(kind=rt), allocatable :: base_state(:,:), base_r(:)

    integer :: ipos
    character (len=500) :: outfile, model_name

    real(kind=rt) :: summ

    integer :: ibegin
    real(kind=rt) :: grav_ener, M_center

    type (eos_t) :: eos_state
    type (burn_t) :: burn_state

    !===========================================================================
    ! Create a 1-d uniform grid that is identical to the mesh that we are
    ! mapping onto, and then we want to force it into HSE on that mesh.
    !===========================================================================

    ! allocate storage
    allocate(xzn_hse(nx))
    allocate(xznl(nx))
    allocate(xznr(nx))
    allocate(model_mesa_hse(nx,nvar))
    allocate(model_hybrid_hse(nx,nvar))
    allocate(M_enclosed(nx))
    allocate(entropy_want(nx))

    ! compute the coordinates of the new gridded function

    delx = (xmax - xmin) / dble(nx)

    do i = 1, nx
       xznl(i) = xmin + (dble(i) - 1.0d0)*delx
       xznr(i) = xmin + (dble(i))*delx
       xzn_hse(i) = 0.5d0*(xznl(i) + xznr(i))
    enddo

    !===========================================================================
    ! read in the MESA model
    !===========================================================================

    call read_file(model_file, base_state, base_r, npts_model)

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

    open (unit=30, file="model.uniform", status="unknown")

1000 format (1x, 30(g26.16, 1x))

    write (30,*) "# initial model just after putting onto a uniform grid"

    do i = 1, nx

       write (30,1000) xzn_hse(i), model_mesa_hse(i,idens), model_mesa_hse(i,itemp), &
            model_mesa_hse(i,ipres), model_mesa_hse(i,iyef), (model_mesa_hse(i,ispec-1+n), n=1,nspec)

    enddo

    close (unit=30)



    !===========================================================================
    ! reset the composition if we are in NSE
    !===========================================================================

    do i = 1, nx

       ! we need to fill a burn_t in order to check for NSE

       eos_state % rho = model_mesa_hse(i,idens)
       eos_state % T = model_mesa_hse(i,idens)
       eos_state % aux(:) = 0.0
       eos_state % aux(iye) = model_mesa_hse(i,iyef)

       do n = 1, nspec
          eos_state % xn(n) = model_mesa_hse(i,ispec-1+n)
       end do


       call eos_to_burn(eos_state, burn_state)

       call in_nse(burn_state, nse_check)

       if (nse_check == 1) then

          ! we are in NSE, so let's call the table to get the correct mass fractions

          call nse_interp(eos_state % T, eos_state % rho, eos_state % aux(iye), &
                          abar_pass, dq_pass, dyedt_pass, eos_state % xn(:))

       else

          ! we are not in NSE, so let's compute a consistent ye from the composition
          eos_state % aux(iye) = sum(eos_state % xn * zion * aion_inv)

       end if

       ! copy the composition variables back

       model_mesa_hse(i,ispec:ispec-1+nspec) = eos_state % xn(:)
       model_mesa_hse(i,iyef) = eos_state % aux(iye)

    end do


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

       eos_state%aux(iye) = model_mesa_hse(ibegin,iyef)

       call set_aux(eos_state)

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

          ye = model_mesa_hse(i,iyef)

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

             eos_state%aux(iye) = ye

             call set_aux(eos_state)

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
             print *, 'integrate down'
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

          eos_state%aux(iye) = ye

          call set_aux(eos_state)

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

       ye = model_mesa_hse(i,iyef)

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

             eos_state%aux(iye) = ye

             call set_aux(eos_state)

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

       eos_state%aux(iye) = ye

       call set_aux(eos_state)

       ! if we were in NSE, then this updated eos_State % xn(:), so copy that over
       model_mesa_hse(i,ispec:ispec-1+nspec) = eos_state % xn(:)

       print *, "output: ", eos_state%T, eos_state%rho, eos_state%aux(iye)

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
