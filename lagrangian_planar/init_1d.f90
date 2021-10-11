!!  The model is placed into HSE by the following differencing:
!!
!!   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
!!
!!  This will be iterated over in tandem with the EOS call,
!!  P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1)
!!

subroutine init_1d() bind(C, name="init_1d")

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use amrex_error_module
  use eos_module, only: eos, eos_init
  use eos_type_module, only: eos_t, eos_input_rt
  use extern_probin_module
  use network, only : nspec, network_species_index, spec_names, network_init
  use fundamental_constants_module, only: Gconst
  use model_module
  use interpolate_module
  use cfmt_module

  implicit none

  integer :: i, n

  character(len=128) :: params_file

  real (kind=rt), DIMENSION(nspec) :: xn

  real (kind=rt), allocatable :: xzn_hse(:), xznl_hse(:), xznr_hse(:)
  real (kind=rt), allocatable :: model_hse(:,:)

  real :: A

  integer :: lun1, lun2

  real (kind=rt) :: dCoord

  real (kind=rt) :: dens_zone, temp_zone, pres_zone, entropy
  real (kind=rt) :: dpd

  real (kind=rt) :: p_want, drho, dtemp, delx

  real (kind=rt) :: g_zone

  real (kind=rt), parameter :: TOL = 1.e-10

  integer, parameter :: MAX_ITER = 250

  integer :: iter

  logical :: converged_hse, fluff

  real (kind=rt) :: max_T

  integer :: index_base

  character (len=256) :: outfile, outfile2
  character (len=8) :: num
  character (len=32) :: dxstr
  character (len=32) :: num_to_unitstring

  real (kind=rt) :: max_hse_error, dpdr, rhog

  type (eos_t) :: eos_state

  integer :: npts_model
  real(kind=rt), allocatable :: model_state(:,:), model_r(:)


  ! start by reading in the Lagrangian initial model
  call read_file(model_file, model_state, model_r, npts_model)

  ! apply the shift
  do i = 1, npts_model
     model_r(i) = model_r(i) - model_shift
  enddo


!-----------------------------------------------------------------------------
! Create a 1-d uniform grid that is identical to the mesh that we are
! mapping onto, and then we want to force it into HSE on that mesh.
!-----------------------------------------------------------------------------

  ! allocate storage
  allocate(xzn_hse(nx))
  allocate(xznl_hse(nx))
  allocate(xznr_hse(nx))
  allocate(model_hse(nx,nvar))


  ! compute the coordinates of the new gridded function
  dCoord = (xmax - xmin) / dble(nx)

  do i = 1, nx
     xznl_hse(i) = xmin + (dble(i) - ONE)*dCoord
     xzn_hse(i)  = xmin + (dble(i) - HALF)*dCoord
     xznr_hse(i) = xmin + (dble(i))*dCoord
  enddo


!-----------------------------------------------------------------------------
! put the model onto our new uniform grid
!-----------------------------------------------------------------------------

  fluff = .false.

  do i = 1, nx
     do n = 1, nvar
        if (n == itemp) then
           model_hse(i,n) = max(temp_cutoff, &
                                interpolate(xzn_hse(i), npts_model, model_r, model_state(:,n)))
        else
           model_hse(i,n) = interpolate(xzn_hse(i), npts_model, model_r, model_state(:,n))
        endif
     enddo

     ! make it all thermodynamically consistent
     eos_state%rho = model_hse(i,idens)
     eos_state%T = model_hse(i,itemp)
     eos_state%xn(:) = model_hse(i,ispec:ispec-1+nspec)

     call eos(eos_input_rt, eos_state)

     model_hse(i,ipres) = eos_state%p
  enddo


  ! find the index to integrate from by looking for the peak temperature
  index_base = -1
  max_T = -1.0_rt

  do i = 1, nx
     if (model_hse(i,itemp) > max_T) then
        index_base = i
        max_T = model_hse(i,itemp)
     endif
  enddo

  if (index_base == -1) then
     call amrex_error('ERROR: invalid base_height')
  endif

  print *, 'index_base = ', index_base

  ! make the base thermodynamics consistent for this base point -- that is
  ! what we will integrate from!
  eos_state%rho = model_hse(index_base,idens)
  eos_state%T = model_hse(index_base,itemp)
  eos_state%xn(:) = model_hse(index_base,ispec:ispec-1+nspec)

  call eos(eos_input_rt, eos_state)

  model_hse(index_base,ipres) = eos_state%p


!-----------------------------------------------------------------------------
! HSE + entropy solve
!-----------------------------------------------------------------------------

  ! the HSE state will be done respecting the interpolated temperature
  ! from the initial model.  When the temperature drops below T_lo,
  ! we floor it.

  !---------------------------------------------------------------------------
  ! integrate up
  !---------------------------------------------------------------------------
  do i = index_base+1, nx

     delx = xzn_hse(i) - xzn_hse(i-1)

     ! compute the gravitation acceleration at the lower edge
     if (do_invsq_grav == 1) then
        g_zone = -Gconst*M_enclosed/xznl_hse(i)**2
     else
        g_zone = g_const
     endif

     ! we've already set initial guesses for density, temperature, and
     ! composition
     dens_zone = model_hse(i,idens)
     temp_zone = max(temp_cutoff, model_hse(i,itemp))
     xn(:) = model_hse(i,ispec:ispec-1+nspec)


     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.

     do iter = 1, MAX_ITER

        ! what pressure does HSE say we want?
        p_want = model_hse(i-1,ipres) + &
             delx*0.5*(dens_zone + model_hse(i-1,idens))*g_zone

        ! (t, rho) -> (p)
        eos_state%T   = temp_zone
        eos_state%rho = dens_zone
        eos_state%xn(:) = xn(:)

        call eos(eos_input_rt, eos_state)

        entropy = eos_state%s
        pres_zone = eos_state%p

        dpd = eos_state%dpdr

        drho = (p_want - pres_zone)/(dpd - 0.5*delx*g_zone)

        dens_zone = max(0.9*dens_zone, &
             min(dens_zone + drho, 1.1*dens_zone))

        if (abs(drho) < TOL*dens_zone) then
           converged_hse = .TRUE.
           exit
        endif

        if (dens_zone < low_density_cutoff) then
           dens_zone = low_density_cutoff
           temp_zone = temp_cutoff
           converged_hse = .TRUE.
           exit
        endif

     enddo


     if (.NOT. converged_hse) then
        print *, 'Error zone', i, ' did not converge in init_1d'
        print *, 'integrate up'
        print *, dens_zone, temp_zone
        print *, p_want
        print *, drho, dtemp
        call amrex_error('Error: HSE non-convergence')
     endif


     ! call the EOS one more time for this zone and then go on to the next
     ! (t, rho) -> (p)
     eos_state%T     = temp_zone
     eos_state%rho   = dens_zone
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     pres_zone = eos_state%p

     ! update the thermodynamics in this zone
     model_hse(i,idens) = dens_zone
     model_hse(i,itemp) = temp_zone
     model_hse(i,ipres) = pres_zone

     ! to make this process converge faster, set the density in the
     ! next zone to the density in this zone
     ! model_hse(i+1,idens) = dens_zone

  enddo


  !---------------------------------------------------------------------------
  ! integrate down -- using the temperature profile defined above
  !---------------------------------------------------------------------------
  do i = index_base-1, 1, -1

     delx = xzn_hse(i+1) - xzn_hse(i)

     ! compute the gravitation acceleration at the upper edge
     if (do_invsq_grav == 1) then
        g_zone = -Gconst*M_enclosed/xznr_hse(i)**2
     else
        g_zone = g_const
     endif

     ! we already set the temperature and composition profiles
     temp_zone = max(temp_cutoff, model_hse(i,itemp))
     xn(:) = model_hse(i,ispec:ispec-1+nspec)

     ! use our previous initial guess for density
     dens_zone = model_hse(i,idens)


     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.

     do iter = 1, MAX_ITER

        ! get the pressure we want from the HSE equation, just the
        ! zone below the current.  Note, we are using an average of
        ! the density of the two zones as an approximation of the
        ! interface value -- this means that we need to iterate for
        ! find the density and pressure that are consistent

        ! HSE differencing
        p_want = model_hse(i+1,ipres) - &
             delx*0.5*(dens_zone + model_hse(i+1,idens))*g_zone


        ! we will take the temperature already defined in model_hse
        ! so we only need to zero:
        !   A = p_want - p(rho)

        ! (t, rho) -> (p)
        eos_state%T     = temp_zone
        eos_state%rho   = dens_zone
        eos_state%xn(:) = xn(:)

        call eos(eos_input_rt, eos_state)

        pres_zone = eos_state%p

        dpd = eos_state%dpdr

        A = p_want - pres_zone

        drho = A/(dpd + 0.5*delx*g_zone)

        dens_zone = max(0.9_rt*dens_zone, &
             min(dens_zone + drho, 1.1_rt*dens_zone))


        if (abs(drho) < TOL*dens_zone) then
           converged_hse = .TRUE.
           exit
        endif


     enddo

     if (.NOT. converged_hse) then

        print *, 'Error zone', i, ' did not converge in init_1d'
        print *, 'integrate down'
        print *, dens_zone, temp_zone
        print *, p_want
        print *, drho
        call amrex_error('Error: HSE non-convergence')

     endif


     ! call the EOS one more time for this zone and then go on to the next
     ! (t, rho) -> (p)
     eos_state%T     = temp_zone
     eos_state%rho   = dens_zone
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     pres_zone = eos_state%p

     ! update the thermodynamics in this zone
     model_hse(i,idens) = dens_zone
     model_hse(i,itemp) = temp_zone
     model_hse(i,ipres) = pres_zone

  enddo


1000 format (1x, 100(g26.16, 1x))
1001 format (a, i5)
1002 format (a)
1003 format (a,a)

  write(num,'(i8)') nx

  dxstr = num_to_unitstring(dCoord)

  outfile = trim(model_prefix) // ".hse" // ".dx_" // trim(adjustl(dxstr))
  outfile2 = trim(outfile) // ".extras"

  open (newunit=lun1, file=outfile, status="unknown")
  open (newunit=lun2, file=outfile2, status="unknown")

  write (lun1,1001) "# npts = ", nx
  write (lun1,1001) "# num of variables = ", 3 + nspec
  write (lun1,1002) "# density"
  write (lun1,1002) "# temperature"
  write (lun1,1002) "# pressure"

  do n = 1, nspec
     write (lun1, 1003) "# ", spec_names(n)
  enddo

  do i = 1, nx
     write (lun1,1000) cfmt(xzn_hse(i)), cfmt(model_hse(i,idens)), &
                       cfmt(model_hse(i,itemp)), cfmt(model_hse(i,ipres)), &
                      (cfmt(model_hse(i,ispec-1+n)), n=1,nspec)
  enddo


  write (lun2,1001) "# npts = ", nx
  write (lun2,1001) "# num of variables = ", 2
  write (lun2,1002) "# entropy"
  write (lun2,1002) "# c_s"

  ! test: bulk EOS call -- Maestro will do this once we are mapped, so make
  ! sure that we are in HSE with updated thermodynamics
  do i = 1, nx
     eos_state%rho = model_hse(i,idens)
     eos_state%T = model_hse(i,itemp)
     eos_state%xn(:) = model_hse(i,ispec:ispec-1+nspec)

     call eos(eos_input_rt, eos_state)

     model_hse(i,ipres) = eos_state%p

     write (lun2,1000) cfmt(xzn_hse(i)), cfmt(eos_state%s), cfmt(eos_state%cs)
  enddo

  ! compute the maximum HSE error
  max_hse_error = -1.d30

  do i = 2, nx-1

     ! compute the gravitation acceleration at the lower edge
     if (do_invsq_grav == 1) then
        g_zone = -Gconst*M_enclosed/xznl_hse(i)**2
     else
        g_zone = g_const
     endif

     dpdr = (model_hse(i,ipres) - model_hse(i-1,ipres))/delx
     rhog = HALF*(model_hse(i,idens) + model_hse(i-1,idens))*g_zone

     if (dpdr /= ZERO .and. model_hse(i+1,idens) > low_density_cutoff) then
        max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))
     endif

  enddo

  print *, 'maximum HSE error = ', max_hse_error
  print *, ' '

  close (unit=lun1)
  close (unit=lun2)

end subroutine init_1d


function num_to_unitstring(value)

  use amrex_fort_module, only: rt => amrex_real
  implicit none

  real (kind=rt) :: value
  character (len=32) :: num_to_unitstring
  character (len=16) :: temp

  if (value > 1.d5) then

     ! work in km
     write(temp,'(f6.3)') value/1.d5
     num_to_unitstring = trim(temp) // "km"
  else

     ! work in cm
     if (value > 1.d3) then
        write(temp,'(f8.3)') value
        num_to_unitstring = trim(temp) // "cm"

     else
        write(temp,'(f6.3)') value
        num_to_unitstring = trim(temp) // "cm"
     endif

  endif

  return
end function num_to_unitstring
