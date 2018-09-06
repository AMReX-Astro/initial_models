!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! These subroutines perform HSE integration and enforce thermodynamic
! consistency on the model
!
! NOTE: We are solving a differential equation of the form y' = f(x,y).  The
! numerical methods involve approximating y_next from y_curr using a
! reconstruction based on several values of f (f_i, f_i-1, etc).  In our case
! f = density * gravitational acceleration, and g is constant (i.e.,
! independent of density), so we can factor g out and reconstruct based on
! values of density.  However our density is not correct; it is an
! approximation based on interpolating the input model, and the density will be
! corrected using a Newton-Raphson iteration that links the density and
! pressure through the EoS.  Thus for density of the current cell (rho_i), we
! are iterating on a guess and can give any value we like for the initial
! guess.  Since we are integrating from the base and correcting as we go, and
! the cells are "small", then the density of the previous zone is a better
! guess than the interpolated density from the initial model.  That is why we
! use rho_i-1 (which has already been corrected) instead of rho_i (which has
! not yet been corrected) for the initial guess of dens_zone.  Then the
! Newton-Raphson iteration corrects this guess, using dens_zone as the value of
! rho_i in the reconstruction.
!
! The arguments are
! - g_type   : in    : the type of gravity
! - out_file : in    : name of file to print final MAESTRO-readable data
! - err_file : in    : name of file to print hydrostatic equilibrium error data

subroutine integrate_HSE(g_type, mass, p_type, temp_fluff, outfile)

   use eos_type_module
   use eos_module
   use network
   use init_1d_variables
   use init_1d_grids
   use bl_types
   use bl_error_module

   implicit none

   !===========================================================================
   ! Declare variables

   ! Arguments ................................................................
   integer,                      intent(in   ) :: g_type
   real(kind=dp_t),              intent(in   ) :: mass
   integer,                      intent(in   ) :: p_type
   real(kind=dp_t),              intent(in   ) :: temp_fluff
   character(len=256),           intent(in   ) :: outfile

   ! Locals ...................................................................
   real(kind=dp_t)              :: g_zone
   real(kind=dp_t)              :: g_const

   ! Temporaries ..............................................................
   integer :: i, n ! loop indices
   real(kind=dp_t) :: dx, r_edge

   ! Functions ................................................................
   real(kind=dp_t) :: m_shell

   type(eos_t) :: eos_state

   !===========================================================================
   ! Enforce thermodynamic consistency at the base
   
   eos_state % T = Ustate(ibase, itemp)
   eos_state % rho  = Ustate(ibase, idens)
   eos_state % xn(:) = Ustate(ibase, ispec:Nvars)

   call eos(eos_input_rt, eos_state)

   Ustate(ibase, ipres) = eos_state % p
   Ustate(ibase, ientr) = eos_state % s

   !===========================================================================
   ! Construct the gravitational constant

   ! construct gravitational constant
   select case(g_type)
   case(GRAVITY_CONST)
      g_const = G * mass * Uradius(ibase)**(-2)
   case(GRAVITY_INVSQ)
      g_const = G * mass
   case(GRAVITY_MENCL)
      ! For enclosed-mass gravity, we assume that the mass supplied on the
      ! command line is equal to the mass of the entire core (but exclusding
      ! the envlope).  We assume that the base zone is the lowest zone of the
      ! envelope.  We define the enclosed mass for a given zone (m_enc(i)) to
      ! be the mass of all shells below that zone:
      !    m_enc(i) = sum(j = 1..i-1)(m_shell(j)),
      ! where m_shell(j) is the mass in shell j:
      !    m_shell(j) = (4/3) * pi * rho * ((r+)**3 - (r-)**3).
      ! Therefor the upward integration would add the shell mass of the (i-1)th
      ! zone to get the current enclosed mass, while the downward integration
      ! would subtract the shell mass of the (i+1)th zone to get the current
      ! enclosed mass.
      g_const = mass
   case default
      call bl_error('ERROR: invalid gravity mode in integrate_HSE subroutine')
   end select

   !===========================================================================
   ! Integrate radially outward from the base

   call integration_loop(ibase+1, NrU, g_type, g_const, p_type, temp_fluff)

   !===========================================================================
   ! Integrate radially inward from the base

   if (g_type == GRAVITY_MENCL) then
      ! Since the command-line mass is the total enclosed mass at r_(ibase-1/2)
      ! and the integration_loop routine assumes the enclosed mass needs to be
      ! updated prior to computing gravity, we must add in the shell mass of
      ! the base zone, as the integration_loop routine will subtract that off
      ! prior to the first iteration.
      g_const = mass + m_shell(ibase)
   end if

   call integration_loop(ibase-1, 1, g_type, g_const, p_type, temp_fluff)

   !===========================================================================
   ! Print

   ! Print MAESTRO-readable file
   open(unit=99, file=outfile)
   write(99,'(a,i5)') "# npts = ", NrU
   write(99,'(a,i5)') "# num of variables = ", Nvars
   write(99,'(a)')    "# density"
   write(99,'(a)')    "# temperature"
   write(99,'(a)')    "# pressure"
   write(99,'(a)')    "# entropy"
   do n = 1, nspec
      write(99,'(a,a)') "# ", spec_names(n)
   end do
   do i = 1, NrU
      write(99,'(1x,30(g17.10,1x))') Uradius(i), &
                                     (Ustate(i, n), n=1,Nvars)
   end do
   close(unit=99)

end subroutine integrate_HSE

!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! Perform HSE integration

subroutine integration_loop(istart, iend, g_type, g_const, p_type, temp_fluff)

   use eos_module
   use eos_type_module
   use init_1d_variables
   use init_1d_grids
   use bl_constants_module
   use bl_types
   use bl_error_module

   implicit none

   !===========================================================================
   ! Declare variables

   ! Arguments ................................................................
   integer,                      intent(in   ) :: istart, iend
   integer,                      intent(in   ) :: g_type
   real(kind=dp_t),              intent(in   ) :: g_const
   integer,                      intent(in   ) :: p_type
   real(kind=dp_t),              intent(in   ) :: temp_fluff

   ! Locals ...................................................................
   integer :: istep
   logical :: fluff
   real(kind=dp_t) :: dx, r_edge
   real(kind=dp_t) :: dens_zone, temp_zone, pres_zone, s_zone, g_zone
   real(kind=dp_t) :: core_entropy
   real(kind=dp_t) :: Mencl

   ! Temporaries ..............................................................
   integer :: i ! loop index

   type(eos_t) :: eos_state

   ! Functions ................................................................
   real(kind=dp_t) :: m_shell

   !===========================================================================
   ! Set up

   ! Set loop direction
   if (istart <= iend) then
      istep = 1
   else
      istep = -1
   end if

   ! Assume we start in the non-fluff region
   fluff = .false.

   ! Set initial enclosed mass
   if (g_type == GRAVITY_MENCL) then
      Mencl = g_const
   end if

   !===========================================================================
   ! Integration loop

   do i = istart, iend, istep

      ! Compute radial step between cell centers, radius of edge between cells
      dx = abs(Uradius(i) - Uradius(i-istep))
      r_edge = HALF * (Uradius(i) + Uradius(i-istep))

      ! Compute gravity
      select case(g_type)
      case(GRAVITY_CONST)
         g_zone = g_const
      case(GRAVITY_INVSQ)
         g_zone = g_const / r_edge**2
      case(GRAVITY_MENCL)
         ! Add previous shell mass on the upward sweep (m_shell(i-1))
         ! Subtract previous shell mass on the downward sweep (m_shell(i+1))
         Mencl = Mencl + istep*m_shell(i-istep)
         if (Mencl .le. 0) then
            call bl_error("ERROR: enclosed mass non-positive")
         end if
         g_zone  = G * Mencl / r_edge**2
      case default
         call bl_error('ERROR: invalid gravity mode in integration_loop &
                            &subroutine')
      end select

      ! set composition data
      ! since xn_eos is in a module, it will still be set correctly in NR loop
      eos_state % xn(:) = Ustate(i, ispec:Nvars)

      ! compute quantities for current zone
      if (fluff) then
         ! if in the fluff region, use fluff values
         dens_zone = dens_fluff_cutoff
         temp_zone = temp_fluff
      else
         ! if not in the fluff, Newton-Raphson loop to correct initial guesses

         ! Always use the density of the last zone as the best guess for the
         ! density of the current zone (because changes from zone to zone are
         ! "small", while deviations from the interpolated model may be "large"
         dens_zone   = Ustate(i-istep, idens)
         ! If we are deep enough into the core that we are using a specific
         ! profile, the isothermal variant always uses the temperature at icore
         ! while the isentropic always uses the temperature of the last zone
         ! as a guess for the current zone.  If we are not using a specific
         ! profile, use the temperature of the current zone to follow the input
         ! model.
         if (i < icore) then
            select case(p_type)
            case(CORE_PROFILE_SECOND_HSE)
               temp_zone = Ustate(i, itemp)
            case(CORE_PROFILE_ISOTHERMAL)
               temp_zone = Ustate(icore, itemp)
            case(CORE_PROFILE_ISENTROPIC)
               temp_zone = Ustate(i-istep, itemp)
            case default
               call bl_error("ERROR: invalid profile type in integration &
                                  &loop subroutine")
            end select
         else
            temp_zone = Ustate(i, itemp)
         end if
         ! either core_entropy is set to the correct value, or the entropy is
         ! not used (isothermal Newton-Raphson), so we can always set the zone
         ! entropy to core_entropy
         s_zone = core_entropy

         if ((p_type == CORE_PROFILE_SECOND_HSE) &
             .or. (p_type == CORE_PROFILE_ISOTHERMAL) .or. (i >= icore)) then
            call HSE_NR_loop_isothermal(dens_zone, Ustate(i-istep, idens),   &
               temp_zone, pres_zone, Ustate(i-istep, ipres), s_zone, g_zone, &
               dx, i, istep, fluff, temp_fluff)
         else ! apply isentropic profile
            call HSE_NR_loop_isentropic(dens_zone, Ustate(i-istep, idens),   &
               temp_zone, pres_zone, Ustate(i-istep, ipres), s_zone, g_zone, &
               dx, i, istep, fluff, temp_fluff)
         end if
      end if

      ! Enforce thermodynamic consistency
      eos_state % T = temp_zone
      eos_state % rho = dens_zone
      call eos(eos_input_rt, eos_state)
      pres_zone = eos_state % p
      s_zone    = eos_state % s

      ! Save entropy at core profile edge once that zone has been corrected.
      ! We will only need the core entropy on the inward integration, which is
      ! also the pass that computes the core entropy since icore is computed in
      ! a way that forces icore < ibase.  Thus we can compute it locally
      ! because that is the only time it is needed.
      if (i == icore) then
         core_entropy = s_zone
      end if

      ! Save results of HSE integration and thermodynamic check
      Ustate(i, idens) = dens_zone
      Ustate(i, itemp) = temp_zone
      Ustate(i, ipres) = pres_zone
      Ustate(i, ientr) = s_zone

   end do
   
end subroutine integration_loop

!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! isothermal Newton-Raphson loop

subroutine HSE_NR_loop_isothermal(dens_curr, dens_prev, temp_curr, pres_curr, &
                                  pres_prev, s_curr, grav, dx, izone, istep,  &
                                  fluff, temp_fluff)

   use eos_module
   use eos_type_module
   use init_1d_variables
   use bl_types
   use bl_constants_module
   use bl_error_module

   implicit none

   !===========================================================================
   ! Declare variables

   ! Arguments ................................................................
   real(kind=dp_t), intent(inout) :: dens_curr
   real(kind=dp_t), intent(in   ) :: dens_prev
   real(kind=dp_t), intent(inout) :: temp_curr
   real(kind=dp_t), intent(  out) :: pres_curr
   real(kind=dp_t), intent(in   ) :: pres_prev
   real(kind=dp_t), intent(inout) :: s_curr
   real(kind=dp_t), intent(in   ) :: grav, dx
   integer,         intent(in   ) :: izone
   integer,         intent(in   ) :: istep
   logical,         intent(  out) :: fluff
   real(kind=dp_t), intent(in   ) :: temp_fluff

   ! Locals ...................................................................
   real(kind=dp_t) :: pres_hse, pres_eos
   real(kind=dp_t) :: dp_hse, dp_eos
   real(kind=dp_t) :: drho
   logical         :: converged

   type(eos_t) :: eos_state

   ! Temporaries ..............................................................
   integer :: i ! loop index

   !===========================================================================
   ! Newton-Raphson loop

   converged = .false.

   do i = 1, MAX_NR_ITER
      
      ! Compute pressure from HSE
      pres_hse = pres_prev + dble(istep)*HALF*(dens_curr+dens_prev)*dx*grav

      ! Compute pressure from EOS
      eos_state % T = temp_curr
      eos_state % rho = dens_curr
      call eos(eos_input_rt, eos_state)
      pres_eos = eos_state % p

      ! Compute change in density based on Taylor expansions of EOS and HSE
      dp_hse = HALF*dx*grav
      dp_eos = eos_state % dpdr
      drho = (pres_hse - pres_eos) / (dp_eos - dp_hse)

      ! Restrict change in density to prevent runaway
      dens_curr = max(0.9*dens_curr, min(dens_curr + drho, 1.1*dens_curr))

      ! Check if convergence achieved
      if (abs(drho / dens_curr) < TOLERANCE) then
         converged = .true.
         exit
      end if

      ! Check if solution wandered into fluff
      if (dens_curr < dens_fluff_cutoff) then
         dens_curr = dens_fluff_cutoff
         temp_curr = temp_fluff
         converged = .true.
         fluff = .true.
         exit
      end if

   end do

   ! Verify the solution converged
   if (.not. converged) then
      write(*,*) "Error zone", izone, " did not converge in init_1d"
      if (izone > ibase) then
         write(*,*) "integrate up"
      else
         write(*,*) "integrate down"
      end if
      write(*,*) "temp = ", temp_curr
      write(*,*) "dens = ", dens_curr
      write(*,*) "drho = ", drho
      write(*,*) "pres (EoS) = ", pres_eos
      write(*,*) "pres (HSE) = ", pres_hse
      call bl_error("ERROR: HSE non-convergence")
   end if

   ! Use the fluff temperature as a floor
   if (temp_curr < temp_fluff) then
      temp_curr = temp_fluff
   end if

   pres_curr = pres_eos

end subroutine HSE_NR_loop_isothermal

!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! isentropic Newton-Raphson loop

subroutine HSE_NR_loop_isentropic(dens_curr, dens_prev, temp_curr, pres_curr, &
                                  pres_prev, s_curr, grav, dx, izone, istep,  &
                                  fluff, temp_fluff)

   use eos_module
   use eos_type_module
   use init_1d_variables
   use bl_types
   use bl_constants_module
   use bl_error_module

   implicit none

   !===========================================================================
   ! Declare variables

   ! Arguments ................................................................
   real(kind=dp_t), intent(inout) :: dens_curr
   real(kind=dp_t), intent(in   ) :: dens_prev
   real(kind=dp_t), intent(inout) :: temp_curr
   real(kind=dp_t), intent(  out) :: pres_curr
   real(kind=dp_t), intent(in   ) :: pres_prev
   real(kind=dp_t), intent(inout) :: s_curr
   real(kind=dp_t), intent(in   ) :: grav, dx
   integer,         intent(in   ) :: izone
   integer,         intent(in   ) :: istep
   logical,         intent(  out) :: fluff
   real(kind=dp_t), intent(in   ) :: temp_fluff

   ! Locals ...................................................................
   real(kind=dp_t) :: pres_hse, pres_eos
   real(kind=dp_t) :: entr_iso, entr_eos
   real(kind=dp_t) :: dpdr, dpdt, dsdr, dsdt, dp_hse
   real(kind=dp_t) :: determinant
   real(kind=dp_t) :: drho, dtmp
   logical         :: converged

   ! Temporaries ..............................................................
   integer         :: i ! loop index
   real(kind=dp_t) :: junk

   type(eos_t) :: eos_state

   !===========================================================================
   ! Newton-Raphson loop

   converged = .false.

   entr_iso = s_curr

   do i = 1, MAX_NR_ITER
      
      ! Compute pressure from HSE
      pres_hse = pres_prev + dble(istep)*HALF*(dens_curr+dens_prev)*dx*grav

      ! Compute pressure from EOS
      eos_state % T = temp_curr
      eos_state % rho = dens_curr
      call eos(eos_input_rt, eos_state)
      pres_eos = eos_state % p
      entr_eos = eos_state % s

      ! Compute change in density based on Taylor expansions of EOS and HSE
      dp_hse = HALF*dx*grav
      dpdr = eos_state % dpdr
      dpdt = eos_state % dpdt
      dsdr = eos_state % dsdr
      dsdt = eos_state % dsdt
      determinant = dsdt*(dpdr - dp_hse) - dpdt*dsdr
      junk = dsdt*(pres_hse - pres_eos) - dpdT*(entr_iso - entr_eos)
      drho = junk / determinant
      junk = (dpdr - dp_hse)*(entr_iso - entr_eos) - dsdr*(pres_hse - pres_eos)
      dtmp = junk / determinant

      ! Restrict change in density to prevent runaway
      dens_curr = max(0.9*dens_curr, min(dens_curr + drho, 1.1*dens_curr))
      temp_curr = max(0.9*temp_curr, min(temp_curr + dtmp, 1.1*temp_curr))

      ! Check if convergence achieved
      if ((abs(drho / dens_curr) < TOLERANCE) .and. &
          (abs(dtmp / temp_curr) < TOLERANCE))then
         converged = .true.
         exit
      end if

      ! Check if solution wandered into fluff
      if (dens_curr < dens_fluff_cutoff) then
         dens_curr = dens_fluff_cutoff
         temp_curr = temp_fluff
         converged = .true.
         fluff = .true.
         exit
      end if

   end do

   ! Verify the solution converged
   if (.not. converged) then
      write(*,*) "Error zone", izone, " did not converge in init_1d"
      if (izone > ibase) then
         write(*,*) "integrate up"
      else
         write(*,*) "integrate down"
      end if
      write(*,*) "temp = ", temp_curr
      write(*,*) "dens = ", dens_curr
      write(*,*) "drho = ", drho
      write(*,*) "dtmp = ", dtmp
      write(*,*) "pres (EoS) = ", pres_eos
      write(*,*) "pres (HSE) = ", pres_hse
      write(*,*) "entr (EoS) = ", entr_eos
      write(*,*) "entr (iso) = ", entr_iso
      call bl_error("ERROR: HSE non-convergence")
   end if

   ! Use the fluff temperature as a floor
   if (temp_curr < temp_fluff) then
      temp_curr = temp_fluff
   end if

   pres_curr = pres_eos
   s_curr    = entr_eos

end subroutine HSE_NR_loop_isentropic

!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! computes the mass contained in the requested shell

function m_shell(i)

   use init_1d_variables
   use init_1d_grids, only: Ustate, Uradius, NrU
   use bl_types
   use bl_constants_module

   implicit none

   integer, intent(in) :: i
   real(kind=dp_t)     :: m_shell
   real(kind=dp_t)     :: r_upper, r_lower
   real(kind=dp_t)     :: dr

   if (i == 1) then
      ! r_upper is average of current cell-centered radius and next
      ! cell-centered radius up; assume cell is symmetric
      r_upper = HALF * (Uradius(2) + Uradius(1))
      dr      = 2.0d0 * (r_upper - Uradius(1))
      r_lower = r_upper - dr
   else if (i == NrU) then
      ! r_lower is average of current cell-centered radius and next
      ! cell-centered radius down; assume cell is symmetric
      r_lower = HALF * (Uradius(NrU) + Uradius(NrU-1))
      dr      = 2.0d0 * (Uradius(NrU) - r_lower)
      r_upper = r_lower + dr
   else
      r_upper = HALF * (Uradius(i+1) + Uradius(i))
      r_lower = HALF * (Uradius(i) + Uradius(i-1))
   end if

   m_shell = M_GEOM * Ustate(i, idens) * (r_upper**3 - r_lower**3)

end function m_shell

!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! This subroutine computes the HSE error and prints the hse_err file
!
! The arguments are
! - g_type  : in    : the type of gravity
! - mass    : in    : the mass of the star
! - outfile : in    : name of file to print hydrostatic equilibrium error data

subroutine compute_HSE_error(g_type, mass, temp_fluff, outfile)

   use eos_module
   use network
   use init_1d_variables
   use init_1d_grids
   use bl_types
   use bl_constants_module

   implicit none

   !===========================================================================
   ! Declare variables

   ! Arguments ................................................................
   integer,                      intent(in   ) :: g_type
   real(kind=dp_t),              intent(in   ) :: mass
   real(kind=dp_t),              intent(in   ) :: temp_fluff
   character(len=256),           intent(in   ) :: outfile

   ! Locals ...................................................................
   real(kind=dp_t), allocatable :: dpdr(:), rhog(:), hse_error(:)
   real(kind=dp_t)              :: max_hse_error, max_hse_radius, g_zone
   real(kind=dp_t)              :: g_const

   ! Temporaries ..............................................................
   integer :: i, n ! loop indices
   real(kind=dp_t) :: dx, r_edge

   ! Functions ................................................................
   real(kind=dp_t) :: m_shell

   !===========================================================================
   ! Set up

   allocate(dpdr(NrU-1))
   allocate(rhog(NrU-1))
   allocate(hse_error(NrU-1))

   !===========================================================================
   ! Construct the gravitational constant

   ! construct gravitational constant
   select case(g_type)
   case(GRAVITY_CONST)
      g_const = G * mass * Uradius(ibase)**(-2)
   case(GRAVITY_INVSQ)
      g_const = G * mass
   case(GRAVITY_MENCL)
      ! "mass" is enclosed mass at i = ibase-1/2; need to subtract off shell
      ! mass for shell ibase-1 to get to enclosed mass for i = ibase-3/2, and
      ! so on to get enclosed mass for i = 3/2, which is the boundary between
      ! the 1st and 2nd cells; thus we subtract shell mass for ibase-1 to 2
      ! However, the HSE loop will add back in the last shell mass prior to
      ! computing gravity for each iteration; thus we have to go down to the
      ! enclosed mass for i = 1/2, because prior to computing gravity at
      ! i = 3/2, the loop will update the enclosed mass by adding in the shell
      ! mass for i = 1
      g_const = mass
      do i = ibase-1, 1, -1
         g_const = g_const - m_shell(i)
      end do
   case default
      call bl_error('ERROR: invalid gravity mode in integrate_HSE subroutine')
   end select

   !===========================================================================
   ! Compute HSE error

   max_hse_error = ZERO
   do i = 1, NrU-1

      ! compute dx and r_mean
      dx = Uradius(i+1) - Uradius(i)
      r_edge = HALF * (Uradius(i+1) + Uradius(i))

      ! set gravity
      select case(g_type)
      case(GRAVITY_CONST)
         g_zone = g_const
      case(GRAVITY_INVSQ)
         g_zone = g_const / r_edge**2
      case(GRAVITY_MENCL)
         g_const = g_const + m_shell(i)
         g_zone = G * g_const / r_edge**2
      ! no 'case default' because that would have been caught in previous
      ! 'select' statement near the beginning of this subroutine
      end select

      ! compute d(pressure)/d(radius) and density*gravity
      dpdr(i) = (Ustate(i+1, ipres) - Ustate(i, ipres)) / dx
      rhog(i) = HALF * (Ustate(i+1, idens) + Ustate(i, idens)) * g_zone

      ! compute HSE error
      if (dpdr(i) /= ZERO) then
         hse_error(i) = abs((dpdr(i) - rhog(i)) / dpdr(i))
         if (Ustate(i, itemp) < 1.1d0*temp_fluff) then
            hse_error(i) = min(ONE, hse_error(i))
         end if
         ! check if we have a new maximum HSE error
         if (hse_error(i) > max_hse_error) then
            max_hse_error  = hse_error(i)
            max_hse_radius = r_edge
         end if
      else
         ! fake it if dpdr is zero
         if (rhog(i) /= ZERO) then
            hse_error(i) = ONE
         else
            hse_error(i) = ZERO
         end if
      end if
   end do

   !===========================================================================
   ! Print

   ! Print maximum HSE error and its location
   write(*,*) "Maximum HSE error is ", max_hse_error
   write(*,*) "Located at a radius of ", max_hse_radius

   ! Print HSE error file
   open(unit=99, file=outfile)
   write(99,'(a)') "# edge-centered radius"
   write(99,'(a)') "# dP/dr"
   write(99,'(a)') "# rho*g"
   write(99,'(a)') "# HSE error"
   do i = 1, NrU-1
      r_edge = HALF * (Uradius(i+1) + Uradius(i))
      write(99,'(1x,30(g17.10,1x))') r_edge, dpdr(i), rhog(i), hse_error(i)
   end do
   close(99)

   !===========================================================================
   ! Clean up

   deallocate(dpdr)
   deallocate(rhog)
   deallocate(hse_error)

end subroutine compute_HSE_error

