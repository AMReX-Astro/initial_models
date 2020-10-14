!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial model builder
! Adapted by Brendan Krueger from pre-existing MAESTRO initial model builder
! 25 March 2011
!
! This code takes an initial model from a Lagrangian code, puts it into a
! uniform grid and enforces hydrostatic equilibrium (HSE) and thermodynamic
! consistency.  It creates several files showing the model at intermediate
! steps, then a final .hse file that can be read directly by MAESTRO and a .err
! file that shows the hydrostatic equilibrium error.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_1d() bind(C, name="init_1d")

  ! modules
  use eos_module
  use network
  use init_1d_variables
  use init_1d_grids
  use extern_probin_module

  ! use implicit none or Doug will segfault your machine with his mind
  implicit none

  ! Declare variables

  character(len=256) :: model_input   ! input model
  character(len=256) :: model_initial ! original model (reprinted to check)
  character(len=256) :: model_interp  ! interpolated model, no smoothing/HSE
  character(len=256) :: model_integ1  ! first HSE integration, no smoothing
  character(len=256) :: model_smooth  ! smoothed model, no final HSE integ
  character(len=256) :: model_maestro ! MAESTRO-readable file (final output)
  character(len=256) :: model_hse_err ! hydrostatic equilibrium error file

  ! temporaries & miscellaneous ..............................................
  integer :: i                  ! loop index
  real(kind=rt) :: temp_fluff ! temperature in fluff region

  ! Setup

  call check_args()

  ! Set parameters
  iH1 = ispec
  Nvars = nspec + ispec - 1

  ! Construct file names
  call construct_file_names(model_initial, model_interp , &
                            model_integ1, model_smooth, model_maestro,  &
                            model_hse_err)

  ! Read initial model

  call read_input_model(model_file, temp_fluff, model_initial)

  ! Interpolate model to uniform grid

  print *, "calling interpolate to uniform"
  call interpolate_to_uniform(model_interp)
  print *, "done"

  ! First HSE pass (must do this before smoothing to account for chaning the
  ! profiles, because smoothing and the changing profiles may interact in
  ! unexpected and detrimental ways); print model_integ1

  call integrate_HSE(p_type, temp_fluff, model_integ1)

  ! Smooth the model; print model_smooth

  call smooth_model(model_smooth)

  ! Second HSE pass (need to correct for the changes made by smoothing; the
  ! question remains: do we do this integration isothermally or
  ! isentropically?); print model_maestro

  call integrate_HSE(CORE_PROFILE_SECOND_HSE, temp_fluff, &
                     model_maestro)

  ! Compute and print hse_err file

  call compute_HSE_error(temp_fluff, model_hse_err)

  ! Cleanup

  deallocate(Ustate)
  deallocate(Uradius)

end subroutine init_1d



subroutine check_args()

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module
  use amrex_error_module
  use init_1d_variables
  use extern_probin_module

   implicit none

   ! Verify command-line arguments

   ! If we are going to smooth the model, the smoothing length must be
   ! positive
   if ((s_type /= SMOOTH_NONE) .and. (s_length .le. ZERO)) then
      call bl_error("ERROR: smoothing length must be greater than zero.")
   end if

   ! The mass of the star must be positive
   if (M .le. ZERO) then
      call bl_error("ERROR: stellar mass must be greater than zero.")
   end if

   ! The minimum radius must be less than the maximum radius; in the
   ! invsq and mencl cases the minimum radius must be positive, in the
   ! constant case the minimum radius must be non-negative
   if (g_type == GRAVITY_CONST) then
      if (rmin .lt. ZERO) then
         call bl_error("ERROR: lower bound must be non-negative when &
                            &using constant gravity.")
      end if
   else
      if (rmin .le. ZERO) then
         call bl_error("ERROR: lower bound must be greater than zero &
                            &when using invsq or mencl gravity.")
      end if
   end if
   if (rmax .le. rmin) then
      call bl_error("ERROR: upper bound must be greater than lower bound.")
   end if

 end subroutine check_args


! This subroutine constructs the output file names
!
! The arguments are:
! - model_initial :   out : the name of the re-printed copy of the initial data
! - model_interp  :   out : the name of the interpolated data
! - model_integ1  :   out : the name of the first-pass HSE-integrated data
! - model_smooth  :   out : the name of the smoothed data
! - model_maestro :   out : the name of the MAESTRO-readable data
! - model_hse_err :   out : the name of the hydrostatic equilibrium error data

subroutine construct_file_names(model_initial, model_interp , &
                                model_integ1, model_smooth, model_maestro,  &
                                model_hse_err)

  use extern_probin_module
   use init_1d_variables

   implicit none

   character(len=256), intent(  out) :: model_initial ! original model
   character(len=256), intent(  out) :: model_interp  ! interpolated model
   character(len=256), intent(  out) :: model_integ1  ! 1st HSE integ model
   character(len=256), intent(  out) :: model_smooth  ! smoothed model
   character(len=256), intent(  out) :: model_maestro ! MAESTRO-readable file
   character(len=256), intent(  out) :: model_hse_err ! HSE error file

   character(len=256) :: model_base    ! base name of model files
   integer            :: ipos

   ! Construct file names
   ipos = index(model_file, '.dat')
   model_base    = model_file(1:ipos-1)
   model_initial = trim(model_base) // '.initial'
   model_interp  = trim(model_base) // '.interp'
   model_integ1  = trim(model_base) // '.integ1'
   model_smooth  = trim(model_base) // '.smooth'
   model_hse_err = trim(model_base) // '.hse_err'

   select case(g_type)
   case(GRAVITY_CONST)
      model_maestro = trim(model_base) // '_Gconst'
   case(GRAVITY_INVSQ)
      model_maestro = trim(model_base) // '_Ginvsq'
   case(GRAVITY_MENCL)
      model_maestro = trim(model_base) // '_Gmencl'
   case default
      model_maestro = trim(model_base) // '_G-----'
   end select

   model_maestro = trim(model_maestro) // '_P'
   select case(p_type)
   case(CORE_PROFILE_ISOTHERMAL)
      model_maestro = trim(model_maestro) // 'th' ! isoTHermal
   case(CORE_PROFILE_ISENTROPIC)
      model_maestro = trim(model_maestro) // 'en' ! isENtropic
   case default
      model_maestro = trim(model_maestro) // '--'
   end select

   select case(c_edge)
   case(CORE_EXTEND)
      model_maestro = trim(model_maestro) // 'ex' ! EXtend
   case(CORE_OVERWR)
      model_maestro = trim(model_maestro) // 'ow' ! OverWrite
   case default
      model_maestro = trim(model_maestro) // '--'
   end select

   select case(s_type)
   case(SMOOTH_NONE)
      model_maestro = trim(model_maestro) // '_Snone'
   case(SMOOTH_KRNL)
      model_maestro = trim(model_maestro) // '_Skrnl'
   case(SMOOTH_TNH0)
      model_maestro = trim(model_maestro) // '_Stnh0'
   case(SMOOTH_TANH)
      model_maestro = trim(model_maestro) // '_Stanh'
   case(SMOOTH_TNH2)
      model_maestro = trim(model_maestro) // '_Stnh2'
   case(SMOOTH_GSSN)
      model_maestro = trim(model_maestro) // '_Sgssn'
   case(SMOOTH_EXPN)
      model_maestro = trim(model_maestro) // '_Sexpn'
   case default
      model_maestro = trim(model_maestro) // '_S----'
   end select

   model_maestro = trim(model_maestro) // '.hse'

end subroutine construct_file_names
