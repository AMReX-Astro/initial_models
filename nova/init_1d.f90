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

  ! use implicit none or Doug will segfault your machine with his mind
  implicit none

  ! Declare variables

  ! files ....................................................................
  character(len=256) :: model_input   ! input model
  character(len=256) :: model_initial ! original model (reprinted to check)
  character(len=256) :: model_interp  ! interpolated model, no smoothing/HSE
  character(len=256) :: model_integ1  ! first HSE integration, no smoothing
  character(len=256) :: model_smooth  ! smoothed model, no final HSE integ
  character(len=256) :: model_maestro ! MAESTRO-readable file (final output)
  character(len=256) :: model_hse_err ! hydrostatic equilibrium error file

  ! gravity ..................................................................
  integer         :: g_type  ! gravity type
  real(kind=rt) :: M       ! mass of underlying star

  ! core profile .............................................................
  integer :: p_type ! core profile type
  integer :: c_edge ! extend from initial model edge or from core edge (ibase)

  ! smoothing ................................................................
  integer         :: s_type     ! smoothing type
  real(kind=rt) :: s_length   !smoothing length

  ! temporaries & miscellaneous ..............................................
  integer :: i                  ! loop index
  real(kind=rt) :: temp_fluff ! temperature in fluff region

  ! Setup

  ! Set defaults and parse command-line arguments
  g_type = GRAVITY_INVSQ           ! default: 1/r^2 gravity
  M = 1.98892d33                   ! default: 1 solar mass
  p_type = CORE_PROFILE_ISOTHERMAL ! default: isothermal profile
  c_edge = CORE_EXTEND             ! default: extend initial model
  s_type = SMOOTH_NONE             ! default: no smoothing
  s_length = 1.0d6                 ! default: 10 km
  NrU = 1024                       ! default: 2^10 zones
  rmin = 6.0d10                    ! default: 60,000 km
  rmax = 7.0d10                    ! default: 70,000 km
  call parse_args(g_type, M, p_type, c_edge, s_type, s_length, NrU, rmin, &
                  rmax, model_input)

  ! Initialize externals
  call eos_init()
  call network_init()

  ! Set parameters
  iH1 = ispec
  Nvars = nspec + ispec - 1

  ! Construct file names
  call construct_file_names(model_input , model_initial, model_interp , &
                            model_integ1, model_smooth, model_maestro,  &
                            model_hse_err, g_type, p_type, c_edge, s_type)

  ! Read initial model

  call read_input_model(model_input, temp_fluff, model_initial)

  ! Interpolate model to uniform grid

  call interpolate_to_uniform(c_edge, model_interp)

  ! First HSE pass (must do this before smoothing to account for chaning the
  ! profiles, because smoothing and the changing profiles may interact in
  ! unexpected and detrimental ways); print model_integ1

  call integrate_HSE(g_type, M, p_type, temp_fluff, model_integ1)

  ! Smooth the model; print model_smooth

  call smooth_model(s_type, s_length, model_smooth)

  ! Second HSE pass (need to correct for the changes made by smoothing; the
  ! question remains: do we do this integration isothermally or
  ! isentropically?); print model_maestro

  call integrate_HSE(g_type, M, CORE_PROFILE_SECOND_HSE, temp_fluff, &
                     model_maestro)

  ! Compute and print hse_err file

  call compute_HSE_error(g_type, M, temp_fluff, model_hse_err)

  ! Cleanup

  deallocate(Ustate)
  deallocate(Uradius)

end subroutine init_1d


!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! This subroutine parses the command-line arguments supplied to the calling
! program.
!
! The arguments are:
! - g_type     : inout : The type of gravity
! - M          : inout : The mass of the star
! - s_type     : inout : The type of smoothing
! - s_length   : inout : The smoothing length
! - nx         : inout : The number of radial zones
! - xmin       : inout : The minimum radius
! - xmax       : inout : The maximum radius
! - model_file :   out : The name of the input file
!
! Note that all parameters but the model_file are inouts to allow for default
! values to be supplied by the calling program.

subroutine parse_args(g_type, M, p_type, c_edge, s_type, s_length, nx, xmin, &
                      xmax, model_file)

  use amrex_fort_module, only: rt => amrex_real
   use amrex_constants_module
   use amrex_error_module
   use init_1d_variables

   implicit none

   !===========================================================================
   ! Declare variables

   ! Arguments ................................................................
   integer,            intent(inout) :: g_type
   real(kind=rt),    intent(inout) :: M
   integer,            intent(inout) :: p_type
   integer,            intent(inout) :: c_edge
   integer,            intent(inout) :: s_type
   real(kind=rt),    intent(inout) :: s_length
   integer,            intent(inout) :: nx
   real(kind=rt),    intent(inout) :: xmin
   real(kind=rt),    intent(inout) :: xmax
   character(len=256), intent(  out) :: model_file

   ! Locals ...................................................................
   integer            :: narg, farg
   character(len=128) :: fname
   character(len=4)   :: s_name
   logical            :: no_model_file = .true.

   !===========================================================================
   ! Read arguments

   narg = command_argument_count()
   farg = 1

   if (narg .eq. 0) then
      write (*,*) "For assistance, use the '--help' option."
      stop
   end if

   do while (farg <= narg)
      call get_command_argument(farg, value=fname)
      select case(fname)

      !------------------------------------------------------------------------
      ! Information options

      ! If the user asks for help, print a summary of command-line arguments
      ! and terminate the program
      case('--help')
         write (*,*) "Options for init_1d are as follows:"
         write (*,*) "  Gravity:"
         write (*,*) "    '-g g_type' to select gravity type; allowed &
                            &values for g_type are:"
         write (*,*) "      const: constant gravity (calibrated for the &
                              &base of the H1 layer)"
         write (*,*) "      invsq: 1/r^2 gravity"
         write (*,*) "      mencl: enclosed-mass gravity (calibrated for the &
                              &base of H1 layer)"
         write (*,*) "    '-m mass' to set the mass of the star (in g)"
         write (*,*) "  Core Profile:"
         write (*,*) "    -p p_type' to select profile type; allowed values &
                            &for p_type are:"
         write (*,*) "      thrm: isothermal profile"
         write (*,*) "      entr: isentropic profile"
         write (*,*) "    -c c_edge' to select core profile edge; &
                            &allowed values for c_edge are:"
         write (*,*) "      extend: extend initial model with profile"
         write (*,*) "      overwr: overwrite profile from base of envelope"
         write (*,*) "  Smoothing:"
         write (*,*) "    '-s s_type' to select smoothing type; allowed &
                            &values for s_type are:"
         write (*,*) "      none: no smoothing"
         write (*,*) "      krnl: Gaussian kernel smoothing"
         write (*,*) "      tnh0: hyperbolic transition across the &
                              &core-envelope boundary."
         write (*,*) "      tanh: hyperbolic transition 0.5 slen below the &
                              &core-envelope boundary"
         write (*,*) "      tnh2: hyperbolic transition 2.0 slen below the &
                              &core-envelope boundary"
         write (*,*) "      gssn: half-Gaussian transition below the &
                              &core-envelope boundary"
         write (*,*) "      expn: decaying-exponential transition below the &
                              &core-envelope boundary"
         write (*,*) "    '-slen s_length' to set smoothing length"
         write (*,*) "  Grid:"
         write (*,*) "    '-nx num_points' to set number of grid points"
         write (*,*) "    '-xmin min_radius' to set minimum radius (in cm)"
         write (*,*) "    '-xmax max_radius' to set maximum radius (in cm)"
         write (*,*) "  Model:"
         write (*,*) "    '--model_file filename' to set name of input model &
                            &file; REQUIRED"
         stop

      !------------------------------------------------------------------------
      ! Gravity options

      ! Set the gravity type: invsq (1/r^2 gravity from a point mass at r = 0),
      ! or constant
      case('-g')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         select case (trim(adjustl(fname)))
         case('const')
            g_type = GRAVITY_CONST
         case('invsq')
            g_type = GRAVITY_INVSQ
         case('mencl')
            g_type = GRAVITY_MENCL
         case default
            write (*,*) "Invalid value for option '-g'."
            write (*,*) "Must be one of 'invsq', 'const'."
            call bl_error("ERROR in argument parsing.")
         end select

      ! Set the value of the mass of the star: the point mass in the
      ! inverse-squared case, or the mass to get the correct gravity at the
      ! core-envelope
      case('-m')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         read(fname,*) M

      !------------------------------------------------------------------------
      ! Core profile options

      ! Set the core profile type: isothermal or isentropic
      case('-p')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         select case(trim(adjustl(fname)))
         case('thrm')
            p_type = CORE_PROFILE_ISOTHERMAL
         case('entr')
            p_type = CORE_PROFILE_ISENTROPIC
         case default
            write (*,*) "Invalid value for option '-p'."
            write (*,*) "Must be one of 'thrm', 'entr'."
            call bl_error("ERROR in argument parsing.")
         end select

      case('-c')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         select case(trim(adjustl(fname)))
         case('extend')
            c_edge = CORE_EXTEND
         case('overwr')
            c_edge = CORE_OVERWR
         case default
            write (*,*) "Invalid value for option '-c'."
            write (*,*) "Must be one of 'extend', 'overwr'."
            call bl_error("ERROR in argument parsing.")
         end select

      !------------------------------------------------------------------------
      ! Smoothing options

      ! Set the smoothing type: none, tanh (hyperbolic tangent smoothing of
      ! core-envelope transition), or krnl (Gaussian kernal smoothing of entire
      ! initial model)
      case('-s')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         fname = adjustl(fname)
         s_name = fname(1:4)
         select case (s_name)
         case('none')
            s_type = SMOOTH_NONE
         case('krnl')
            s_type = SMOOTH_KRNL
         case('tnh0')
            s_type = SMOOTH_TNH0
         case('tanh')
            s_type = SMOOTH_TANH
         case('tnh2')
            s_type = SMOOTH_TNH2
         case('gssn')
            s_type = SMOOTH_GSSN
         case('expn')
            s_type = SMOOTH_EXPN
         case default
            write (*,*) "Invalid value for option '-s'."
            write (*,*) "Must be one of 'none', 'krnl, 'tnh0', 'tanh', &
                        &'tnh2', 'gssn', 'expn'."
            call bl_error("ERROR in argument parsing.")
         end select

      ! Set the smoothing length to be used by tanh or krnl smoothing
      case('-slen')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         read(fname,*) s_length

      !------------------------------------------------------------------------
      ! Grid options

      ! Set the number of radial zones in the uniform grid
      case('-nx')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         read(fname,*) nx

      ! Set the minimum radius of the domain
      case('-xmin')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         read(fname,*) xmin

      ! Set the maximum radius of the domain
      case('-xmax')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         read(fname,*) xmax

      !------------------------------------------------------------------------
      ! Filename options

      ! Set the name of the input model file
      case('--model_file')
         farg = farg+1
         call get_command_argument(farg, value=fname)
         model_file = trim(fname)
         no_model_file = .false.

      case default
         write (*,*) "UNKNOWN option = ", fname
         call bl_error("ERROR in argument parsing.")

      end select

      farg = farg+1
   enddo

   !===========================================================================
   ! Verify command-line arguments

   ! If we are going to smooth the model, the smoothing length must be positive
   if ((s_type /= SMOOTH_NONE) .and. (s_length .le. ZERO)) then
      call bl_error("ERROR: smoothing length must be greater than zero.")
   end if

   ! A model file must be supplied
   if (no_model_file) then
      call bl_error("ERROR: no model file supplied.  Use &
                         &'--model_file option.")
   end if

   ! The mass of the star must be positive
   if (M .le. ZERO) then
      call bl_error("ERROR: stellar mass must be greater than zero.")
   end if

   ! The minimum radius must be less than the maximum radius; in the invsq and
   ! mencl cases the minimum radius must be positive, in the constant case the
   ! minimum radius must be non-negative
   if (g_type == GRAVITY_CONST) then
      if (xmin .lt. ZERO) then
         call bl_error("ERROR: lower bound must be non-negative when &
                            &using constant gravity.")
      end if
   else
      if (xmin .le. ZERO) then
         call bl_error("ERROR: lower bound must be greater than zero &
                            &when using invsq or mencl gravity.")
      end if
   end if
   if (xmax .le. xmin) then
      call bl_error("ERROR: upper bound must be greater than lower bound.")
   end if

end subroutine parse_args

!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! This subroutine constructs the output file names
!
! The arguments are:
! - model_input   : in    : the name of the input file
! - model_initial :   out : the name of the re-printed copy of the initial data
! - model_interp  :   out : the name of the interpolated data
! - model_integ1  :   out : the name of the first-pass HSE-integrated data
! - model_smooth  :   out : the name of the smoothed data
! - model_maestro :   out : the name of the MAESTRO-readable data
! - model_hse_err :   out : the name of the hydrostatic equilibrium error data
! - g_type        : in    : indicator for gravity type
! - s_type        : in    : indicator for smoothing type
! - p_type        : in    : indicator for core profile type
! - c_edge        : in    : indicator for core edge type

subroutine construct_file_names(model_input , model_initial, model_interp , &
                                model_integ1, model_smooth, model_maestro,  &
                                model_hse_err, g_type, p_type, c_edge, s_type)

   use init_1d_variables

   implicit none

   character(len=256), intent(in   ) :: model_input   ! input model
   character(len=256), intent(  out) :: model_initial ! original model
   character(len=256), intent(  out) :: model_interp  ! interpolated model
   character(len=256), intent(  out) :: model_integ1  ! 1st HSE integ model
   character(len=256), intent(  out) :: model_smooth  ! smoothed model
   character(len=256), intent(  out) :: model_maestro ! MAESTRO-readable file
   character(len=256), intent(  out) :: model_hse_err ! HSE error file
   integer,            intent(in   ) :: g_type
   integer,            intent(in   ) :: s_type
   integer,            intent(in   ) :: p_type
   integer,            intent(in   ) :: c_edge

   character(len=256) :: model_base    ! base name of model files
   integer            :: ipos

   ! Construct file names
   ipos = index(model_input, '.dat')
   model_base    = model_input(1:ipos-1)
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
