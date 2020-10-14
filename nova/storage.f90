module init_1d_variables

  use amrex_fort_module, only: rt => amrex_real

   implicit none

   ! geometry factor for computer Mencl = (4/3) pi
   real(kind=rt), parameter :: M_GEOM  = 4.0d0/3.0d0 * 4.0d0*datan(1.0d0)

   ! character string constants
   integer, parameter :: MAX_VARNAME_LENGTH=80

   ! Newton-Raphson loop
   integer, parameter         :: MAX_NR_ITER = 250
   real(kind=rt), parameter :: TOLERANCE = 1.0d-9

   ! state storage indices and sizes
   integer, parameter :: idens = 1
   integer, parameter :: itemp = 2
   integer, parameter :: ipres = 3
   integer, parameter :: ientr = 4
   integer, parameter :: ispec = 5
   integer            :: iH1           ! index of H1 (generally iH1 == ispec)
   integer            :: ibase, iTpeak ! radial indices of H1 base, T_peak
   integer            :: icore         ! index below which core profile is used
   integer            :: Nvars

   ! thermodynamic quantities constants
   real(kind=rt), parameter :: smallX = 1.0e-6
   real(kind=rt), parameter :: H1_base = 6.5d-1
   real(kind=rt), parameter :: dens_fluff_cutoff = 1.0d-6

   ! smoothing constants
   integer, parameter :: SMOOTH_NONE = 0
   integer, parameter :: SMOOTH_KRNL = 1
   integer, parameter :: SMOOTH_TNH0 = 2
   integer, parameter :: SMOOTH_TANH = 3
   integer, parameter :: SMOOTH_TNH2 = 4
   integer, parameter :: SMOOTH_GSSN = 5
   integer, parameter :: SMOOTH_EXPN = 6

   ! core profile constants
   integer, parameter :: CORE_PROFILE_SECOND_HSE = 0  ! for the second HSE pass
   integer, parameter :: CORE_PROFILE_ISOTHERMAL = 1
   integer, parameter :: CORE_PROFILE_ISENTROPIC = 2

   integer, parameter :: CORE_EXTEND = 0
   integer, parameter :: CORE_OVERWR = 1

   ! gravity constants
   real(kind=rt), parameter :: G = -6.67428d-8 ! Newton's grav constant
   integer,         parameter :: GRAVITY_CONST = 0
   integer,         parameter :: GRAVITY_INVSQ = 1
   integer,         parameter :: GRAVITY_MENCL = 2
end module init_1d_variables

module init_1d_grids

   use init_1d_variables
   use amrex_fort_module, only: rt => amrex_real


   implicit none

   real(kind=rt), allocatable :: Istate(:,:)   ! input grid
   real(kind=rt), allocatable :: Iradius(:)    ! cell-centered radii
   real(kind=rt), allocatable :: Ustate(:,:)   ! uniform grid
   real(kind=rt), allocatable :: Uradius(:)    ! cell-centered radii
   real(kind=rt)              :: dr            ! spacing of grid
   integer                      :: NrI      ! number of grid zones 
end module init_1d_grids
