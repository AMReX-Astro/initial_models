!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! This subroutine reads the input model from a file and sorts it into the
! appropriate data structures for later use
!
! The arguments are:
! - infile  : in    : the name of the input file to read from
! - outfile : in    : the name of the output file to re-print the input data
!                     for purposes of verification

subroutine read_input_model(infile, temp_fluff, outfile)

   use network
   use init_1d_variables
   use init_1d_grids
   use bl_types
   use bl_constants_module

   implicit none

   !===========================================================================
   ! Declare variables

   ! Arguments ................................................................
   character(len=256),           intent(in   ) :: infile
   real(kind=dp_t),              intent(  out) :: temp_fluff
   character(len=256),           intent(in   ) :: outfile

   ! Locals ...................................................................
   integer                                        :: Nvars_inp ! # inputs vars
   real(kind=dp_t), allocatable                   :: vars(:)   ! zone values
   character(len=MAX_VARNAME_LENGTH), allocatable :: varnames(:)
   logical :: found

   ! Temporaries ..............................................................
   character(len=256) :: line ! a line read from the input file
   integer            :: ipos ! a string position for trimming
   integer            :: i, j, n ! loop indices
   real(kind=dp_t)    :: Xsum ! for normalizing abundances
   real(kind=dp_t), dimension(nspec) :: species ! for normalizing abundances

   !===========================================================================
   ! Get number of input variables and number of zones; then allocate

   ! Open input file
   open(99, file=infile, status="old")

   ! Get number of data zones from first line
   read(99, '(a256)') line
   ipos = index(line, '=') + 1
   read(line(ipos:),*) NrI
   write(*,*) NrI, " points found in the initial model file."

   ! Get number of variables (less radial coordinate) from second line
   read(99, '(a256)') line
   ipos = index(line, '=') + 1
   read(line(ipos:),*) Nvars_inp
   write(*,*) Nvars_inp, " variables found in the initial model file."

   ! Allocate temporary arrays
   allocate(vars(Nvars_inp))
   allocate(varnames(Nvars_inp))


   ! Allocate storage arrays
   allocate(Iradius(NrI))
   allocate(Istate(NrI, Nvars))

   !===========================================================================
   ! Read headers

   ! Extract and store column headers
   do i = 1, Nvars_inp
      read(99, '(a256)') line
      ipos = index(line, '#') + 1
      varnames(i) = trim(adjustl(line(ipos:)))
   end do

   !===========================================================================
   ! Read data

   do i = 1, NrI

      ! Read a line of data
      read(99,*) Iradius(i), (vars(j), j = 1, Nvars_inp)

      ! Clear state for this line
      Istate(i,:) = ZERO


      ! Loop over all variables in data line
      do j = 1, Nvars_inp

         found = .false.

         select case(trim(varnames(j)))
         case("dens") ! Save density
            Istate(i, idens) = vars(j)
            found = .true.
         case("temp") ! Save temperature
            Istate(i, itemp) = vars(j)
            found = .true.
         case("pres") ! Save pressure
            Istate(i, ipres) = vars(j)
            found = .true.
         case("vely") ! Throw away radial velocity
            found = .true.
         case default ! Save species
            do n = 1, Nspec
               if (varnames(j) == spec_names(n)) then
                  Istate(i, ispec-1+n) = min(ONE, max(vars(j), smallX))
                  found = .true.
                  exit
               end if
            end do
         end select

         ! Error if we can't find this variable
         if (.not. found) then
            write(*,*) "ERROR: variable not found: ", varnames(j)
         end if

      end do

      ! Normalize species
      Xsum = ZERO
      species(:) = Istate(i, ispec:Nvars)
      do n = 1, Nspec
         Xsum = Xsum + species(n)
      end do
      Istate(i, ispec:Nvars) = species(:) / Xsum

   end do

   ! Set fluff temperature
   temp_fluff = minval(Istate(:,itemp))

   !===========================================================================
   ! Clean up

   ! Close inputs file
   close(unit=99)

   ! Deallocate local arrays
   deallocate(vars)
   deallocate(varnames)

   !===========================================================================
   ! Re-print
   open(99, file=outfile)
   write(99,*) "# initial model read in"
   do i = 1, NrI
      write(99,'(1x,30(g17.10,1x))') Iradius(i), (Istate(i,j), j = 1, Nvars)
   end do
   close(unit=99)

end subroutine read_input_model

