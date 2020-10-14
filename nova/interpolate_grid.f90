!------------------------------------------------------------------------------
!==============================================================================
!------------------------------------------------------------------------------
! This subroutine interpolates the input model to a uniform grid
!
! The arguments are
! - outfile : in    : name of file to print interpolated uniform grid to

subroutine interpolate_to_uniform(outfile)

   use network
   use init_1d_variables
   use init_1d_grids
   use amrex_fort_module, only: rt => amrex_real
   use amrex_constants_module
   use extern_probin_module

   implicit none

   !===========================================================================
   ! Declare variables

   ! Arguments ................................................................
   character(len=256), intent(in   ) :: outfile

   ! Temporaries ..............................................................
   integer :: i, j, n               ! loop indices
   logical :: at_H1_base = .false.  ! have we found base of the H1 layer yet?
   logical :: at_Tpeak = .false.    ! have we found peak temperature yet?
   real(kind=rt), dimension(nspec) :: composition
   real(kind=rt) :: Xsum ! accumulator
   real(kind=rt) :: Tmin ! minimum core temperature

   ! Functions ................................................................
   real(kind=rt) :: interpolate   ! interpolate a quantity

   !===========================================================================
   ! Allocate uniform grid and compute radii

   ! construct Uradius and Ustate
   allocate(Ustate(NrU,Nvars))
   allocate(Uradius(NrU))
   dr = (rmax - rmin) / dble(NrU)
   do i = 1, NrU
      Uradius(i) = rmin + (dble(i) - HALF)*dr
   end do

   Tmin = Istate(1, itemp)

   !===========================================================================
   ! Loop over all grid zones, all variables

   do i = 1, NrU        ! Loop over uniform grid
      
      do n = 1, Nvars   ! Loop over all variables

         ! Don't interpolate entropy; that will be calculated later
         if (n == ientr) then
            cycle
         end if

         !---------------------------------------------------------------------
         ! Extend uniform grid below input grid
         if (Uradius(i) < Iradius(1)) then
            Ustate(i,n) = Istate(1,n)
            ! last time in this branch will set icore to lowest index that
            ! is inside the initial model
            icore = i + 1

         !---------------------------------------------------------------------
         ! Interpolate input grid to uniform grid
         else if (Uradius(i) < Iradius(NrI)) then
            Ustate(i,n) = interpolate(Uradius(i), NrI, Iradius, Istate(:,n))

            if (.not. at_H1_base) then

               ! Find the index of the minimum temperature below the H-1 base
               if ((c_edge == CORE_OVERWR) .and. (n == itemp) &
                  .and. (Ustate(i, itemp) <= Tmin)) then
                  Tmin = Ustate(i, itemp)
                  icore = i
               end if

               ! Find the index of H-1 base on the uniform grid
               if ((n == iH1) .and. (Ustate(i,iH1) >= H1_base)) then
                  at_H1_base = .true.
                  ibase = i
                  write (*,*) 'core profile edge:'
                  write (*,*) '   i     = ', icore
                  write (*,*) '   X[i]  = ', Uradius(icore)
                  write (*,*) '   T[i]  = ', Ustate(icore, itemp)
                  write (*,*) 'H1 base:'
                  write (*,*) '   i     = ', ibase
                  write (*,*) '   X[i]  = ', Uradius(ibase)
                  write (*,*) '   H1[i] = ', Ustate(ibase,iH1)
               endif
            end if

            ! Find the first temperature peak after the H-1 base
            if ((n == itemp) .and. at_H1_base .and. (.not. at_Tpeak) &
                             .and.  (Ustate(i,itemp) < Ustate(i-1,itemp))) then
               iTpeak = i-1
               at_Tpeak = .true.
               write (*,*) 'peak temperature:'
               write (*,*) '   i    = ', iTpeak
               write (*,*) '   X[i] = ', Uradius(iTpeak)
               write (*,*) '   T[i] = ', Ustate(iTpeak,itemp)
            endif

         !---------------------------------------------------------------------
         ! Extend uniform grid above input grid
         else
            Ustate(i,n) = Istate(NrI,n)
         endif

      enddo    ! end loop over variables

      !------------------------------------------------------------------------
      ! Enforce composition limits and renormalize
      composition(:) = Ustate(i,ispec:Nvars)
      Xsum = 0.0
      do n = 1, nspec
         composition(n) = min(max(smallX, composition(n)), ONE)
      enddo
      do n = 1, nspec
         Xsum = Xsum + composition(n)
      enddo
      do n = 1, nspec
         composition(n) = composition(n)/Xsum
      enddo
      Ustate(i,ispec:Nvars) = composition(:)

   enddo        ! end loop over zones

   !===========================================================================
   ! Deallocate unneccesary storage
   deallocate(Istate)
   deallocate(Iradius)

   !===========================================================================
   ! Print
   open(99, file=outfile)
   write(99,*) "# initial model just after putting onto a uniform grid"
   do i = 1, NrU
      write(99,'(1x,30(g17.10,1x))') Uradius(i), (Ustate(i,j), j = 1, Nvars)
   end do
   close(unit=99)

end subroutine interpolate_to_uniform

! Interpolation function taken from old init_1d

function interpolate(r, npts, model_r, model_var)

   use init_1d_variables
   use amrex_fort_module, only: rt => amrex_real


   implicit none


   ! given the array of model coordinates (model_r), and variable (model_var),
   ! find the value of model_var at point r using linear interpolation.
   ! Eventually, we can do something fancier here.

   real(kind=rt) :: interpolate
   real(kind=rt), intent(in) :: r
   integer :: npts
   real(kind=rt), dimension(npts) :: model_r, model_var
   real(kind=rt) :: val, slope
   real(kind=rt) :: minvar, maxvar

   logical :: found

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
      minvar = min(model_var(id+1),model_var(id))
      maxvar = max(model_var(id+1),model_var(id))
      interpolate = max(interpolate,minvar)
      interpolate = min(interpolate,maxvar)

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

