module interpolate_module

  implicit none

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

end module interpolate_module
