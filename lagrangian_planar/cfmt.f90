module cfmt_module

  implicit none

contains
  function cfmt(x) result (safe_x)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real (kind=rt), intent(in) :: x
    real (kind=rt) :: safe_x

    safe_x = x
    if (abs(x) < 1.e-98_rt) then
       safe_x = 1.e-98_rt
    end if

  end function cfmt
end module cfmt_module

