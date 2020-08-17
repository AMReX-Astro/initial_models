module initialization_module

contains

  subroutine do_initialization(inputs_name, inputs_len) bind(C)

    use microphysics_module

    implicit none

    integer, value :: inputs_len
    integer :: inputs_name(inputs_len)
    integer :: i, narg

    call runtime_init(inputs_name, inputs_len)

    call microphysics_init()

  end subroutine do_initialization

end module initialization_module
