#include "fabm_driver.h"

module imr_norwecom

    use fabm_types

    implicit none

    private

    type, extends(type_base_model), public :: type_imr_norwecom
    contains
        procedure :: initialize
        procedure :: do_pelagic
        procedure :: get_vertical_movement
    end type

contains

    subroutine initialize(self, configunit)
        class(type_imr_norwecom), intent(inout), target :: self
        integer, intent(in) :: configunit
        
    end subroutine initialize

    subroutine do_pelagic(self, _ARGUMENTS_DO_)
        class(type_imr_norwecom), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        _LOOP_BEGIN_

        _LOOP_END_
    end subroutine do_pelagic

    subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
        class(type_imr_norwecom), intent(in) :: self
        _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

        _LOOP_BEGIN_

        _LOOP_END_
    end subroutine get_vertical_movement

end module imr_norwecom