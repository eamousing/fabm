#:include 'setup_num.fypp'
#include "fabm_driver.h"

module imr_norwecom_num

    use fabm_types

    implicit none

    private

    type, extends(type_base_model), public :: type_imr_norwecom_num
    contains
        procedure :: initialize
    end type type_imr_norwecom_num

contains

    subroutine initialize(self, configunit)
        class(type_imr_norwecom_num), intent(inout), target :: self
        integer, intent(in) :: configunit
    end subroutine initialize

end module imr_norwecom_num