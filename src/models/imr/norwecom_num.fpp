#:include 'setup_num.fypp'
#include "fabm_driver.h"

module imr_norwecom_num

    use fabm_types
    use globals, only: dp, rhoCN
    use NUMmodel, only: setupGeneralistsOnly, calcDerivatives, idxN, idxDOC, idxB, nGrid, printRates, group

    implicit none

    private

    real(rk), parameter :: daysec = 86400.0_rk

    type, extends(type_base_model), public :: type_imr_norwecom_num
        ! State variables
        #:for p in range(NGROUPS)
        type(type_state_variable_id) :: id_p${p+1}$
        #:endfor
        type(type_state_variable_id) :: id_no3
        type(type_state_variable_id) :: id_doc
        type(type_state_variable_id) :: id_poc

        ! Dependencies
        type(type_dependency_id) :: id_temp
        type(type_dependency_id) :: id_par
    contains
        procedure :: initialize
        procedure :: do
        procedure :: get_vertical_movement
    end type type_imr_norwecom_num

    real(dp), allocatable :: u(:)
    real(dp), allocatable :: dudt(:)
    character(len=20) :: errorstr
    logical(1) :: errorio = .false.

contains

    subroutine initialize(self, configunit)
        class(type_imr_norwecom_num), intent(inout), target :: self
        integer, intent(in) :: configunit
        ! State variables
        #:for p in range(NGROUPS)
        call self%register_state_variable(self%id_p${p+1}$, "p${p+1}$", "ugC l-1", "Generalist size group ${p+1}$", minimum = 0.0_rk, initial_value = 0.005_rk + 0.01_rk * ${p+1}$)
        #:endfor
        call self%register_state_variable(self%id_no3, "no3", "ugN l-1", "Nitrate concentration", minimum = 0.0_rk, initial_value = 150.0_rk)
        call self%register_state_variable(self%id_doc, "doc", "ugC l-1", "Dissolved organic carbon", minimum = 0.0_rk, initial_value = 100.0_rk)
        call self%register_state_variable(self%id_poc, "poc", "ugC l-1", "Particulate organic matter", minimum = 0.0_rk, initial_value = 10.0_rk)

        ! Dependencies
        call self%register_dependency(self%id_temp, standard_variables%temperature)
        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

        ! Initialize unicellular size spectrum
        call setupGeneralistsOnly(${NGROUPS}$, errorio, errorstr)
        if (errorio .eqv. .true.) then
            print *, "Error reading parameter ", errorstr
            stop
        end if
        allocate(u(nGrid))
        allocate(dudt(nGrid))
    end subroutine initialize

    subroutine do(self, _ARGUMENTS_DO_)
        class(type_imr_norwecom_num), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        real(rk) :: par
        real(rk) :: temp

        _LOOP_BEGIN_ 

        u = 0.0_dp

        ! Get local concentrations
        #:for p in range(NGROUPS)
        _GET_(self%id_p${p+1}$, u(idxB+${p}$))
        #:endfor
        _GET_(self%id_no3, u(idxN))
        _GET_(self%id_doc, u(idxDOC))
        _GET_(self%id_temp, temp)
        _GET_(self%id_par, par)

        ! Convert photosynthetic radiative flux
        par = par / 0.217_rk ! W m-2 -> uE m-2 s-1

        call calcDerivatives(u, par, temp, 1.0_dp, dudt, .false.)

        #:for p in range(NGROUPS)
        _ADD_SOURCE_(self%id_p${p+1}$, (real(dudt(idxB+${p}$), kind=rk)/daysec))
        #:endfor
        _ADD_SOURCE_(self%id_no3, (real(dudt(idxN), kind=rk)/daysec))
        _ADD_SOURCE_(self%id_doc, (real(dudt(idxDOC), kind=rk)/daysec))

        _LOOP_END_
    end subroutine do

    subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
        !! Sets the vertical movement of state variables
        !!
        !! Status: velocity for generalists is currently initialized to 0, so no sinking
        !!   takes place. Use Stokes law?
        class(type_imr_norwecom_num), intent(in) :: self
        _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

        _LOOP_BEGIN_

        #:for p in range(NGROUPS)
        _ADD_VERTICAL_VELOCITY_(self%id_p${p+1}$, (real(group(1)%spec%velocity(${p+1}$), kind=rk)/daysec))
        #:endfor

        _LOOP_END_
    end subroutine get_vertical_movement

end module imr_norwecom_num
