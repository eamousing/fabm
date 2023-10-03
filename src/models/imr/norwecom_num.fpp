#:include 'setup_num.fypp'
#include "fabm_driver.h"

module imr_norwecom_num

    use fabm_types
    use globals, only: dp, rhoCN
    use NUMmodel
    use POM, only: spectrumPOM

    implicit none

    private

    real(rk), parameter :: daysec = 86400.0_rk

    type, extends(type_base_model), public :: type_imr_norwecom_num
        ! State variables
        #:for p in range(NGENERALISTS)
        type(type_state_variable_id) :: id_p${p+1}$
        #:endfor
        #:for p in range(NPOM)
        type(type_state_variable_id) :: id_pom${p+1}$
        #:endfor
        type(type_state_variable_id) :: id_no3
        type(type_state_variable_id) :: id_doc
        type(type_state_variable_id) :: id_poc

        ! Dependencies
        type(type_dependency_id) :: id_temp
        type(type_dependency_id) :: id_par

        ! Diagnostic variables
        type(type_diagnostic_variable_id) :: id_n2p1
    contains
        procedure :: initialize
        procedure :: do
        procedure :: get_vertical_movement
    end type type_imr_norwecom_num

    real(dp), allocatable :: u(:)
    real(dp), allocatable :: dudt(:)
    real(dp), allocatable :: mort_mat(:,:)
    character(len=20) :: errorstr
    logical(1) :: errorio = .false.

    integer, save :: count = 0

contains

    subroutine initialize(self, configunit)
        class(type_imr_norwecom_num), intent(inout), target :: self
        integer, intent(in) :: configunit
        ! State variables
        #:for p in range(NGENERALISTS)
        call self%register_state_variable(self%id_p${p+1}$, "p${p+1}$", "ugC l-1", "Generalist size group ${p+1}$", minimum = 0.0_rk, initial_value = 0.005_rk + 0.01_rk * ${p+1}$)
        #:endfor
        #:for p in range(NPOM)
        call self%register_state_variable(self%id_pom${p+1}$, "pom${p+1}$", "ugC l-1", "Particulate organic matter size group ${p+1}$", minimum = 0.0_rk, initial_value = 10.0_rk)
        #:endfor
        call self%register_state_variable(self%id_no3, "no3", "ugN l-1", "Nitrate concentration", minimum = 0.0_rk, initial_value = 150.0_rk)
        call self%register_state_variable(self%id_doc, "doc", "ugC l-1", "Dissolved organic carbon", minimum = 0.0_rk, initial_value = 150.0_rk)
        call self%register_state_variable(self%id_poc, "poc", "ugC l-1", "Particulate organic matter", minimum = 0.0_rk, initial_value = 10.0_rk)

        ! Dependencies
        call self%register_dependency(self%id_temp, standard_variables%temperature)
        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

        ! Diagnostics
        call self%register_diagnostic_variable(self%id_n2p1, "n2p1", "ugN", "nitrogen to p1") ! for testing

        ! Initialize unicellular size spectrum
        call setupGeneralistsPOM(${NGENERALISTS}$, ${NPOM}$, errorio, errorstr)
        if (errorio .eqv. .true.) then
            print *, "Error reading parameter ", errorstr
            stop
        end if
        allocate(u(nGrid))
        allocate(dudt(nGrid))
        allocate(mort_mat(nGrid,nGrid))
    end subroutine initialize

    subroutine do(self, _ARGUMENTS_DO_)
        class(type_imr_norwecom_num), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        real(rk) :: par
        real(rk) :: temp

        ! count = count + 1
        ! if (mod(count*600, 86400) == 0) print *, count

        _LOOP_BEGIN_ 

        u = 0.0_dp

        ! Get local concentrations
        #:for p in range(NGENERALISTS)
        _GET_(self%id_p${p+1}$, u(idxB+${p}$))
        #:endfor
        #:for p in range(NPOM)
        _GET_(self%id_pom${p+1}$, u(ixStart(idxPOM)+${p}$))
        #:endfor
        _GET_(self%id_no3, u(idxN))
        _GET_(self%id_doc, u(idxDOC))
        _GET_(self%id_temp, temp)
        _GET_(self%id_par, par)

        ! Convert photosynthetic radiative flux
        par = par / 0.217_rk ! W m-2 -> uE m-2 s-1

        mort_mat = 0.0_dp
        call calcDerivatives(u, par, temp, 1.0_dp, dudt, .false., mort_mat)

        _SET_DIAGNOSTIC_(self%id_n2p1, dudt(idxN)/rhoCN) ! for testing

        #:for p in range(NGENERALISTS)
        _ADD_SOURCE_(self%id_p${p+1}$, (real(dudt(idxB+${p}$), kind=rk)/daysec))
        #:endfor
        #:for p in range(NPOM)
        _ADD_SOURCE_(self%id_pom${p+1}$, (real(dudt(ixStart(idxPOM)+${p}$), kind=rk)/daysec))
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

        integer :: iGroup

        _LOOP_BEGIN_

        do iGroup = 1, nGroups
            select type(spec => group(iGroup)%spec)
            type is(spectrumGeneralists)
                #:for p in range(NGENERALISTS)
                _ADD_VERTICAL_VELOCITY_(self%id_p${p+1}$, -(real(spec%velocity(${p+1}$), kind=rk)/daysec))
                #:endfor
            type is(spectrumPOM)
                #:for p in range(NPOM)
                _ADD_VERTICAL_VELOCITY_(self%id_pom${p+1}$, -(real(spec%velocity(${p+1}$), kind=rk)/daysec))
                #:endfor
            end select
        end do

        _LOOP_END_
    end subroutine get_vertical_movement

end module imr_norwecom_num
