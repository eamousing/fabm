#include "fabm_driver.h"

module imr_norwecom

    use fabm_types

    implicit none

    private

    type, extends(type_base_model), public :: type_imr_norwecom
        ! Register parameter identifiers
        real(rk) :: AA1 !! Diatom production maximum at 0 degC (s-1)
        real(rk) :: AA2 !! Diatom temperature dependent pmax (degC-1)
        real(rk) :: AA3 !! Flagellate production maximum at 0 degC (s-1)
        real(rk) :: AA4 !! Flagellate temperature dependent pmax (degC-1)
        real(rk) :: AA5 !! Phytoplankton respiration rate at 0 degC (s-1)
        real(rk) :: AA6 !! Phytoplankton respiration rate temperature dependence (degC-1)
        real(rk) :: CC3 !! Phytoplankton death rate (s-1)
        real(rk) :: CC4 !! Detritus decomposition rate (s-1)
        real(rk) :: SCC4 !! Biogenic silica decomposition rate (s-1)
        real(rk) :: KN_DIA !! Diatom affinity for nitrate ((mmol N)-1 s-1)
        real(rk) :: KN_FLA !! Flagellate affinity for nitrate ((mmol N)-1 s-1)
        real(rk) :: KP_DIA !! Diatom affinity for phosphate ((mmol P)-1 s-1)
        real(rk) :: KP_FLA !! Flagellate affinity for phosphate ((mmol P)-1 s-1)
        real(rk) :: KR_DIA !! Diatom affinity for radition (m2 (W s)-1)
        real(rk) :: KR_FLA !! Flagellate affinity for radition (m2 (W s)-1)
        real(rk) :: KS_DIA !! Diatom affinity for silicate ((mmol Si)-1 s-1)
        real(rk) :: SR_DET !! Detritus sinking rate (m s-1)
        real(rk) :: SR_DIAMAX !! Diatom maximum sinking rate (m s-1)
        real(rk) :: SR_DIAMIN !! Diatom minimum sinking rate (m s-1)
        real(rk) :: SR_FLA !! Flagellate sinking rate (m s-1)
        real(rk) :: SR_SIS !! Biogenic silica sinking rate (m s-1)       
        real(rk) :: SIB !! Concentration of SIL when max sinking speed of DIA (mmol Si m-3)
        real(rk) :: P_N_RATIO !! Phosphorus to nitrogen ratio (mmol P (mmol N)-1)
        real(rk) :: SI_N_RATIO !! Silicon to nitrogen ratio (mmol Si (mmol N)-1)

        ! Register variable identifiers
        type(type_state_variable_id) :: id_nit !! Nitrate
        type(type_state_variable_id) :: id_pho !! Phosphate
        type(type_state_variable_id) :: id_sil !! Silicate
        type(type_state_variable_id) :: id_sis !! Biogenic silica
        type(type_state_variable_id) :: id_det !! Nitrogen detritus
        type(type_state_variable_id) :: id_detp !! Phosphorus detritus
        type(type_state_variable_id) :: id_oxy !! Oxygen
        type(type_state_variable_id) :: id_fla !! Flagellates
        type(type_state_variable_id) :: id_dia !! Diatoms

        ! Register dependencies identifiers
        type(type_dependency_id) :: id_par !! Photoactive radiation
        type(type_dependency_id) :: id_temp !! Temperature

    contains

        procedure :: initialize
        procedure :: do_pelagic
        procedure :: get_vertical_movement
    end type

contains

    subroutine initialize(self, configunit)
        class(type_imr_norwecom), intent(inout), target :: self
        integer, intent(in) :: configunit

        ! Initialize model parameters
        call self%get_parameter(self%AA1, "AA1", "s-1", "Diatom production maximum at 0 degC", default=1.505e-5_rk)
        call self%get_parameter(self%AA2, "AA2", "degC-1", "Diatom temperature dependent pmax", default=0.063_rk)
        call self%get_parameter(self%AA3, "AA3", "s-1", "Flagellate production maximum at 0 degC", default=1.042e-5_rk)
        call self%get_parameter(self%AA4, "AA4", "degC-1", "Flagellate temperature dependent pmax", default=0.063_rk)
        call self%get_parameter(self%AA5, "AA5", "s-1", "Phytoplankton respiration rate at 0 degC", default=8.102e-8_rk)
        call self%get_parameter(self%AA6, "AA6", "degC-1", "Phytoplankton respiration rate temperature dependence", default=0.07_rk)
        call self%get_parameter(self%CC3, "CC3", "s-1", "Phytoplankton death rate", default=1.620e-6_rk)
        call self%get_parameter(self%CC4, "CC4", "s-1", "Detritus decomposition rate", default=1.620e-6_rk)
        call self%get_parameter(self%SCC4, "SCC4", "s-1", "Biogenic silica decomposition rate", default=1.505e-8_rk)
        call self%get_parameter(self%KN_DIA, "KN_DIA", "(mmol N)-1 s-1", "Diatom affinity for nitrate", default=1.736e-5_rk)
        call self%get_parameter(self%KN_FLA, "KN_FLA", "(mmol N)-1 s-1", "Flagellate affinity for nitrate", default=1.505e-5_rk)
        call self%get_parameter(self%KP_DIA, "KP_DIA", "(mmol P)-1 s-1", "Diatom affinity for phosphate", default=2.697e-4_rk)
        call self%get_parameter(self%KP_FLA, "KP_FLA", "(mmol P)-1 s-1", "Flagellate affinity for phosphate", default=2.5e-4_rk)
        call self%get_parameter(self%KR_DIA, "KR_DIA", "m2 (W s)-1", "Diatom affinity for radition", default=1.5e-6_rk)
        call self%get_parameter(self%KR_FLA, "KR_FLA", "m2 (W s)-1", "Flagellate affinity for radition", default=4.6e-7_rk)
        call self%get_parameter(self%KS_DIA, "KS_DIA", "(mmol Si)-1 s-1", "Diatom affinity for silicate", default=2.546e-5_rk)
        call self%get_parameter(self%SR_DET, "SR_DET", "m s-1", "Detritus sinking rate", default=-3.472e-5_rk)
        call self%get_parameter(self%SR_DIAMAX, "SR_DIAMAX", "m s-1", "Diatom maximum sinking rate", default=-3.472e-5_rk)
        call self%get_parameter(self%SR_DIAMIN, "SR_DIAMIN", "m s-1", "Diatom minimum sinking rate", default=-3.472e-6_rk)
        call self%get_parameter(self%SR_FLA, "SR_FLA", "m s-1", "Flagellate sinking rate", default=-2.894e-6_rk)
        call self%get_parameter(self%SR_SIS, "SR_SIS", "m s-1", "Biogenic silica sinking rate", default=-3.472e-5_rk)
        call self%get_parameter(self%SIB, "SIB", "mmol Si m-3", "Concentration of SIL when max sinking speed of DIA", default=1.0_rk)
        call self%get_parameter(self%P_N_RATIO, "P_N_RATIO", "mmol P (mmol N)-1", "Phosphorus to nitrogen ratio", default=0.0625_rk)
        call self%get_parameter(self%SI_N_RATIO, "SI_N_RATIO", "mmol Si (mmol N)-1", "Silicon to nitrogen ratio", default=0.875_rk)
        
        ! Initialize state variables
        call self%register_state_variable(self%id_nit, "NIT", "mmol N m-3", "Nitrate", minimum=0.0_rk, initial_value=10.0_rk)
        call self%register_state_variable(self%id_pho, "PHO", "mmol P m-3", "Phosphate", minimum=0.0_rk, initial_value=0.7_rk)
        call self%register_state_variable(self%id_sil, "SIL", "mmol Si m-3", "Silicate", minimum=0.0_rk, initial_value=9.0_rk)
        call self%register_state_variable(self%id_sis, "SIS", "mmol Si m-3", "Biogenic silica", minimum=0.0_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_det, "DET", "mmol N m-3", "Nitrogen detritus", minimum=0.0_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_detp, "DETP", "mmol P m-3", "Phosphorus detritus", minimum=0.0_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_oxy, "OXY", "ml l-1", "Oxygen", minimum=0.0_rk, initial_value=8.0_rk)
        call self%register_state_variable(self%id_fla, "FLA", "mmol N m-3", "Flagellates", minimum=0.01_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_dia, "DIA", "mmol N m-3", "Diatoms", minimum=0.01_rk, initial_value=0.1_rk)

        ! Initialize model dependencies
        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_temp, standard_variables%temperature)
    end subroutine initialize

    subroutine do_pelagic(self, _ARGUMENTS_DO_)
        class(type_imr_norwecom), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        ! Local state variables
        real(rk) :: temp, par
        real(rk) :: nit, pho, sil, sis, det, detp
        real(rk) :: dia, fla

        ! Temporary arrays
        real(rk) :: pmax, nit_lim, pho_lim, sil_lim, rad_lim, growth_lim
        real(rk) :: growth_dia, mort_dia, resp_dia
        real(rk) :: growth_fla, mort_fla, resp_fla
        real(rk) :: nit_dia, nit_fla, dia_nit, fla_nit, dia_det, fla_det, det_nit
        real(rk) :: pho_dia, pho_fla, dia_pho, fla_pho, dia_detp, fla_detp, detp_pho
        real(rk) :: sil_dia, dia_sis, sis_sil
        real(rk) :: dnit, dpho, dsil, dsis, ddet, ddetp, ddia, dfla

        _LOOP_BEGIN_

            ! Get local values of state variables
            _GET_(self%id_nit, nit)
            _GET_(self%id_pho, pho)
            _GET_(self%id_sil, sil)
            _GET_(self%id_sis, sis)
            _GET_(self%id_det, det)
            _GET_(self%id_detp, detp)
            _GET_(self%id_dia, dia)
            _GET_(self%id_fla, fla)
            _GET_(self%id_temp, temp)
            _GET_(self%id_par, par)

            ! Phytoplankton

            ! Diatoms
            pmax = self%AA1 * exp(self%AA2 * temp)
            nit_lim = slim(self%KN_DIA, pmax, nit)
            pho_lim = slim(self%KP_DIA, pmax, pho)
            sil_lim = slim(self%KS_DIA, pmax, sil)
            rad_lim = slim(self%KR_DIA, pmax, par)
            growth_lim = min(nit_lim, pho_lim, sil_lim, rad_lim)
            growth_dia = pmax * growth_lim * dia
            mort_dia = self%CC3 * dia
            resp_dia = self%AA5 * dia * exp(self%AA6 * temp)

            ! Flagellates
            pmax = self%AA3 * exp(self%AA4 * temp)
            nit_lim = slim(self%KN_FLA, pmax, nit)
            pho_lim = slim(self%KP_FLA, pmax, pho)
            rad_lim = slim(self%KR_FLA, pmax, par)
            growth_lim = min(nit_lim, pho_lim, rad_lim)
            growth_fla = pmax * growth_lim * fla
            mort_fla = self%CC3 * fla
            resp_fla = self%AA5 * fla * exp(self%AA6 * temp)

            ! Sinks and sources terms

            ! Nitrogen flows
            nit_dia = growth_dia
            nit_fla = growth_fla
            dia_nit = resp_dia + 0.1_rk * mort_dia
            fla_nit = resp_fla + 0.1_rk * mort_fla
            dia_det = 0.9_rk * mort_dia
            fla_det = 0.9_rk * mort_fla
            det_nit = self%CC4 * det

            ! Phosphorus flows
            pho_dia = self%P_N_RATIO * nit_dia
            pho_fla = self%P_N_RATIO * nit_fla
            dia_pho = self%P_N_RATIO * (resp_dia + 0.25_rk * mort_dia)
            fla_pho = self%P_N_RATIO * (resp_fla + 0.25_rk * mort_fla)
            dia_detp = self%P_N_RATIO * 0.75_rk * mort_dia
            fla_detp = self%P_N_RATIO * 0.75_rk * mort_fla
            detp_pho = 1.3_rk * self%CC4 * detp

            ! Silicon flows
            sil_dia = self%SI_N_RATIO * nit_dia
            dia_sis = self%SI_N_RATIO * (resp_dia + mort_dia)
            sis_sil = self%SCC4 * sis

            ! Rates of change
            dnit = dia_nit + fla_nit + det_nit - (nit_dia + nit_fla)
            dpho = dia_pho + fla_pho + detp_pho - (pho_dia + pho_fla)
            dsil = sis_sil - sil_dia
            dsis = dia_sis - sis_sil
            ddet = dia_det + fla_det - det_nit
            ddetp = dia_detp + fla_detp - detp_pho
            ddia = nit_dia - dia_nit
            dfla = nit_fla - fla_nit

            ! Update state variables
            _SET_ODE_(self%id_nit, dnit)
            _SET_ODE_(self%id_pho, dpho)
            _SET_ODE_(self%id_sil, dsil)
            _SET_ODE_(self%id_sis, dsis)
            _SET_ODE_(self%id_det, ddet)
            _SET_ODE_(self%id_detp, ddetp)
            _SET_ODE_(self%id_dia, ddia)
            _SET_ODE_(self%id_fla, dfla)
        
        _LOOP_END_

    contains
        
        function slim(alpha, pmax, s) result(lim)
            !! Calculates the Michaelis-Menten limitation term
            real(rk), intent(in) :: alpha !! Substrate affinity
            real(rk), intent(in) :: pmax !! Production maximum
            real(rk), intent(in) :: s !! Substrate concentration
            
            real(rk) :: k, lim

            k = pmax / alpha
            if (s .lt. 0.0_rk) then
                lim = 0.0_rk
            else
                lim = s / (k + s)
            end if
        end function slim

    end subroutine do_pelagic

    subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
        !! Calculates sinking speed for particular matter
        class(type_imr_norwecom), intent(in) :: self
        _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

        real(rk) :: sil
        real(rk) :: sr_dia

        _LOOP_BEGIN_

            _GET_(self%id_sil, sil)

            if (sil .lt. self%SIB) then
                sr_dia = self%SR_DIAMAX
            else
                sr_dia = self%SR_DIAMIN + ((self%SR_DIAMAX - self%SR_DIAMIN) / sil)
            end if

            _SET_VERTICAL_MOVEMENT_(self%id_dia, -sr_dia)
            _SET_VERTICAL_MOVEMENT_(self%id_det, -self%SR_DET)
            _SET_VERTICAL_MOVEMENT_(self%id_fla, -self%SR_FLA)
            _SET_VERTICAL_MOVEMENT_(self%id_sis, -self%SR_SIS)

        _LOOP_END_
    end subroutine get_vertical_movement

end module imr_norwecom