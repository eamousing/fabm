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
        real(rk) :: C_N_RATIO !! Carbon to nitrogen ratio (mmol C (mmol N)-1)
        real(rk) :: O_N_RATIO !! Oxygen to nitrogen ratio (mg O (mmol N)-1)
        real(rk) :: MIN_PHYTO !! Minimum phytoplankton concentration (mmol N m-3)
        real(rk) :: PI_MES_DIA !! Mesozooplankton preference for diatoms (0-1)
        real(rk) :: PI_MES_MIC !! Mesozooplankton preference for microzooplankton (0-1)
        real(rk) :: PI_MES_DET !! Mesozooplankton preference for detritus (0-1)
        real(rk) :: PI_MIC_FLA !! Microzooplankton preference for flagellates (0-1)
        real(rk) :: PI_MIC_DET !! Microzooplankton preference for detritus (0-1)
        real(rk) :: K3 !! Half-saturation constant for zooplankton ingestion (mmol m-3)
        real(rk) :: K6 !! Half-saturation constant for zooplankton loss (mmol m-3)
        real(rk) :: BETA !! Zooplankton assimilation efficiency (0-1)
        real(rk) :: DELTA !! Fraction of loss from zooplankton to detritus (0-1)
        real(rk) :: MJU2 !! Microzooplankton maximum loss rate (s-1)
        real(rk) :: MJU3 !! Mesozooplankton maximum loss rate (s-1)
        real(rk) :: G_MIC !! Microzooplankton maximum growth rate (s-1)
        real(rk) :: G_MES !! Mesozooplankton maximum growth rate (s-1)
        real(rk) :: Q10 !! Temperature coefficient for zooplankton growth

        ! Register state variable identifiers
        type(type_state_variable_id) :: id_nit !! Nitrate (mmol N m-3)
        type(type_state_variable_id) :: id_pho !! Phosphate (mmol P m-3)
        type(type_state_variable_id) :: id_sil !! Silicate (mmol Si m-3)
        type(type_state_variable_id) :: id_sis !! Biogenic silica (mmol Si m-3)
        type(type_state_variable_id) :: id_det !! Nitrogen detritus (mmol N m-3)
        type(type_state_variable_id) :: id_detp !! Phosphorus detritus (mmol P m-3)
        type(type_state_variable_id) :: id_oxy !! Oxygen (ml l-1)
        type(type_state_variable_id) :: id_fla !! Flagellates (mmol N m-3)
        type(type_state_variable_id) :: id_dia !! Diatoms (mmol N m-3)
        type(type_state_variable_id) :: id_mic !! Microzooplankton (mmol N m-3)
        type(type_state_variable_id) :: id_mes !! Mesozooplankton (mmol N m-3)

        ! Register diagnostic variable identifiers
        type(type_diagnostic_variable_id) :: id_dia_gpp !! Diatom gross primary production (gC m-3 s-1)
        type(type_diagnostic_variable_id) :: id_dia_npp !! Diatom net primary production (gC m-3 s-1)
        type(type_diagnostic_variable_id) :: id_fla_gpp !! Flagellate gross primary production (gC m-3 s-1)
        type(type_diagnostic_variable_id) :: id_fla_npp !! Flagellate net primary production (gC m-3 s-1)
        type(type_diagnostic_variable_id) :: id_gpp !! Total gross primary production (gC m-3 s-1)
        type(type_diagnostic_variable_id) :: id_npp !! Total net primary production (gC m-3 s-1)

        ! Register dependencies identifiers
        type(type_dependency_id) :: id_par !! Photoactive radiation (m2 (W s)-1)
        type(type_dependency_id) :: id_temp !! Temperature (degC)
        type(type_dependency_id) :: id_sal !! Salinity
        type(type_dependency_id) :: id_dens !! Density (kg m-3)

    contains

        procedure :: initialize
        procedure :: do_surface
        procedure :: do
        procedure :: get_vertical_movement
    end type

contains

    subroutine initialize(self, configunit)
        class(type_imr_norwecom), intent(inout), target :: self
        integer, intent(in) :: configunit

        ! Initialize model parameters
        call self%get_parameter(self%AA1, "AA1", "s-1", "Diatom production maximum at 0 degC", default=1.53e-5_rk) ! Checked
        call self%get_parameter(self%AA2, "AA2", "degC-1", "Diatom temperature dependent pmax", default=0.063_rk) ! Checked
        call self%get_parameter(self%AA3, "AA3", "s-1", "Flagellate production maximum at 0 degC", default=1.02e-5_rk) ! Checked
        call self%get_parameter(self%AA4, "AA4", "degC-1", "Flagellate temperature dependent pmax", default=0.063_rk) ! Checked
        call self%get_parameter(self%AA5, "AA5", "s-1", "Phytoplankton respiration rate at 0 degC", default=8.05e-8_rk) ! Checked
        call self%get_parameter(self%AA6, "AA6", "degC-1", "Phytoplankton respiration rate temperature dependence", default=0.07_rk) ! Checked
        call self%get_parameter(self%CC3, "CC3", "s-1", "Phytoplankton death rate", default=1.6e-6_rk) ! Checked
        call self%get_parameter(self%CC4, "CC4", "s-1", "Detritus decomposition rate", default=1.52e-7_rk) ! Checked
        call self%get_parameter(self%SCC4, "SCC4", "s-1", "Biogenic silica decomposition rate", default=6.43e-8_rk) ! Checked
        call self%get_parameter(self%KN_DIA, "KN_DIA", "(mmol N)-1 s-1", "Diatom affinity for nitrate", default=1.736e-5_rk)
        call self%get_parameter(self%KN_FLA, "KN_FLA", "(mmol N)-1 s-1", "Flagellate affinity for nitrate", default=1.505e-5_rk)
        call self%get_parameter(self%KP_DIA, "KP_DIA", "(mmol P)-1 s-1", "Diatom affinity for phosphate", default=2.697e-4_rk)
        call self%get_parameter(self%KP_FLA, "KP_FLA", "(mmol P)-1 s-1", "Flagellate affinity for phosphate", default=2.5e-4_rk)
        call self%get_parameter(self%KR_DIA, "KR_DIA", "m2 (W s)-1", "Diatom affinity for radition", default=1.5e-6_rk)
        call self%get_parameter(self%KR_FLA, "KR_FLA", "m2 (W s)-1", "Flagellate affinity for radition", default=4.6e-7_rk)
        call self%get_parameter(self%KS_DIA, "KS_DIA", "(mmol Si)-1 s-1", "Diatom affinity for silicate", default=2.546e-5_rk)
        call self%get_parameter(self%SR_DET, "SR_DET", "m s-1", "Detritus sinking rate", default=-3.472e-5_rk) ! Checked
        call self%get_parameter(self%SR_DIAMAX, "SR_DIAMAX", "m s-1", "Diatom maximum sinking rate", default=-3.472e-5_rk) ! Checked
        call self%get_parameter(self%SR_DIAMIN, "SR_DIAMIN", "m s-1", "Diatom minimum sinking rate", default=-3.472e-6_rk) ! Checked
        call self%get_parameter(self%SR_FLA, "SR_FLA", "m s-1", "Flagellate sinking rate", default=-2.894e-6_rk) ! Checked
        call self%get_parameter(self%SR_SIS, "SR_SIS", "m s-1", "Biogenic silica sinking rate", default=-3.472e-5_rk) !Checked
        call self%get_parameter(self%SIB, "SIB", "mmol Si m-3", "Concentration of SIL when max sinking speed of DIA", default=1.0_rk) ! Checked
        call self%get_parameter(self%P_N_RATIO, "P_N_RATIO", "mmol P (mmol N)-1", "Phosphorus to nitrogen ratio", default=0.0625_rk)
        call self%get_parameter(self%SI_N_RATIO, "SI_N_RATIO", "mmol Si (mmol N)-1", "Silicon to nitrogen ratio", default=0.875_rk)
        call self%get_parameter(self%C_N_RATIO, "C_N_RATIO", "mmol C (mmol N)-1", "Carbon to nitrogen ratio", default=6.625_rk)
        call self%get_parameter(self%O_N_RATIO, "O_N_RATIO", "ml O (mmol N)-1", "Oxygen to nitrogen ratio", default=0.2_rk) ! TODO: check this!
        call self%get_parameter(self%MIN_PHYTO, "MIN_PHYTO", "mmol N m-3", "Minimum phytoplankton concentration", default=0.1_rk)
        call self%get_parameter(self%PI_MES_DIA, "PI_MES_DIA", "[0:1]", "Mesozooplankton preference for diatoms", default=0.450_rk)
        call self%get_parameter(self%PI_MES_MIC, "PI_MES_MIC", "[0:1]", "Mesozooplankton preference for microzooplankton", default=0.275_rk)
        call self%get_parameter(self%PI_MES_DET, "PI_MES_DET", "[0:1]", "Mesozooplankton preference for detritus", default=0.275_rk)
        call self%get_parameter(self%PI_MIC_FLA, "PI_MIC_FLA", "[0:1]", "Microzooplankton preference for flagellates", default=0.633_rk)
        call self%get_parameter(self%PI_MIC_DET, "PI_MIC_DET", "[0:1]", "Microzooplankton preference for detritus", default=0.367_rk)
        call self%get_parameter(self%K3, "K3", "mmol m-3", "Half-saturation constant for zooplankton ingestion", default=1.0_rk)
        call self%get_parameter(self%K6, "K6", "mmol m-3", "Half-saturation constant for zooplankton loss", default=0.2_rk)
        call self%get_parameter(self%BETA, "BETA", "[0:1]", "Zooplankton assimilation efficiency", default=0.75_rk)
        call self%get_parameter(self%DELTA, "DELTA", "[0:1]", "Fraction of loss from zooplankton to detritus", default=0.5_rk)
        call self%get_parameter(self%MJU2, "MJU2", "s-1", "Microzooplankton maximum loss rate", default=2.315e-6_rk)
        call self%get_parameter(self%MJU3, "MJU3", "s-1", "Mesozooplankton maximum loss rate", default=2.315e-6_rk)
        call self%get_parameter(self%G_MIC, "G_MIC", "s-1", "Microzooplankton maximum growth rate", default=4.63e-6_rk)
        call self%get_parameter(self%G_MES, "G_MES", "s-1", "Mesozooplankton maximum growth rate", default=4.63e-6_rk)
        call self%get_parameter(self%Q10, "Q10", " ", "Temperature coefficient for zooplankton growth", default=1.5_rk)

        ! Initialize state variables
        call self%register_state_variable(self%id_nit, "NIT", "mmol N m-3", "Nitrate", minimum=0.0_rk, initial_value=10.0_rk)
        call self%register_state_variable(self%id_pho, "PHO", "mmol P m-3", "Phosphate", minimum=0.0_rk, initial_value=0.7_rk)
        call self%register_state_variable(self%id_sil, "SIL", "mmol Si m-3", "Silicate", minimum=0.0_rk, initial_value=9.0_rk)
        call self%register_state_variable(self%id_sis, "SIS", "mmol Si m-3", "Biogenic silica", minimum=0.0_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_det, "DET", "mmol N m-3", "Nitrogen detritus", minimum=0.0_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_detp, "DETP", "mmol P m-3", "Phosphorus detritus", minimum=0.0_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_oxy, "OXY", "ml l-1", "Oxygen", minimum=0.0_rk, initial_value=8.0_rk)
        call self%register_state_variable(self%id_fla, "FLA", "mmol N m-3", "Flagellates", minimum=0.0001_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_dia, "DIA", "mmol N m-3", "Diatoms", minimum=0.0001_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_mic, "MIC", "mmol N m-3", "Microzooplankton", minimum=0.0001_rk, initial_value=0.1_rk)
        call self%register_state_variable(self%id_mes, "MES", "mmol N m-3", "Mesozooplankton", minimum=0.0001_rk, initial_value=0.1_rk)

        ! Initialze diagnostic variables
        call self%register_diagnostic_variable(self%id_gpp, "GPP", "gC m-3 s-1", "Gross primary production rate", output=output_time_step_averaged)
        call self%register_diagnostic_variable(self%id_npp, "NPP", "gC m-3 s-1", "Net primary production rate", output=output_time_step_averaged)
        call self%register_diagnostic_variable(self%id_dia_gpp, "DIA_GPP", "gC m-3 s-1", "Diatom gross primary production rate", output=output_time_step_averaged)
        call self%register_diagnostic_variable(self%id_dia_npp, "DIA_NPP", "gC m-3 s-1", "Diatom net primary production rate", output=output_time_step_averaged)
        call self%register_diagnostic_variable(self%id_fla_gpp, "FLA_GPP", "gC m-3 s-1", "Flagellate gross primary production rate", output=output_time_step_averaged)
        call self%register_diagnostic_variable(self%id_fla_npp, "FLA_NPP", "gC m-3 s-1", "Flagellate net primary production rate", output=output_time_step_averaged)

        ! Initialize model dependencies
        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_temp, standard_variables%temperature)
        call self%register_dependency(self%id_sal, standard_variables%practical_salinity)
        call self%register_dependency(self%id_dens, standard_variables%density)
    end subroutine initialize

    subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
        class(type_imr_norwecom), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_SURFACE_

        real(rk), parameter :: a1 = -173.4292_rk
        real(rk), parameter :: a2 = 249.6339_rk
        real(rk), parameter :: a3 = 143.3483_rk
        real(rk), parameter :: a4 = -21.8492_rk
        real(rk), parameter :: b1 = -0.033096_rk
        real(rk), parameter :: b2 = 0.014259_rk
        real(rk), parameter :: b3 = -0.001700_rk
        real(rk), parameter :: kelvin = 273.16_rk
        real(rk), parameter :: pvel = 5.0_rk / 86400.0_rk
        real(rk) :: osat_weiss, oxy_flux, tk
        real(rk) :: temp, sal, oxy
        
        _SURFACE_LOOP_BEGIN_
        
            _GET_(self%id_temp, temp)
            _GET_(self%id_sal, sal)
            _GET_(self%id_oxy, oxy)
            
            ! Oxygen saturation state           
            tk = temp + kelvin
            osat_weiss = a1 + a2*100.0_rk/tk + a3*log(0.01_rk*tk) + a4*0.01_rk*tk
            osat_weiss = osat_weiss + sal*(b1 + b2*0.01_rk*tk + b3*(0.01_rk*tk)**2)
            osat_weiss = exp(osat_weiss)
            
            ! Surface oxygen flux
            oxy_flux = pvel*(osat_weiss - oxy)

            ! Update state variables
            _SET_SURFACE_EXCHANGE_(self%id_oxy, oxy_flux)

        _SURFACE_LOOP_END_
    end subroutine do_surface

    subroutine do(self, _ARGUMENTS_DO_)
        class(type_imr_norwecom), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        ! Local state variables
        real(rk) :: temp, par
        real(rk) :: nit, pho, sil, sis, det, detp
        real(rk) :: dia, fla
        real(rk) :: mic, mes
        real(rk) :: gpp, npp, gpp_dia, npp_dia, gpp_fla, npp_fla

        ! Temporary arrays
        real(rk) :: pmax, nit_lim, pho_lim, sil_lim, rad_lim, growth_lim
        real(rk) :: growth_dia, mort_dia, resp_dia
        real(rk) :: growth_fla, mort_fla, resp_fla
        real(rk) :: nit_dia, nit_fla, dia_nit, fla_nit, dia_det, fla_det, det_nit
        real(rk) :: pho_dia, pho_fla, dia_pho, fla_pho, dia_detp, fla_detp, detp_pho
        real(rk) :: sil_dia, dia_sis, sis_sil
        real(rk) :: oxy_dia, oxy_fla, dia_oxy, fla_oxy, oxy_det
        real(rk) :: dnit, dpho, dsil, dsis, ddet, ddetp, ddia, dfla, doxy
        real(rk) :: denum, tmp
        real(rk) :: p_fla_mic, p_det_mic, g_fla_mic, g_det_mic, mort_mic
        real(rk) :: p_dia_mes, p_mic_mes, p_det_mes, g_dia_mes, g_mic_mes, g_det_mes, mort_mes
        real(rk) :: mic_mes, mic_det, fla_mic, det_mic, mic_nit 
        real(rk) :: mes_det, mes_nit, dia_mes, det_mes
        real(rk) :: mic_pho, mic_detp, mes_pho, mes_detp
        real(rk) :: mes_sis
        real(rk) :: oxy_mic, oxy_mes

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
            _GET_(self%id_mic, mic)
            _GET_(self%id_mes, mes)

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
            if (dia .le. self%MIN_PHYTO) then
                mort_dia = 0.0_rk
                resp_dia = 0.0_rk
            end if

            ! Flagellates
            pmax = self%AA3 * exp(self%AA4 * temp)
            nit_lim = slim(self%KN_FLA, pmax, nit)
            pho_lim = slim(self%KP_FLA, pmax, pho)
            rad_lim = slim(self%KR_FLA, pmax, par)
            growth_lim = min(nit_lim, pho_lim, rad_lim)
            growth_fla = pmax * growth_lim * fla
            mort_fla = self%CC3 * fla
            resp_fla = self%AA5 * fla * exp(self%AA6 * temp)
            if (fla .le. self%MIN_PHYTO) then
                mort_fla = 0.0_rk
                resp_fla = 0.0_rk
            end if

            ! Gross and net primary production
            gpp_dia = growth_dia*self%C_N_RATIO*12.01_rk*1.0e-3_rk
            npp_dia = (growth_dia - resp_dia)*self%C_N_RATIO*12.01_rk*1.0e-3_rk
            gpp_fla = growth_fla*self%C_N_RATIO*12.01_rk*1.0e-3_rk
            npp_fla = (growth_fla - resp_fla)*self%C_N_RATIO*12.01_rk*1.0e-3_rk
            gpp = gpp_dia + gpp_fla
            npp = npp_dia + npp_fla

            ! Zooplankton terms

            ! Microzooplankton
            denum = self%PI_MIC_FLA*fla + self%PI_MIC_DET*det + 1.0e-8_rk
            p_fla_mic = (self%PI_MIC_FLA*fla)/denum
            p_det_mic = (self%PI_MIC_DET*det)/denum
            tmp = self%Q10**((temp-10.0_rk)/10.0_rk)*self%G_MIC/(self%K3 + p_fla_mic*fla + p_det_mic*det)
            g_fla_mic = tmp*p_fla_mic*fla*mic
            g_det_mic = tmp*p_det_mic*det*mic
            mort_mic = self%MJU2*(mic/(self%K6 + mic))*mic

            ! Mezozooplankton
            denum = self%PI_MES_DIA*dia + self%PI_MES_MIC*mic + self%PI_MES_DET*det + 1.0e-8_rk
            p_dia_mes = (self%PI_MES_DIA*dia)/denum
            p_mic_mes = (self%PI_MES_MIC*mic)/denum
            p_det_mic = (self%PI_MES_DET*det)/denum
            tmp = self%Q10**((temp-10.0_rk)/10.0_rk)*self%G_MES/(self%K3 + p_dia_mes*dia + p_mic_mes*mic + p_det_mes*det)
            g_dia_mes = tmp*p_dia_mes*dia
            g_mic_mes = tmp*p_mic_mes*mic
            g_det_mes = tmp*p_det_mes*det
            mort_mes = self%MJU3*(mes/(self%K6 + mes))*mes

            ! Sinks and sources terms

            ! Nitrogen flows
            nit_dia = growth_dia
            nit_fla = growth_fla
            dia_nit = resp_dia + 0.1_rk * mort_dia
            fla_nit = resp_fla + 0.1_rk * mort_fla
            dia_det = 0.9_rk * mort_dia
            fla_det = 0.9_rk * mort_fla
            det_nit = self%CC4 * det
            dia_mes = g_dia_mes
            mic_mes = g_mic_mes
            det_mes = g_det_mes
            fla_mic = g_fla_mic
            det_mic = g_det_mic
            mic_det = self%DELTA*mort_mic + (1.0_rk - self%BETA)*(g_fla_mic + g_det_mic)
            mic_nit = (1.0_rk - self%DELTA)*mort_mic
            mes_det = self%DELTA*mort_mes + (1.0_rk - self%BETA)*(dia_mes + mic_mes + det_mes)
            mes_nit = (1.0_rk - self%DELTA)*mort_mes

            ! Phosphorus flows
            pho_dia = self%P_N_RATIO * nit_dia
            pho_fla = self%P_N_RATIO * nit_fla
            dia_pho = self%P_N_RATIO * (resp_dia + 0.25_rk * mort_dia)
            fla_pho = self%P_N_RATIO * (resp_fla + 0.25_rk * mort_fla)
            dia_detp = self%P_N_RATIO * 0.75_rk * mort_dia
            fla_detp = self%P_N_RATIO * 0.75_rk * mort_fla
            mic_pho = self%P_N_RATIO * mic_nit
            mic_detp = self%P_N_RATIO * mic_det
            mes_pho = self%P_N_RATIO * mes_nit
            mes_detp = self%P_N_RATIO * mes_det
            detp_pho = 1.3_rk * self%CC4 * detp

            ! Silicon flows
            sil_dia = self%SI_N_RATIO*nit_dia
            dia_sis = self%SI_N_RATIO*(resp_dia + mort_dia)
            mes_sis = self%SI_N_RATIO*mes_det
            sis_sil = self%SCC4*sis

            ! Oxygen flows
            dia_oxy = nit_dia * self%O_N_RATIO
            fla_oxy = nit_fla * self%O_N_RATIO
            oxy_dia = dia_nit * self%O_N_RATIO
            oxy_fla = fla_nit * self%O_N_RATIO
            oxy_det = det_nit * self%O_N_RATIO
            oxy_mic = mic_nit * self%O_N_RATIO
            oxy_mes = mes_nit * self%O_N_RATIO

            ! Rates of change
            dnit = dia_nit + fla_nit + det_nit + mic_nit + mes_nit - (nit_dia + nit_fla)
            dpho = dia_pho + fla_pho + detp_pho + mic_pho + mes_pho - (pho_dia + pho_fla)
            dsil = sis_sil - sil_dia
            dsis = dia_sis + mes_sis - sis_sil
            ddet = dia_det + fla_det + mic_det + mes_det - det_nit
            ddetp = dia_detp + fla_detp + mic_detp + mes_detp - detp_pho
            ddia = nit_dia - (dia_nit + dia_det + dia_mes)
            dfla = nit_fla - (fla_nit + fla_det + fla_mic)
            doxy = dia_oxy + fla_oxy - (oxy_dia + oxy_fla + oxy_det + oxy_mic + oxy_mes)

            ! Update state variables
            _SET_ODE_(self%id_nit, dnit)
            _SET_ODE_(self%id_pho, dpho)
            _SET_ODE_(self%id_sil, dsil)
            _SET_ODE_(self%id_sis, dsis)
            _SET_ODE_(self%id_det, ddet)
            _SET_ODE_(self%id_detp, ddetp)
            _SET_ODE_(self%id_dia, ddia)
            _SET_ODE_(self%id_fla, dfla)
            _SET_ODE_(self%id_oxy, doxy)

            ! Update diagnostic variables
            _SET_DIAGNOSTIC_(self%id_dia_gpp, gpp_dia)
            _SET_DIAGNOSTIC_(self%id_dia_npp, npp_dia)
            _SET_DIAGNOSTIC_(self%id_fla_gpp, gpp_fla)
            _SET_DIAGNOSTIC_(self%id_fla_npp, npp_fla)
            _SET_DIAGNOSTIC_(self%id_gpp, gpp)
            _SET_DIAGNOSTIC_(self%id_npp, npp)
        
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

    end subroutine do

    subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
        !! Calculates sinking speed for particular matter
        class(type_imr_norwecom), intent(in) :: self
        _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

        real(rk) :: sil, dia, det, fla, detp, sis
        real(rk) :: sr_dia, sr_det, sr_fla, sr_detp, sr_sis

        _LOOP_BEGIN_

            _GET_(self%id_sil, sil)
            _GET_(self%id_sis, sis)
            _GET_(self%id_dia, dia)
            _GET_(self%id_fla, fla)
            _GET_(self%id_det, det)
            _GET_(self%id_detp, detp)

            if (sil .lt. self%SIB) then
                sr_dia = self%SR_DIAMAX
            else
                sr_dia = self%SR_DIAMIN + ((self%SR_DIAMAX - self%SR_DIAMIN) / sil)
            end if

            sr_det = self%SR_DET
            sr_fla = self%SR_FLA
            sr_detp = self%SR_DET
            sr_sis = self%SR_SIS

            if (dia .lt. 0.1_rk) sr_dia = 0.0_rk
            if (fla .lt. 0.1_rk) sr_fla = 0.0_rk
            if (det .lt. 0.1_rk) sr_det = 0.0_rk 
            if (detp .lt. 0.1_rk * self%P_N_RATIO) sr_detp = 0.0_rk
            if (sis .lt. 0.1_rk * self%SI_N_RATIO) sr_sis = 0.0_rk

            _SET_VERTICAL_MOVEMENT_(self%id_dia, sr_dia)
            _SET_VERTICAL_MOVEMENT_(self%id_det, sr_det)
            _SET_VERTICAL_MOVEMENT_(self%id_fla, sr_fla)
            _SET_VERTICAL_MOVEMENT_(self%id_sis, sr_sis)
            _SET_VERTICAL_MOVEMENT_(self%id_detp, sr_detp)

        _LOOP_END_
    end subroutine get_vertical_movement

end module imr_norwecom