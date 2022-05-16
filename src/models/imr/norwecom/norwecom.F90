#include "fabm_driver.h"

module imr_norwecom

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_imr_norwecom
      ! Parameters
      real(rk) :: a(6)
      real(rk) :: ad(4)
      real(rk) :: af(3)
      real(rk) :: cc(4)
      real(rk) :: diamin
      real(rk) :: flamin
      real(rk) :: srdia_min
      real(rk) :: srdia_max
      real(rk) :: sib
      real(rk) :: scc(7)

      ! State variables
      type(type_state_variable_id) :: id_nit
      type(type_state_variable_id) :: id_pho
      type(type_state_variable_id) :: id_sil
      type(type_state_variable_id) :: id_det
      type(type_state_variable_id) :: id_sis
      type(type_state_variable_id) :: id_dia
      type(type_state_variable_id) :: id_fla
      type(type_state_variable_id) :: id_oxy

      ! Diagnostic variables

      ! Dependencies
      type(type_dependency_id) :: id_temp
      type(type_dependency_id) :: id_salt
      type(type_dependency_id) :: id_dens
      type(type_dependency_id) :: id_par
   contains
      procedure :: initialize
      procedure :: do_surface
      procedure :: do
      ! procedure :: do_bottom
      procedure :: get_vertical_movement
   end type

contains

   subroutine initialize(self, configunit)
      class(type_imr_norwecom), intent(inout), target :: self
      integer, intent(in) :: configunit
      
      ! Initialize parameters
      call self%get_parameter(self%a(1), "a(1)", "s-1", "Diatoms production maximum at 0 degC", default=1.53e-5_rk)
      call self%get_parameter(self%a(2), "a(2)", "degC-1", "Diatoms temp dependent Pmax", default=0.063_rk)
      call self%get_parameter(self%a(3), "a(3)", "s-1", "Flagellates production max at 0 degC", default=1.02e-5_rk)
      call self%get_parameter(self%a(4), "a(4)", "degC-1", "Flagellate temp dependent Pmax", default=0.063_rk)
      call self%get_parameter(self%a(5), "a(5)", "s-1", "Phyto metabolic loss rate at 0 degC", default=8.05e-7_rk)
      call self%get_parameter(self%a(6), "a(6)", "degC-1", "Phyto metabolic loss rate temp dependence", default=0.07_rk)
      call self%get_parameter(self%ad(1), "ad(1)", "m2 uE-1", "Diatoms growth affinity for irradiance", default=3.6e-7_rk)
      call self%get_parameter(self%ad(2), "ad(2)", "s-1 uM-1", "Diatoms growth affinity for nitrate", default=1.7e-5_rk)
      call self%get_parameter(self%ad(3), "ad(3)", "s-1 uM-1", "Diatoms growth affinity for phosphate", default=2.7e-4_rk)
      call self%get_parameter(self%ad(4), "ad(4)", "s-1 uM-1", "Diatoms growth affinity for silicate", default=2.5e-5_rk)
      call self%get_parameter(self%af(1), "af(1)", "m2 uE-1", "Flagellates growth affinity for irradience", default=1.1e-7_rk)
      call self%get_parameter(self%af(2), "af(2)", "s-1 uM-1", "Flagellates growth affinity for nitrate", default=1.5e-5_rk)
      call self%get_parameter(self%af(3), "af(3)", "s-1 uM-1", "Flagellates growth affinity for phosphate", default=2.5e-4_rk)
      call self%get_parameter(self%cc(1), "cc(1)", "mgP mgN-1", "Intercellular P/N ratio", default=0.138_rk)
      call self%get_parameter(self%cc(2), "cc(2)", "mgSi mgN-1", "Intercellular Si/N ratio", default=1.75_rk)
      call self%get_parameter(self%cc(3), "cc(3)", "s-1", "Phyto death rate", default=1.6e-6_rk)
      call self%get_parameter(self%cc(4), "cc(4)", "s-1", "Detritus decomposition rate", default=1.52e-7_rk)
      call self%get_parameter(self%diamin, "diamin", "mgN m-3", "Minimum diatoms concentration", default=0.1_rk)
      call self%get_parameter(self%flamin, "flamin", "mgN m-3", "Minimum flagellates concentration", default=0.1_rk)
      call self%get_parameter(self%srdia_min, "srdia_min", "m s-1", "Diatoms minimum sinking rate", default=3.47e-6_rk)
      call self%get_parameter(self%srdia_max, "srdia_max", "m s-1", "Diatoms maximum sinking rate", default=3.47e-5_rk)
      call self%get_parameter(self%sib, "sib", "uM", "Conc. of silicate where diatom max sinking speed", default=1.0_rk)
      call self%get_parameter(self%scc(1), "scc(1)", "mgO mgN-1", "O/N consumption ratio", default=19.71_rk)
      call self%get_parameter(self%scc(4), "scc(4)", "s-1", "Biogenic silica decomposition rate", default=1.45e-8_rk)

      ! Initialize state variables
      call self%register_state_variable(self%id_nit, "nit", "mgN m-3", "Nitrate concentration", &
         minimum=0.0_rk, initial_value=10.0_rk)
      call self%register_state_variable(self%id_pho, "pho", "mgP m-3", "Phosphate concentration", &
         minimum=0.0_rk, initial_value=2.0_rk)
      call self%register_state_variable(self%id_sil, "sil", "mgSi m-3", "Silicate concentrations", &
         minimum=0.0_rk, initial_value=10.0_rk)
      call self%register_state_variable(self%id_det, "det", "mgN m-3", "Detritus concentration", &
         minimum=0.0_rk, initial_value=0.1_rk, vertical_movement=-3.47e-5_rk)
      call self%register_state_variable(self%id_sis, "sis", "mgSi m-3", "Biogenic silica concentration", &
         minimum=0.0_rk, initial_value=0.1_rk, vertical_movement=-3.47e-5_rk)
      call self%register_state_variable(self%id_dia, "dia", "mgN m-3", "Diatoms concentration", &
         minimum=0.0001_rk, initial_value=0.1_rk, vertical_movement=-2.89e-6_rk) ! -3.47e-5_rk)
      call self%register_state_variable(self%id_fla, "fla", "mgN m-3", "Flagellates concentration", &
         minimum=0.0001_rk, initial_value=0.1_rk, vertical_movement=-2.89e-6_rk)
      call self%register_state_variable(self%id_oxy, "oxy", "mgO l-1", "Oxygen concentration", &
         minimum=0.0_rk, initial_value=10.0_rk)

      ! Initialize dependencies
      call self%register_dependency(self%id_temp, standard_variables%temperature)
      call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
      call self%register_dependency(self%id_dens, standard_variables%density)
      call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)

      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_dia)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_fla)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_det)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_nit)
   end subroutine initialize

   subroutine do_surface(self, _ARGUMENTS_DO_SURFACE_)
      !! Oxygen air-water flux
      !!
      class(type_imr_norwecom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_SURFACE_

      real(rk) :: temp, salt, dens, oxy
      real(rk) :: tempk, pvel, osat, doxy

      pvel = 5.0

      _SURFACE_LOOP_BEGIN_

      _GET_(self%id_temp, temp)
      _GET_(self%id_salt, salt)
      _GET_(self%id_dens, dens)
      _GET_(self%id_oxy, oxy)

      tempk = (temp + 273.15) / 100.0
      osat = exp(-173.9894 + 255.5907 / tempk + 146.4813 * log(tempk) - 22.2040 * tempk + &
         salt * (-0.037376 + 0.016504 * tempk - 0.0020564 * tempk * tempk)) ! mol kg-1
      osat = osat * dens ! mol m-3
      osat = osat * 32.0 * 1e-6 ! mg m-3
      doxy = (pvel / 86400.0) * (osat - oxy)

      _ADD_SURFACE_FLUX_(self%id_oxy, doxy)
      
      _SURFACE_LOOP_END_
   end subroutine do_surface

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_imr_norwecom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: temp, dia, fla, nit, pho, sil, det, sis, oxy, par, rad
      real(rk) :: umax, rad_lim, nit_lim, pho_lim, sil_lim
      real(rk) :: prod_dia, resp_dia, mort_dia, prod_fla, resp_fla, mort_fla
      real(rk) :: dnit, dpho, dsil, ddet, dsis, ddia, dfla, doxy

      _LOOP_BEGIN_

      _GET_(self%id_temp, temp)
      _GET_(self%id_par, par)
      _GET_(self%id_dia, dia)
      _GET_(self%id_fla, fla)
      _GET_(self%id_nit, nit)
      _GET_(self%id_pho, pho)
      _GET_(self%id_sil, sil)
      _GET_(self%id_det, det)
      _GET_(self%id_sis, sis)
      _GET_(self%id_oxy, oxy)

      rad = par / 0.217 ! W m-2 to uE m-2 s-1

      ! Diatoms
      umax = self%a(1) * exp(self%a(2) * temp)
      rad_lim = slim(self%ad(1), umax, rad)
      nit_lim = slim(self%ad(2), umax, nit)
      pho_lim = slim(self%ad(3), umax, pho) 
      sil_lim = slim(self%ad(4), umax, sil)
      prod_dia = umax * rad_lim * min(nit_lim, pho_lim, sil_lim) * dia
      resp_dia = self%a(5) * dia * exp(self%a(6) * temp)
      mort_dia = self%cc(3) * dia

      if (dia < self%diamin) then
         resp_dia = 0.0
         mort_dia = 0.0
      end if

      ! Flagellates
      umax = self%a(3) * exp(self%a(4) * temp)
      rad_lim = slim(self%af(1), umax, rad)
      nit_lim = slim(self%af(2), umax, nit) 
      pho_lim = slim(self%ad(3), umax, pho)
      prod_fla = umax * rad_lim * min(nit_lim, pho_lim) * fla
      resp_fla = self%a(5) * fla * exp(self%a(6) * temp)
      mort_fla = self%cc(3) * fla

      if (fla < self%flamin) then
         resp_fla = 0.0
         mort_fla = 0.0
      end if

      ! Fluxes
      dnit = resp_dia + resp_fla + self%cc(4) * det - (prod_dia + prod_fla)
      dpho = self%cc(1) * dnit
      dsil = self%scc(4) * sis - self%cc(2) * prod_dia
      ddet = mort_dia + mort_fla - self%cc(4) * det
      dsis = self%cc(2) * (resp_dia + mort_dia) - self%scc(4) * sis
      ddia = prod_dia - (resp_dia + mort_dia)
      dfla = prod_fla - (resp_fla + mort_fla)
      doxy = self%scc(1) * (prod_dia + prod_fla - (resp_dia + resp_fla + self%cc(4) * det))

      ! Update FABM
      _ADD_SOURCE_(self%id_nit, dnit)
      _ADD_SOURCE_(self%id_pho, dpho)
      _ADD_SOURCE_(self%id_sil, dsil)
      _ADD_SOURCE_(self%id_det, ddet)
      _ADD_SOURCE_(self%id_sis, dsis)
      _ADD_SOURCE_(self%id_dia, ddia)
      _ADD_SOURCE_(self%id_fla, dfla)
      _ADD_SOURCE_(self%id_oxy, doxy*1e-3) ! mg m-3 -> mg l-1

      _LOOP_END_
   
   contains

      real(rk) function slim(alpha, pmax, s)
         real(rk), intent(in) :: alpha
         real(rk), intent(in) :: pmax
         real(rk), intent(in) :: s

         if (s < 0.0) then
            slim = 0.0
         else
            slim = s / (s + (pmax / alpha))
         end if
      end function

   end subroutine do

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class(type_imr_norwecom), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_


   end subroutine do_bottom

   subroutine get_vertical_movement(self, _ARGUMENTS_GET_VERTICAL_MOVEMENT_)
      class(type_imr_norwecom), intent(in) :: self
      _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

      real(rk) :: dia, sil, vdia, csil

      csil = 28.09_rk

      _LOOP_BEGIN_

      _GET_(self%id_dia, dia)
      _GET_(self%id_sil, sil)

      if (sil / csil < self%sib) then
         vdia = self%srdia_max
      else
         vdia = self%srdia_min + (self%srdia_max - self%srdia_min) / (sil / csil)
      end if

      ! Stop sinking if concentration is too low
      if (dia < self%diamin) vdia = 0.0_rk

      _ADD_VERTICAL_VELOCITY_(self%id_dia, -1.0 * vdia)

      _LOOP_END_
   end subroutine get_vertical_movement

end module imr_norwecom
