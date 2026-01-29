MODULE circadian_mod
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CIRCADIAN_MOD'

PRIVATE

PUBLIC tdq10_factor, update_resp_nocturnal

CONTAINS

! *********************************************************************
! *********************************************************************
SUBROUTINE update_resp_nocturnal(ipar, tstar, tstar_ref, resp_ref,             &
                                 resp_fac, resp_p)

USE timestep_mod, ONLY : timestep

IMPLICIT NONE

REAL(KIND=real_jlslsm), PARAMETER :: ipar_min = 1.0E-6
!                           ! PAR threshold used to determine day/night.
REAL(KIND=real_jlslsm), PARAMETER :: resp_fac_min = 0.55
!                           ! Lower limit for respiration supression.
!                           ! Equivalent to 24h of Bruhn et al reduction.
REAL(KIND=real_jlslsm), PARAMETER :: resp_fac_max = 0.99
!                           ! Upper limit for respiration supression.  Must be
!                           ! less than 1.0 for ODE solver to work.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 ipar                                                                          &
!                           ! IN Incident PAR (W/m2).
,tstar
!                           ! IN Surface temperature (K).

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                       &
 tstar_ref                                                                     &
!                           ! INOUT Reference temperature (K).
,resp_ref                                                                      &
!                           ! INOUT Reference plant respiration rate
!                           ! (kg C/m2/sec).
,resp_fac                                                                      &
!                           ! INOUT Plant respiration nocturnal state
!                           ! (dimensionless).
,resp_p
                            ! INOUT Plant respiration rate
!                           ! (kg C/m2/sec).

REAL(KIND=real_jlslsm) :: q10_fac     ! WORK Respiration temperature adjustment.
REAL(KIND=real_jlslsm) :: k1, k2, k3  ! WORK ODE solver factors.


IF ( ipar < ipar_min ) THEN
  ! Night-time.
  ! Update the respiration multiplcation factor using Huen's ODE method.
  k1 = bruhn_ode(resp_fac)
  k2 = bruhn_ode(resp_fac + timestep*k1)
  k3 = bruhn_ode(resp_fac + timestep*(0.5*k1+0.5*k2))
  resp_fac = resp_fac + timestep*k3

  q10_fac = tdq10_factor(tstar, tstar_ref)

  ! Recalculate plant respiration relative to the reference respiration
  ! rate, which is the last value before sunset.
  resp_p = resp_ref * resp_fac * q10_fac
ELSE
  ! Day-time.
  ! Set (or reset) the multiplication factor and update the reference
  ! respiration rate.
  tstar_ref = tstar
  resp_fac = resp_fac_max
  resp_ref = resp_p
END IF

resp_fac = MIN(MAX(resp_fac, resp_fac_min), resp_fac_max)

END SUBROUTINE update_resp_nocturnal

! *********************************************************************
! *********************************************************************
FUNCTION bruhn_ode(resp_fac) RESULT(factor)

IMPLICIT NONE

REAL(KIND=real_jlslsm), PARAMETER :: bruhn_b = 0.54
REAL(KIND=real_jlslsm), PARAMETER :: bruhn_a = 0.08 / (3600.0**bruhn_b)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 resp_fac
!                           ! IN Plant respiration nocturnal state
!                           ! (dimensionless).

REAL(KIND=real_jlslsm) :: power
REAL(KIND=real_jlslsm) :: factor

power = (bruhn_b - 1.0) / bruhn_b
factor = -1.0 * bruhn_a * bruhn_b * ((1.0-resp_fac)/bruhn_a)**power

END FUNCTION bruhn_ode

! *********************************************************************
! *********************************************************************
FUNCTION tdq10_factor(tstar, tref) RESULT(factor)

IMPLICIT NONE

REAL(KIND=real_jlslsm), PARAMETER :: tjoel_a = 3.22
REAL(KIND=real_jlslsm), PARAMETER :: tjoel_b = 0.046

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 tstar,                                                                        &
!                           ! IN Temperature (K).
 tref
!                           ! IN Reference temperature (K).

REAL(KIND=real_jlslsm) ::                                                      &
 q10,                                                                          &
!                           ! WORK Temperature-dependent Q10 (dimensionless).
 power
!                           ! WORK Power (dimensionless).

REAL(KIND=real_jlslsm) :: factor
!                           ! OUT Temperature adjustment (dimensionless).

q10 = tjoel_a - tjoel_b * (0.5*(tstar + tref) - 273.15)
power = 0.1 * (tstar - tref)

factor = q10 ** power

END FUNCTION tdq10_factor

END MODULE circadian_mod
