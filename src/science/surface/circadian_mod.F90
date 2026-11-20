MODULE circadian_mod
USE um_types, ONLY : real_jlslsm

IMPLICIT NONE

CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='CIRCADIAN_MOD'

REAL(KIND=real_jlslsm), PARAMETER :: ipar_min = 1.0E-6
!                           ! PAR threshold used to determine day/night.

PRIVATE

PUBLIC ipar_min, rd_faber

CONTAINS

! *********************************************************************
! Return the circadian multiplication factor from Eqn 8 of Faber et al
! (2022).  The coefficients come from the caption to Fig S15 of that
! paper.
! *********************************************************************
FUNCTION rd_faber(model_time, longitude, latitude) RESULT(factor)

USE conversions_mod, ONLY : rsec_per_day, rsec_per_hour, rhour_per_day

IMPLICIT NONE

INTEGER, INTENT(in) :: model_time
REAL(KIND=real_jlslsm), INTENT(in) :: longitude, latitude

REAL(KIND=real_jlslsm), PARAMETER :: a_aus = 0.6487044
REAL(KIND=real_jlslsm), PARAMETER :: b_aus = 0.0936338
REAL(KIND=real_jlslsm), PARAMETER :: b1_aus = -0.0044895

REAL(KIND=real_jlslsm), PARAMETER :: a_den = 0.3795814
REAL(KIND=real_jlslsm), PARAMETER :: b_den = 0.1168153
REAL(KIND=real_jlslsm), PARAMETER :: b1_den = -0.0046665

REAL(KIND=real_jlslsm), PARAMETER :: base_hour = 8.0  ! Local hour of day.

REAL(KIND=real_jlslsm) :: local_hour
REAL(KIND=real_jlslsm) :: factor

local_hour = (model_time + (longitude / 360.0) * rsec_per_day) / rsec_per_hour
local_hour = MODULO(local_hour, rhour_per_day)

factor = (a_aus + b_aus * local_hour + b1_aus * local_hour * local_hour) /     &
         (a_aus + b_aus * base_hour + b1_aus * base_hour * base_hour)

END FUNCTION rd_faber

END MODULE circadian_mod
