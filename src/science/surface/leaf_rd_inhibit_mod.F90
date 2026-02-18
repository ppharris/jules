MODULE leaf_rd_inhibit_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LEAF_RD_INHIBIT_MOD'

PRIVATE
PUBLIC leaf_rd_inhibit

CONTAINS

! *********************************************************************
!
! Reduces leaf dark respiration to account for its inhibition in
! daylight.
!
! *********************************************************************
SUBROUTINE leaf_rd_inhibit(land_pts, veg_index, open_pts, open_index, acr      &
,                          faparv_layer, rd)

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
 land_pts                                                                      &
!                           ! IN Number of land points to be processed.
,veg_index(land_pts)                                                           &
!                           ! IN Index of vegetated points on the land grid.
,open_pts                                                                      &
!                           ! IN Number of land points with open stomata.
,open_index(land_pts)
!                           ! IN Index of land points with open stomata.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 acr(land_pts)                                                                 &
!                           ! IN Absorbed PAR (mol photons/m2/s).
,faparv_layer(land_pts)
                            ! IN Fraction of absorbed PAR

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                       &
 rd(land_pts)
!                           ! INOUT Dark respiration, including any effect of
!                           ! light inhibition (mol CO2/m2/s).

INTEGER ::                                                                     &
 l, m                       ! WORK loop counters.

REAL(KIND=real_jlslsm) ::                                                      &
 apar
!                           ! WORK Absorbed PAR in layer (W/m2).

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
 apar_min = 10.0                                                               &
!                           ! Threshold above which respiration is inhibited. 
,convapar = 1.0e6
!                           ! Conversion factor.

!$OMP PARALLEL DO IF(open_pts > 1) DEFAULT(NONE) PRIVATE(apar, l, m)           &
!$OMP             SHARED(acr, apar_min, convapar, faparv_layer,                &
!$OMP                    open_index, open_pts, rd, veg_index) SCHEDULE(STATIC)
DO m = 1,open_pts
  l = veg_index(open_index(m))
  apar = acr(l) * convapar * faparv_layer(l)
  IF ( apar > apar_min ) THEN
    rd(l) = ( 0.5 - 0.05 * LOG(apar) ) * rd(l)
  END IF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE leaf_rd_inhibit
END MODULE leaf_rd_inhibit_mod

