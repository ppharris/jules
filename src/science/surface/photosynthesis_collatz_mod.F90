! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE photosynthesis_collatz_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHOTOSYNTHESIS_COLLATZ_MOD'

CONTAINS

SUBROUTINE prep_collatz(ft, land_field, veg_pts, veg_index, oa, tstar,         &
                        denom, qtenf_term, ccp, kc, ko)

USE conversions_mod, ONLY: zerodegc

USE pftparm, ONLY:                                                             &
! imported arrays that are not changed
  c3, q10_leaf, tlow, tupp

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with intent(in).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
 ft                                                                            &
                            ! Plant functional type.
,land_field                                                                    &
                            ! Total number of land points.
,veg_pts                                                                       &
                            ! Number of vegetated points.
,veg_index(land_field)
                            ! Index of vegetated points
                            ! on the land grid.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 oa(land_field)                                                                &
                            ! Atmospheric O2 pressure (Pa).
,tstar(land_field)
                            ! Surface temperature (K).

!-----------------------------------------------------------------------------
! Arguments with INTENT(out).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 denom(land_field)                                                             &
    ! Denominator in equation for Vcmax with the Collatz model.
,qtenf_term(land_field)                                                        &
   ! Q10 temperature term used for Vcmax with the Collatz model.
,ccp(land_field)                                                               &
   ! Photorespiratory compensatory point (Pa). This is zero for C4 plants.
,kc(land_field)                                                                &
   ! Michaelis-Menten constant for CO2 (Pa).
,ko(land_field)
   ! Michaelis-Menten constant for O2 (Pa).

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 l,m

REAL(KIND=real_jlslsm) ::                                                      &
 power                                                                         &
   ! Exponent used in Q10 term.
,tau                                                                           &
   ! Rubisco specificty for CO2 relative to O2.
,tdegc
   ! Temperature (deg C).

  ! Use the Collatz model (for C3 or C4 plants).
!$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l,m,power,tau,tdegc)   &
!$OMP SHARED(c3, veg_pts, veg_index, ccp, denom, ft, kc, ko, oa, tlow,         &
!$OMP        q10_leaf,  qtenf_term, tstar, tupp) SCHEDULE(STATIC)
  DO m  = 1,veg_pts
    l = veg_index(m)
    tdegc         = tstar(l) - zerodegc
    power         = 0.1 * (tdegc- 25.0)
    denom(l)      = (1.0 + EXP (0.3 * (tdegc - tupp(ft)))) *                   &
                    (1.0 + EXP (0.3 * (tlow(ft) - tdegc)))
    qtenf_term(l) = q10_leaf(ft)** power

    IF ( c3(ft) == 1 ) THEN
      ! Calculate terms that are only needed for C3 plants.
      ! Although oa, kc and ko are always used together we keep them separate
      ! to maintain bit comparability.
      tau    = 2600.0  * (0.57 ** power)
      ccp(l) = 0.5 * oa(l) / tau
      kc(l)  = 30.0    * (2.1 ** power)
      ko(l)  = 30000.0 * (1.2 ** power)
    END IF

  END DO
!$OMP END PARALLEL DO

END SUBROUTINE prep_collatz



! *********************************************************************
! Purpose:
! Calculates leaf internal CO2 pressure using either:
!       (i) Jacobs (1994) CI/CA closure.
!   or (ii) Ci/Ca from the Medlyn et al. (2011) conductance model.
!
! Calculates leaf-level gross photosynthesis using:
!       (i) Collatz et al. (1992) model for C3 plants
!           and Collatz et al. (1991) model for C4 plants
!
! References:
! Collatz et al., 1991, Agr. Forest Meteorol., 54: 107--136,
!   https://doi.org/10.1016/0168-1923(91)90002-8
! Collatz et al., 1992, Aust. J. Plant Physiol., 19: 519--538,
!   https://doi.org/10.1071/PP992051.
! Jacobs, 1994, Ph.D. thesis, Wageningen Agricultural University.
! Medlyn et al., 2011, Global Change Biology, 17: 2134--2144,
!   https://doi.org/10.1111/j.1365-2486.2010.02375.x
! *********************************************************************
SUBROUTINE leaf_limits_collatz(ft, land_field, veg_pts, veg_index              &
,                      acr, apar, ca, ccp, dq, fsmc, kc, ko, oa                &
,                      pstar, vcmax                                            &
,                      clos_pts, open_pts, clos_index, open_index              &
,                      ci, wcarb, wexpt, wlite )

USE pftparm, ONLY: alpha, c3, dqcrit, f0, g1_stomata
USE planet_constants_mod, ONLY: repsilon
USE jules_surface_mod, ONLY: fwe_c3, fwe_c4
USE jules_vegetation_mod, ONLY:                                                &
! imported parameters
    stomata_jacobs,                                                            &
! imported scalars that are not changed
    stomata_model

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with intent(in).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
 ft                                                                            &
                            ! Plant functional type.
,land_field                                                                    &
                            ! Total number of land points.
,veg_pts                                                                       &
                            ! Number of vegetated points.
,veg_index(land_field)
                            ! Index of vegetated points
                            ! on the land grid.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 acr(land_field)                                                               &
                            ! Absorbed PAR (mol photons/m2/s).
,apar(land_field)                                                              &
                            ! Absorbed PAR (W m-2).
,ca(land_field)                                                                &
                            ! Canopy CO2 pressure (Pa).
,ccp(land_field)                                                               &
                            ! Photorespiratory compensatory point (Pa).
,dq(land_field)                                                                &
                            ! Canopy level specific humidity deficit
                            ! (kg H2O/kg air).
,fsmc(land_field)                                                              &
                            ! Soil water factor.
,kc(land_field)                                                                &
                            ! Michaelis-Menten constant for CO2 (Pa).
,ko(land_field)                                                                &
                            ! Michaelis-Menten constant for O2 (Pa).
,oa(land_field)                                                                &
                            ! Atmospheric O2 pressure (Pa).
,pstar(land_field)                                                             &
                            ! Atmospheric pressure (Pa).
,vcmax(land_field)
                            ! Maximum rate of carboxylation of Rubisco
                            ! (mol CO2/m2/s).

INTEGER, INTENT(OUT) ::                                                        &
 clos_pts                                                                      &
                            ! Number of land points with closed stomata.
,open_pts                                                                      &
                            ! Number of land points with open stomata.
,clos_index(land_field)                                                        &
                            ! Index of land points with closed stomata.
,open_index(land_field)
                            ! Index of land points with open stomata.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 ci(land_field)                                                                &
                            ! Internal CO2 pressure (Pa).
,wcarb(land_field)                                                             &
                            ! Carboxylation-limited gross photosynthetic
!                           ! rate (mol CO2/m2/s).
,wexpt(land_field)                                                             &
                            ! Export-limited gross photosynthetic rate
!                           ! (mol CO2/m2/s). Not used with Farquhar model.
,wlite(land_field)
                            ! Light-limited gross photosynthetic rate
!                           ! (mol CO2/m2/s).

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  errcode                                                                      &
                            ! Error code to pass to ereport.
 ,j,l                       ! Loop counters.

REAL(KIND=real_jlslsm) ::                                                      &
  vpd_factor
                            ! Factor used in the calculation of humidity
                            ! deficit (kPa) from the deficit expressed in
                            ! terms of specific humidity.

LOGICAL ::                                                                     &
  l_closed(land_field)      ! Logical to mark closed points to help
                            ! parallel performance.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_LIMITS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Calculate a constant for the Medlyn model.
!-----------------------------------------------------------------------------
vpd_factor = 1.0 / ( repsilon * 1.0e3 )

!-----------------------------------------------------------------------------
! Flag open and closed points, and calculate the internal CO2 pressure.
!-----------------------------------------------------------------------------
!$OMP PARALLEL DO IF(veg_pts > 1)                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,l)                                                             &
!$OMP SHARED(veg_pts,veg_index,ft,                                             &
!$OMP        ccp,vcmax,ci,ca,f0,dq,dqcrit,l_closed,fsmc,apar,g1_stomata,       &
!$OMP        stomata_model,pstar,vpd_factor)
DO j = 1,veg_pts
  l = veg_index(j)

  ! Calculate the internal CO2 pressure.
  ! Although this is only required at points with open stomata we calculate
  ! at all points to retain bit comparability.
  IF ( stomata_model == stomata_jacobs ) THEN

    ci(l) = (ca(l) - ccp(l)) * f0(ft) * (1.0 - dq(l) / dqcrit(ft)) + ccp(l)

    ! Identify points at which the stomata are closed.
    ! Note that we test apar rather than acr (which is apar but in different
    ! units) to retain bit comparability with older versions.
    IF (fsmc(l) == 0.0 .OR. dq(l) >= dqcrit(ft) .OR. apar(l) == 0.0) THEN
      l_closed(l) = .TRUE.
    ELSE
      l_closed(l) = .FALSE.
    END IF

  ELSE

    ! stomata_model == stomata_medlyn

    ! Calculate the internal CO2 pressure.
    ! This is Eqn.13 of Medlyn et al. (2012),
    ! doi: 10.1111/j.1365-2486.2012.02790.x, also converting specific humidity
    ! deficit to vapour pressure deficit.
    ci(l) = ca(l) * g1_stomata(ft)                                             &
              / ( g1_stomata(ft) + SQRT( dq(l) * pstar(l) * vpd_factor ) )

    ! Flag where the stomata are closed.
    IF (fsmc(l) == 0.0 .OR. apar(l) == 0.0) THEN
      l_closed(l) = .TRUE.
    ELSE
      l_closed(l) = .FALSE.
    END IF

  END IF  !  stomata_model

END DO
!$OMP END PARALLEL DO

!Isolate the piece of work that won't go easily into OpenMP
clos_pts = 0
open_pts = 0
DO j = 1,veg_pts
  l = veg_index(j)
  IF ( l_closed(l) ) THEN
    clos_pts = clos_pts + 1
    clos_index(clos_pts) = j
  ELSE
    open_pts = open_pts + 1
    open_index(open_pts) = j
  END IF
END DO

!-----------------------------------------------------------------------------
! Use the Collatz models (for C3 or C4 plants).
! Calculate the gross photosynthesis for RuBP-Carboxylase-, light- and
! export-limited photosynthesis.
!-----------------------------------------------------------------------------
  IF (c3(ft) == 1) THEN

!$OMP PARALLEL DO IF(open_pts > 1)                                             &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,j)                                                             &
!$OMP SHARED(open_pts,veg_index,open_index,wcarb,vcmax,ci,ccp,kc,oa,ko,wlite,  &
!$OMP        ft,wexpt,fwe_c3,alpha,acr)
    DO j = 1,open_pts
      l = veg_index(open_index(j))

      ! The numbers in these equations are from Cox, HCTN 24,
      ! "Description ... Vegetation Model", equations 54 and 55.
      wcarb(l) = vcmax(l) * (ci(l) - ccp(l))                                   &
                 / (ci(l) + kc(l) * (1.0 + oa(l) / ko(l)))
      wlite(l) = alpha(ft) * acr(l) * (ci(l) - ccp(l)) / (ci(l) + 2.0 * ccp(l))
      wlite(l) = MAX(wlite(l), TINY(1.0e0))
      wexpt(l) = fwe_c3 * vcmax(l)
    END DO
!$OMP END PARALLEL DO

  ELSE
    !  C4
!$OMP PARALLEL DO IF(open_pts > 1)                                             &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,j)                                                             &
!$OMP SHARED(open_pts,veg_index,open_index,wcarb,vcmax,wlite,ft,wexpt,pstar,   &
!$OMP        alpha,fwe_c4,ci,acr)
    DO j = 1,open_pts
      l = veg_index(open_index(j))
      wcarb(l) = vcmax(l)
      wlite(l) = alpha(ft) * acr(l)
      wlite(l) = MAX(wlite(l), TINY(1.0e0))
      wexpt(l) = fwe_c4 * vcmax(l) * ci(l) / pstar(l)
    END DO
!$OMP END PARALLEL DO

  END IF  !  c3

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE leaf_limits_collatz

!#############################################################################
!#############################################################################

SUBROUTINE calc_photo_collatz( ft, land_pts, veg_pts, veg_index,               &
                               denom, nleaf, qtenf_term, vcmax_temp,           &
                               rd_dark, vcmax )

! Calculate the maximum rates of carboxylation of Rubisco and dark
! respiration without light inhibition.

USE jules_vegetation_mod, ONLY:                                                &
! imported scalars that are not changed
    l_trait_phys

USE pftparm, ONLY: fd, neff, vint, vsl

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
INTEGER,INTENT(IN) ::                                                          &
  ft,                                                                          &
    ! Index of plant functional type.
  land_pts,                                                                    &
    ! Number of land points.
  veg_pts,                                                                     &
    ! Number of vegetated points.
  veg_index(land_pts)
    ! Index of vegetated points on the land grid.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  denom(land_pts),                                                             &
    ! Denominator in equation for Vcmax with the Collatz model.
  nleaf(land_pts),                                                             &
    ! Leaf nitrogen concentration.
    ! If l_trait_phys = (kg N m-2),  else = (kgN [kgC]-1).
  qtenf_term(land_pts),                                                        &
   ! Q10 temperature term used for Vcmax with the Collatz model.
  vcmax_temp(land_pts)
    ! Factor expressing the effect of temperature on Vcmax.
    ! Only used with the Farquhar model.

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  rd_dark(land_pts),                                                           &
    ! Dark respiration before light inhibition (mol CO2/m2/s).
  vcmax(land_pts)
    ! Maximum rate of carboxylation of Rubisco (mol CO2/m2/s).

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  l, m
    ! Indices.

!-----------------------------------------------------------------------------
! Local array variables.
!-----------------------------------------------------------------------------
REAL ::                                                                        &
  vcmax_ref(land_pts)
    ! Maximum rate of carboxylation of Rubisco at the reference temperature,
    ! ignoring the effects of acclimation and N allocation (mol CO2/m2/s).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_PHOTO_PARAMETERS'

!-----------------------------------------------------------------------------
!end of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Calculate Vcmax at the reference temperature, without any acclimation.
!-----------------------------------------------------------------------------
!$OMP PARALLEL IF(veg_pts > 1)  DEFAULT(NONE)                                  &
!$OMP PRIVATE(l, m, n_total)                                                   &
!$OMP SHARED(ft, veg_index, veg_pts, denom, fd, neff, nleaf, qtenf_term,       &
!$OMP        rd_dark, vcmax, vcmax_ref, vint, vsl, l_trait_phys )

!$OMP DO SCHEDULE(STATIC)
DO m = 1,veg_pts

  l = veg_index(m)

  IF (l_trait_phys) THEN
    vcmax_ref(l) = (vsl(ft) * nleaf(l) + vint(ft)) * 1.0e-6  ! Kattge 2009
  ELSE
    vcmax_ref(l) = neff(ft) * nleaf(l)
  END IF

END DO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------------
! Calculate Vcmax.
!-----------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  ! Use the Collatz model.
  !---------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO m = 1,veg_pts
    l = veg_index(m)
    ! Using brackets here to recreate existing results.
    vcmax(l) = ( vcmax_ref(l) * qtenf_term(l) ) / denom(l)
  END DO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------------
! Calculate dark respiration. Any effect of light inhibition is added later.
!-----------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO m = 1,veg_pts
  l = veg_index(m)
  rd_dark(l) = fd(ft) * vcmax(l)
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE calc_photo_collatz

END MODULE photosynthesis_collatz_mod
