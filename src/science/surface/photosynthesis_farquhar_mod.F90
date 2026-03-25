! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE photosynthesis_farquhar_mod
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHOTOSYNTHESIS_FARQUHAR_MOD'

CONTAINS
! *********************************************************************
! Purpose:
! Calculates leaf internal CO2 pressure using either:
!       (i) Jacobs (1994) CI/CA closure.
!   or (ii) Ci/Ca from the Medlyn et al. (2011) conductance model.
!
! Calculates leaf-level gross photosynthesis using:
!           Farquhar et al. (1980) model for C3 plants.
!
! References:
! Farquhar et al., 1980, Planta, 149: 78--90,
!   https://doi.org/10.1007/BF0038623
! Jacobs, 1994, Ph.D. thesis, Wageningen Agricultural University.
! Medlyn et al., 2011, Global Change Biology, 17: 2134--2144,
!   https://doi.org/10.1111/j.1365-2486.2010.02375.x
! *********************************************************************
SUBROUTINE leaf_limits_farquhar(ft, land_field, veg_pts, veg_index             &
,                      apar, ca, ccp, dq, fsmc, je, km                         &
,                      pstar, vcmax                                            &
,                      clos_pts, open_pts, clos_index, open_index              &
,                      ci, wcarb, wlite )

USE pftparm, ONLY: dqcrit, f0, g1_stomata
USE planet_constants_mod, ONLY: repsilon
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
 apar(land_field)                                                              &
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
,je(land_field)                                                                &
                            ! Electron transport rate (mol m-2 s-1)
,km(land_field)                                                                &
                            ! A combination of Michaelis-Menten and other
                            ! terms.
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
! Use the Farquhar model (for C3 plants only).
! Calculate the gross photosynthesis for RuBP-Carboxylase- and light-limited
! photosynthesis.
!-----------------------------------------------------------------------------

!$OMP PARALLEL DO IF(open_pts > 1)                                             &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,j)                                                             &
!$OMP SHARED(open_pts,veg_index,open_index,wcarb,vcmax,wlite,ci,ccp,km,je)
  DO j = 1,open_pts
    l = veg_index(open_index(j))
    wcarb(l) = vcmax(l) * ( ci(l) - ccp(l) ) / ( ci(l) + km(l) )
    wlite(l) = je(l) / 4.0 * ( ci(l) - ccp(l) ) / ( ci(l) + 2.0 * ccp(l) )
    wlite(l) = MAX(wlite(l), TINY(1.0e0))
  END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE leaf_limits_farquhar



!#############################################################################
!#############################################################################

SUBROUTINE calc_photo_farquhar( ft, land_pts, veg_pts, veg_index,              &
                                jmax_temp, jv25, nleaf, vcmax_temp,            &
                                jmax, rd_dark, vcmax )

! Calculate the maximum rates of carboxylation of Rubisco and electron
! transport, and dark respiration without light inhibition.

USE jules_vegetation_mod, ONLY:                                                &
! imported parameters
    jv_ntotal, jv_scale,                                                       &
! imported scalars that are not changed
    n_alloc_jmax, n_alloc_vcmax, l_trait_phys, photo_jv_model

USE pftparm, ONLY: fd, jv25_ratio, neff, vint, vsl

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
  jmax_temp(land_pts),                                                         &
    ! Factor expressing the effect of temperature on Jmax.
    ! Only used with the Farquhar model.
  jv25(land_pts),                                                              &
    ! Ratio of Jmax to Vcmax at 25 degC, including any acclimation.
    ! Only used with the Farquhar model.
  nleaf(land_pts),                                                             &
    ! Leaf nitrogen concentration.
    ! If l_trait_phys = (kg N m-2),  else = (kgN [kgC]-1).
  vcmax_temp(land_pts)
    ! Factor expressing the effect of temperature on Vcmax.
    ! Only used with the Farquhar model.

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  jmax(land_pts),                                                              &
    ! Maximum rate of electron transport (mol CO2 m-2 s-1).
    ! Only calculated with the Farquhar model.
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

REAL(KIND=real_jlslsm) ::                                                      &
  n_total,                                                                     &
    ! Total N allocated to photosynthetic components (kg m-2).
  recip_j,                                                                     &
    ! Reciprocal of n_alloc_jmax (kg m-2 of N [mol CO2 m-2 s-1]).
  recip_v
    ! Reciprocal of n_alloc_vcmax (kg m-2 of N [mol CO2 m-2 s-1]).

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
! Calculate some constants.
!-----------------------------------------------------------------------------
IF ( photo_jv_model == jv_ntotal ) THEN
  recip_j  = 1.0 / n_alloc_jmax
  recip_v  = 1.0 / n_alloc_vcmax
END IF

!-----------------------------------------------------------------------------
! Calculate Vcmax at the reference temperature, without any acclimation.
!-----------------------------------------------------------------------------
!$OMP PARALLEL IF(veg_pts > 1)  DEFAULT(NONE)                                  &
!$OMP PRIVATE(l, m, n_total)                                                   &
!$OMP SHARED(ft, photo_jv_model, veg_index, veg_pts,                           &
!$OMP        fd, jmax, jmax_temp, jv25, jv25_ratio, neff, nleaf,               &
!$OMP        rd_dark, recip_j, recip_v, vcmax, vcmax_ref,                      &
!$OMP        vcmax_temp, vint, vsl, l_trait_phys )

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
! Calculate Vcmax and Jmax.
!-----------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Use the Farquhar model (for C3 plants).
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Calculate values at the reference temperature, including any acclimation
  ! of jv25 (but excluding other acclimation terms).
  !---------------------------------------------------------------------------
  SELECT CASE ( photo_jv_model )

  CASE ( jv_scale )
    ! Find J25 by scaling V25.
!$OMP DO SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      vcmax(l) = vcmax_ref(l)
      jmax(l)  = vcmax_ref(l) * jv25(l)
    END DO
!$OMP END DO NOWAIT

  CASE ( jv_ntotal )
    ! Assume the total N allocated to photosynthetic capacity is constant.
!$OMP DO SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      ! Calculate total N allocated to photosynthetic capacity, using the
      ! prescribed parameters at the reference temperature.
      ! This is Eq.5 of Mercado et al. (2018).
      n_total = vcmax_ref(l) * recip_v                                         &
                + vcmax_ref(l) * jv25_ratio(ft) * recip_j
      ! Calculate Vcmax and Jmax at 25degC, including temperature acclimation
      ! of J:V.
      vcmax(l) = n_total / ( recip_v + jv25(l) * recip_j )
      jmax(l)  = n_total / ( recip_v / jv25(l) + recip_j )
    END DO
!$OMP END DO NOWAIT

  END SELECT  !  photo_jv_model

  !---------------------------------------------------------------------------
  ! Calculate final values, including temperature effect.
  !---------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO m = 1,veg_pts
    l = veg_index(m)

    ! Calculate rates according to acclimated ratio and N allocation to
    ! photosynthesis, and including temperature term.
    ! At present neither acclimation nor N allocation are represented.
    vcmax(l) = vcmax(l) * vcmax_temp(l)
    jmax(l)  = jmax(l)  * jmax_temp(l)

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

END SUBROUTINE calc_photo_farquhar


!#############################################################################
!#############################################################################

SUBROUTINE calc_electron_flux( land_pts, veg_pts, veg_index, i2, jmax, je )

! Calculate the electron flux for the Farquhar model.

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
INTEGER,INTENT(IN) ::                                                          &
  land_pts,                                                                    &
    ! Number of land points.
  veg_pts,                                                                     &
    ! Number of vegetated points.
  veg_index(land_pts)
    ! Index of vegetated points on the land grid.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  i2(land_pts),                                                                &
    ! Radiation that goes to Photosystem II, expressed as an electron flux
    ! (mol m-2 s-1).
  jmax(land_pts)
    ! Maximum rate of electron transport (mol CO2 m-2 s-1).

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  je(land_pts)
    ! Electron transport rate (mol m-2 s-1).

!-----------------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  light_curvature = 0.90
    ! Curvature of the light response function. Used with Farquhar model of
    ! photosynthesis. See Eq.4 of Medlyn et al. (2002).

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  l, m
    ! Indices.

REAL(KIND=real_jlslsm) ::                                                      &
 recip_denom
   ! The reciprocal of the denominator.


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_ELECTRON_FLUX'

!-----------------------------------------------------------------------------
!end of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Calculate a constant.
!-----------------------------------------------------------------------------
recip_denom = 1.0 / ( 2.0 * light_curvature )

!-----------------------------------------------------------------------------
! Calculate electron flux by finding a root of a quadratic equation.
! This is the solution of Eq.4 of Medlyn et al. (2002).
!-----------------------------------------------------------------------------
DO m = 1,veg_pts
  l = veg_index(m)
  je(l)  = ( i2(l) + jmax(l)                                                   &
                    - SQRT( ( i2(l) + jmax(l) )**2                             &
                            - 4.0 * light_curvature * i2(l) * jmax(l) )        &
           ) * recip_denom
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE calc_electron_flux

END MODULE photosynthesis_farquhar_mod
