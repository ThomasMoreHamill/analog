SUBROUTINE interpolate_training_data(nxf, nyf, nxa, nya, ntrainsamps, nsupp, nsuppmax, ixa, jya, &
     xlocation, ylocation, rlons_fcst, rlats_fcst, rlons_anal, rlats_anal, inearest, jnearest, &
     forecast, forecast_interpolated)

! simple bilinear interpolation of the forecast to the analyzed location for each of the 
! supplemental data locations.  called from analog_forecast_supp_v2.f90.  Bilinear interpolation
! taken from Press et al. 1992 Numerical Recipes text.

INTEGER, INTENT(IN) :: nxf, nyf, nxa, nya, ntrainsamps, nsupp, nsuppmax, ixa, jya
INTEGER, INTENT(IN), DIMENSION(nxa,nya,nsuppmax) :: xlocation ! supplemental data locations
INTEGER, INTENT(IN), DIMENSION(nxa,nya,nsuppmax) :: ylocation
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlons_fcst
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlats_fcst
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlons_anal
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlats_anal
INTEGER, INTENT(IN), DIMENSION(nxa,nya) :: inearest
INTEGER, INTENT(IN), DIMENSION(nxa,nya) :: jnearest
REAL, INTENT(IN), DIMENSION(nxf, nyf, ntrainsamps)  :: forecast ! precipitation fcst mean
REAL, DIMENSION(nsupp,ntrainsamps) :: forecast_interpolated

DO isupp = 1, nsupp
   ilocna = xlocation(ixa,jya,isupp)
   jlocna = ylocation(ixa,jya,isupp)
   ixlow  = inearest(ilocna,jlocna)
   jylow  = jnearest(ilocna,jlocna)
   ixhigh = ixlow + 1
   jyhigh = jylow + 1
   weightxlow = 1. - (rlons_anal(ilocna,jlocna)-rlons_fcst(ixlow,1))/ &
      (rlons_fcst(ixhigh,1)-rlons_fcst(ixlow,1))
   weightxhigh = 1.-weightxlow
   weightylow = 1. - (rlats_anal(ilocna,jlocna)-rlats_fcst(1,jylow))/ &
      (rlats_fcst(1,jyhigh)-rlats_fcst(1,jylow))
   weightyhigh = 1.-weightylow
   DO i = 1, ntrainsamps
      forecast_interpolated(isupp,i) = &
         weightxlow *weightylow *forecast(ixlow,jylow,i) + &
         weightxlow *weightyhigh*forecast(ixlow,jyhigh,i) + &
         weightxhigh*weightylow *forecast(ixhigh,jylow,i) + &
         weightxhigh*weightyhigh*forecast(ixhigh,jyhigh,i)
   END DO
END DO
RETURN
END SUBROUTINE interpolate_training_data
