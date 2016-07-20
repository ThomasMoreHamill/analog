SUBROUTINE interpolate_weights(nxf, nyf, nxa, nya, rlons_fcst, rlats_fcst, &
   rlons_anal, rlats_anal, raw_weightf, raw_weight)

! simple bilinear interpolation

INTEGER, INTENT(IN) :: nxf, nyf, nxa, nya
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlons_fcst
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlats_fcst
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlons_anal
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlats_anal
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: raw_weightf
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: raw_weight

DO ixa = 1, nxa
   !PRINT *,'interpolate_weights: ixa = ',ixa

   ! ---- find the relative bounding indices of this analysis longitude with respect to forecast longitude
   ixlow = 1
   ixhigh = ixlow+1
   DO ixf = 1, nxf-1
      IF (rlons_anal(ixa,1) .gt. rlons_fcst(ixf,1) .and. rlons_anal(ixa,1) .le. rlons_fcst(ixf+1,1)) THEN
         ixlow = ixf
         ixhigh = ixf+1
         GOTO 7777
      ENDIF
   END DO
   7777 CONTINUE

   DO jya = 1, nya
      !PRINT *, 'jya= ',jya

      ! ---- find the relative bounding indices of this analysis latitude with respect to the forecast latitude

      jylow = 1
      jyhigh = jylow+1
      DO jyf = 1, nyf-1
         IF (rlats_anal(1,jya) .gt. rlats_fcst(1,jyf) .and. rlats_anal(1,jya) .le. rlats_fcst(1,jyf+1)) THEN
            jylow  = jyf
            jyhigh = jyf+1
            GOTO 8888
         ENDIF
      END DO
      8888 CONTINUE

      !print *,ixlow, ixhigh, jylow, jyhigh

      ! ---- compute weights and interpolate

      weightxlow = 1. - (rlons_anal(ixa,1)-rlons_fcst(ixlow,1))/ &
           (rlons_fcst(ixhigh,1)-rlons_fcst(ixlow,1))
      weightxhigh = 1.-weightxlow
      weightylow = 1. - (rlats_anal(1,jya)-rlats_fcst(1,jylow))/ &
           (rlats_fcst(1,jyhigh)-rlats_fcst(1,jylow))
      weightyhigh = 1.-weightylow
      !print *,weightxlow, weightxhigh, weightylow, weightyhigh
      !print *,raw_weightf(ixlow,jylow), raw_weightf(ixlow,jyhigh), raw_weightf(ixhigh,jylow), raw_weightf(ixhigh,jyhigh)

      !print *,'ixa,jya = ',ixa,jya, nxa,nya
      raw_weight(ixa,jya) = &
         weightxlow *weightylow *raw_weightf(ixlow,jylow) + &
         weightxlow *weightyhigh*raw_weightf(ixlow,jyhigh) + &
         weightxhigh*weightylow *raw_weightf(ixhigh,jylow) + &
         weightxhigh*weightyhigh*raw_weightf(ixhigh,jyhigh)

   END DO ! jya

END DO
RETURN
END SUBROUTINE interpolate_weights
