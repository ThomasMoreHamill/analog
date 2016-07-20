SUBROUTINE compute_all_climatologies_ccpa2 (nxa, nya, ntrainsamps, &
   apcp_anal, iccpa_mask, thresh_95, thresh_uppertercile, thresh_upperquintile, &
   thresh_upperdecile, thresh_median, climo_pop, climo_1mm, climo_2p5mm, &
   climo_5mm, climo_10mm, climo_25mm, climo_50mm, climo_median, &
   climo_uppertercile, climo_upperquintile, climo_upperdecile, climo_95)

! compute the climatological probability of a given event.  Note that 
! in this case, the climatology is not cross-validated, e.g., when we will
! be computing the skill scores for 2002 forecasts with respect to climatology,
! 2002 data will be included in computing the climatology.  In this respect,
! the final skill scores for the post-processed data, computed with respect to 
! climatology, will be somewhat diminished, because the climatological
! reference will be more skillful than it ought.  This is a simplification
! to save memory.
!
! coded by Tom Hamill, tom.hamill@noaa.gov, (303) 497-3060

INTEGER, INTENT(IN) :: nxa, nya, ntrainsamps
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: apcp_anal
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: iccpa_mask
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_95
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_uppertercile
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_upperquintile
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_upperdecile
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_median
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_pop
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_1mm
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_2p5mm
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_5mm
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_10mm
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_25mm
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_50mm
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_median
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_uppertercile
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_upperquintile
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_upperdecile
REAL, INTENT(OUT), DIMENSION(nxa,nya) ::  climo_95

DO ix = 1, nxa
   DO jy = 1, nya
      IF (iccpa_mask(ix,jy) .eq. 1) THEN
         ktr_pop = 0
         ktr_1mm = 0
         ktr_2p5mm = 0
         ktr_5mm = 0
         ktr_10mm = 0
         ktr_25mm = 0
         ktr_50mm = 0
         ktr_median = 0
         ktr_uppertercile = 0
         ktr_upperquintile = 0
         ktr_upperdecile = 0
         ktr_95 = 0
         ktr = 0
         DO isamp = 1, ntrainsamps
            IF (apcp_anal(ix,jy,isamp) .ge. 0.1) ktr_pop   = ktr_pop + 1
            IF (apcp_anal(ix,jy,isamp) .ge. 1.0) ktr_1mm   = ktr_1mm + 1
            IF (apcp_anal(ix,jy,isamp) .ge. 2.5) ktr_2p5mm = ktr_2p5mm + 1
            IF (apcp_anal(ix,jy,isamp) .ge. 5.0) ktr_5mm   = ktr_5mm + 1
            IF (apcp_anal(ix,jy,isamp) .ge. 10.) ktr_10mm  = ktr_10mm + 1
            IF (apcp_anal(ix,jy,isamp) .ge. 25.) ktr_25mm  = ktr_25mm + 1
            IF (apcp_anal(ix,jy,isamp) .ge. 50.) ktr_50mm  = ktr_50mm + 1
            IF (apcp_anal(ix,jy,isamp) .ge. thresh_median(ix,jy)) &
                 ktr_median = ktr_median + 1
            IF (apcp_anal(ix,jy,isamp) .ge. thresh_uppertercile(ix,jy)) &
                 ktr_uppertercile = ktr_uppertercile + 1
            IF (apcp_anal(ix,jy,isamp) .ge. thresh_upperquintile(ix,jy)) &
                 ktr_upperquintile = ktr_upperquintile + 1
            IF (apcp_anal(ix,jy,isamp) .ge. thresh_upperdecile(ix,jy)) &
                 ktr_upperdecile = ktr_upperdecile + 1
            IF (apcp_anal(ix,jy,isamp) .ge. thresh_95(ix,jy)) ktr_95 = ktr_95 + 1
            IF (apcp_anal(ix,jy,isamp) .ge. 0.0) ktr = ktr+1
         END DO
         rktr = REAL(ktr)

         IF (rktr .ge. 1) THEN
            climo_pop(ix,jy) = REAL(ktr_pop) / rktr
            climo_1mm(ix,jy) = REAL(ktr_1mm) / rktr 
            climo_2p5mm(ix,jy) = REAL(ktr_2p5mm) / rktr 
            climo_5mm(ix,jy) = REAL(ktr_5mm) / rktr 
            climo_10mm(ix,jy) = REAL(ktr_10mm) / rktr 
            climo_25mm(ix,jy) = REAL(ktr_25mm) / rktr 
            climo_50mm(ix,jy) = REAL(ktr_50mm) / rktr 
            climo_median(ix,jy) = REAL(ktr_median) / rktr 
            climo_uppertercile(ix,jy) = REAL(ktr_uppertercile) /  rktr
            climo_upperquintile(ix,jy) = REAL(ktr_upperquintile) / rktr 
            climo_upperdecile(ix,jy) = REAL(ktr_upperdecile) / rktr 
            climo_95(ix,jy) = REAL(ktr_95) / rktr 
         ELSE
            climo_pop(ix,jy) = -99.99
            climo_1mm(ix,jy) = -99.99
            climo_2p5mm(ix,jy) = -99.99
            climo_5mm(ix,jy) = -99.99
            climo_10mm(ix,jy) = -99.99
            climo_25mm(ix,jy) = -99.99
            climo_50mm(ix,jy) = -99.99
            climo_median(ix,jy) = -99.99
            climo_uppertercile(ix,jy) = -99.99
            climo_upperquintile(ix,jy) = -99.99
            climo_upperdecile(ix,jy) = -99.99
            climo_95(ix,jy) = -99.99
         ENDIF
      ELSE
         climo_pop(ix,jy) = -99.99
         climo_1mm(ix,jy) = -99.99 
         climo_2p5mm(ix,jy) = -99.99 
         climo_5mm(ix,jy) = -99.99 
         climo_10mm(ix,jy) = -99.99 
         climo_25mm(ix,jy) = -99.99 
         climo_50mm(ix,jy) = -99.99 
         climo_median(ix,jy) = -99.99 
         climo_uppertercile(ix,jy) = -99.99 
         climo_upperquintile(ix,jy) = -99.99 
         climo_upperdecile(ix,jy) = -99.99
         climo_95(ix,jy) = -99.99 
      END IF

   END DO ! jy
END DO    ! ix

RETURN
END SUBROUTINE compute_all_climatologies_ccpa2
