SUBROUTINE compute_quantiles_ccpa2 (nx, ny, ntrainsamps, verifppn, iccpa_mask, &
  thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, thresh_median, &
  thresh_95)

! compute the precipitation amount thresholds associated with various
! quantiles of the CCPA precipitation analysis for this month of the year.
! Note the lack of cross validation; information includes precipitation forecasts
! from the year currently being processed.

! coded by Tom Hamill, tom.hamill@noaa.gov, (303) 497-3060

INTEGER, INTENT(IN) :: nx, ny, ntrainsamps
INTEGER*2, INTENT(IN), DIMENSION(nx,ny) :: iccpa_mask
REAL, INTENT(IN), DIMENSION(nx,ny,ntrainsamps) :: verifppn
REAL, INTENT(OUT), DIMENSION(nx,ny) :: thresh_uppertercile
REAL, INTENT(OUT), DIMENSION(nx,ny) :: thresh_upperquintile
REAL, INTENT(OUT), DIMENSION(nx,ny) :: thresh_upperdecile
REAL, INTENT(OUT), DIMENSION(nx,ny) :: thresh_median
REAL, INTENT(OUT), DIMENSION(nx,ny) :: thresh_95

REAL, DIMENSION(ntrainsamps) :: xsort

thresh_uppertercile = -99.99
thresh_upperquintile = -99.99
thresh_upperdecile = -99.99
thresh_median = -99.99
thresh_95 = -99.99

ns = ntrainsamps
nsampsgood = ns
ibottom = 1
itop = nsampsgood
DO ix = 1, nx
   DO jy = 1, ny
      IF (iccpa_mask(ix,jy) .eq. 1) THEN ! over CONUS

         ! ---- sort the data

         xsort(:) = verifppn(ix,jy,:)
         DO i = 1, ntrainsamps
            IF (xsort(i) .gt. 1000) xsort(i) = -99.99
         END DO
         CALL heapsort(ns,xsort)
         DO isamp = 1, ns  ! may be some missing values w. -99.99 at beginning
            IF (xsort(isamp) .GE. 0.0) THEN
               ibottom = isamp
               itop = ntrainsamps
               nsampsgood = itop - ibottom + 1
               GOTO 6000
            END IF
         END DO
6000     CONTINUE
         ifirstnonz = 1   ! first nonzero
         DO isamp = 1, ns  
            IF (xsort(isamp) .GT. 0.0) THEN
               ifirstnonz = isamp
               goto 7000
            ENDIF
         END DO
7000     continue
         rfirstquant = REAL(ifirstnonz-ibottom) / REAL(nsampsgood)

         ! ---- we will have some situations where the precip amount associated with 
         !      a given quantile of the distribution is zero.  In this case (following our
         !      procedure in Hamill et al. June 2004 MWR, bottom p. 1438) reassign the 
         !      precipitation threshold associated with that quantile to the first 
         !      nonzero precipitation value.

         IF (rfirstquant .GT. 0.667) THEN
            thresh_uppertercile(ix,jy) = MAX(xsort(ifirstnonz),0.0001)
         ELSE
            i = ibottom + NINT(0.667*REAL(nsampsgood))
            thresh_uppertercile(ix,jy) = xsort(i)
         END IF

         IF (rfirstquant .GT. 0.8) THEN
			thresh_upperquintile(ix,jy) = MAX(xsort(ifirstnonz),0.0001)
         ELSE
            i = ibottom + NINT(0.8*REAL(nsampsgood))
            thresh_upperquintile(ix,jy) = xsort(i)
         ENDIF

         IF (rfirstquant .GT. 0.9) THEN
			thresh_upperdecile(ix,jy) = MAX(xsort(ifirstnonz),0.0001)
         ELSE
            i = ibottom + NINT(0.9*REAL(nsampsgood))
            thresh_upperdecile(ix,jy) = xsort(i)
         ENDIF

         IF (rfirstquant .GT. 0.5) THEN
			thresh_median(ix,jy) = MAX(xsort(ifirstnonz),0.0001)
         ELSE
            i = ibottom + NINT(0.5*REAL(nsampsgood))
            thresh_median(ix,jy) = xsort(i)
         ENDIF

         IF (rfirstquant .GT. 0.95) THEN
			thresh_95(ix,jy) = MAX(xsort(ifirstnonz),0.0001)
         ELSE
            i = ibottom + NINT(0.95*REAL(nsampsgood))
            thresh_95(ix,jy) = xsort(i)
         ENDIF

      ENDIF ! inside conus
   END DO  ! jy
END DO     ! ix

RETURN
END

