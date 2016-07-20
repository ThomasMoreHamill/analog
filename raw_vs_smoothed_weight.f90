! ==============================================================================

SUBROUTINE raw_vs_smoothed_weight(nxa, nya, nxf, nyf, rlons_fcst, rlats_fcst, &
   rlons_anal, rlats_anal, input_data_directory, raw_weight)

! ---- we want a linear combination of the analog ensemble (raw) probability and the 
!      Savitzky-Golay smoothed probabilities, weighting more toward the former
!      in complex terrain and latter in smooth terrain.  Determine the weight
!      to apply to the raw analog ensemble.

! coded by Tom Hamill, tom.hamill@noaa.gov, +1 (303) 497-3060

INTEGER, INTENT(IN) :: nxa, nya, nxf, nyf
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlons_fcst, rlats_fcst
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlons_anal, rlats_anal
CHARACTER*(*), INTENT(IN) :: input_data_directory
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: raw_weight

CHARACTER*80 infile
REAL work(9), wmean, workdev(9)
REAL, DIMENSION(nxf,nyf) :: raw_weightf
REAL, DIMENSION(nxf,nyf) :: terrain

! --- read in the terrain data

infile = TRIM(input_data_directory) // '/conus_orography_t254.dat'
print *,'reading from ', infile
OPEN (unit=1,file=infile,status='old',form='unformatted')
READ (1) nxfin, nyfin
READ (1) terrain
CLOSE (1)

PRINT *,'max(terrain) = ',maxval(terrain)

! ---- determine the standard deviation of the terrain about its mean and then
!      the weight.  See Figs. A4 and A5 from
!      www.esrl.noaa.gov/psd/people/tom.hamill/Analog-CCPA-MWRexpedited-Hamill-AppA-v2.pdf

workdev(:) = 0.0
DO i = 1, nxf
   imin = MAX(1,i-1)
   imax = MIN(nxf,i+1)
   DO j = 1, nyf
      jmin = MAX(1,j-1)
      jmax = MIN(nyf,j+1)
      ktr = 0
	  work(:) = 0.0
      DO i2 = imin, imax
         DO j2 = jmin, jmax
            ktr = ktr+1
            work(ktr) = terrain(i2,j2)
         END DO
      END DO
	  IF (ktr .lt. 1) ktr = 1
      wmean = SUM(work(1:ktr))/REAL(ktr)
      workdev(1:ktr) = work(1:ktr) - wmean
      var = 0.
      DO k = 1, ktr
         var = var + workdev(k)**2
      END DO
	  IF (ktr .gt. 1) THEN 
		  stdev = SQRT(var/REAL(ktr-1))
	  ELSE
		  stdev = 0.
	  ENDIF
      sqrt_stdev = SQRT(stdev)
      weight = 0.2 + (sqrt_stdev-8.)/16.6667
      IF (weight .lt. 0.2) weight = 0.2
      IF (weight .gt. 0.8) weight = 0.8
      raw_weightf(i,j) = weight
   END DO
END DO

! ---- now interpolate the weights to the analysis grid

print *,'from raw_vs_smoothed_weight, calling interpolate_weights'
CALL interpolate_weights(nxf, nyf, nxa, nya, rlons_fcst, rlats_fcst, &
   rlons_anal, rlats_anal, raw_weightf, raw_weight)

RETURN
END

 
