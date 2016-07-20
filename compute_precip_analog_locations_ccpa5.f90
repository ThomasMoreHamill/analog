PROGRAM compute_precip_analog_locations_ccpa5

! purpose:  For every (i,j) in CONUS, determine a set of "supplemental" data locations.
!    Base the supplemental locations on the similarity of analyzed CDFs as well as the
!    local terrain heights and gradients.  Do this separately for every month of the year,
!    as precipitation climatology may change with the seasons.  
!    We want grid points that will have somewhat independent data, so make sure that 
!    the set of analog points are not too close together (minseparation, below).
! 
! coded by Tom Hamill, most recently touched in Mar 2016.  tom.hamill@noaa.gov
! 
USE netcdf

PARAMETER (nyears        = 14)   ! 2002-2015 here for the CCPA data
PARAMETER (ndaysub       = 91)   ! number of days of training data for each yr, 15th of month +/- 45 days
PARAMETER (nxa           = 464)  ! CCPA 1/8-degree precip analysis grid dimensions, surrounding CONUS
PARAMETER (nya           = 224)  !
PARAMETER (nmonths       = 12)   ! # months of the year
PARAMETER (nmaxlist      = 16)   ! max number of CCPA grid pts associated with any forecast grid pt.
PARAMETER (nanalogs      = 20)   ! max number of CCPA grid analogs locations to be stored for a given forecast grid box ! prev 20
PARAMETER (minseparation = 8)    ! user-defined minimum separation between analysis grid locations for analog (in grid pts)
PARAMETER (maxseparation = 160)  ! user-defined minimum separation between analysis grid locations for analog (in grid pts) ! prev 120
PARAMETER (npct          = 107)  ! number of quantiles of CDF to compute at
PARAMETER (paweight      = 0.4)  ! relative weight to give to analysis CDF differences
PARAMETER (ter_weight    = 0.2)  ! relative weight to give to terrain differences
PARAMETER (gradx_weight  = 0.2)  ! relative weight to give to terrain gradx
PARAMETER (grady_weight  = 0.2)  ! relative weight to give to terrain grady
PARAMETER (distpenalty   = .01)  ! penalty function constant for distance ! prev 0.02
PARAMETER (dx         = 52083.)  ! forecast grid spacing at equator

INTEGER*2, DIMENSION(nxa,nya)        :: conusmask  ! 1/0 mask for whether precip analysis is in CONUS
INTEGER*2, DIMENSION(nxa,nya)        :: maskout    ! array used to block out nearby points from further consideration
INTEGER, DIMENSION(nxa,nya,nanalogs) :: xlocationa ! for each forecast point, a list of which other forecast points
INTEGER, DIMENSION(nxa,nya,nanalogs) :: ylocationa !    have the closest climatologies and forecast-obs relationship
INTEGER minx
INTEGER miny

REAL, DIMENSION(nxa,nya)      :: gradx_mean  ! mean terrain height gradient in E-W direction
REAL, DIMENSION(nxa,nya)      :: gradx_stddev  ! local standard deviation of E-W terrain height gradient
REAL, DIMENSION(nxa,nya)      :: grady_mean  ! as above but for N-S direction
REAL, DIMENSION(nxa,nya)      :: grady_stddev
REAL, DIMENSION(nxa,nya)      :: ter_mean ! mean terrain height in region surrounding a particular (i,j)
REAL, DIMENSION(nxa,nya)      :: ter_stddev ! local standard deviation of terrain height
REAL, DIMENSION(nxa,nya)      :: difference
REAL, DIMENSION(nxa,nya)      :: lonsa  ! analysis grid lons and lats
REAL, DIMENSION(nxa,nya)      :: latsa
REAL, DIMENSION(nxa,nya,npct) :: panal_quantiles  ! CDF information for analyzed precip on CCPA grid
REAL, DIMENSION(npct)         :: pctv
REAL, DIMENSION(npct)         :: pa_here  ! CDF at current location
REAL, DIMENSION(npct)         :: pa_there ! CDF at alternative location
REAL, DIMENSION(nxa,nya)      :: pa_stddev  ! local standard deviation of differences in CDFs
REAL, DIMENSION(nxa,nya)      :: pa_mean  ! local mean of differences in CDFs
REAL, DIMENSION(nxa,nya)      :: terrain  ! terrain height data

REAL diffmin

CHARACTER*120 outfile, output_data_directory, input_data_directory
CHARACTER*3 cleadb
CHARACTER*3 cleade
CHARACTER*3, DIMENSION(12) :: cmonth
CHARACTER*3 csupp

DATA cmonth /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/
DATA input_data_directory /'/Projects/Reforecast2/netcdf'/
DATA output_data_directory /'/Projects/Reforecast2/netcdf'/


! --- hardcoded here to use 24-h accumulated precipitation CDFs as part of determining
!     similarity between locations.  Here we'll pluck this 24-h analyzed precipitation
!     CDFs from one of the netCDF files we've previously computed, here a file for
!     the 24- to 48-hour forecasts.

cleadb = '024'
cleade = '048'
!READ (cleadb,'(i3)') ileadb
!READ (cleade,'(i3)') ileade
IF (nanalogs .ge. 100) THEN
   WRITE (csupp,'(i3)') nanalogs
ELSE
   WRITE (csupp,'(i2)') nanalogs
   csupp = '0' // csupp(1:2)
ENDIF
print *,'csupp = ',csupp

! ---- for each month, load the quantile information and rank correlation previously
!      computed, and then determine analogs

xlocationa = -99
ylocationa = -99
DO imonth = 1, 12

   ! ---- load the CDF data for those 24- to 48-h forecasts

   CALL load_precipquantiles_ccpa_analonly (imonth, nxa, nya, npct, &
   	cleadb, cleade, input_data_directory, &
	conusmask, lonsa, latsa, pctv, panal_quantiles, terrain)

   ! ---- set terrain values for non-CONUS to zero
   
   terrain = terrain * REAL(conusmask)

   ! ---- loop thru grid points 

   PRINT *,'----- conusmask -----'
   DO j = nya/2+10,nya/2-10,-1
      PRINT 122, j,(conusmask(i,j),i=nxa/2-10,nxa/2+10)
      122 format(i3,1x,21(i1))
   END DO

   PRINT *,'preprocessing to determine standard devs'
   DO ixa = 1, nxa
      IF (MOD(ixa,10) .eq. 0) PRINT *,'processing  ixa = ',ixa,' of ',nxa
      ! --- only bother searching in a box +/- 10 grid points to limit computations.
      ixmin = MAX(1,ixa-10)
      ixmax = MIN(nxa,ixa+10)
      DO jya = 1, nya
         !PRINT *,'ixa, jya = ',ixa,jya
         jymin = MAX(1,jya-10)
         jymax = MIN(nya,jya+10)
         sum_pa2 = 0.0
         sum_pa  = 0.0
         sum_t   = 0.0
         sum_t2  = 0.0
         sum_gradx  = 0.0
         sum_gradx2 = 0.0
         sum_grady  = 0.0
         sum_grady2 = 0.0

         nsamps  = 0
         ter_here = 0.
         gradx_here = 0.
         grady_here = 0.
         IF (conusmask(ixa,jya) .gt. 0) THEN ! grid point "inside" US
            pa_here(:) = panal_quantiles(ixa, jya,:)
            ter_here = terrain(ixa,jya)
            IF (ixa .ne. 1 .and. ixa .ne. nxa .and. jya .ne. 1 .and. jya .ne. nya) THEN
               gradx_here = terrain(ixa+1,jya)-terrain(ixa-1,jya)
               grady_here = terrain(ixa,jya+1)-terrain(ixa,jya-1)
            ENDIF
            ter_there = 0.
            gradx_there = 0.
            grady_there = 0
            DO ix2 = ixmin, ixmax
               DO jy2 = jymin, jymax
                  IF (conusmask(ix2,jy2).gt. 0) THEN
                     pa_there(:) = panal_quantiles(ix2, jy2,:)
                     ter_there = terrain(ix2,jy2)
                     IF (ix2 .ne. 1 .and. ix2 .ne. nxa .and. jy2 .ne. 1 .and. jy2 .ne. nya) THEN
                        gradx_there = terrain(ix2+1,jy2)-terrain(ix2-1,jy2)
                        grady_there = terrain(ix2,jy2+1)-terrain(ix2,jy2-1)
                     ENDIF
                     ! --- quantify the difference between forecast distributions at 
                     !     the two grid points
                     CALL quantify_difference (npct, pctv, pa_here, pa_there, padiff)
                     sum_pa2 = sum_pa2 + padiff**2
                     sum_pa  = sum_pa  + padiff
                     sum_t2  = sum_t2  + ter_there**2
                     sum_t   = sum_t   + ter_there
                     sum_gradx  = sum_gradx + gradx_there
                     sum_gradx2 = sum_gradx2 + gradx_there**2
                     sum_grady  = sum_grady + grady_there
                     sum_grady2 = sum_grady2 + grady_there**2
                     nsamps  = nsamps  + 1
                  ENDIF
               END DO
            END DO

            ! --- now calculate std deviation by shortcut (Wilks 2006, eq. 3.20)
            pa_mean(ixa,jya)   = sum_pa / REAL(nsamps)
            pa_stddev(ixa,jya) = SQRT((sum_pa2 - REAL(nsamps)*pa_mean(ixa,jya)**2) / REAL(nsamps-1))
			pa_stddev(ixa,jya) = MAX(pa_stddev(ixa,jya),0.01)
            ter_mean(ixa,jya)    = sum_t  / REAL(nsamps)
            ter_stddev(ixa,jya)  = SQRT((sum_t2  - REAL(nsamps)*ter_mean(ixa,jya)**2) / REAL(nsamps-1))
			ter_stddev(ixa,jya) = MAX(ter_stddev(ixa,jya),1.0)
            gradx_mean(ixa,jya)    = sum_gradx  / REAL(nsamps)
            gradx_stddev(ixa,jya)  = SQRT((sum_gradx2  - REAL(nsamps)*gradx_mean(ixa,jya)**2) / REAL(nsamps-1))
			gradx_stddev(ixa,jya)  = MAX(0.0001,gradx_stddev(ixa,jya))
            grady_mean(ixa,jya)    = sum_grady  / REAL(nsamps)
            grady_stddev(ixa,jya)  = SQRT((sum_grady2  - REAL(nsamps)*grady_mean(ixa,jya)**2) / REAL(nsamps-1))
			grady_stddev(ixa,jya)  = MAX(0.0001,grady_stddev(ixa,jya))
         ELSE
            pa_mean(ixa,jya)   = -99.99
            pa_stddev(ixa,jya) = -99.99
            ter_mean(ixa,jya)  = -99.99
            ter_stddev(ixa,jya) = -99.99
            gradx_mean(ixa,jya) = -99.99
            gradx_stddev(ixa,jya) = -99.99
            grady_mean(ixa,jya) = -99.99
            grady_stddev(ixa,jya) = -99.99
         ENDIF

      END DO ! jya
   END DO    ! ixa

   PRINT *,' ------  pa_mean ------'
   DO j = nya/2+10,nya/2-10,-1
      PRINT 123, j,(pa_mean(i,j),i=nxa/2-10,nxa/2+10)
      123 format(i3,1x,21(f6.2,1x))
   END DO

   PRINT *,' ------  pa_stddev ------'
   DO j = nya/2+10,nya/2-10,-1
      PRINT 123, j,(pa_stddev(i,j),i=nxa/2-10,nxa/2+10)
   END DO

   PRINT *,' ------  ter_mean ------'
   DO j = nya/2+10,nya/2-10,-1
      PRINT 123, j,(ter_mean(i,j),i=nxa/2-10,nxa/2+10)
   END DO

   PRINT *,' ------  ter_stddev ------'
   DO j = nya/2+10,nya/2-10,-1
      PRINT 123, j,(ter_stddev(i,j),i=nxa/2-10,nxa/2+10)
   END DO

   ! ---- find locations of supplemental locations that are the best fit for each grid point

   PRINT *, 'finding supplemental locations'
   DO ixa = 1, nxa
      IF (MOD(ixa,10) .eq. 0) PRINT *,'processing ixa = ',ixa,' of ',nxa
      DO jya = 1, nya

         maskout(:,:) = 0
         difference(:,:) = 99999999.
         IF (conusmask(ixa,jya) .gt. 0) THEN ! grid point inside US

            ter_here = 0.
            gradx_here = 0.
            grady_here = 0.
            IF (ixa .eq. 1 .or. ixa .eq. nxa .or. jya .eq. 1 .or. jya .eq. nya) THEN
               tweight = 0.0
               gxweight = 0.0
               gyweight = 0.0
            ELSE 
               tweight = ter_weight
               gxweight = gradx_weight
               gyweight = grady_weight
               ter_here = terrain(ixa,jya)
               gradx_here = terrain(ixa+1,jya) - terrain(ixa-1,jya)
               grady_here = terrain(ixa,jya+1) - terrain(ixa,jya-1)
            ENDIF
            pa_here(:) = panal_quantiles(ixa, jya,:)

            imin = max(1,ixa-maxseparation)
            imax = min(nxa,ixa+maxseparation)
            jmin = max(1,jya-maxseparation)
            jmax = min(nya,jya+maxseparation)

            DO ix2 = imin, imax
               DO jy2 = jmin, jmax

                  ! ---- only consider this grid point if inside CONUS, and is less than
                  !      max separation, computed in 1/8-degree grid points.

                  d2 = REAL((ixa-ix2)**2 + (jya-jy2)**2)
                  dist = SQRT(d2)


                  IF (conusmask(ix2,jy2) .gt. 0 .and. maskout(ix2,jy2) .eq. 0 .and. &
                  dist .le. maxseparation) THEN

                     pa_there(:) = panal_quantiles(ix2, jy2,:)

                     ter_there = 0.0
                     gradx_there = 0.0
                     grady_there = 0.0

                     IF (ix2 .ne. 1 .and. ix2 .ne. nxa .and. jy2 .ne. 1 .and. jy2 .ne. nya) THEN
                        ter_there = terrain(ix2,jy2)
                        gradx_there = terrain(ix2+1,jy2) - terrain(ix2-1,jy2)
                        grady_there = terrain(ix2,jy2+1) - terrain(ix2,jy2-1)
                     ENDIF

                     ! --- quantify the difference between forecast distributions at 
                     !     the two grid points 
                     CALL quantify_difference (npct, pctv, pa_here, pa_there, padiff)

                     ! --- using information on the difference between forecast CDFs
                     !     (pfdiff), and F/O regression relationship information
                     !     (b0_mean, b1_mean), determine the single number that 
                     !     quantifies the overall strength of relationship between
                     !     forecast characteristics at these two grid points.
                     CALL normify (nxa, nya, ixa, jya, dist, distpenalty, &
                          pa_stddev, padiff,  paweight, &
                          ter_stddev, ter_here, ter_there, tweight, &
                          gradx_stddev, gradx_here, gradx_there, gxweight, &
                          grady_stddev, grady_here, grady_there, gyweight, &
                          difference(ix2,jy2))

                  ENDIF
               END DO ! jy2
            END DO ! ix2
            
            DO iana = 1, nanalogs  ! number of supplemental locations

               IF (iana .eq. 1) THEN ! the first supplemental location is simply that orginal grid point
                  xlocationa(ixa,jya,iana) = ixa
                  ylocationa(ixa,jya,iana) = jya
                  ! ---- don't consider points right around grid pt of interest,
                  !      too strong a correlation (want quasi-independent samples)
                  CALL mask_around_thisgridpt(ixa, jya, nxa, nya, minseparation, maskout)
               ELSE
                  ! ---- now find the analysis grid point with the next closest similarity. Set this 
                  !      as supplemental location point and then mask around this to eliminate nearby
                  !      points from any future consideration

                  minx = 0
                  miny = 0
                  diffmin = 999999.
                  DO ix2 = 1, nxa
                     DO jy2 = 1, nya
                        IF (difference(ix2,jy2) .lt. diffmin .and. maskout(ix2,jy2) .eq. 0 &
                        .and. conusmask(ix2,jy2) .gt. 0) THEN
                           diffmin = difference(ix2,jy2)
                           minx = ix2
                           miny = jy2
                        ENDIF
                     END DO ! jy2
                  END DO  ! ix2
               
                  ! ---- finally, define the supplemental location on the forecast grid to be the 
                  !      grid point that had this closest fit as defined above.
                  !
                  xlocationa(ixa,jya,iana) = minx
                  ylocationa(ixa,jya,iana) = miny
               
                  ! ---- make sure no other grid points very close to this new supplemental 
                  !      location are considered as another potential supplemental location.

                  CALL mask_around_thisgridpt(minx, miny, nxa, nya, minseparation, maskout)
               ENDIF ! iana > 1
            END DO ! iana
			
			lmin = MINVAL(xlocationa(ixa,jya,:))
			IF (lmin .lt. 0) THEN
				PRINT *,'For location', ixa, jya,' there is a missing value'
				PRINT *,'xlocationa = ',xlocationa(ixa,jya,:)
				PRINT *,'ylocationa = ',ylocationa(ixa,jya,:)
				STOP
			ENDIF

         ENDIF   ! conusmask
      END DO  ! jy
   END DO ! ix
   
   ! --- save out data for this month

   outfile = TRIM(output_data_directory)//'/ppn_analog_locns_ccpa5_'//cmonth(imonth)//'_'// &
        cleadb//'_to_'//cleade//'.dat'
   PRINT *, 'writing to ', TRIM(outfile)
   OPEN (UNIT=1, FILE=outfile, STATUS='replace', FORM='unformatted')
   WRITE (1) nxf, nyf, nxa, nya, nanalogs
   WRITE (1) xlocationa
   WRITE (1) ylocationa
   WRITE (1) conusmask
   WRITE (1) lonsa
   WRITE (1) latsa
   CLOSE (1)

END DO  ! imonth
PRINT *,'Done'
CONTAINS                        

! ====================================================================================

SUBROUTINE mask_around_thisgridpt(minx, miny, nx, ny, minseparation, maskout)
! set maskout array around grid point of interest to 1 to indicate nonconsideration
!    of this in the future.
INTEGER, INTENT(IN) :: minx, miny, nx, ny, minseparation
INTEGER*2, INTENT(INOUT), DIMENSION(nx,ny) :: maskout
imin = MAX(1,minx-minseparation)
imax = MIN(nx,minx+minseparation)
jmin = MAX(1,miny-minseparation)
jmax = MIN(ny,miny+minseparation)
DO i = imin, imax
   DO j = jmin, jmax
      dist2 = REAL((minx-i)**2 + (miny-j)**2)
      dist = SQRT(dist2)
      IF (dist .le. minseparation) maskout(i,j) = 1
   END DO
END DO
RETURN
END SUBROUTINE mask_around_thisgridpt

! ====================================================================================

SUBROUTINE quantify_difference (npct, pctv, pa_here, pa_there, padiff )
! quantify the differences in the precipitation fcst amounts at two grid points 
! associated with pre-specified quantiles of the distributions.  
INTEGER, INTENT(IN) :: npct
REAL, INTENT(IN), DIMENSION(npct) :: pa_here, pa_there, pctv
REAL, INTENT(OUT) :: padiff

imin = 1
imax = npct
DO i = 1,npct
   IF (pa_here(i) .gt. 0.) THEN
      imin = i
      GOTO 123
   ENDIF
END DO

123 CONTINUE
DO i = npct,1,-1
   IF (pctv(i) .le. 0.95) THEN
      imax = i
      GOTO 124
   ENDIF
END DO

124 CONTINUE

IF (imax .le. imin) THEN
   imax = npct
   imin = 1
ENDIF

IF (imin .lt. 1 .or. imax .gt. npct) THEN
   PRINT *,'imax, imin = ',imax, imin
   PRINT *,'pa_here = ',pa_here
   PRINT *,'pctv = ',pctv
ENDIF

padiff = SUM(ABS(pa_here(imin:imax) - pa_there(imin:imax))) / REAL(imax-imin+1)
RETURN
END SUBROUTINE quantify_difference

! ====================================================================================

SUBROUTINE  normify (nxa, nya, ixa, jya, dist, distpenalty, &
    pa_stddev, padiff,  paweight, &
    ter_stddev, ter_here, ter_there, ter_weight, &
    gradx_stddev, gradx_here, gradx_there, gradx_weight, &
    grady_stddev, grady_here, grady_there, grady_weight, &
    difference) 

! determine a single number that quantifies the similarity between F, A information
! at two grid points.  Effectively this is sort of like a Mahalanobis norm together with a
! user-defined weighting.

INTEGER, INTENT(IN) :: nxa, nya, ixa, jya
REAL, INTENT(IN) :: distpenalty, dist
REAL, INTENT(IN), DIMENSION(nxa,nya) :: pa_stddev
REAL, INTENT(IN) :: padiff, paweight
REAL, INTENT(IN), DIMENSION(nxa,nya) :: ter_stddev
REAL, INTENT(IN) :: ter_here, ter_there, ter_weight
REAL, INTENT(IN), DIMENSION(nxa,nya) :: gradx_stddev
REAL, INTENT(IN) :: gradx_here, gradx_there, gradx_weight
REAL, INTENT(IN), DIMENSION(nxa,nya) :: grady_stddev
REAL, INTENT(IN) :: grady_here, grady_there, grady_weight

REAL, INTENT(OUT) :: difference

ter_Z   = ABS(ter_here - ter_there) / ter_stddev(ixa,jya)
gradx_Z = ABS(gradx_here - gradx_there) / gradx_stddev(ixa,jya)
grady_Z = ABS(grady_here - grady_there) / grady_stddev(ixa,jya)

difference = (paweight*padiff / pa_stddev(ixa,jya) + &
     ter_weight*ter_Z + gradx_weight*gradx_Z + grady_weight*grady_Z) / &
     (paweight + ter_weight + gradx_weight + grady_weight) + dist*distpenalty

RETURN
END SUBROUTINE normify

! ====================================================================================

END PROGRAM compute_precip_analog_locations_ccpa5
