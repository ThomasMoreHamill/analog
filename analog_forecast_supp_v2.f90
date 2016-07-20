SUBROUTINE analog_forecast_supp_v2 (nxf, nyf, nxa, nya, ntrainsamps, nanalogs, &
   nsupp, nsuppmax, imonth, window_size, order, ileade, apcp_fcst, pwat_fcst_mean, apcp_anal, &
   iccpa_mask, date_anal, inearest, jnearest, rlons_fcst, rlats_fcst, &
   thresh_median, thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, &
   thresh_95, xlocation, ylocation, rlons_anal, rlats_anal, input_data_directory, &
   prob_pop, prob_1mm, prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, &
   prob_50mm, prob_median, prob_uppertercile, &
   prob_upperquintile, prob_upperdecile, prob_95, fmean_analog)

! ---- conduct the analog forecast procedure to generate post-processed probabilistic
!      precipitation forecasts.  For each model forecast grid point, 
!      extract the forecast at that point.  Compare "past" 
!      (with the cross validation procedure, it could actually be future as well) forecasts 
!      to these forecasts.  Determine which past forecasts have the closest 
!      match to today's forecast, and form an ensemble from the analyzed data 
!      associated with the dates of the past matching forecasts.  In this application,
!      the analyzed data is on a higher resolution grid, so for a given forecast 
!      point, we find all the analyzed grid points which are closest to that forecast 
!      point and construct the analyzed ensemble using the closest dates associated 
!      with that nearby forecast grid point.
!
!      This version uses extra training data from supplemental data locations, so that
!      we have more past forecasts to choose from that may resemble today's forecast.
!
!      calling tree:
!         cpu_time
!         interpolate_training_data
!         R_rnkpar
!         probfcst_relfreq
!         sgolay_2d_weights
!             mpinv
!         raw_vs_smoothed_weight
!         sgolay_smooth

!      coded by Tom Hamill, tom.hamill@noaa.gov, (303) 497-3060.



INTEGER, INTENT(IN) :: nxf, nyf, nxa, nya, ntrainsamps, window_size, order, imonth
INTEGER, INTENT(IN) :: nanalogs, nsupp, nsuppmax, ileade
INTEGER, INTENT(IN), DIMENSION(ntrainsamps)      :: date_anal 
    ! in yyyymmddhh format, initialization time
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya)        :: iccpa_mask   ! a mask for US.
INTEGER, INTENT(IN), DIMENSION(nxa,nya) :: inearest ! which fcst (i,j) is closest to this
INTEGER, INTENT(IN), DIMENSION(nxa,nya) :: jnearest ! analysis (i,j)
INTEGER, INTENT(IN), DIMENSION(nxa,nya,nsuppmax) :: xlocation ! supplemental data locations
INTEGER, INTENT(IN), DIMENSION(nxa,nya,nsuppmax) :: ylocation

CHARACTER*(*), INTENT(IN) :: input_data_directory

REAL, INTENT(IN), DIMENSION(nxa, nya, ntrainsamps)  :: apcp_anal ! precip anal on 0.125-deg grid
REAL, INTENT(IN), DIMENSION(nxf, nyf, ntrainsamps)  :: apcp_fcst ! precipitation fcst mean
REAL, INTENT(IN), DIMENSION(nxf, nyf, ntrainsamps)  :: pwat_fcst_mean ! precipitable water fcst mean

REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)   :: fmean_analog ! mean of the analog forecasts
REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_median  ! 50th %ile
REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_uppertercile ! 67th %ile
REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_upperquintile ! 80th %ile
REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_upperdecile  ! 90th %ile
REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_95 ! 95th %ile
REAL, INTENT(IN), DIMENSION(nxf, nyf) :: rlons_fcst
REAL, INTENT(IN), DIMENSION(nxf, nyf) :: rlats_fcst
REAL, INTENT(IN), DIMENSION(nxa, nya) :: rlons_anal
REAL, INTENT(IN), DIMENSION(nxa, nya) :: rlats_anal
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_1mm ! prob forecast for > 1mm event
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_2p5mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_5mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_10mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_25mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_50mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_median 
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_uppertercile
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_upperquintile
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_upperdecile
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_95
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps)  :: prob_pop

! --- local variables 

INTEGER, DIMENSION(ntrainsamps*nsupp) :: indexsort    ! vector indices to be sorted alongside rank differences
INTEGER, DIMENSION(ntrainsamps*nsupp) :: indexsort2    ! vector indices to be sorted alongside rank differences
INTEGER, DIMENSION(ntrainsamps*nsupp) :: isuppsort
INTEGER, DIMENSION(ntrainsamps*nsupp) :: isuppsort2

REAL, DIMENSION(ntrainsamps*nsupp)     :: diffsort ! sorted array of forecast differences (see also indexsort)
REAL, DIMENSION(nanalogs)              :: ensfcst  ! 1-d vector produced by this calibration algorithm
REAL, DIMENSION(nsupp,ntrainsamps)     :: precip_interp
REAL, DIMENSION(nsupp,ntrainsamps)     :: pwat_interp
REAL, DIMENSION(nxa,nya)               :: raw_weight
REAL, DIMENSION(window_size,window_size) :: weights ! used in Savitzky-Golay
REAL, DIMENSION(nxa,nya)               :: work ! used in Savitzky-Golay stuff

REAL mtfact

! ---- initialize output forecast probabilities to missing data value

prob_1mm           = -99.99 
prob_2p5mm         = -99.99
prob_5mm           = -99.99
prob_10mm          = -99.99
prob_25mm          = -99.99
prob_50mm          = -99.99
prob_median        = -99.99
prob_uppertercile  = -99.99
prob_upperquintile = -99.99
prob_upperdecile   = -99.99
prob_pop           = -99.99
prob_95            = -99.99
fmean_analog       = -99.99
iseed = -12345 ! seed for random number generator

! ---- initialize the Savitzky-Golay smoother weights that are used to smooth the 
!      probability forecast

istat = 0
PRINT *,'analog_forecast_supp_v2 calling sgolay_2d_weights', window_size, order
CALL sgolay_2d_weights(window_size, order, istat, weights)

! ---- determine how much we will weight the Savitzky-Golay smoothed fields vs. 
!      the raw input.  It makes sense to weight the former more in regions
!      where there is not much topographic variation, and the latter where
!      there is.

PRINT *, 'analog_forecast_supp_v2 calling raw_vs_smoothed_weight'
CALL raw_vs_smoothed_weight(nxa, nya, nxf, nyf, rlons_fcst, rlats_fcst, &
    rlons_anal, rlats_anal, input_data_directory, raw_weight)
	
! ---- determine the number of analogs to use for estimating probabilities.  Previously 
!      (Hamill, T. M., J. S. Whitaker, and S. L. Mullen, 2006: Reforecasts, an important 
!      dataset for improving weather predictions. Bull. Amer. Meteor. Soc., 87,33-46.)
!      it was shown that the optimal number of analogs depended on the rarity of the 
!      event as well as the forecast lead time.  Here we only address lead-time dependence.
!
!      Also note that the number of analogs here was tuned somewhat by 
!      trial and error to give a reasonably result when 10 supplemental locations are
!      used.  Fewer analogs should be used if fewer or no supplemental locations
!      are used.

if (ileade .le. 24) mtfact = 1.0
if (ileade .gt. 24 .and. ileade .le. 48)  mtfact = 1.5
if (ileade .gt. 48 .and. ileade .le. 96)  mtfact = 2.0
if (ileade .gt. 96 .and. ileade .le. 120) mtfact = 3.0
if (ileade .gt. 120) mtfact = 4.0
nana = NINT(mtfact*50.)
IF (nana .gt. ntrainsamps) nana = ntrainsamps

! ---- process all grid points on the analysis grid

DO ixa = 1, nxa
   CALL cpu_time(t1)
   ihigh = MAX(1,ixa-1)
   ipart = SUM(INT(iccpa_mask(1:ihigh,:)))
   ifull = SUM(INT(iccpa_mask))
   IF (MOD(ixa,10) .eq. 0) PRINT *,'ixa = ',ixa,' of ',nxa,'. CPU time = ',&
   	  t1,'. Fraction done = ', REAL(ipart) / REAL(ifull)
   DO jya = 1, nya

      ! ---- only bother running the analog machinery for this analysis grid point 
      !      is inside conus

      IF (iccpa_mask(ixa,jya) .eq. 1) THEN
               
         CALL cpu_time(t1)
		 
		 ! ---- interpolate the mean forecast to the analysis location for this particular analysis 
		 !      location and for the supplemental locations associated with this analysis grid pt.

		 !PRINT *,'ixa,jya, minval(inearest) = ', ixa, jya, minval(inearest)

         CALL interpolate_training_data(nxf, nyf, nxa, nya, ntrainsamps, nsupp, nsuppmax, ixa, jya, &
              xlocation, ylocation, rlons_fcst, rlats_fcst, rlons_anal, rlats_anal, &
              inearest, jnearest, apcp_fcst, precip_interp)

         CALL interpolate_training_data(nxf, nyf, nxa, nya, ntrainsamps, nsupp, nsuppmax, ixa, jya, &
              xlocation, ylocation, rlons_fcst, rlats_fcst, rlons_anal, rlats_anal, &
              inearest, jnearest, pwat_fcst_mean, pwat_interp)

         ! ---- loop thru all the days in the data set and determine how different every
         !      other day's forecast is from today's forecast

         DO iday = 1, ntrainsamps
	     
            iyyyy     = date_anal(iday) / 1000000
            immddhh   = date_anal(iday) - iyyyy*1000000
            imm       = immddhh/10000
            IF (apcp_anal(ixa,jya,iday) .ge. 0.0 .and. imm .eq. imonth) THEN

               ! ---- compare other forecast training samples to this one, at least 
			   !      if they aren't from the same year of data.

               indexsortktr = 1
               diffsort(:) = 9999999999.
               ktrday = 1
               DO jday = 1, ntrainsamps

                  ! ---- extract year, month for samples to compare

                  jyyyy     = date_anal(jday) / 1000000
                  jmmddhh   = date_anal(jday) - jyyyy*1000000
                  jmm       = jmmddhh/10000

                  ! ---- only bother doing the comparison for samples that aren't 
                  !      in the same calendar year, and if Dec or Jan, aren't in
                  !      the month from the subsequent or prior year. 

                  iok = 1  ! 1=yes, consider; 0 = no, don't consider
                  IF (iyyyy .eq. jyyyy) iok = 0 ! same year
                  IF (iyyyy .lt. 1985 .or. jyyyy .lt. 1985) iok = 0
                  IF (jyyyy .eq. iyyyy-1 .and. jmm .eq. 12 .and. imm .eq. 1)  iok = 0
                  IF (jyyyy .eq. iyyyy+1 .and. jmm .eq. 1  .and. imm .eq. 12) iok = 0
                  IF (imm .lt. 1 .or. imm .gt. 12 .or. jmm .lt. 1 .or. jmm .gt. 12) iok = 0
                  IF (iok .eq. 0) THEN

                     ! ---- for cross validation purposes, flag these forecasts as having
                     !      extremely high error so that they effectively eliminated from
                     !      consideration as analog dates

                     diffsort(ktrday) = 9999999999.
					 
                  ELSE
                     ! ---- here we determine how different the jday's forecast is from
                     !      the iday's forecast.  At this point determining how different the 
                     !      forecasts are from each other is done only on the difference in
                     !      interpolated mean precip and pwat.  Use data from the original 
                     !      forecast location (isupp=1) and from nsupp-1 supplementary locations.
					 
                     DO isupp = 1, nsupp
                        ixl = xlocation(ixa,jya,isupp)
                        jyl = ylocation(ixa,jya,isupp)
                        IF (ixl .gt. 0 .and. jyl .gt. 0) THEN
                           diffsort(ktrday) = &
                              0.7*ABS(precip_interp(1,iday)- precip_interp(isupp,jday)) + &
                              0.3*ABS(pwat_interp(1,iday)- pwat_interp(isupp,jday))
                           indexsort(ktrday) = jday
                           isuppsort(ktrday) = isupp

                           ! ---- add small random number to deal with the possibility of ties
                        
                           IF (diffsort(ktrday) .GT. 10) THEN
                              diffsort(ktrday) = diffsort(ktrday) + (ran3(iseed)-0.5)*0.05
                           ELSE IF (diffsort(ktrday) .LE. 10 .and. diffsort(ktrday) .GT. 2) THEN
                              diffsort(ktrday) = diffsort(ktrday) + (ran3(iseed)-0.5)*0.02
                           ELSE IF (diffsort(ktrday) .LE. 2 .AND. diffsort(ktrday) .GT. 0.5) THEN
                              diffsort(ktrday) = diffsort(ktrday) + (ran3(iseed)-0.5)*0.006
                           ELSE IF (diffsort(ktrday) .LE. 0.5 .AND. diffsort(ktrday) .GT. 0.1) THEN
                              diffsort(ktrday) = diffsort(ktrday) + (ran3(iseed)-0.5)*0.003
                           ELSE
                              diffsort(ktrday) = diffsort(ktrday) + (ran3(iseed)-0.5)*0.001
                           ENDIF
                           ktrday = ktrday + 1
                        END IF ! ixl > 0, jxl > 0
                     END DO ! isupp
                  END IF ! iok
               END DO  ! jday
               ktrday = ktrday-1

               ! ---- now sort lowest to highest, and use associated sorted indices to indicate
               !      which elements of arrays to extract later.  

               CALL R_rnkpar(diffsort, indexsort2, ktrday, nana)
               DO i = 1, nana
                  isuppsort2(i) = isuppsort(indexsort2(i))
               END DO

               ! ---- form an ensemble of analogs of the observed on the dates of other 
               !      forecasts that have the closest pattern match to today's forecast

               nmax = 0
               DO iana = 1, nana

                  ! ---- get the supplemental index number, and from that the associated 
                  !      analysis grid point indices

                  isno = isuppsort2(iana)
                  ixsupp = xlocation(ixa,jya,isno)
                  jysupp = ylocation(ixa,jya,isno)
                     
                  ! ---- fetch the analysis associated with this particular date that 
                  !      had a close forecast (the date defined by indexsort(indexsort2(iana))

                  pa = apcp_anal(ixsupp, jysupp, indexsort(indexsort2(iana)))
                  IF (ixsupp .ge. 1 .and. ixsupp .le. nxa .and. jysupp .ge. 1 &
                       .and. jya .le. nya .and. pa .ge. 0) THEN
                     nmax = nmax+1
                     ensfcst(nmax) = pa
                  END IF
               END DO

               ! ---- now simply set the ensemble probabilities from the relative frequency of
               !      the analogs of observed

               IF (nmax .ge. 10) THEN

                  CALL probfcst_relfreq (nmax, 0.1,  ensfcst,  prob_pop  (ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, 1.0,  ensfcst,  prob_1mm  (ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, 2.5,  ensfcst,  prob_2p5mm(ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, 5.0,  ensfcst,  prob_5mm  (ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, 10.0, ensfcst,  prob_10mm (ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, 25.0, ensfcst,  prob_25mm (ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, 50.0, ensfcst,  prob_50mm (ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, thresh_median(ixa,jya), ensfcst, &
                    prob_median(ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, thresh_uppertercile(ixa,jya), ensfcst, &
                    prob_uppertercile(ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, thresh_upperquintile(ixa,jya), ensfcst, &
                    prob_upperquintile(ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, thresh_upperdecile(ixa,jya), ensfcst, &
                    prob_upperdecile (ixa,jya,iday))
                  CALL probfcst_relfreq (nmax, thresh_95(ixa,jya), ensfcst, prob_95(ixa,jya,iday))
                  fmean_analog(ixa,jya,iday) = SUM(ensfcst(1:nana))/REAL(nana)
               ELSE  ! not enough analogs; set to missing.
                  prob_pop(ixa,jya,iday) = -99.99
                  prob_1mm(ixa,jya,iday) = -99.99
                  prob_2p5mm(ixa,jya,iday) = -99.99
                  prob_5mm(ixa,jya,iday) = -99.99
                  prob_10mm(ixa,jya,iday) = -99.99
                  prob_25mm(ixa,jya,iday) = -99.99
                  prob_50mm(ixa,jya,iday) = -99.99
                  prob_median(ixa,jya,iday) = -99.99
                  prob_uppertercile(ixa,jya,iday) = -99.99
                  prob_upperquintile(ixa,jya,iday) = -99.99
                  prob_upperdecile(ixa,jya,iday) = -99.99
                  prob_95(ixa,jya,iday) = -99.99
               ENDIF
            ELSE  ! apcp_anal(ixa,jya,iday) .lt. 0.0 .and. imm .eq. imonth
               prob_pop(ixa,jya,iday) = -99.99
               prob_1mm(ixa,jya,iday) = -99.99
               prob_2p5mm(ixa,jya,iday) = -99.99
               prob_5mm(ixa,jya,iday) = -99.99
               prob_10mm(ixa,jya,iday) = -99.99
               prob_25mm(ixa,jya,iday) = -99.99
               prob_50mm(ixa,jya,iday) = -99.99
               prob_median(ixa,jya,iday) = -99.99
               prob_uppertercile(ixa,jya,iday) = -99.99
               prob_upperquintile(ixa,jya,iday) = -99.99
               prob_upperdecile(ixa,jya,iday) = -99.99
               prob_95(ixa,jya,iday) = -99.99
               fmean_analog(ixa,jya,iday) = -99.99
            END IF
         END DO      ! iday = 1, ntrainsamps
      END IF         ! IF (iccpa_mask(ixa,jya) .eq. 1) 
   END DO            ! jya = 1, nya
END DO               ! ixa = 1, nxa

! ---- now produce each forecast day's final probability forecasts from a 
!      linear combination of the smoothed and unsmoothed, using more unsmoothed
!      in mountainous regions.  There we judge that the small-scale variability
!      is more likely to be signal related to terrain than noise.

DO iday = 1, ntrainsamps
	
   IF (MOD(ixa,50) .eq. 0) PRINT *,'smoothing date_anal(iday)= ',&
   	   date_anal(iday), iday, ntrainsamps

   iyyyy     = date_anal(iday) / 1000000
   immddhh   = date_anal(iday) - iyyyy*1000000
   imm       = immddhh/10000
   
   ! ---- use a sample point in the middle of the domain to check to see if this day's 
   !      data is missing for this day.  Produce final linear combination of
   !      raw and smoothed probabilities if data is available for this day.
    
   IF (apcp_anal(nxa/2,nya/2,iday) .ge. 0.0 .and. imm .eq. imonth) THEN

	  CALL cpu_time(t1)

      ! ---- now do the Savitzky-Golay smoothing, and the blending of raw and smoothed

      work(:,:) = prob_pop(:,:,iday)
	  
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_pop(:,:,iday) = raw_weight(:,:)*work(:,:) + (1.-raw_weight(:,:))*prob_pop(:,:,iday)

      work(:,:) = prob_1mm(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_1mm(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_1mm(:,:,iday)

      work(:,:) = prob_2p5mm(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_2p5mm(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_2p5mm(:,:,iday)

      work(:,:) = prob_5mm(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_5mm(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_5mm(:,:,iday)

      work(:,:) = prob_10mm(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_10mm(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_10mm(:,:,iday)

      work(:,:) = prob_25mm(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_25mm(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_25mm(:,:,iday)
      
      work(:,:) = prob_50mm(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_50mm(:,:,iday)  = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_50mm(:,:,iday) 

      work(:,:) = prob_median(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_median(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_median(:,:,iday)

      work(:,:) = prob_uppertercile(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_uppertercile(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_uppertercile(:,:,iday)

      work(:,:) = prob_upperquintile(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_upperquintile(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_upperquintile(:,:,iday)

      work(:,:) = prob_upperdecile(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_upperdecile(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	  	raw_weight(:,:)*prob_upperdecile(:,:,iday)

      work(:,:) = prob_95(:,:,iday)
      CALL sgolay_smooth(nxa, nya, work, iccpa_mask, weights, &
	  	window_size, order, istat)
      prob_95(:,:,iday) = (1.-raw_weight(:,:))*work(:,:) + &
	    raw_weight(:,:)*prob_95(:,:,iday)

   ENDIF

END DO

PRINT *,'finished in analog_forecast_supp_v2'

RETURN
END SUBROUTINE analog_forecast_supp_v2
