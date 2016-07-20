SUBROUTINE direct_forecast (nxf, nyf, nxa, nya, ntrainsamps, nmembers, nmaxlist, &
   apcp_fcst_ens, apcp_anal, iccpa_ktr, iccpa_mask, imonth, date_anal, &
   iccpa_xlist, iccpa_ylist, thresh_median, thresh_uppertercile, &
   thresh_upperquintile, thresh_upperdecile, thresh_95, prob_pop, prob_1mm, &
   prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, prob_50mm, prob_median, &
   prob_uppertercile, prob_upperquintile, prob_upperdecile, prob_95)

! ---- conduct the direct ensemble forecast procedure.  In this application,
!      the analyzed data is on a higher resolution grid, so for a given forecast 
!      point, we find all the analyzed grid points which are closest to that forecast 
!      point and construct the forecast ensemble using the closest forecast data
!
!      coded by Tom Hamill, Sep 2013, tom.hamill@noaa.gov, (303) 497-3060.

INTEGER, INTENT(IN) :: nxf, nyf, nxa, nya, ntrainsamps, nmembers, nmaxlist, imonth
INTEGER, INTENT(IN), DIMENSION(nxf,nyf,nmaxlist) :: iccpa_xlist  
    ! list of which ccpa grid points are closest
INTEGER, INTENT(IN), DIMENSION(nxf,nyf,nmaxlist) :: iccpa_ylist  
    ! to this particular forecast lat/lon point
INTEGER, INTENT(IN), DIMENSION(nxf,nyf)          :: iccpa_ktr    
    ! how many analysis pts are closest to this fcst point
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya)        :: iccpa_mask   ! a mask for US.
INTEGER, INTENT(IN), DIMENSION(ntrainsamps)      :: date_anal

REAL, INTENT(IN), DIMENSION(nxa, nya, ntrainsamps)  :: apcp_anal ! precip anal on 0.125-deg grid
REAL, INTENT(IN), DIMENSION(nxf, nyf, nmembers, ntrainsamps)  :: apcp_fcst_ens ! precipitation fcst mean

REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_median  ! 50th %ile
REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_uppertercile ! 67th %ile
REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_upperquintile ! 80th %ile
REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_upperdecile  ! 90th %ile
REAL, INTENT(IN), DIMENSION(nxa, nya) :: thresh_95 ! 95th %ile
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

REAL, DIMENSION(nmembers)            :: ensfcst  ! 1-d vector produced by this calibration algorithm

! ---- initialize to missing data value

prob_1mm       = -99.99 
prob_2p5mm     = -99.99
prob_5mm       = -99.99
prob_10mm      = -99.99
prob_25mm      = -99.99
prob_50mm      = -99.99
prob_median        = -99.99
prob_uppertercile  = -99.99
prob_upperquintile = -99.99
prob_upperdecile   = -99.99
prob_pop           = -99.99
prob_95            = -99.99

! ---- process all grid points on the forecast grid

!DO ixf = 3,3
DO ixf = 1, nxf
   CALL cpu_time(t1)
   PRINT *,'ixf = ',ixf,' of ',nxf,'. CPU time = ',t1
   istart = ixf - nboxsize 
   iend   = ixf + nboxsize
   IF (istart .lt. 1) istart = 1
   IF (iend .gt. nxf) iend = nxf
   nxs = iend - istart + 1
   DO jyf = 1, nyf
   !DO jyf = 50,50

      ! ---- see if there are any CONUS analysis points contained within this forecast box by 
      !      count of analysis points nearest to this fcst grid pt.  Flag igood = 1 if there are
      !      any such points over CONUS

      !print *,'ixf, jyf = ',ixf, jyf
      igood = 0
      IF (iccpa_ktr(ixf,jyf) .GT. 0) THEN
         DO ik = 1, iccpa_ktr(ixf,jyf)
            IF (iccpa_mask(iccpa_xlist(ixf,jyf,ik),iccpa_ylist(ixf,jyf,ik)) .gt. 0) igood=1
         END DO
      END IF

      ! ---- only bother running the analog machinery for this forecast grid point 
      !      if there are indeed analysis grid points nearby it.

      !print *,'igood = ',igood
      IF (igood .eq. 1) THEN

         ! ---- define indices for the grid box of data that we will seek pattern matches for.
               
         CALL cpu_time(t1)

         ! ---- loop thru all the days in the data set and determine how different every
         !      other day's forecast is from today's forecast

         DO iday = 1, ntrainsamps
            !print *,'date_anal(iday) = ',date_anal(iday)
            iyyyy     = date_anal(iday) / 1000000
            immddhh   = date_anal(iday) - iyyyy*1000000
            imm       = immddhh/10000
            !print *,'iday = ', iday, 'iyyyy immddhh, imm = ',iyyyy,immddhh, imm
            IF (imm .eq. imonth) THEN

               ! ---- produce raw ensemble forecast

               !print *,'iccpa_ktr(ixf,jyf)  =',iccpa_ktr(ixf,jyf)
               DO ianal = 1, iccpa_ktr(ixf,jyf)
                  ixa = iccpa_xlist(ixf,jyf,ianal)
                  jya = iccpa_ylist(ixf,jyf,ianal)
                  !print *,'ixa, jya= ',ixa,jya
                  IF (iccpa_mask(ixa, jya) .gt. 0) THEN
                     DO iana = 1, nmembers
                        ensfcst(iana) = apcp_fcst_ens(ixf,jyf,iana,iday)
                     END DO

                     ! ---- now simply set the ensemble probabilities from the relative frequency of
                     !      the analogs of observed

                     CALL probfcst_relfreq (nmembers, 0.1,  ensfcst,  prob_pop  (ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, 1.0,  ensfcst,  prob_1mm  (ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, 2.5,  ensfcst,  prob_2p5mm(ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, 5.0,  ensfcst,  prob_5mm  (ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, 10.0, ensfcst,  prob_10mm (ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, 25.0, ensfcst,  prob_25mm (ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, 50.0, ensfcst,  prob_50mm (ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, thresh_median(ixa,jya), ensfcst, &
                          prob_median(ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, thresh_uppertercile(ixa,jya), ensfcst, &
                          prob_uppertercile(ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, thresh_upperquintile(ixa,jya), ensfcst, &
                          prob_upperquintile(ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, thresh_upperdecile(ixa,jya), ensfcst, &
                          prob_upperdecile (ixa,jya,iday))
                     CALL probfcst_relfreq (nmembers, thresh_95(ixa,jya), ensfcst, prob_95(ixa,jya,iday))
                  ELSE
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
                  END IF
               END DO   ! ianal = 1, iccpa_ktr(ix,jy)
            ELSE     ! imm .ne. imonth
               DO ianal = 1, iccpa_ktr(ixf,jyf)
                  ixa = iccpa_xlist(ixf,jyf,ianal)
                  jya = iccpa_ylist(ixf,jyf,ianal)
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
               END DO
            ENDIF
         END DO      ! iday = 1, ntrainsamps
      END IF         ! igood = 1
   END DO            ! jy = 1, nyf
END DO               ! ix = 1, nxf
RETURN
END SUBROUTINE direct_forecast
