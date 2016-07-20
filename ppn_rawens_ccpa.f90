PROGRAM ppn_rawens_ccpa

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! raw ensemble forecasts validated against 0.125-degree CCPA grid to 
! produce calibrated probabilisitic total precipitation amount forecasts.   
!
! written by Tom Hamill, March 2014 (tom.hamill@noaa.gov), 303.497.3060

! usage:  > ppn_rawens_ccpa.x cleadb cleade ccompute
!    where cleadb = beginning lead time in hours (e.g., 012)
!    cleade = end lead time in hours (e.g., 024)
!    ccompute = Y when we're computing from scratch, N when reading
!    previously computed output from file
!
! calling tree:
!	load_precipfcst_and_ccpa
!   list_to_nearest_fcstgridpt
!   compute_quantiles_ccpa2
!   compute_all_climatologies_ccpa2
!   direct_forecast 
!   write_ccpa_analogs_tofile_direct2 
!   read_ccpa_analogs_fromfile_direct2 
!   calculate_reliabilities_rawens 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!USE netcdf

PARAMETER (nyears      = 14)  ! 2002-2015 here for the CCPA data
PARAMETER (ndaysub     = 91)  ! number of days of training data for each yr, 15th of month +/- 45 days
PARAMETER (nxa         = 464) ! CCPA grid dimensions, surrounding CONUS
PARAMETER (nya         = 224) !
PARAMETER (nxf         = 128) ! Grid dimensions of reforecast data, here on ~ 1/2 deg Gaussian grid 
PARAMETER (nyf         = 61)  !    surrounding CONUS
PARAMETER (ntrainsamps = ndaysub*nyears) ! number of training data samples loaded up
PARAMETER (nboxsize    = 5)   ! a box +/- nboxsize analysis grid points in each dir used for finding analogs
PARAMETER (nbs         = 2*nboxsize+1) ! grid pts in each direction for box
PARAMETER (nanalogs    = 100) ! max number of analogs used
PARAMETER (nmembers    = 11)  ! number of ensemble members in reforecast
PARAMETER (nmaxlist    = 16)  ! presumed maximum number of CCPA points that are closest to any forecast pt

INTEGER, DIMENSION(ntrainsamps)      :: date_anal ! in yyyymmddhh format, initialization time
INTEGER, DIMENSION(ntrainsamps)      :: date_fcst ! in yyyymmddhh format, end of fcst period
INTEGER, DIMENSION(nxf,nyf,nmaxlist) :: iccpa_xlist  ! list of which ccpa grid points are closest
INTEGER, DIMENSION(nxf,nyf,nmaxlist) :: iccpa_ylist  ! to this particular forecast lat/lon point
INTEGER, DIMENSION(nxf,nyf)          :: iccpa_ktr    ! how many analysis pts are closest to this fcst point
INTEGER*2, DIMENSION(nxa,nya)        :: iccpa_mask   ! a mask for US.

REAL, DIMENSION(nxa, nya, ntrainsamps)  :: apcp_anal ! precipitation analysis on 0.125-deg grid
REAL, DIMENSION(nxf, nyf, nmembers, ntrainsamps) :: apcp_fcst_ens ! precip fcst ensemble
REAL, DIMENSION(nxf, nyf, ntrainsamps)  :: apcp_fcst ! precip fcst mean
REAL, DIMENSION(nxa,nya)                :: climo_pop ! climatological probabilities for various
REAL, DIMENSION(nxa,nya)                :: climo_1mm ! thresholds
REAL, DIMENSION(nxa,nya)                :: climo_2p5mm
REAL, DIMENSION(nxa,nya)                :: climo_5mm
REAL, DIMENSION(nxa,nya)                :: climo_10mm
REAL, DIMENSION(nxa,nya)                :: climo_25mm
REAL, DIMENSION(nxa,nya)                :: climo_50mm
REAL, DIMENSION(nxa,nya)                :: climo_median
REAL, DIMENSION(nxa,nya)                :: climo_uppertercile
REAL, DIMENSION(nxa,nya)                :: climo_upperquintile
REAL, DIMENSION(nxa,nya)                :: climo_upperdecile
REAL, DIMENSION(nxa,nya)                :: climo_95
REAL, DIMENSION(nxa,nya)                :: climop
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_1mm ! output probability forecasts
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_2p5mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_5mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_10mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_25mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_50mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_median 
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_uppertercile
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_upperquintile
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_upperdecile
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_95
REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: prob_pop
REAL, DIMENSION(nxf,nyf,ntrainsamps)    :: pwat_fcst_mean ! total-col precip. water mean fcst
REAL, DIMENSION(nxf,nyf)                :: rlons_fcst
REAL, DIMENSION(nxf,nyf)                :: rlats_fcst ! forecast lat/lon grids
REAL, DIMENSION(nxa,nya)                :: rlons_anal 
REAL, DIMENSION(nxa,nya)                :: rlats_anal ! 0.125-degree CCPA analysis grid
REAL, DIMENSION(nxa,nya)                :: thresh_95 ! thresholds for various quantiles of climo dist
REAL, DIMENSION(nxa,nya)                :: thresh_uppertercile
REAL, DIMENSION(nxa,nya)                :: thresh_upperquintile
REAL, DIMENSION(nxa,nya)                :: thresh_upperdecile 
REAL, DIMENSION(nxa,nya)                :: thresh_median
REAL, DIMENSION(nxa,nya)                :: threshfield

CHARACTER*120 outfile, infile, outfile_dates, outfile_data
CHARACTER*120 outfile_determ, outfile_quantiles, infile_dates, infile_data
CHARACTER*80 input_data_directory, output_data_directory
CHARACTER*3 cleadb, cleade
CHARACTER*3, DIMENSION(12) :: cmonth
CHARACTER*1 ccompute

ALLOCATE(prob_1mm(nxa,nya,ntrainsamps),prob_2p5mm(nxa,nya,ntrainsamps),prob_5mm(nxa,nya,ntrainsamps),&
   prob_10mm(nxa,nya,ntrainsamps), prob_25mm(nxa,nya,ntrainsamps),prob_50mm(nxa,nya,ntrainsamps),&
   prob_median(nxa,nya,ntrainsamps), prob_uppertercile(nxa,nya,ntrainsamps),&
   prob_upperquintile(nxa,nya,ntrainsamps), prob_upperdecile(nxa,nya,ntrainsamps),&
   prob_95(nxa,nya,ntrainsamps), prob_pop(nxa,nya,ntrainsamps))

DATA input_data_directory /'/Projects/Reforecast2/netcdf'/
DATA output_data_directory /'/Users/thamill/refcst2/test/git'/
   
DATA cmonth /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

! ---- query user for forecast lead day (i1)

CALL getarg(1,cleadb)    ! lead time in hours, start of window
CALL getarg(2,cleade)    ! lead time in hours, end of window
CALL getarg(3,ccompute)  ! Y or N for computing.  If N, read previously computed data from file
PRINT *,'running ppn_rawense_ccpa.x for lead times of ',cleadb,' to ',cleade,' hours.'
IF (ccompute .eq. 'Y' .or. ccompute .eq. 'y') THEN
   PRINT *,'Will generate raw forecasts here as opposed to reading them in from file.'
ELSE
   PRINT *,'Will read in from file previously generated raw forecast probabilities'
ENDIF

! ---- loop thru months and load up all the data for this month, as
!      well as the prior and subsequent 30 days of data.  Loads the 
!      nyears worth of data

DO imonth = 1, 12
   
   CALL cpu_time(t1)
   PRINT *,'----- processing month  = ', imonth, ' cpu time = ',t1
	 
     ! ---- read in the forecast and analyzed data for this month and the surrounding 
     !      two months, over all the nyears worth of data
   
   PRINT *,'calling load_precipfcst_and_ccpa'
   CALL load_precipfcst_and_ccpa (imonth, nxf, nyf, nxa, nya, ntrainsamps, &
     nmembers, cleadb, cleade, input_data_directory, &
  	 rlons_fcst, rlats_fcst, rlons_anal, rlats_anal, &
     apcp_anal, apcp_fcst_ens, pwat_fcst_mean, date_anal, &
	 date_fcst, iccpa_mask)
	 
   ! --- ensemble mean precip. Calc this.

   apcp_fcst = SUM(apcp_fcst_ens, DIM=3) / REAL(nmembers)

   ! ---- the analog procedure will work by finding the dates of past forecasts
   !      that have the closest pattern match to current forecasts.  Then we form an
   !      ensemble from the analyzed precipitation on those dates.  The analyzed 
   !      precip is on a different, higher-resolution grid, so we will need to 
   !      determine, for an analyzed grid point in the conus, which is the closest 
   !      forecast grid point to it, and then use the set of dates for that forecast
   !      grid point.  Call a procedure here that makes a list of all the analyzed
   !      grid points that are closest to a given forecast grid point.

   PRINT *,'calling list_to_nearest_fcstgridpt'
   CALL list_to_nearest_fcstgridpt( nxa, nya, nxf, nyf, nmaxlist, rlons_anal, rlats_anal, &
        rlons_fcst, rlats_fcst, iccpa_mask, iccpa_xlist, iccpa_ylist, iccpa_ktr )

   ! --- calculate precipitation amounts associated with certain quantiles of the analyzed data

   PRINT *,'calling compute_quantiles_ccpa'
   CALL compute_quantiles_ccpa2 (nxa, nya, ntrainsamps, apcp_anal, iccpa_mask, &
        thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, thresh_median, &
        thresh_95)

   ! --- compute climatologies for analyzed precipitation at fixed thresholds and for
   !     thresholds associated with quantiles of the analyzed distribution

   PRINT *,'calling compute_all_climatologies_ccpa2'
   CALL compute_all_climatologies_ccpa2 (nxa, nya, ntrainsamps, &
        apcp_anal, iccpa_mask, thresh_95, thresh_uppertercile, thresh_upperquintile, &
        thresh_upperdecile, thresh_median, climo_pop, climo_1mm, climo_2p5mm, &
        climo_5mm, climo_10mm, climo_25mm, climo_50mm, climo_median, &
        climo_uppertercile, climo_upperquintile, climo_upperdecile, climo_95)

   ! ---- only do the actual analog computations and save them to file if 
   !      ccompute = Y or y (read in on command line).  

   PRINT *,'ccompute = ', ccompute
   IF (ccompute .eq. 'Y' .or. ccompute .eq. 'y') THEN

      ! ---- process all grid points on the forecast grid with the analog technique,
      !      returning probabilities. 
      ! **** NOTE TO DEVELOPERS OF OTHER FORECAST METHODS ****
      !      HERE IS WHERE YOU WOULD INSERT YOUR OWN ALERNATIVE 
      !      POST-PROCESSING ROUTINE

      CALL direct_forecast (nxf, nyf, nxa, nya, ntrainsamps, nmembers, nmaxlist, &
         apcp_fcst_ens, apcp_anal, iccpa_ktr, iccpa_mask, imonth, date_anal, &
         iccpa_xlist, iccpa_ylist, thresh_median, thresh_uppertercile, &
         thresh_upperquintile, thresh_upperdecile, thresh_95, prob_pop, prob_1mm, &
         prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, prob_50mm, prob_median, &
         prob_uppertercile, prob_upperquintile, prob_upperdecile, prob_95)

      ! ---- write this month's post-processed analog forecast data to file
      
      outfile_dates = TRIM(output_data_directory)//'/ppn_rawens_ccpa_dates_'//&
           cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
      outfile_data = TRIM(output_data_directory)//'/ppn_rawens_ccpa_data_'//&
           cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
	  outfile_determ = TRIM(output_data_directory)//'/ensmean_'//&
	       cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
	  outfile_quantiles = TRIM(output_data_directory)//'/quantiles_'//cmonth(imonth)//'_'//&
	 	   cleadb//'_to_'//cleade//'.dat'
      CALL write_ccpa_analogs_tofile_direct2 (outfile_dates, outfile_data, outfile_determ, &
	       outfile_quantiles, nxa, nya, nxf, nyf, ntrainsamps, iccpa_mask, prob_pop, prob_1mm, &
	 	   prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, prob_50mm, prob_median, prob_uppertercile, &
	 	   prob_upperquintile,prob_upperdecile, prob_95, apcp_anal, apcp_fcst, apcp_fcst, &
	       thresh_median, thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, &
	       thresh_95, rlons_anal, rlats_anal, rlons_fcst, rlats_fcst, date_anal, date_fcst)

   ELSE ! ccompute = "N" or "n"
	  
      infile_dates = TRIM(output_data_directory)//'/ppn_rawens_ccpa_dates_'//&
           cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
      infile_data = TRIM(output_data_directory)//'/ppn_rawens_ccpa_data_'//&
           cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
      CALL read_ccpa_analogs_fromfile_direct2 (infile_dates, infile_data, &
           nxa, nya, nxf, nyf, ntrainsamps, iccpa_mask, prob_pop, prob_1mm, prob_2p5mm, &
           prob_5mm, prob_10mm, prob_25mm, prob_50mm, prob_median, prob_uppertercile, &
           prob_upperquintile, prob_upperdecile, prob_95, nxfin, nyfin, apcp_anal, &
           apcp_fcst, rlons_anal, rlats_anal, rlons_fcst, rlats_fcst, date_anal, date_fcst)

   ENDIF ! ccompute

   ! ---- calculate brier skill scores, reliability information, and store
   !      that information to file.

   PRINT *, 'calling calculate_reliabilities_rawens'
   CALL calculate_reliabilities_rawens (cmonth(imonth), cleadb, cleade, &
        output_data_directory, nxa, nya, ntrainsamps, nyears, iccpa_mask, &
		rlons_anal, rlats_anal, apcp_anal, date_anal, imonth, thresh_95, &
        thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, thresh_median, climo_pop, &
        climo_1mm, climo_2p5mm, climo_5mm, climo_10mm, climo_25mm, climo_50mm, &
        climo_median, climo_uppertercile, climo_upperquintile,  climo_upperdecile, &
        climo_95, prob_pop, prob_1mm, prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, &
        prob_50mm, prob_median, prob_uppertercile, prob_upperquintile, prob_upperdecile, prob_95)

END DO ! imonth = 1, nmonths

DEALLOCATE(prob_1mm,prob_2p5mm,prob_5mm,prob_10mm, prob_25mm,prob_50mm,&
   prob_median, prob_uppertercile,prob_upperquintile, prob_upperdecile,&
   prob_95, prob_pop)

END PROGRAM ppn_rawens_ccpa

