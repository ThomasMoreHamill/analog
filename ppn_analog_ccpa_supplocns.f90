PROGRAM ppn_analog_ccpa_supplocns

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! driver for rank analog technique to produce calibrated probabilisitic total 
! precipitation amount forecasts  on 0.125-degree CCPA grid .  Uses 
! "supplemental locations" to increase sample size. 
!
! written by Tom Hamill, Sep 2013-Mar 2016 (tom.hamill@noaa.gov), 303.497.3060
! and generally follows algorithm from Hamill, T. M., M. Scheuerer,
! and G. T. Bates, 2015: Analog probabilistic precipitation forecasts using GEFS
! Reforecasts and Climatology-Calibrated Precipitation Analyses. 
! Mon. Wea. Rev., 143, 3300-3309.  Also: online appendix A and appendix B.
!
! Calling tree:
! ppn_analog_ccpa_supplocns
!    getarg (multiple times, a system call)
!    cpu_time
!    load_precipfcst_and_ccpa
!       (check.f90 and various netCDF library I/O routines)  
!    read_precip_analog_locations_ccpa5
!    determine_nearest_nint
!    compute_quantiles_ccpa2
!    compute_all_climatologies_ccpa2
!    analog_forecast_supp_v2
!       (calls several routines therein)
!    write_ccpa_analogs_tofile_direct2
!    read_ccpa_analogs_fromfile_direct2
!    calculate_reliabilities_supp
!
!  Notes: (1) in routines write_ccpa_analogs_tofile_direct2.f90 and 
!  read_ccpa_analogs_fromfile_direct2.f90, I believe record length is 
!  measured in bytes for gfortran, words for ifort compiler.  Change there
!  as needed. (2) the computational speed of the algorithm is going
!  to be slower with more supplemental locations, but the algorithm
!  will also produce more accurate forecasts, especially for extreme 
!  large forecast events.
!
! Usage: Compile with > make -f make_ppn_analog_ccpa_supplocns
! > ppn_analog_ccpa_supplocns.x cleadb cleade csupp ccompute
!    where cleadb is the beginning lead time in hours
!    cleade is the ending lead time in hours
!    csupp is the number of supplemental locations (01 to 20)
!       (with supplemental locations previously determined via 
!        compute_precip_analog_locations_ccpa5.x)
!    ccompute = Y if computations to be performed, N if previously
!       computed forecasts are to be read from file and validated
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PARAMETER (nyears      = 14)  ! 2002-2015 here for the CCPA data
PARAMETER (ndaysub     = 91)  ! number of days of training data for each yr, 15th of month +/- 45 days
PARAMETER (nxa         = 464) ! CCPA grid dimensions, surrounding CONUS
PARAMETER (nya         = 224) !
PARAMETER (nxf         = 128) ! Grid dimensions of reforecast data, here on ~ 1/2 deg Gaussian grid 
PARAMETER (nyf         = 61)  !    surrounding CONUS
PARAMETER (ntrainsamps = ndaysub*nyears) ! number of training data samples loaded up
PARAMETER (nanalogs    = 300) ! max number of analogs used
PARAMETER (nmembers    = 11)  ! number of ensemble members in reforecast
PARAMETER (nmaxlist    = 16)  ! presumed maximum number of 1/8-deg CCPA points that 
	                          ! are closest to a particular 1/2-deg forecast pt
PARAMETER (nsuppmax    = 20)  ! maximum number of supplemental locations

INTEGER, PARAMETER :: window_size = 9 ! (+/- 4 grid pts)  for Savitzky-Golay smoother
INTEGER, PARAMETER :: order = 3 ! for Savitzky-Golay smoother (see Numerical Recipes)

INTEGER, DIMENSION(ntrainsamps)      :: date_anal ! in yyyymmddhh format, initialization time
INTEGER, DIMENSION(ntrainsamps)      :: date_fcst ! in yyyymmddhh format, end of fcst period
INTEGER*2, DIMENSION(nxa,nya)        :: iccpa_mask   ! a mask for US.
INTEGER, DIMENSION(nxa,nya)          :: inearest ! i location on fcst grid nearest (below) this analysis grid location
INTEGER, DIMENSION(nxa,nya)          :: jnearest ! j location on fcst grid nearest (below) this analysis grid location

INTEGER, DIMENSION(nxa,nya,nsuppmax) :: xlocation  ! these will contain the grid indices of the supplemental
INTEGER, DIMENSION(nxa,nya,nsuppmax) :: ylocation  ! training locations for each CCPA grid point

REAL, ALLOCATABLE, DIMENSION(:,:,:)  :: apcp_anal ! precipitation analysis on 0.125-deg grid
REAL, DIMENSION(nxf, nyf, nmembers, ntrainsamps) :: apcp_fcst_ens ! precip fcst ensemble
REAL, DIMENSION(nxf, nyf, ntrainsamps) :: apcp_fcst ! precip fcst mean
REAL, DIMENSION(nxa,nya)                  :: climo_pop  ! climatological probabilities
REAL, DIMENSION(nxa,nya)                  :: climo_1mm
REAL, DIMENSION(nxa,nya)                  :: climo_2p5mm
REAL, DIMENSION(nxa,nya)                  :: climo_5mm
REAL, DIMENSION(nxa,nya)                  :: climo_10mm
REAL, DIMENSION(nxa,nya)                  :: climo_25mm
REAL, DIMENSION(nxa,nya)                  :: climo_50mm
REAL, DIMENSION(nxa,nya)                  :: climo_median
REAL, DIMENSION(nxa,nya)                  :: climo_uppertercile ! set to 1/3
REAL, DIMENSION(nxa,nya)                  :: climo_upperquintile ! set to 1/5
REAL, DIMENSION(nxa,nya)                  :: climo_upperdecile ! set to 1/10
REAL, DIMENSION(nxa,nya)                  :: climo_95 ! set to 1/20th

REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: fmean_analog
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_1mm ! output arrays of post-processed, here for P(obs>1mm)
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_2p5mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_5mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_10mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_25mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_50mm
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_median 
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_uppertercile
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_upperquintile
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_upperdecile
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_95
REAL, ALLOCATABLE, DIMENSION(:,:,:)       :: prob_pop
REAL, DIMENSION(nxf,nyf,ntrainsamps)      :: pwat_fcst_mean ! total-col precip. water mean fcst
REAL, DIMENSION(nxf,nyf)                  :: rlons_fcst
REAL, DIMENSION(nxf,nyf)                  :: rlats_fcst ! forecast lat/lon grids
REAL, DIMENSION(nxa,nya)                  :: rlons_anal 
REAL, DIMENSION(nxa,nya)                  :: rlats_anal ! 0.125-degree CCPA analysis grid
REAL, DIMENSION(nxa,nya)                  :: thresh_95 ! precipitation amount threshold assoc'd with 95th %ile
REAL, DIMENSION(nxa,nya)                  :: thresh_uppertercile
REAL, DIMENSION(nxa,nya)                  :: thresh_upperquintile
REAL, DIMENSION(nxa,nya)                  :: thresh_upperdecile 
REAL, DIMENSION(nxa,nya)                  :: thresh_median

CHARACTER*120 outfile_dates, outfile_data, outfile_determ, outfile_quantiles, infile_dates, infile_data
CHARACTER*80 input_data_directory, output_data_directory
CHARACTER*3 cleadb, cleade
CHARACTER*3, DIMENSION(12) :: cmonth
CHARACTER*1 ccompute
CHARACTER*2 csupp

DATA cmonth /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/
DATA input_data_directory /'/Projects/Reforecast2/netcdf'/
DATA output_data_directory /'/Users/thamill/refcst2/test/git'/

! --- use dynamic memory allocation for large arrays

ALLOCATE (apcp_anal(nxa,nya,ntrainsamps), prob_1mm(nxa,nya,ntrainsamps), &
   prob_2p5mm(nxa,nya,ntrainsamps), prob_5mm(nxa,nya,ntrainsamps), &
   prob_10mm(nxa,nya,ntrainsamps), prob_25mm(nxa,nya,ntrainsamps), &
   prob_50mm(nxa,nya,ntrainsamps), prob_median(nxa,nya,ntrainsamps), &
   prob_uppertercile(nxa,nya,ntrainsamps), prob_upperquintile(nxa,nya,ntrainsamps), &
   prob_upperdecile(nxa,nya,ntrainsamps), prob_95(nxa,nya,ntrainsamps), &
   prob_pop(nxa,nya,ntrainsamps), fmean_analog(nxa,nya,ntrainsamps))

! ---- query user for forecast lead day (i1)

CALL getarg (1,cleadb)    ! lead time in hours, start of window
CALL getarg (2,cleade)    ! lead time in hours, end of window
CALL getarg (3,csupp)     ! number of supplemental locations to use, 01 to 20
READ (csupp,'(i2)') nsupp
READ (cleade, '(i3)') ileade
CALL getarg (4,ccompute)  ! Y or N for computing.  If N, read previously computed data from file
PRINT *,'running ppn_analog_rank_ccpa.x for lead times of ',cleadb,' to ',cleade,' hours.'
IF (ccompute .eq. 'Y' .or. ccompute .eq. 'y') THEN
   PRINT *,'Will generate analog forecasts here as opposed to reading them in from file.'
ELSE
   PRINT *,'Will read in from file previously generated analog forecast probabilities'
ENDIF

! ---- loop thru months and process each 

DO imonth = 1, 12
   
   CALL cpu_time(t1)
   PRINT *,'----- processing month  = ', imonth, ' cpu time = ',t1

   ! ---- read in the forecast and analyzed data for this month and the surrounding 
   !      two months, over all the nyears worth of data
   
   PRINT *,'calling load_precipfcst_and_ccpa'
   CALL load_precipfcst_and_ccpa (imonth, nxf, nyf, nxa, nya, ntrainsamps, &
     nmembers, cleadb, cleade, input_data_directory, &
	 rlons_fcst, rlats_fcst, rlons_anal, rlats_anal, &
     apcp_anal, apcp_fcst_ens, pwat_fcst_mean, date_anal, date_fcst, iccpa_mask)

   ! --- this analog approach here will be based on the ensemble mean precip. Calc this.

   apcp_fcst = SUM(apcp_fcst_ens, DIM=3) / REAL(nmembers)

   !  --- read in the supplemental data locations for each 1/8-degree output point.  See 
   !  www.esrl.noaa.gov/psd/people/tom.hamill/Analog-CCPA-MWRexpedited-Hamill-AppA-v2.pdf

   PRINT *, 'calling read_precip_analog_locations_ccpa5'
   CALL read_precip_analog_locations_ccpa5 (imonth, nxa, nya, &
       nsuppmax, input_data_directory, xlocation, ylocation)

   ! ---- for each analysis (i,j) location, identify the (i,j) index of the 
   !      forecast grid point just below and to the left of the analysis location.  

   PRINT *, 'calling determine_nearest_nint'
   CALL determine_nearest_nint(nxa, nya, nxf, nyf, rlons_anal, rlats_anal, &
      rlons_fcst, rlats_fcst, inearest, jnearest)

   ! --- calculate precipitation amounts associated with certain quantiles of the analyzed data

   PRINT *,'calling compute_quantiles_ccpa2'
   CALL compute_quantiles_ccpa2 (nxa, nya, ntrainsamps, apcp_anal, iccpa_mask, &
        thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, thresh_median, &
        thresh_95)

   ! --- compute climatologies for analyzed precipitation at fixed thresholds and for
   !     thresholds associated with quantiles of the analyzed distribution

   PRINT *,'calling compute_all_climatologies_ccpa2'
   CALL  compute_all_climatologies_ccpa2 (nxa, nya, ntrainsamps, &
        apcp_anal, iccpa_mask, thresh_95, thresh_uppertercile, thresh_upperquintile, &
        thresh_upperdecile, thresh_median, climo_pop, climo_1mm, climo_2p5mm, &
        climo_5mm, climo_10mm, climo_25mm, climo_50mm, climo_median, &
        climo_uppertercile, climo_upperquintile, climo_upperdecile, climo_95)

   ! ---- only do the actual analog computations and save them to file if 
   !      ccompute = Y or y (read in on command line).  

   PRINT *,'ccompute = ', ccompute
   IF (ccompute .eq. 'Y' .or. ccompute .eq. 'y') THEN

      ! ---- process all grid points on the forecast grid with the analog technique,
      !      returning probabilities.   This version potentially uses data from
      !      other supplemental training locations.
      ! **** NOTE TO DEVELOPERS OF OTHER FORECAST METHODS ****
      !      HERE IS WHERE YOU WOULD INSERT YOUR OWN ALERNATIVE 
      !      POST-PROCESSING ROUTINE

      PRINT *, 'calling analog_forecast_supp_v2'
      CALL analog_forecast_supp_v2 (nxf, nyf, nxa, nya, ntrainsamps, nanalogs, &
         nsupp, nsuppmax, imonth, window_size, order, ileade, apcp_fcst, pwat_fcst_mean, apcp_anal, &
         iccpa_mask, date_anal, inearest, jnearest, rlons_fcst, rlats_fcst, &
         thresh_median, thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, &
         thresh_95, xlocation, ylocation, rlons_anal, rlats_anal, input_data_directory, &
		 prob_pop, prob_1mm, prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, &
		 prob_50mm, prob_median, prob_uppertercile, prob_upperquintile, &
         prob_upperdecile, prob_95, fmean_analog)

      ! ---- write this month's post-processed analog forecast data to file
      
      outfile_dates = TRIM(output_data_directory)//'/ppn_pwat_ccpa_supp'//csupp//'_analog_dates_'//&
           cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
      outfile_data = TRIM(output_data_directory)//'/ppn_pwat_ccpa_supp'//csupp//'_analog_data_'//&
           cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
      outfile_determ = TRIM(output_data_directory)//'/ppn_pwat_ccpa_supp'//csupp//'_analog_determ_'//&
           cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
	  outfile_quantiles = TRIM(output_data_directory)//'/quantiles_'//cmonth(imonth)//'_'//&
	  	   cleadb//'_to_'//cleade//'.dat'
      CALL write_ccpa_analogs_tofile_direct2 (outfile_dates, outfile_data, outfile_determ, &
           outfile_quantiles, nxa, nya, nxf, nyf, ntrainsamps, iccpa_mask, prob_pop, prob_1mm, &
		   prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, prob_50mm, prob_median, prob_uppertercile, &
		   prob_upperquintile,prob_upperdecile, prob_95, apcp_anal, apcp_fcst, fmean_analog, &
           thresh_median, thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, &
           thresh_95, rlons_anal, rlats_anal, rlons_fcst, rlats_fcst, date_anal, date_fcst)

   ELSE ! ccompute = "N" or "n" - then read in previously computed forecast probabilities

      infile_dates = TRIM(output_data_directory)//'/ppn_pwat_ccpa_supp'//csupp//'_analog_dates_'//&
           cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
      infile_data = TRIM(output_data_directory)//'/ppn_pwat_ccpa_supp'//csupp//'_analog_data_'//&
           cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'
      CALL read_ccpa_analogs_fromfile_direct2 (infile_dates, infile_data, &
           nxa, nya, nxf, nyf, ntrainsamps, iccpa_mask, prob_pop, prob_1mm, prob_2p5mm, &
           prob_5mm, prob_10mm, prob_25mm, prob_50mm, prob_median, prob_uppertercile, &
           prob_upperquintile, prob_upperdecile, prob_95, nxfin, nyfin, apcp_anal, &
           apcp_fcst, rlons_anal, rlats_anal, rlons_fcst, rlats_fcst, date_anal, date_fcst)

   ENDIF ! ccompute

   ! ---- calculate brier skill scores, reliability information, and store
   !      that information to file.

   PRINT *, 'calling calculate_reliabilities_supp, cleadb, cleade = ',cleadb, cleade
   CALL calculate_reliabilities_supp (cmonth(imonth), cleadb, cleade, csupp, &
        nxa, nya, ntrainsamps, nyears, date_anal, imonth, output_data_directory, &
        iccpa_mask, rlons_anal, rlats_anal, apcp_anal, thresh_95, &
        thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, thresh_median, &
        climo_pop, climo_1mm, climo_2p5mm, climo_5mm, climo_10mm, climo_25mm, climo_50mm, &
        climo_median, climo_uppertercile, climo_upperquintile, climo_upperdecile, &
        climo_95, prob_pop, prob_1mm, prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, &
        prob_50mm, prob_median, prob_uppertercile, prob_upperquintile, prob_upperdecile, &
        prob_95)

END DO ! imonth = 1, nmonths

DEALLOCATE (apcp_anal, prob_1mm, prob_2p5mm, prob_5mm, prob_10mm, prob_25mm,prob_50mm, prob_median, &
   prob_uppertercile, prob_upperquintile, prob_upperdecile, prob_95, prob_pop, fmean_analog)

END PROGRAM ppn_analog_ccpa_supplocns

