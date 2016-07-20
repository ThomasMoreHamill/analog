SUBROUTINE load_precipfcst_and_ccpa (imonth, nxf, nyf, nxa, nya, ntrainsamps, &
	nmembers, cleadb, cleade, input_data_directory, &
	rlons_fcst, rlats_fcst, rlons_anal, rlats_anal, &
	apcp_anal, apcp_fcst_ens, pwat_fcst_mean, date_anal, date_fcst, iccpa_mask)
!
! Purpose: for the month of interest, read in the multi-yearly reforecast and 
! associated CCPA precipitation data.  Reforecast data is on the native ~0.5 degree
! grid and the CCPA data is on its 0.125 degree grid.  Read in data
! not only the month of interest but the two surrounding months, e.g., read
! Dec-Jan-Feb data when Jan is the month of interest.
!
! called by ppn_analog_ccpa_supplocns
! 
USE netcdf

INTEGER, INTENT(IN) :: imonth    ! month of interest to read in
INTEGER, INTENT(IN) :: nxf, nyf  ! forecast grid x and y dimensions
INTEGER, INTENT(IN) :: nxa, nya  ! precipitation analysis grid x and y dimensions
INTEGER, INTENT(IN) :: ntrainsamps   ! expected max number of training samples (91*nears)
   ! since we'll read in a given month and the surrounding 2 months of data.
INTEGER, INTENT(IN) :: nmembers  ! number of ensemble members
CHARACTER*3, INTENT(IN) :: cleadb, cleade  ! beginning and ending lead times in hours, e.g., '012'
CHARACTER*(*), INTENT(IN) :: input_data_directory ! location of input data

REAL, INTENT(OUT), DIMENSION(nxf, nyf) :: rlons_fcst, rlats_fcst ! forecast lat/lon grids
REAL, INTENT(OUT), DIMENSION(nxa, nya) :: rlons_anal, rlats_anal ! 0.125-degree CCPA analysis grid
REAL, INTENT(OUT), DIMENSION(nxa, nya, ntrainsamps) :: apcp_anal ! precipitation analysis on 0.125-deg grid
REAL, INTENT(OUT), DIMENSION(nxf, nyf, nmembers, ntrainsamps) :: apcp_fcst_ens ! precip fcst ensemble
REAL, INTENT(OUT), DIMENSION(nxf, nyf, ntrainsamps) :: pwat_fcst_mean ! total-col precip. water mean fcst
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: iccpa_mask

INTEGER, DIMENSION(ntrainsamps), INTENT(OUT) :: date_anal ! in yyyymmddhh format, initialization time
INTEGER, DIMENSION(ntrainsamps), INTENT(OUT) :: date_fcst ! in yyyymmddhh format, end of fcst period

REAL, DIMENSION (nxa,nya) :: pverif
REAL, DIMENSION (nxf,nyf) :: pfcst
REAL, DIMENSION (nxf,nyf,nmembers) :: pfcst_ens

CHARACTER*120 infile
CHARACTER*20 cfield, dimname

INTEGER, ALLOCATABLE, DIMENSION(:) :: iyyyymmddhh_init 
INTEGER, ALLOCATABLE, DIMENSION(:) :: iyyyymmddhh_fcst
INTEGER, DIMENSION(366) :: idxuse
INTEGER, DIMENSION(12)  :: imid ! julian day of middle of month
INTEGER, DIMENSION(12)  :: ibegin_noleap ! julian day of begin of month, non-leap year
INTEGER, DIMENSION(12)  :: ibegin_leap   ! julian day of middle of month
INTEGER ntimes

DATA imid /15,46,74,105,135,166,196,227,258,288,319,349/ ! middle julian day of each month
DATA ibegin_noleap /1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/
DATA ibegin_leap   /1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336/

! ---- initialize data to missing

apcp_anal      = -99.99
apcp_fcst_ens  = -99.99
pwat_fcst_mean = -99.99

! ---- set an index for which days of the year to include in output sample for a given month

idxuse(:) = 0
IF (imonth .eq. 1) THEN
   idxuse(1:60) = 1
   idxuse(335:365) = 1
ELSE IF (imonth .gt. 1 .and. imonth .lt. 12) THEN
   idxuse(imid(imonth)-45:imid(imonth)+45) = 1
ELSE
   idxuse(305:365) = 1
   idxuse(1:30) = 1
ENDIF   
date_anal(:) = -99  ! initialize to missing
date_fcst(:) = -99 

! ---- open the netcdf file

infile = TRIM(input_data_directory)//'/refcstv2_precip_ccpav3_' // cleadb // &
  '_to_' // cleade // '.nc'
PRINT *,'netid, reading from ',netid, TRIM(infile)
CALL check (nf90_open(infile,NF90_NOWRITE,netid))

! ---- read in the latitude and longitude arrays, forecast and analysis, as
!      well as the conus mask for the ccpa analysis grid.

CALL check (nf90_open(infile,NF90_NOWRITE,netid))
cfield = 'conusmask'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,iccpa_mask,start=(/1,1/),count=(/nxa,nya/)))

cfield = 'lats_anal'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,rlats_anal,start=(/1,1/),count=(/nxa,nya/)))

cfield = 'lons_anal'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,rlons_anal,start=(/1,1/),count=(/nxa,nya/)))

cfield = 'lats_fcst'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,rlats_fcst,start=(/1,1/),count=(/nxf,nyf/)))

cfield = 'lons_fcst'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,rlons_fcst,start=(/1,1/),count=(/nxf,nyf/)))
   
! ---- load in the list of initial and forecast times

cfield = 'time'
CALL check (nf90_inq_dimid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_Inquire_Dimension(netid,ivar,dimname,ntimes))

ALLOCATE (iyyyymmddhh_init(ntimes))
ALLOCATE (iyyyymmddhh_fcst(ntimes))

! ---- read in the initial time of the forecast

cfield = 'yyyymmddhh_init'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,iyyyymmddhh_init,start=(/1/),count=(/ntimes/)))

! ---- read in the time for the end of the forecast period

cfield = 'yyyymmddhh_fcste'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,iyyyymmddhh_fcst,start=(/1/),count=(/ntimes/)))

! ---- for a given year and (input) month, loop thru the days of the year, and if
!      this day is flagged for being appropriate to read in as training data, do so

ktrday = 1
DO itime = 1, ntimes
   iyear   = iyyyymmddhh_init(itime) / 1000000
   immddhh = iyyyymmddhh_init(itime) - iyear*1000000
   imo     = immddhh / 10000
   iddhh   = immddhh - imo*10000
   iday    = iddhh / 100
   IF (mod(iyear,4) .eq. 0) THEN
      ijulday = ibegin_leap(imo) + iday - 1
   ELSE
      ijulday = ibegin_noleap(imo) + iday -1
   ENDIF
   IF (idxuse(ijulday) .eq. 1) THEN
 
      IF (iddhh .eq. 100 .and. imo .eq. imonth) PRINT *, 'Reading in ',iyear,' reforecast and CCPA precip data.'

      ! --- read in precipitation analysis information

      pverif = -99.99
      cfield = 'apcp_anal'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pverif,&
         start=(/1,1,itime/),count=(/nxa,nya,1/)))
      apcp_anal(:,:,ktrday) = pverif(:,:)

      ! ---- read in precip ensemble forecast

      pfcst = -99.99
      cfield ='apcp_fcst_ens'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid,ivar,pfcst_ens,&
         start = (/1,1,1,itime/), count = (/nxf,nyf,nmembers,1/)))
      apcp_fcst_ens(:,:,:,ktrday) = pfcst_ens(:,:,:)

      ! ---- read in precipitable water mean forecast

      pfcst = -99.99
      cfield ='pwat_fcst_mean'
      CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
      CALL check (nf90_get_var(netid, ivar,pfcst,&
          start = (/1,1,itime/), count = (/nxf,nyf,1/)))
      pwat_fcst_mean(:,:,ktrday) = pfcst(:,:)

      ! ---- set the initialization and forecast date for this sample
      
      date_anal(ktrday) = iyyyymmddhh_init(itime)
      date_fcst(ktrday) = iyyyymmddhh_fcst(itime)

      ! ---- increment counter

      ktrday = ktrday+1

   ENDIF ! (idxuse(iday) .eq. 1) 

END DO  ! itime = 1, ntimes
ktrday = ktrday-1

! --- close netcdf file.

CALL check(nf90_close(netid))

DEALLOCATE (iyyyymmddhh_init)
DEALLOCATE (iyyyymmddhh_fcst)

! --- flag points with excessively high analyzed data as bad.

isum = 0
DO jy = nya, 1, -1
   DO ix = 1, nxa
      DO isamp = 1,ktrday
         IF (apcp_anal(ix,jy,isamp) .gt. 1000.) apcp_anal(ix,jy,isamp) = -99.99
      END DO
   END DO
END DO
PRINT *,'returning from load_precipfcst_and_ccpa'

RETURN
END subroutine load_precipfcst_and_ccpa
