SUBROUTINE load_precipquantiles_ccpa_analonly (imonth, nxa, nya, npct, &
     cleadb, cleade, input_data_directory, conusmask, lonsa, latsa, &
	 pctv, panal_quantiles, terrain)
!
! purpose: read in the data set for this month of the precipitation quantiles, 
!    both forecast and analyzed, as well as the rank correlation between 
!    forecast and analysis.
! 
USE netcdf

INTEGER, INTENT(IN) :: imonth, nxa, nya, npct
CHARACTER*3, INTENT(IN) :: cleadb, cleade
CHARACTER*(*), INTENT(IN) :: input_data_directory
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(OUT), DIMENSION(nxa,nya,npct) :: panal_quantiles
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: lonsa, latsa
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: terrain
REAL, INTENT(OUT), DIMENSION(npct) :: pctv

CHARACTER*80 infile
CHARACTER*20 cfield
CHARACTER*3, DIMENSION(12) :: cmonths

DATA cmonths /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

! ---- open the file

infile = TRIM(input_data_directory)//'/refcstv2_apcp_CCPAgrid_CDF_fhour' // &
  cleadb//'_to_'//cleade//'_'//cmonths(imonth)//'.nc'
PRINT *,'netid, reading from ',netid, TRIM(infile)
CALL check (nf90_open(infile,NF90_NOWRITE,netid))

! ---- read in the conus mask (1 if inside conus) and lats and lons 

cfield ='conusmask'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,conusmask,&
     start=(/1,1/),count=(/nxa,nya/)))

!print *,'conusmask(nxa/2,:) = ',conusmask(nxa/2,:)
!print *,'conusmask(:,nya/2) = ',conusmask(:,nya/2)

cfield ='lonsa'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,lonsa,&
  start=(/1,1/),count=(/nxa,nya/)))

cfield ='latsa'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,latsa,&
  start=(/1,1/),count=(/nxa,nya/)))

cfield = 'pct'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,pctv,start=(/1/),count=(/npct/)))

cfield = 'panal_quantiles'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,panal_quantiles,start=(/1,1,1/),count=(/nxa,nya,npct/)))

! --- close file.

CALL check(nf90_close(netid))

! ---- open the terrain file 

terrain(:,:) = 0.0
infile = TRIM(input_data_directory)//'/elevation_ccpa.nc'
PRINT *,'netid, reading from ',netid, TRIM(infile)
CALL check (nf90_open(infile,NF90_NOWRITE,netid))
cfield ='terrain_height'
CALL check (nf90_inq_varid(netid,trim(adjustl(cfield)),ivar))
CALL check (nf90_get_var(netid,ivar,terrain,&
  start=(/1,1/),count=(/nxa,nya/)))

CALL check(nf90_close(netid))

RETURN
END subroutine load_precipquantiles_ccpa_analonly
