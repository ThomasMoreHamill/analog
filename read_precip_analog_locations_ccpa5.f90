SUBROUTINE read_precip_analog_locations_ccpa5 (imonth, nxa, nya, &
   nsupp, input_data_directory, xlocation, ylocation)

! this routine reads in the data containing the supplemental locations for each 
! 1/8-degree analysis grid point.  The analog post-processing procedure is very
! sensitive for extreme forecasts to a lack of data, so the supplemental locations
! are used to boost the training sample size.  Effectively, the analog procedure
! finds past forecasts that match, either at the particular grid point in question,
! or at any of the supplemental locations.  The analogs are formed from the 
! analyzed data at the same time as the forecasts.

INTEGER, INTENT(IN) :: imonth ! 1-12
INTEGER, INTENT(IN) :: nxa, nya, nsupp
CHARACTER*(*), INTENT(IN) :: input_data_directory
INTEGER, DIMENSION(nxa,nya,nsupp), INTENT(OUT) :: xlocation
INTEGER, DIMENSION(nxa,nya,nsupp), INTENT(OUT) :: ylocation

CHARACTER*120 infile
CHARACTER*3, DIMENSION(12) :: cmonth

DATA cmonth /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

! ---- simplifying assumption here: all forecast leads will use the supplemental locations that
!      were used for the 024 to 048-h period.  The assumption is that there is no reason to
!      expect that a short-lead or a long-lead forecast would have different errors and different
!      relationships to terrain characteristics such as to require supplemental locations to 
!      be different for different lead times.  This could be an unwarranted assumption, we admit.

infile = TRIM(input_data_directory) // '/ppn_analog_locns_ccpa5_'// &
   cmonth(imonth)//'_024_to_048.dat'
PRINT *, 'reading from ', TRIM(infile)
OPEN (UNIT=1, FILE=infile, STATUS='old', FORM='unformatted')
READ (1) nxfin, nyfin, nxain, nyain, nanalogsin
IF (nxain .ne. nxa .and. nyain .ne. nya .and. nanalogsin .ne. nsupp) THEN
   PRINT *, 'Error in subroutine read_precip_analog_locations_ccpa3'
   PRINT *, '   Dimensions do not match'
   PRINT *, '   nxa, nya, nsupp = ',nxa, nya, nsupp
   PRINT *, '   nxain, nyain, nanalogsin = ',nxain, nyain, nanalogsin
   PRINT *, '   Stopping.'
   STOP
ENDIF
READ (1) xlocation
READ (1) ylocation
CLOSE (1)


RETURN
END SUBROUTINE read_precip_analog_locations_ccpa5
