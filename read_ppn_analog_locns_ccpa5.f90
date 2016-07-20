SUBROUTINE read_ppn_analog_locns_ccpa5(cmonth,cleadb, cleade, nxa, nya, &
    nanalogs, xlocationa, ylocationa, conusmaska, rlona, rlata)

! compile with f2py -c -m read_ppn_analog_locns_ccpa5 read_ppn_analog_locns_ccpa5.f90

CHARACTER*3, INTENT(IN) :: cmonth, cleadb, cleade

INTEGER, INTENT(IN) :: nxa, nya, nanalogs

INTEGER, INTENT(OUT), DIMENSION(nxa,nya,nanalogs) :: xlocationa, ylocationa
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: conusmaska
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlona, rlata

! f2py intent(in) cmonth, cleadb, cleade, nya, nxa, nanalogs
! f2py intent(out) xlocationa, ylocationa, conusmaska, rlona, rlata
! f2py depend(nxa,nya,nanalogs) xlocationa, ylocationa
! f2py depend(nxa,nya) conusmaska, rlona, rlata

CHARACTER*120 outfile

outfile = '/Projects/Reforecast2/netcdf/ppn_analog_locns_ccpa5_'//cmonth//'_'// &
   cleadb//'_to_'//cleade//'.dat'
PRINT *, 'reading from ', TRIM(outfile)
OPEN (UNIT=1, FILE=outfile, STATUS='old', FORM='unformatted')
READ (1) nxfin, nyfin, nxain, nyain, nanalogsin
print *,nxain,nyain
print *,nxa,nya
READ (1) xlocationa
READ (1) ylocationa
READ (1) conusmaska
READ (1) rlona
READ (1) rlata
CLOSE (1)

RETURN
END SUBROUTINE read_ppn_analog_locns_ccpa5
