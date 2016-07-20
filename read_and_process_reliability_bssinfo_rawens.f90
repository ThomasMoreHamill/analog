! FILE: read_and_process_reliability_bssinfo_rawens.f90
!  to compile this so that python can call the fortran subroutines herein, 
!  the user will need to install the f2py module (here for python version 2.7).  
!  Then to run this, one would type on the command line:  
!
! f2py -c -m read_and_process_reliability_bssinfo_rawens read_and_process_reliability_bssinfo_rawens.f90
!
! from read_and_process_reliability_bssinfo_rawens import read_and_process_reliability_bssinfo_rawens
!
! see documentation of f2py for more.  The python syntax of the call is different
! than the fortran syntax; suppose you have a routine x with a fortran interface like
! SUBROUTINE x(indata,outdata)
!
! The call to this from python would be outdata = x(indata).

SUBROUTINE read_and_process_reliability_bssinfo_rawens(cleadb, cleade, cthresh, &
	output_data_directory, nxa, nya, bssyearly, bssmonthly, bssmap, &
	relia, frequse, rlons_anal, rlats_anal, iccpa_mask)

! this subroutine reads in reliability diagram information related to the raw ensemble.

INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*3, INTENT(IN) :: cleadb, cleade
CHARACTER*(*), INTENT(IN) :: cthresh
CHARACTER*(*), INTENT(IN) :: output_data_directory
REAL, INTENT(OUT) :: bssyearly
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: bssmap
REAL, INTENT(OUT), DIMENSION(12) :: bssmonthly(12)
REAL, INTENT(OUT), DIMENSION(0:20) :: frequse, relia
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlons_anal, rlats_anal
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: iccpa_mask

! below are f2py compiler information .... don't erase.
!f2py intent(in) cleadb, cleade, cthresh, output_data_directory, nxa, nya
!f2py intent(out) bssyearly, bssmonthly, bssmap, relia, frequse
!f2py intent(out) rlons_anal, rlats_anal, iccpa_mask
!f2py depend(12) bssmonthly
!f2py depend(nxa,nya) bssmap, rlons_anal, rlats_anal, iccpa_mask
!f2py depend(0:20) relia, frequse

REAL*8 bsclimo ! accumulating lots of small numbers.  To keep roundoff from being an issue
REAL*8 bs     ! use real*8's not real*4's
REAL bss_mo
REAL*8 bsclimo_mo ! accumulating lots of small numbers.  To keep roundoff from being an issue
REAL*8 bs_mo     ! use real*8's not real*4's
REAL, DIMENSION(nxa,nya)   :: bssmap_mo
REAL*8, DIMENSION(nxa,nya) :: bsmap_mo
REAL*8, DIMENSION(nxa,nya) :: bsmap      ! map of brier scores
REAL*8, DIMENSION(nxa,nya) :: bsclimomap ! map of climatological forecast brier scores
REAL*8, DIMENSION(nxa,nya) :: bsclimomap_mo ! map of climatological forecast brier scores
REAL*8, DIMENSION(0:20,2)  :: contab
REAL*8, DIMENSION(0:20,2)  :: contab_mo
REAL*8                     :: ctot
REAL, DIMENSION(0:20)      :: frequse_mo
REAL, DIMENSION(0:20)      :: relia_mo

CHARACTER*120 infile
CHARACTER*3, DIMENSION(12) :: cmonth
DATA cmonth /'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'/

! ---- initialize

bs = 0.
bsclimo = 0.
contab = 0.
bsmap = 0.
bsclimomap = 0.
bssmonthly = 0.
rlons_anal = -99.99
rlats_anal = -99.99
iccpa_mask = -99

! ---- loop thru months

DO imonth = 1, 12

   infile = TRIM(output_data_directory) // 'reliability_ppn_rawens_'//&
	      cthresh//'_'//cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'

   ! ---- read reliability, BSS information to file
 
   PRINT *,'reading reliability information from ',TRIM(infile)
   OPEN (unit=1, file=infile, status='old', form='unformatted')
   READ (1) bss_mo, bs_mo, bsclimo_mo
   READ (1) frequse_mo
   READ (1) relia_mo
   READ (1) contab_mo
   READ (1) bssmap_mo
   READ (1) bsmap_mo
   READ (1) bsclimomap_mo
   !READ (1) rlons_anal
   !READ (1) rlats_anal
   !READ (1) iccpa_mask
   CLOSE (1)
   PRINT *, 'max, min of bsmap_mo = ',maxval(bsmap_mo), minval(bsmap_mo)
   PRINT *, 'max, min of bsclimomap_mo = ',maxval(bsclimomap_mo), minval(bsclimomap_mo)
   PRINT *, 'bsclimomap_mo(nxa/2,nya/2) = ',bsclimomap_mo(nxa/2,nya/2)

   ! --- tally up statistics over months, etc.

   bs = bs + bs_mo
   bsclimo = bsclimo + bsclimo_mo
   contab = contab + contab_mo
   bsmap = bsmap + bsmap_mo
   bsclimomap = bsclimomap + bsclimomap_mo
   bssmonthly(imonth) = bss_mo

END DO

! ---- with tallied contingency tables, now set reliability and frequency of use

ctot      = SUM(contab)
bssyearly = 1. - bs / bsclimo
relia(:)  = -99.99
DO icat = 0,20
  frequse(icat) = (contab(icat,2)  + contab(icat,1)) / ctot
  IF ((contab(icat,1) + contab(icat,2)) .gt. 0) &
     relia(icat) = contab(icat,2) / (contab(icat,1) + contab(icat,2))
  !PRINT 203,float(icat*5),relia(icat)*100.,frequse(icat)
  !203 format(f5.1,3x,2(f8.3,3x))
END DO  !icat

! --- compute the skill score map as a function of location.

DO i = 1, nxa
  DO j = 1, nya
     IF (bsclimomap(i,j) .GT. 0.0) THEN
        bssmap(i,j) = 1. - bsmap(i,j) / bsclimomap(i,j)
     ELSE
        bssmap(i,j) = -999.99
     ENDIF
  END DO
END DO
!print *,'bssmap(:,nya/2) = ',bssmap(:,nya/2)

RETURN
END SUBROUTINE read_and_process_reliability_bssinfo_rawens