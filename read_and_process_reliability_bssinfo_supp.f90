!
! routine called from python script bss_fmonth_ccpa_compare2.py
!
! compiled with 
! f2py -c -m read_and_process_reliability_bssinfo_supp read_and_process_reliability_bssinfo_supp.f90
!
! which assumes you have f2py on your machine.
!
SUBROUTINE read_and_process_reliability_bssinfo_supp(cleadb, cleade, cthresh, &
	csupp, output_data_directory, nxa, nya, bssyearly, bssmonthly, bssmap, &
	relia, frequse, rlons_anal, rlats_anal, iccpa_mask)

INTEGER, INTENT(IN) :: nxa, nya
CHARACTER*3, INTENT(IN) :: cleadb, cleade
CHARACTER*(*), INTENT(IN) :: cthresh, output_data_directory
CHARACTER*2, INTENT(IN) :: csupp
REAL, INTENT(OUT) :: bssyearly
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: bssmap
REAL, INTENT(OUT), DIMENSION(12) :: bssmonthly(12)
REAL, INTENT(OUT), DIMENSION(0:20) :: frequse, relia
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: rlons_anal, rlats_anal
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: iccpa_mask

! below are f2py compiler information .... don't erase.
!f2py intent(in) cleadb, cleade, cthresh, csupp, output_data_directory, nxa, nya
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

   infile = TRIM(output_data_directory)//'/reliability_ppn_pwat_ccpa_supp'//csupp//&
      '_analog_'//cthresh//'_'//cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'

   !infile = '/Rf2_tests/ccpa/reliability_ppn_pwat_ccpa_rankanalog_'//cthresh//'_'//&
   !   cmonth(imonth)//'_'//cleadb//'_to_'//cleade//'.dat'

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
END SUBROUTINE read_and_process_reliability_bssinfo_supp
