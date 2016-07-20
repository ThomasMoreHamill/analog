! FILE: read_analog_probabilities.f90
!  to compile this so that python can call the fortran subroutines herein, 
!  the user will need to install the f2py module (here for python version 2.7).  
!  Then to run this, one would type on the command line:  
!
! f2py -c -m read_analog_probabilities read_analog_probabilities.f90
!
! in python one then has at the beginning of their
! code:
!
! from read_analog_probabilities import read_analog_probabilities
!
! see documentation of f2py for more.  The python syntax of the call is different
! than the fortran syntax; suppose you have a routine x with a fortran interface like
! SUBROUTINE x(indata,outdata)
!
! The call to this from python would be outdata = x(indata).

! ================================================================================================

SUBROUTINE read_analog_probabilities(infile_dates, infile_data, date_anal_in, nxa, nya, nxf, nyf, &
   cthresh, probfcst, iccpa_mask)

CHARACTER*(*), INTENT(IN) :: infile_dates, infile_data, cthresh
INTEGER, INTENT(IN) :: date_anal_in, nxa, nya, nxf, nyf
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: probfcst
INTEGER*2, INTENT(OUT), DIMENSION(nxa,nya) :: iccpa_mask

!f2py intent(in) infile_dates, infile_data, chtresh, date_anal_in, nxa, nya
!f2py intent(out) probfcst
!f2py intent(out) iccpa_mask
!f2py depend(nxa,nya) probfcst
!f2py depend(nxa,nya) iccpa_mask

REAL, DIMENSION(nxa,nya) :: rlons_anal, rlats_anal
REAL, DIMENSION(nxa,nya) :: apcp_anal
REAL, DIMENSION(nxa,nya) :: prob_pop
REAL, DIMENSION(nxa,nya) :: prob_1mm
REAL, DIMENSION(nxa,nya) :: prob_2p5mm
REAL, DIMENSION(nxa,nya) :: prob_5mm
REAL, DIMENSION(nxa,nya) :: prob_10mm
REAL, DIMENSION(nxa,nya) :: prob_25mm
REAL, DIMENSION(nxa,nya) :: prob_50mm
REAL, DIMENSION(nxa,nya) :: prob_median
REAL, DIMENSION(nxa,nya) :: prob_uppertercile
REAL, DIMENSION(nxa,nya) :: prob_upperdecile
REAL, DIMENSION(nxa,nya) :: prob_upperquintile
REAL, DIMENSION(nxa,nya) :: prob_95
REAL, DIMENSION(nxf,nyf) :: apcp_fcst
REAL, DIMENSION(nxf,nyf) :: rlons_fcst
REAL, DIMENSION(nxf,nyf) :: rlats_fcst

INTEGER, DIMENSION(:), ALLOCATABLE :: date_anal, date_fcst

! ----- read the array of dates associated with each of the ntrainsamps samples

PRINT *,'reading from ',TRIM(infile_dates)
OPEN (unit=1, file=infile_dates, status='old',form='unformatted')
READ (1) ntrainsamps
ALLOCATE (date_anal(ntrainsamps), date_fcst(ntrainsamps))
READ (1) date_anal  ! initial time of forecast
READ (1) date_fcst ! verif time of forecast
CLOSE (1)

! ---- find the index of the date that matches the input date

DO idate = 1, ntrainsamps
   IF (date_anal(idate) .eq. date_anal_in) GOTO 1234
ENDDO
PRINT *, 'Did not find the date ',date_anal_in, ' in the list of valid dates.  Stopping'
STOP
1234 CONTINUE

DEALLOCATE (date_anal, date_fcst)

! ---- read this particular record...

irecl = 4*(3 + (nxa*nya)/2 + 2*nxa*nya + 13*nxa*nya + 2 + 3*nxf*nyf + 2)
PRINT *,'reading from ',TRIM(infile_data)
OPEN (unit=1, file=infile_data, status='old',access='direct',recl=irecl)
READ (1,rec=idate) nxain, nyain, ntrainsampsin, iccpa_mask,rlons_anal, &
	rlats_anal, apcp_anal, prob_pop, prob_1mm, prob_2p5mm, &
    prob_5mm, prob_10mm, prob_25mm, prob_50mm, &
    prob_median, prob_uppertercile, prob_upperquintile, &
	prob_upperdecile, prob_95, nxfin, nyfin, &
    apcp_fcst, rlons_fcst, rlats_fcst, idate_analin, &
    idate_fcstin
CLOSE (1)

probfcst = -99.99
IF (trim(cthresh) .eq. 'POP') probfcst(:,:) = prob_pop(:,:)
IF (trim(cthresh) .eq. '1mm') probfcst(:,:) = prob_1mm(:,:)
IF (trim(cthresh) .eq. '2p5mm') probfcst(:,:) = prob_2p5mm(:,:)
IF (trim(cthresh) .eq. '5mm') probfcst(:,:) = prob_5mm(:,:)
IF (trim(cthresh) .eq. '10mm') probfcst(:,:) = prob_10mm(:,:)
IF (trim(cthresh) .eq. '25mm') probfcst(:,:) = prob_25mm(:,:)
IF (trim(cthresh) .eq. '50mm') probfcst(:,:) = prob_50mm(:,:)
IF (trim(cthresh) .eq. 'q50') probfcst(:,:) = prob_median(:,:)
IF (trim(cthresh) .eq. 'q67') probfcst(:,:) = prob_uppertercile(:,:)
IF (trim(cthresh) .eq. 'q80') probfcst(:,:) = prob_upperquintile(:,:)
IF (trim(cthresh) .eq. 'q90') probfcst(:,:) = prob_upperdecile(:,:)
IF (trim(cthresh) .eq. 'q95') probfcst(:,:) = prob_95(:,:)

RETURN
END subroutine  read_analog_probabilities