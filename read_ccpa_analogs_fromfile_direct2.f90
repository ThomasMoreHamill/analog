SUBROUTINE read_ccpa_analogs_fromfile_direct2 (infile_dates, infile_data, &
   nxa, nya, nxf, nyf, ntrainsamps, iccpa_mask, prob_pop, prob_1mm, prob_2p5mm, &
   prob_5mm, prob_10mm, prob_25mm, prob_50mm, prob_median, prob_uppertercile, &
   prob_upperquintile, prob_upperdecile, prob_95, nxfin, nyfin, apcp_anal, &
   apcp_fcst, rlons_anal, rlats_anal, rlons_fcst, rlats_fcst, date_anal, date_fcst)

! ---- assuming that the analog information for this time has already been computed,
!      read it in from a file rather than processing again.
!      coded by Tom Hamill, (303) 497-3060, tom.hamill@noaa.gov

CHARACTER*(*), INTENT(IN) :: infile_dates, infile_data
INTEGER, INTENT(IN) :: nxa, nya, ntrainsamps, nxf, nyf
INTEGER*2, DIMENSION(nxa,nya), INTENT(OUT) :: iccpa_mask
REAL, DIMENSION(nxa,nya), INTENT(OUT) :: rlons_anal, rlats_anal
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: apcp_anal
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_pop
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_1mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_2p5mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_5mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_10mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_25mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_50mm
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_median
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_uppertercile
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_upperdecile
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_upperquintile
REAL, INTENT(OUT), DIMENSION(nxa,nya,ntrainsamps) :: prob_95
REAL, INTENT(OUT), DIMENSION(nxf,nyf,ntrainsamps) :: apcp_fcst
REAL, INTENT(OUT), DIMENSION(nxf,nyf) :: rlons_fcst
REAL, INTENT(OUT), DIMENSION(nxf,nyf) :: rlats_fcst

INTEGER, INTENT(OUT), DIMENSION(ntrainsamps) :: date_anal, date_fcst

! ----- write out an array of dates associated with each of the ntrainsamps samples

PRINT *,'reading from ',TRIM(infile_dates)
OPEN (unit=2, file=infile_dates, status='old',form='unformatted')
READ (2) date_anal  ! initial time of forecast
READ (2) date_fcst ! verif time of forecast
CLOSE (2)

! ---- write out a direct access file, one record for each training sample

irecl = 3 + (nxa*nya)/2 + 2*nxa*nya + 13*nxa*nya + 2 + 3*nxf*nyf + 2
PRINT *,'reading from ',TRIM(infile_data)

! i believe record length is in bytes for gfortran, words for ifort compiler

OPEN (unit=3, file=infile_data, status='old',access='direct',recl=irecl*4)
DO i = 1, ntrainsamps
   READ (3,rec=i) nxain, nyain, ntrainsampsin,iccpa_mask,rlons_anal,rlats_anal,&
        apcp_anal(:,:,i), prob_pop(:,:,i), prob_1mm(:,:,i), prob_2p5mm(:,:,i), &
        prob_5mm(:,:,i), prob_10mm(:,:,i), prob_25mm(:,:,i), prob_50mm(:,:,i), &
        prob_median(:,:,i), prob_uppertercile(:,:,i), prob_upperquintile(:,:,i),  &
        prob_upperdecile(:,:,i), prob_95(:,:,i), nxfin, nyfin, apcp_fcst(:,:,i), &
        rlons_fcst(:,:), rlats_fcst(:,:), date_anal(i), date_fcst(i)
END DO
CLOSE (3)
PRINT *,'Done reading from files.'



RETURN
END
