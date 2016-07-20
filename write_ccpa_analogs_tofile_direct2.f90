SUBROUTINE write_ccpa_analogs_tofile_direct2 (outfile_dates, outfile_data, outfile_determ, &
   outfile_quantiles,nxa, nya, nxf, nyf, ntrainsamps, iccpa_mask, prob_pop, prob_1mm, &
   prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, prob_50mm, prob_median, prob_uppertercile, &
   prob_upperquintile, prob_upperdecile, prob_95, apcp_anal, apcp_fcst, fmean_analog, &
   thresh_median, thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, &
   thresh_95, rlons_anal, rlats_anal, rlons_fcst, rlats_fcst, date_anal, date_fcst)

! ---- write analog forecast information to file
!      coded by Tom Hamill, (303) 497-3060, tom.hamill@noaa.gov

CHARACTER*120, INTENT(IN) :: outfile_dates, outfile_data, outfile_determ, outfile_quantiles
INTEGER, INTENT(IN) :: nxa, nya, nxf, nyf, ntrainsamps
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: iccpa_mask
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlons_anal, rlats_anal
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: apcp_anal
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: fmean_analog
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_pop
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_1mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_2p5mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_5mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_10mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_25mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_50mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_median
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_uppertercile
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_upperdecile
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_upperquintile
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob_95
REAL, INTENT(IN), DIMENSION(nxf,nyf,ntrainsamps) :: apcp_fcst
REAL, INTENT(IN), DIMENSION(nxa,nya) :: thresh_median
REAL, INTENT(IN), DIMENSION(nxa,nya) :: thresh_uppertercile
REAL, INTENT(IN), DIMENSION(nxa,nya) :: thresh_upperquintile
REAL, INTENT(IN), DIMENSION(nxa,nya) :: thresh_upperdecile
REAL, INTENT(IN), DIMENSION(nxa,nya) :: thresh_95


REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlons_fcst
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlats_fcst

INTEGER, INTENT(IN), DIMENSION(ntrainsamps) :: date_anal, date_fcst

! ----- write out an array of dates associated with each of the ntrainsamps samples

PRINT *,'in subroutine write_ccpa_analogs_tofile_direct2'
PRINT *,'writing to ',TRIM(outfile_dates)
OPEN (unit=1, file=outfile_dates, status='replace',form='unformatted')
WRITE (1) ntrainsamps
WRITE (1) date_anal  ! initial time of forecast
WRITE (1) date_fcst ! verif time of forecast
CLOSE (1)

! ---- write out a direct access file, one record for each training sample

irecl = 3 + (nxa*nya)/2 + 2*nxa*nya + 13*nxa*nya + 2 + 3*nxf*nyf + 2

! ---- Record length above is in words, appropriate for ifort compiler.
!      If using gfortran, multiply by 4.  Instead use words for ifort compiler

irecl = irecl*4  ! for gfortran

PRINT *,'writing to ',TRIM(outfile_data)
OPEN (unit=1, file=outfile_data, status='replace',access='direct',recl=irecl)
DO i = 1, ntrainsamps
   WRITE (1,rec=i) nxa, nya, ntrainsamps,iccpa_mask,rlons_anal,rlats_anal, &
        apcp_anal(:,:,i), prob_pop(:,:,i), prob_1mm(:,:,i), prob_2p5mm(:,:,i), &
        prob_5mm(:,:,i), prob_10mm(:,:,i), prob_25mm(:,:,i), prob_50mm(:,:,i), &
        prob_median(:,:,i), prob_uppertercile(:,:,i), prob_upperquintile(:,:,i), &
        prob_upperdecile(:,:,i), prob_95(:,:,i),nxf, nyf, apcp_fcst(:,:,i), &
        rlons_fcst(:,:), rlats_fcst(:,:), date_anal(i), date_fcst(i)
END DO
CLOSE (1)

PRINT *,'writing to ',TRIM(outfile_determ)
OPEN (unit=1, file=outfile_determ, status='replace',form='unformatted')
WRITE (1) nxa, nya, ntrainsamps   !max days = ndaysm
WRITE (1) apcp_anal
WRITE (1) fmean_analog
CLOSE (1) 
PRINT *,'Done writing to files.'
PRINT *, 'Done with subroutine write_ccpa_analogs_tofile_direct2.'

PRINT *,'writing to ',TRIM(outfile_quantiles)
OPEN (unit=1, file=outfile_quantiles, status='replace',form='unformatted')
WRITE (1) nxa, nya
WRITE (1) thresh_median
WRITE (1) thresh_uppertercile
WRITE (1) thresh_upperquintile
WRITE (1) thresh_upperdecile
WRITE (1) thresh_95
CLOSE (1) 

RETURN
END
