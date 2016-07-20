SUBROUTINE calculate_reliabilities_rawens (cmonth, cleadb, cleade, &
	output_data_directory, nxa, nya, ntrainsamps, &
    nyears, iccpa_mask, rlons_anal, rlats_anal, apcp_anal, date_anal, imonth, thresh_95, &
    thresh_uppertercile, thresh_upperquintile, thresh_upperdecile, thresh_median, climo_pop, &
    climo_1mm, climo_2p5mm, climo_5mm, climo_10mm, climo_25mm, climo_50mm, &
    climo_median, climo_uppertercile, climo_upperquintile, climo_upperdecile, &
    climo_95, prob_pop, prob_1mm, prob_2p5mm, prob_5mm, prob_10mm, prob_25mm, &
    prob_50mm, prob_median, prob_uppertercile, prob_upperquintile, prob_upperdecile, prob_95)

! ---- calculate brier skill scores, reliability information, and store that information to file.
!      coded by Tom Hamill, Sep 2013, (303) 497-3060, tom.hamill@noaa.gov

CHARACTER*3, INTENT(IN) :: cmonth, cleadb, cleade
CHARACTER*(*), INTENT(IN) :: output_data_directory

INTEGER, INTENT(IN) :: nxa, nya, ntrainsamps, nyears
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: iccpa_mask
INTEGER, INTENT(IN) :: imonth
INTEGER, INTENT(IN), DIMENSION(ntrainsamps) :: date_anal

REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlons_anal, rlats_anal
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: apcp_anal
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_95
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_uppertercile
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_upperquintile
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_upperdecile
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  thresh_median
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_pop
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_1mm
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_2p5mm
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_5mm
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_10mm
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_25mm
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_50mm
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_median
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_uppertercile
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_upperquintile
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_upperdecile
REAL, INTENT(IN), DIMENSION(nxa,nya)  ::  climo_95
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_95
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_pop
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_1mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_2p5mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_5mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_10mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_25mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_50mm
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_median
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_uppertercile
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_upperquintile
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps)  ::  prob_upperdecile

REAL, DIMENSION(nxa,nya) :: threshfield
REAL, DIMENSION(nxa,nya) :: climop

CHARACTER*120 outfile

threshfield(:,:) = 0.1 ! POP
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_POP_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask, imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_pop, climo_pop)

threshfield = 1.0
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_1mm_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask, imonth, nyears, &
   date_anal,  rlons_anal, rlats_anal, apcp_anal, threshfield, prob_1mm, climo_1mm )

threshfield = 2.5
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_2p5mm_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_2p5mm, climo_2p5mm )

threshfield = 5.0
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_5mm_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth,  nyears,&
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_5mm, climo_5mm )

threshfield = 10.0
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_10mm_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_10mm, climo_10mm )

threshfield = 25.0
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_25mm_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_25mm, climo_25mm )

threshfield = 50.0
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_50mm_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth,  nyears,&
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_50mm, climo_50mm)

! normally below you'd think that for the quantiles of the distribution the climatological
! probability would simply be the quantile, e.g., for the 5th percentile, the climo
! probability would be 0.05.  However, with precip, the 5th as well as 10th, 20th, ... 
! may all be 0.0.  So previously if the quantile q was actually associated with the value
! 0, the quantile threshold was reset to the first nonzero precipitation value, and 
! q was set to the climatological probability of this event, i.e., the probability of 
! nonzero precipitation.

threshfield = thresh_median
climop = climo_median
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_q50_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_median, climop )

threshfield = thresh_uppertercile
climop = climo_uppertercile
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_q67_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_uppertercile, climop )

threshfield = thresh_upperquintile
climop = climo_upperquintile
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_q80_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_upperquintile, climop )

threshfield = thresh_upperdecile
climop = climo_upperdecile
PRINT *, '-------------------------------------------------------------------------'
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_q90_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_upperdecile, climop )

threshfield = thresh_95
climop = climo_95
outfile = TRIM(output_data_directory) // 'reliability_ppn_rawens_q95_'//&
   cmonth//'_'//cleadb//'_to_'//cleade//'.dat'
PRINT *, '-------------------------------------------------------------------------'
print *,'writing reliability information to ',TRIM(outfile)
CALL relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask,  imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob_95, climop )

RETURN
END
