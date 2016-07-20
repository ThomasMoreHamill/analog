SUBROUTINE determine_nearest_nint (nxa, nya, nxf, nyf, lonsa, latsa, lonsf, latsf, &
   inearest, jnearest)

! purpose:  for each analysis (i,j) location, identify the (i,j) index of the
!    forecast grid point lower than or equal to the lat/lon of the analyzed.
!
! called by ppn_analog_ccpa_supplocns.f90

INTEGER, INTENT(IN) :: nxa, nya, nxf, nyf
REAL, INTENT(IN), DIMENSION(nxa,nya) :: lonsa, latsa
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: lonsf, latsf
INTEGER, INTENT(OUT), DIMENSION(nxa, nya) :: inearest, jnearest

inearest(:,:) = -99
jnearest(:,:) = -99

! --- for each analysis x-grid point index, find the forecast index closest
!     to it but with a lower value

DO ixa = 1, nxa
   imin = ixf-1
   rlona = lonsa(ixa,1) ! presume regular lat/lon grid
   DO ixf = 1, nxf-1
      IF (rlona .ge. lonsf(ixf,1) .and. rlona .lt. lonsf(ixf+1,1)) THEN
         imin = ixf
      ENDIF
   END DO
   inearest(ixa,:) = imin
END DO

! ---- same but for latitudes.

DO jya = 1, nya
   jmin = nyf-1
   rlata = latsa(1,jya) ! presume regular lat/lon grid
   DO jyf = 1, nyf-1
      IF (rlata .ge. latsf(1,jyf) .and. rlata .lt. latsf(1,jyf+1)) THEN
         jmin = jyf
      ENDIF
   END DO
   jnearest(:,jya) = jmin
END DO

RETURN
END SUBROUTINE determine_nearest_nint
