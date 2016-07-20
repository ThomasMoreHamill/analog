SUBROUTINE list_to_nearest_fcstgridpt ( nxa, nya, nxf, nyf, nmaxlist, rlons_anal, &
   rlats_anal, rlons_fcst, rlats_fcst, iccpa_mask, iccpa_xlist, iccpa_ylist, iccpa_ktr )

! Purpose: determine, for each analysis grid point, which forecast grid point is nearest
! to it.  Make a list for each forecast grid point of the indices of the 
! analysis grid points that are closest to it.  Return that list 
! (iccpa_xlist, iccpa_ylist) and how many are closest for each grid point (iccpa_ktr)
!
! coded by Tom Hamill, August 2013.  tom.hamill@noaa.gov, (303) 497-3060

INTEGER, INTENT(IN) :: nxa, nya, nxf, nyf, nmaxlist
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: iccpa_mask
REAL, INTENT(IN), DIMENSION(nxa,nya) :: rlons_anal, rlats_anal
REAL, INTENT(IN), DIMENSION(nxf,nyf) :: rlons_fcst, rlats_fcst
INTEGER, INTENT(OUT), DIMENSION(nxf,nyf,nmaxlist) :: iccpa_xlist, iccpa_ylist 
INTEGER, INTENT(OUT), DIMENSION(nxf,nyf) :: iccpa_ktr

INTEGER, DIMENSION(nxa) :: imin
INTEGER, DIMENSION(nya) :: jmin

!PRINT *,'subroutine list_to_nearest_gridpoint'

! ---- find the index of the nearest forecast longitude to the longitude of 
!      this analysis point.  This will be repeatable for every analysis latitude
!      since it's on a regular lat-lon grid.

DO ixa = 1, nxa
   distmin_lon = 99999999.
   imin(ixa) = 9999
   DO ixf = 1, nxf
      dist2 = (rlons_anal(ixa,1)-rlons_fcst(ixf,1))**2
      IF (dist2 .lt. distmin_lon) THEN
         distmin_lon = dist2
         imin(ixa) = ixf
      ENDIF
   END DO
   !print *,'ixa, distmin_lon, imin = ',ixa,distmin_lon,imin(ixa)
END DO

! ---- find the index of the nearest forecast latitude to the latitude of 
!      this analysis point.  This will be repeatable for every analysis longitude
!      since it's on a regular lat-lon grid.

DO jya = 1, nya
   distmin_lat = 99999999.
   jmin(jya) = 9999
   DO jyf = 1, nyf
      dist2 = (rlats_anal(1,jya)-rlats_fcst(1,jyf))**2
      IF (dist2 .lt. distmin_lat) THEN
         distmin_lat = dist2
         jmin(jya) = jyf
      ENDIF
   END DO
END DO

! ---- Now, finally populate the lists for each forecast grid point with 
!      the indices of the analysis grid points that are closest to it.
!      Keep track also for each forecast grid point of just how many 
!      analysis points are closest.

iccpa_ktr(:,:) = 0
DO ixa = 1, nxa
   DO jya = 1, nya
      IF (iccpa_mask(ixa,jya) .gt. 0) THEN
         iccpa_ktr(imin(ixa),jmin(jya)) = iccpa_ktr(imin(ixa),jmin(jya)) + 1
         iccpa_xlist(imin(ixa),jmin(jya),iccpa_ktr(imin(ixa),jmin(jya))) = ixa
         iccpa_ylist(imin(ixa),jmin(jya),iccpa_ktr(imin(ixa),jmin(jya))) = jya
      END IF
   END DO
END DO
!PRINT *,'iccpa_xlist(nxf/2,nyf/2,:) = ',iccpa_xlist(nxf/2,nyf/2,:) 

!INTEGER, INTENT(OUT), DIMENSION(nxf,nyf,nmaxlist) :: iccpa_xlist, iccpa_ylist 
!INTEGER, INTENT(OUT), DIMENSION(nxf,nyf) :: iccpa_ktr
!OPEN (unit=1,file='HMT_SE/list_to_nearest_gridpt.dat',status='replace',&
!   form='unformatted')
!WRITE (1) nxf, nyf, nmaxlist
!WRITE (1) iccpa_xlist 
!WRITE (1) iccpa_ylist 
!WRITE (1) iccpa_ktr
!CLOSE (1) 
!STOP

RETURN
END
