! FILE: get_ceedeef_precip.f90
!  f2py -c -m get_ceedeef_precip get_ceedeef_precip.f90
!
! this routine contains fortran code, compiled via f2py, into a python interface,
! that gets the CDF for precipitation.
! -------


SUBROUTINE get_ceedeef_precip(nthresh, nj, ni, thresh, apcpf, CDFwork, istat)
INTEGER, INTENT(IN) :: nthresh, nj, ni
REAL, INTENT(IN), DIMENSION(nthresh) :: thresh
!REAL, INTENT(IN), DIMENSION(ni,nj) :: apcpf
REAL, INTENT(IN), DIMENSION(nj,ni) :: apcpf
REAL*8, INTENT(OUT), DIMENSION(nthresh,nj,ni) :: CDFwork
INTEGER, INTENT(OUT) :: istat

!f2py intent(in) nthresh, nj, ni, thresh, apcpf
!f2py depend(nthresh) thresh
!f2py depend(nj,ni) apcpf
!f2py intent(out) CDFwork, istat
!f2py depend(nthresh, nj, ni) CDFwork

CDFwork = 0.
apmax = MAXVAL(apcpf)
IF (apmax .ge. 0.0 .and. apmax .le. 10000.) THEN
   istat = 1
   DO jy = 1, nj
      DO ix = 1, ni
         DO ithr = 1, nthresh
            IF (apcpf(jy,ix) .LE. thresh(ithr)) &
                CDFwork(ithr,jy,ix) = CDFwork(ithr,jy,ix) + 1.0
         END DO
      END DO
   END DO
ELSE
   istat = -1
ENDIF

RETURN
END SUBROUTINE get_ceedeef_precip








