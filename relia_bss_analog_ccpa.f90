SUBROUTINE relia_bss_analog_ccpa ( outfile, nxa, nya, ntrainsamps, iccpa_mask, imonth, nyears, &
     date_anal, rlons_anal, rlats_anal, apcp_anal, threshfield, prob, climop )

! ------------------------------------------------------------------------------
! Calculate reliability and brier skill score of analog forecasts.
! Save information to file.
!
! Coded by Tom Hamill, Sep 2013. tom.hamill@noaa.gov, (303) 497-3060.
! ------------------------------------------------------------------------------

CHARACTER*120, INTENT(IN) :: outfile ! name of file where output should be stored
INTEGER, INTENT(IN) :: nxa, nya      ! grid dimensions
INTEGER, INTENT(IN) :: ntrainsamps   ! number of training samples
INTEGER, INTENT(IN) :: nyears        ! number of years
INTEGER, INTENT(IN) :: imonth        ! month of year, 1-12
INTEGER, INTENT(IN), DIMENSION(ntrainsamps) :: date_anal
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya)   :: iccpa_mask  ! mask, 1 if inside conus, 0 if not
REAL, INTENT(IN), DIMENSION(nxa,nya)        :: rlons_anal, rlats_anal ! longitudes, latitudes
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: apcp_anal ! analyzed precip
REAL, INTENT(IN), DIMENSION(nxa,nya)        :: threshfield ! array of threshold amounts
REAL, INTENT(IN), DIMENSION(nxa,nya,ntrainsamps) :: prob ! forecast probability of event
REAL, INTENT(IN), DIMENSION(nxa,nya)        :: climop ! climatological probability

REAL*8 bsclimo ! accumulating lots of small numbers.  To keep roundoff from being an issue
REAL*8 bs      ! use real*8's not real*4's
REAL bss
REAL, DIMENSION(nxa,nya)   :: bssmap     ! map of brier skill scores
REAL*8, DIMENSION(nxa,nya) :: bsmap      ! map of brier scores
REAL*8, DIMENSION(nxa,nya) :: bsclimomap ! map of climatological forecast brier scores
REAL*8, DIMENSION (0:20,2) :: contab     ! monthly contingency table data, (1mm event)
REAL*8                     :: ctot       ! sum of contingency table
REAL, DIMENSION(0:20)      :: frequse    ! frequency of usage 0,5 ..., 100%
REAL, DIMENSION(0:20)      :: relia      ! reliability, 0,5, ... , 100%
REAL, DIMENSION(nyears)    :: bs_yearly
REAL, DIMENSION(nyears)    :: bsclimo_yearly

! ---- initialize

contab      = 0.
bs          = 0.
bsclimo     = 0.
bsmap       = 0.
bsclimomap  = 0.
bs_yearly(:) = 0.
bsclimo_yearly(:) = 0.
pid180      = 3.1415926 / 180.

! ---- loop thru grid points and dates.  Increment contingency tables and scores
!      with this sample's information.

DO iday = 1, ntrainsamps

   ! ---- determine the year and month

   iyear = date_anal(iday) / 1000000
   immddhh = date_anal(iday) - iyear*1000000
   imm = immddhh / 10000
   iyy = iyear - 2001 ! data starts in 2002

   DO i = 1, nxa
      DO j = 1, nya
         IF (iccpa_mask(i,j) .gt. 0 .and. prob(i,j,iday) .GE. 0.0 .and. &
         prob(i,j,iday) .le. 1.0 .and. apcp_anal(i,j,iday) .GE. 0.0 .and. imm .eq. imonth) THEN
            cfac = cos(rlats_anal(i,j)*pid180)  ! cos of latitude to acct for unequal grid box size
            pclimo = climop(i,j) 
            p      = prob(i,j,iday)
            ipcat  = nint(p*20)
            v      = apcp_anal(i,j,iday)
            IF (v .GE. threshfield(i,j)) THEN
               contab(ipcat,2) = contab(ipcat,2) + cfac
               bs              = bs + cfac*(1.-p)**2
               bsclimo         = bsclimo + cfac*(1.-pclimo)**2
               bsmap(i,j)      = bsmap(i,j) + cfac*(1.-p)**2
               bsclimomap(i,j) = bsclimomap(i,j) + cfac*(1.-pclimo)**2
               bs_yearly(iyy)  = bs_yearly(iyy) + cfac*(1.-p)**2
               bsclimo_yearly(iyy) = bsclimo_yearly(iyy) + cfac*(1.-pclimo)**2
            ELSE
               contab(ipcat,1) = contab(ipcat,1) + cfac
               bs              = bs + cfac * p**2
               bsclimo         = bsclimo + cfac*pclimo**2
               bsmap(i,j)      = bsmap(i,j) + cfac*p**2
               bsclimomap(i,j) = bsclimomap(i,j) + cfac*pclimo**2
               bs_yearly(iyy)  = bs_yearly(iyy) + cfac*p**2
               bsclimo_yearly(iyy) = bsclimo_yearly(iyy) + cfac*pclimo**2
            ENDIF
         ENDIF 
      END DO   ! j = 1, ny
   END DO      ! i = 1, nx
END DO         ! iday

! ---- with tallied contingency tables, now set reliability and frequency of use

ctot     = SUM(contab)
bss      = 1. - bs / bsclimo
relia(:) = -99.9999
PRINT *,'  bss = ',bss
PRINT *,'  p   reliability   freq of usage'
DO icat = 0,20
  frequse(icat) = (contab(icat,2)  + contab(icat,1)) / ctot
  IF ((contab(icat,1) + contab(icat,2)) .gt. 0) THEN
    relia(icat) = contab(icat,2) / (contab(icat,1) + contab(icat,2))
    PRINT 203,float(icat*5),relia(icat)*100.,frequse(icat)
  ELSE
    PRINT 203,float(icat*5),relia(icat),frequse(icat)
  ENDIF
  203 format(f5.1,3x,2(f8.3,3x))
END DO  !icat

DO i = 1, nxa
  DO j = 1, nya
     IF (bsclimomap(i,j) .GT. 0.0) THEN
        bssmap(i,j) = 1. - bsmap(i,j) / bsclimomap(i,j)
     ELSE
        bssmap(i,j) = -999.99
     ENDIF
  END DO
END DO

! ---- write reliability, BSS information to file
 
PRINT *,'writing reliability information to ',TRIM(outfile)
OPEN (unit=1, file=outfile, status='replace', form='unformatted')
WRITE (1) bss, bs, bsclimo
WRITE (1) frequse
WRITE (1) relia
WRITE (1) contab
WRITE (1) bssmap
WRITE (1) bsmap
PRINT *,'max, min of bsclimomap = ',maxval(bsclimomap), minval(bsclimomap)
WRITE (1) bsclimomap
WRITE (1) rlons_anal
WRITE (1) rlats_anal
WRITE (1) iccpa_mask
WRITE (1) nyears
WRITE (1) bs_yearly
WRITE (1) bsclimo_yearly
CLOSE (1)

RETURN
END SUBROUTINE relia_bss_analog_ccpa
