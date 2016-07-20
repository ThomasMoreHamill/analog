SUBROUTINE barnes_like(conusmask, precip_in, nxa, nya, precip_out)

! ---- perform a Barnes-like objective analysis of points over water, i.e., a weighted sum
!      of the values of nearby land points, with the weights being an exponential function of distance.

INTEGER, INTENT(IN) :: nxa, nya
INTEGER*2, INTENT(IN), DIMENSION(nxa,nya) :: conusmask
REAL, INTENT(INOUT), DIMENSION(nxa,nya) :: precip_in
REAL, INTENT(OUT), DIMENSION(nxa,nya) :: precip_out

REAL*8, DIMENSION(nxa,nya) :: numer, denom
REAL*8 weight

pavg = SUM(precip_in*REAL(conusmask))/SUM(real(conusmask))

! ---- now do something akin to a Barnes filtering, weighting the data by exp(-dist**2/lengthscale**2)
!      this loop tallies up the numerator and denominator needed for weighted sum calculation

icut = 10
rlenscale = 3.
rlenscale2 = rlenscale**2
numer = 0.
denom = 0.
DO ixa = 1, nxa
   !IF (mod(ixa,10) .eq. 0) PRINT *,'ixa = ',ixa,' of ',nxa
   imin = MAX(1,ixa-icut)
   imax = MIN(nxa,ixa+icut)
   DO jya = 1, nya
      jmin = MAX(1,jya-icut)
      jmax = MIN(nya,jya+icut)
      IF (conusmask(ixa,jya) .eq. 0) THEN
         DO iloop = imin, imax
            DO jloop = jmin, jmax
               dist = SQRT(REAL((ixa-iloop)**2 + (jya-jloop)**2))
               weight = exp(-dist**2/rlenscale2)
               IF (dist .le. icut .and. conusmask(iloop,jloop) .eq. 1 .and. precip_in(iloop,jloop) .ge. 0.0) THEN
                  numer(ixa,jya) = numer(ixa,jya)  + weight*precip_in(iloop,jloop)
                  denom(ixa,jya) = denom(ixa,jya)  + weight
               ENDIF
            END DO
         END DO
      END IF
   END DO ! jya
END DO ! ixa

! --- set the output over water to the weighted sum.  over land use original values.  If
!     at a point that is too far from land to have any points counted in weight, just set the 
!     value to the domain-averaged land probability.

DO ixa = 1, nxa
   DO jya = 1, nya
      IF (denom(ixa,jya) .gt. 0.0 .and. conusmask(ixa,jya) .eq. 0) THEN 
         precip_out(ixa,jya) = numer(ixa, jya) /denom (ixa, jya) 
      ELSE
         precip_out(ixa,jya) = precip_in(ixa,jya)
      ENDIF
      IF (precip_out(ixa,jya) .lt. 0.0 .or. precip_out(ixa,jya) .gt. 1.0) precip_out(ixa,jya) = pavg
   END DO
END DO

RETURN
END SUBROUTINE barnes_like
