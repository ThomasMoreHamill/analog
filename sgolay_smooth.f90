SUBROUTINE sgolay_smooth(nx, ny, z, conusmask, weights, window_size, order, istat)

! --- a conversion to Fortran90 of the sgolay2d algorithm for python found on the web
!     (http://wiki.scipy.org/Cookbook/SavitzkyGolay).  The original python code is
!     included too, in the comments.  This subroutine inputs pre-calculated weights
!     (computed in sgolady_2d_weights) and smooths the array z.  For the time being,
!     the smoothing is done in physical (grid point) space as opposed to an FFT 
!     convolution.
!
!     Coded (kludged/adapted) by Tom Hamill, NOAA, tom.hamill@noaa.gov
!  

INTEGER, INTENT(IN) :: nx ! row dimension of z
INTEGER, INTENT(IN) :: ny ! column dimension of z
INTEGER, INTENT(IN) :: window_size ! as defined in Press et al. Numerical Recipes
INTEGER, INTENT(IN) :: order ! as defined in Press et al. Numerical Recipes
INTEGER*2, INTENT(IN), DIMENSION(nx,ny) :: conusmask
REAL*4, INTENT(INOUT), DIMENSION(nx,ny) :: z
REAL*4, INTENT(IN), DIMENSION(window_size, window_size) :: weights
INTEGER, INTENT(OUT) :: istat

INTEGER half_size
REAL n_terms
REAL*4, DIMENSION(ny) :: bandtopnbot
REAL*4, DIMENSION(nx) :: bandrandl

REAL*4, ALLOCATABLE, DIMENSION(:) :: bandcorner
REAL*4, ALLOCATABLE, DIMENSION(:,:) :: corner, corner_flipped, trcorner, &
    trcorner_flipped, blcorner, blcorner_flipped, ZZ, z2

ALLOCATE (z2(nx,ny)) 

! ---- input data will have zero or missing data at water points.  If Savitzky-
!      Golay smoother is run directly on this data, it will thus taper to zero
!      precip at coasts.  Hence, fill in ocean points with a smoothed field
!      reflecting the values at nearby land points.

CALL barnes_like(conusmask, z, nx, ny, z2)
z = z2

! ---- # number of terms in the polynomial expression
!      n_terms = ( order + 1 ) * ( order + 2)  / 2.0
!    if  window_size % 2 == 0:
!        raise ValueError('window_size must be odd')
!
!    if window_size**2 < n_terms:
!        raise ValueError('order is too high for the window size')
!
!    half_size = window_size // 2

n_terms = ( order + 1 ) * ( order + 2)  / 2.0
istat = 0
IF (MOD(window_size,2) .eq. 0) THEN
   PRINT *,'window_size must be odd'
   istat = -1
   RETURN
ELSE IF (window_size**2 < n_terms) THEN
   PRINT *,'order is too high for the window size'
   istat = -2
   RETURN
ENDIF
half_size = window_size / 2
ALLOCATE(corner(half_size,half_size), corner_flipped(half_size,half_size))

! ----  pad input array with appropriate values at the four borders

nx2 = nx + 2*half_size
ny2 = ny + 2*half_size
ALLOCATE (ZZ(nx2,ny2))
ZZ(:,:) = 0.0

! ----# left band (in Fortran)
!     band = z[0, :]
!     Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )

bandtopnbot(:) = z(1,:)
DO j = half_size+1, ny2-half_size
   jz = j - half_size
   ibegin = 1
   iend = half_size
   idiff = iend-ibegin
   ibeginz = 2
   iendz = 1+half_size
   DO i = 0,idiff
       ZZ(ibegin+i,j) = bandtopnbot(jz) - ABS(z(iendz-i,jz) -bandtopnbot(jz))
   END DO
END DO

!---- # right band
!     band = z[-1, :]
!     Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band )
!
bandtopnbot = z(nx,:)
DO j = half_size+1, ny2-half_size
   jz = j - half_size
   ibegin = nx2 - half_size + 1
   iend = nx2
   ibeginz = nx - half_size 
   iendz = nx - 1
   idiff = iend - ibegin
   DO i = 0, idiff
       ZZ(ibegin+i,j) = bandtopnbot(jz) + ABS(z(iendz-i,jz) -bandtopnbot(jz))
   END DO
END DO

!---- # bottom band
!     band = np.tile( z[:,0].reshape(-1,1), [1,half_size]) # -1 refers to the inferred dimension
!     Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )

bandrandl = z(:,1)
jbegin = 1
jend = half_size
jdiff = jend - jbegin
jbeginz = 2
jendz = 1+half_size
DO j = 0, jdiff
   DO i = half_size+1, nx2-half_size
      iz = i - half_size
      ZZ(i,jbegin+j) = bandrandl(iz) - ABS (z(iz,jendz-j)-bandrandl(iz))
   END DO
END DO

! ---- # top band
!      band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
!      Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )

bandrandl = z(:,ny)
jbegin = ny2-half_size+1
jend = ny2
jdiff = jend - jbegin
jbeginz = ny - half_size-1
jendz = ny - 1
DO j = 0, jdiff
   DO i = half_size+1, nx2-half_size
      iz = i - half_size
      ZZ(i,jbegin+j) = bandrandl(iz) + ABS (z(iz,jendz-j)-bandrandl(iz))
   END DO
END DO

! ---- # central band
!      Z[half_size:-half_size, half_size:-half_size] = z
!
ZZ(half_size+1:nx2-half_size, half_size+1:ny2-half_size) = z(:,:)

! ---- # bottom left corner
!      band = z[0,0]
!      Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )

band = z(1,1)
corner(:,:) = z(2:half_size+1,2:half_size+1)
DO i = 1, half_size
   iud = half_size+1-i
   DO j = 1, half_size
       jud = half_size+1-j
       corner_flipped(i,j) = corner(iud,jud)
   END DO
END DO
ZZ(1:half_size,1:half_size) = band - ABS(corner_flipped(:,:)-band)
!
! -----  # top right corner
!        band = z[-1,-1]
!        Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band )

band = z(nx,ny)
corner(:,:) = z(nx-half_size:nx-1,ny-half_size:ny-1)
DO i = 1, half_size
   iud = half_size+1-i
   DO j = 1, half_size
       jud = half_size+1-j
       corner_flipped(i,j) = corner(iud,jud)
   END DO
END DO
ZZ(nx2-half_size+1:nx2,ny2-half_size+1:ny2) = band + ABS(corner_flipped(:,:)-band)

! ---- # top left corner
!      band = Z[half_size,-half_size:]
!      Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band )
!
ALLOCATE(bandcorner(half_size))
bandcorner(:) = ZZ(half_size+1,ny2-half_size+1:ny2)

ALLOCATE(trcorner(half_size,half_size), trcorner_flipped(half_size,half_size))
trcorner(:,:) = ZZ(half_size+2:2*half_size+1,ny2-half_size+1:ny2)

DO i = 1, half_size
   trcorner_flipped(i,:) = trcorner(half_size+1-i,:)
END DO

DO i = 1, half_size
   ZZ(i,ny2-half_size+1:ny2) = bandcorner(:) - ABS(trcorner_flipped(i,:) - bandcorner(:))
END DO

! ---- # bottom right corner
!      band = Z[-half_size:,half_size].reshape(-1,1)
!      Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )
!
bandcorner(:) = ZZ(nx2-half_size+1:nx2,half_size+1)
ALLOCATE(blcorner(half_size,half_size), blcorner_flipped(half_size,half_size))
blcorner(:,:) = ZZ(nx2-half_size+1:nx2,half_size+2:2*half_size+1)
DO j = 1, half_size
   blcorner_flipped(:,j) = blcorner(:,half_size+1-j)
END DO

ktr=1
DO j = 1, half_size
   ZZ(nx2-half_size+1:nx2,j) = bandcorner(:) - ABS(blcorner_flipped(:,j) - bandcorner(:))
END DO

! --- smooth the data

DO i = 1 + half_size, nx + half_size
	iorig = i - half_size
	DO j = 1 + half_size, ny + half_size
		jorig = j - half_size
		z(iorig, jorig) = SUM(weights(:,:) * zz(i-half_size:i+half_size, j-half_size:j+half_size))
	END DO
END DO

DEALLOCATE (blcorner, blcorner_flipped, trcorner, trcorner_flipped, bandcorner,&
   corner, corner_flipped, ZZ, z2)

! ---- last step:  re-set the smoothed values to missing outside of CONUS

DO i = 1, nx
   DO j = 1, ny
      IF (conusmask(i,j) .ne. 1)  z(i,j) = -99.99
   END DO
END DO

RETURN
END SUBROUTINE sgolay_smooth






