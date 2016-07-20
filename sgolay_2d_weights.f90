SUBROUTINE sgolay_2d_weights(window_size, order, istat, weights)

! --- a conversion to Fortran90 of the sgolay2d algorithm for python found on the web
!     (http://wiki.scipy.org/Cookbook/SavitzkyGolay).  The original python code is
!     included too, in the comments.  This version only returns the weight for the
!     smoothing kernel.  The actual smoothing is done by a different subroutine.
!
!     Coded (kludged/adapted) by Tom Hamill, NOAA, tom.hamill@noaa.gov, June 2014
!  
!     Uses numerical recipes routines for SVD decomposition, 
!     INTEGER, PARAMETER :: window_size = 9 (+/- 4 grid pts)
!     INTEGER, PARAMETER :: order = 3

INTEGER, INTENT(IN) :: window_size ! as defined in Press et al. Numerical Recipes
INTEGER, INTENT(IN) :: order ! as defined in Press et al. Numerical Recipes

REAL*4, INTENT(OUT), DIMENSION(window_size,window_size) :: weights

INTEGER half_size
REAL n_terms
REAL*8, DIMENSION((order+2)**2,2):: exps ! more elements than we need

REAL*8, ALLOCATABLE, DIMENSION(:) :: E
REAL*8, ALLOCATABLE, DIMENSION(:) :: ind, dx, dy, W, Aplus_row, work
REAL*8, ALLOCATABLE, DIMENSION(:,:) :: A, U, V, Aplus

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

!# exponents of the polynomial.
!    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ...
!    # this line gives a list of two item tuple. Each tuple contains
!    # the exponents of the k-th term. First element of tuple is for x
!    # second element for y.
!    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
!    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]
ktr = 1
DO k = 0,order
   DO n = 0,k
      !print *,'k, n = ',k,n
      exps(ktr,1) = k-n
      exps(ktr,2) = n
      ktr = ktr+1
   END DO
END DO
ktr = ktr-1

! ---- initialize coordinates of points
!      ind = np.arange(-half_size, half_size+1, dtype=np.float64)
!      dx = np.repeat( ind, window_size )
!      dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

lenind = 2*half_size+1
ALLOCATE (ind(lenind))
DO i = -half_size, half_size
   ind(i+half_size+1) = DBLE(i)
END DO
lendxdy = window_size*lenind
ALLOCATE (dx(lendxdy), dy(lendxdy))

DO i = 1, window_size
   istart = (i-1)*lenind+1
   iend = istart+lenind-1
   dx(istart:iend) = ind(i)
   dy(istart:iend) = ind(1:lenind)
END DO

! ---- build matrix of system of equation
!      A = np.empty( (window_size**2, len(exps)) )
!      for i, exp in enumerate( exps ):
!          A[:,i] = (dx**exp[0]) * (dy**exp[1])

nrows = window_size**2
ncols = ((order+2)*(order+1))/2
ALLOCATE(A(nrows,ncols), U(nrows,nrows), W(ncols), E(ncols), V(nrows,ncols), &
   Aplus(nrows,nrows), work(nrows), Aplus_row(nrows))

A(:,:)  = 0.0
E(:)    = 0.0
work(:) = 0.0
DO i = 1, ncols
   A(:,i) = (dx(:)**exps(i,1))*(dy(:)**exps(i,2))
END DO
!
!    # solve system and convolve
!        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
!        return scipy.signal.fftconvolve(Z, m, mode='valid')

! ---- calculate pseudo inverse.  Given SVD decomposition A = U W V^T, then pseudo-inverse is
!      A+ = V W^(-1) U^T.  Use Moore-Penrose pseudo-inverse, gotten from netlib

print *,'nrows, ncols = ', nrows, ncols
CALL mpinv(nrows, nrows, ncols, ncols, A, Aplus, W, E, U, V, work, irank, ierr)

Aplus_row(:) = Aplus(1,:)
ktr = 0
DO i = 1, window_size
   DO j = 1, window_size
      ktr = ktr+1
      weights(i,j) = Aplus_row(ktr)
   END DO
END DO

PRINT *,'Savitzky-Golay smoother weights = '
DO i = 1, 9
   print 368, (weights(i,j),j=1,9)
   368 format(9(f7.4,1x))
END DO

DEALLOCATE (E, ind, dx, dy, W, Aplus_row, work, A, U, V, Aplus)

RETURN
END SUBROUTINE sgolay_2d_weights

