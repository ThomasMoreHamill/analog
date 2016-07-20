        SUBROUTINE MPINV (SIZE,N,M,MINNM,A,AINV,S,E,U,V,
     1                     WORK,IRANK,IERR)
c	 CALL mpinv(nrows, nrows, ncols, ncols, A, Aplus, W, E, U, V, work, irank, ierr)
C	 ALLOCATE(A(nrows,ncols), U(nrows,nrows), W(ncols), E(ncols), V(nrows,ncols), &
C	    Aplus(nrows,nrows), work(nrows), Aplus_row(nrows))
C	 
C
C       FUNCTION:
CF
CF      Calculates the Moore-Penrose pseudo-inverse of the N by M
CF      matrix A, and stores the result in AINV.  The singular
CF      values of A are returned in S.  The left singular vectors
CF      are returned in U, and the right singular vectors are
CF      returned in V.
CF
C       USAGE:
CU
CU      CALL MPINV (SIZE,N,M,MINNM,A,AINV,S,E,U,V,WORK,IERR)
CU
C       INPUTS:
CI
CI      SIZE    integer - leading dimension of all 2-dim. arrays.
CI              SIZE must be greater than or equal to max(N,M).
CI
CI      N       integer - first dimension of matrix A.
CI
CI      M       integer - second dimension of matrix A.
CI
CI      MINNM   integer - MINNM = min(N,M).
CI
CI      A       double precision - array of dimension SIZE by M
CI              containing the N by M matrix A.
CI
C       OUTPUTS:
CO     
CO      AINV    double precision - array of dimension SIZE by N 
CO              containing the M by N matrix AINV which is the 
CO              pseudo-inverse of A.
CO
CO      S       double precision - vector of dimension MINNM containing
CO              the singular values of A if IERR = 0; otherwise, see
CO              the description of LINPACK routine DSVDC.
CO
CO      E       double precision - vector of dimension MINNM used by
CO              LINPACK routine DSVDC; normally zero.
CO
CO      U       double precision - array of dimension SIZE by N
CO              containing the left singular vectors of A.
CO
CO      V       double precision - array of dimension size by M
CO              containing the right singular vectors of A.
CO
CO      WORK    double precision - vector of dimension N used by
CO              LINPACK routine DSVDC.
CO
CO      IRANK   integer - the calculated rank of the matrix A used
CO              in the construction of the pseudo-inverse.
co
CO      IERR    integer - error condition returned from dsvdc.  If
CO              IERR is not zero, then AINV is not correct.
CO
C       ALGORITHM:
CA
CA      Computes the singular value decomposition of A and forms
CA      the sum of the outer products of the left and right singular
CA      vectors of the matrix A divided by the corresponding singular
CA      value, for all non-zero singular values.
CA
C       MACHINE DEPENDENCIES:
CM
CM       none
CM
C       HISTORY:
CH
CH      written by:             Birdwell & Laub
CH      date:                   May 17, 1985
CH      current version:        2.0
CH      modifications:          modifications so that it works 
CH                              under standard F77:bb:8-86.
CH                              add minnm argument:bb:8-86.
C       ROUTINES CALLED:
CC
CC      LINPACK: DSVDC, DROT, DAXPY, DDOT, DSCAL, DSWAP, DNRM2, DROTG
CC
C       COMMON MEMORY USED:
CM
CM       none
CM
C----------------------------------------------------------------------
C       written for:    The CASCADE Project
C                       Oak Ridge National Laboratory
C                       U.S. Department of Energy
C                       Contract number DE-AC05-840R21400
C                       Subcontract number 37B-7685 S13
C                       Organization:   The University of Tennessee
C----------------------------------------------------------------------
C       THIS SOFTWARE IS IN THE PUBLIC DOMAIN
C       NO RESTRICTIONS ON ITS USE ARE IMPLIED
C----------------------------------------------------------------------
C
C--global variables:
C
        INTEGER                 SIZE
        INTEGER                 N
        INTEGER                 M
        INTEGER                 IRANK
        INTEGER                 IERR
C
        DOUBLE PRECISION        A(SIZE,M)
        DOUBLE PRECISION        AINV(SIZE,N)
        DOUBLE PRECISION        S(MINNM)
        DOUBLE PRECISION        E(MINNM)
        DOUBLE PRECISION        U(SIZE,N)
        DOUBLE PRECISION        V(SIZE,M)
        DOUBLE PRECISION        WORK(N)
C
C--local variables:
C
        INTEGER                 I
        INTEGER                 J
C
        DOUBLE PRECISION        EPS
C
C--code:
C
C  Call LINPACK routine to compute SVD.
C
		print *, 'in mpinv, size, n, m = ',size,n, m
		i21 = 21
		CALL DSVDC (A,SIZE,N,M,MINNM,S,E,U,SIZE,V,SIZE,WORK,i21,IERR)
        IF (IERR .NE. 0) RETURN

C        subroutine dsvdc(x, ldx,n,p,MINNM, s,e,u, ldu, v, ldv,work, job,info)
C        integer ldx,n,p,ldu,ldv,job,info, MINNM
C        double precision x(ldx,p),s(MINNM),e(MINNM),u(ldu,n),v(ldv,p),work(n)


C
C  Compute machine epsilon.
C
        EPS = 1.0D0
 10     IF (EPS + 1.D0 .NE. 1.D0) THEN
          EPS = EPS / 2.D0
          GO TO 10
        END IF
        EPS = EPS * 2.D0
C
        EPS = EPS * S(1)
C
        DO 20, J = 1, N
          DO 20, I = 1, M
            AINV(I,J) = 0.D0
 20     CONTINUE
C
        IRANK = 1
 30     IF (IRANK .LE. MIN(N,M)) THEN
 			IF (S(IRANK) .GT. EPS) THEN 
          		DO 41, J = 1, N
            		DO 40, I = 1, M
              		  	AINV(I,J) = AINV(I,J) +
     1          		V(I,IRANK)*U(J,IRANK)/S(IRANK)
 40       			CONTINUE
 41				CONTINUE
          		IRANK = IRANK + 1
          	  	GO TO 30
			END IF
        END IF
C
        IRANK = IRANK - 1
C
        RETURN
        END


