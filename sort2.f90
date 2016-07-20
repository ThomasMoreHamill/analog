SUBROUTINE SORT2 (N,Ra,Rb)
IMPLICIT NONE
!
! Dummy arguments
!
INTEGER :: N
REAL , DIMENSION(N) :: Ra , Rb
INTENT (IN) N
INTENT (INOUT) Ra , Rb 
!
! Local variables
!
INTEGER :: i , ir , j , l
REAL :: rra, rrb
!
!-----------------------------------------------------------------------
!
!     argrument list.
!
!     name      type              description.
!     n         integer           length of arrays ra, rb, rc & rd.
!     (unchanged on exit).
!
!     ra        array of real     array of length n to be sorted into
!     ascending numerical order.
!     (contains result on exit).
!
!     rb        array of real    array of length n to be sorted into
!     an order corresponding to that of ra.
!     (contains result on exit).
!
!-----------------------------------------------------------------------
!
!     looks at array ra and sorts it into ascending order, as well as
!     making the corresponding rearrangements to rb. the
!     heapsort algorithm is used.
!
!     based on numerical recipes, the art of scientific computing,
!     section 8.2, by press, flannery, teukolsky & vetterling,
!     cambridge university press, 1987.
!
!     modified by: david lary
!     ----------
!
!     date started : 28/9/1991
!
!     last modified: 28/9/1991
!
!-----------------------------------------------------------------------
!
l = N/2 + 1
ir = N
DO WHILE ( .TRUE. )
  IF ( l>1 ) THEN
    l = l - 1
    rra = Ra(l)
    rrb = Rb(l)
  ELSE
    rra = Ra(ir)
    rrb = Rb(ir)
    Ra(ir) = Ra(1)
    Rb(ir) = Rb(1)
    ir = ir - 1
    IF ( ir<=1 ) THEN
      Ra(1) = rra
      Rb(1) = rrb
      EXIT
    ENDIF
  ENDIF
  i = l
  j = l + l
  DO WHILE ( .TRUE. )
    IF ( .NOT..TRUE. ) THEN
      RETURN
    ELSE IF ( j<=ir ) THEN
      IF ( j<ir ) THEN
        IF ( Ra(j)<Ra(j+1) ) j = j + 1
      ENDIF
      IF ( rra<Ra(j) ) THEN
        Ra(i) = Ra(j)
        Rb(i) = Rb(j)
        i = j
        j = j + j
      ELSE
        j = ir + 1
      ENDIF
      CYCLE
    ENDIF
    Ra(i) = rra
    Rb(i) = rrb
    GOTO 100
  ENDDO
  EXIT
100 ENDDO
END SUBROUTINE SORT2
