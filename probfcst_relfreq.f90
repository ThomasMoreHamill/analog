SUBROUTINE probfcst_relfreq (nsize, thresh,  ensfcst,  prob)

INTEGER, INTENT(IN) :: nsize
REAL, INTENT(IN):: thresh
REAL, INTENT(IN), DIMENSION(nsize) :: ensfcst
REAL, INTENT(OUT) :: prob

prob = 0.
probtot = 0.
DO i = 1, nsize
   IF (ensfcst(i) .GE. thresh) prob = prob+1.
   IF (ensfcst(i) .GE. 0.0)    probtot = probtot+1.
END DO
prob = prob / probtot
RETURN
END subroutine probfcst_relfreq
