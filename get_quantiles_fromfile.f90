SUBROUTINE get_quantiles_fromfile(infile_quantiles, nxain, nyain, &
	thresh_median, thresh_uppertercile, thresh_upperquintile, &
	thresh_upperdecile, thresh_95)
!  get_quantiles_fromfile.f90
!  compile with > f2py -c -m get_quantiles_fromfile get_quantiles_fromfile.f90
!  assuming you have got f2py working on your system (it's a fortran compiler
!  made to interface with python; google f2py on the web for more).

CHARACTER*(*), INTENT(IN) :: infile_quantiles
INTEGER, INTENT(IN) :: nxain, nyain
REAL, INTENT(OUT), DIMENSION(nxain, nyain) :: thresh_median
REAL, INTENT(OUT), DIMENSION(nxain, nyain) :: thresh_uppertercile
REAL, INTENT(OUT), DIMENSION(nxain, nyain) :: thresh_upperquintile
REAL, INTENT(OUT), DIMENSION(nxain, nyain) :: thresh_upperdecile
REAL, INTENT(OUT), DIMENSION(nxain, nyain) :: thresh_95

! directives for python:
! f2py intent(in) infile_quantiles, nxain, nyain
! f2py intent(out) thresh_median, thresh_uppertercile 
! f2py intent(out) thresh_upperquintile
! f2py intent(out) thresh_upperdecile, thresh_95
! f2py depend(nxain, nyain) thresh_median, thresh_uppertercile 
! f2py depend(nxain, nyain) thresh_upperquintile
! f2py depend(nxain, nyain) thresh_upperdecile, thresh_95

PRINT *,'Fortran subroutine get_quantiles_fromfile'
PRINT *,'reading from ',TRIM(infile_quantiles)
OPEN (unit=1, file=infile_quantiles, status='old',form='unformatted')
READ (1) nxa, nya
READ (1) thresh_median
READ (1) thresh_uppertercile
READ (1) thresh_upperquintile
READ (1) thresh_upperdecile
READ (1) thresh_95
CLOSE (1) 

END SUBROUTINE get_quantiles_fromfile

