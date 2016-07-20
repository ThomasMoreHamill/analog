Subroutine R_inspar (XDONT, iX2, nsize, NORD)
!  Sorts partially XDONT, bringing the NORD lowest values at the
!  begining of the array
!  Also sorts iX2 to the same order as XDONT
! __________________________________________________________
!  This subroutine uses insertion sort, limiting insertion
!  to the first NORD values. It does not use any work array
!  and is faster when NORD is very small (2-5), but worst case
!  behavior can happen fairly probably (initially inverse sorted)
!  In many cases, the refined quicksort method is faster.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! _________________________________________________________
      Integer, Intent (In) :: NORD, nsize
      Real, Dimension (nsize), Intent (InOut) :: XDONT
      Integer, DIMENSION (nsize) :: iX2
! __________________________________________________________
      Real    :: XWRK, XWRK1
      Integer iXWRK2, iXWRK3
!
      Integer :: ICRS, IDCR
!
      Do ICRS = 2, NORD
         !print *,'icrs = ',icrs
         XWRK  = XDONT (ICRS)
         iXWRK2 = iX2 (ICRS)
         Do IDCR = ICRS - 1, 1, - 1
            If (XWRK >= XDONT(IDCR)) Exit
            XDONT (IDCR+1) = XDONT (IDCR)
            iX2 (IDCR+1) = iX2 (IDCR)
         End Do
         XDONT (IDCR+1) = XWRK
         iX2 (IDCR+1)  = iXWRK2
      End Do
!
      XWRK1 = XDONT (NORD)
      iXWRK2 = iX2 (NORD)
      Do ICRS = NORD + 1, nsize
         !print *,'icrs = ',icrs
         If (XDONT(ICRS) < XWRK1) Then
            XWRK  = XDONT (ICRS)
            iXWRK3 = iX2 (ICRS)
            XDONT (ICRS) = XWRK1
            iX2 (ICRS) = iXWRK2
            Do IDCR = NORD - 1, 1, - 1
               If (XWRK >= XDONT(IDCR)) Exit
               XDONT (IDCR+1) = XDONT (IDCR)
               iX2 (IDCR+1) = iX2 (IDCR)
            End Do
            XDONT (IDCR+1) = XWRK
            iX2 (IDCR+1) = iXWRK3
            XWRK1 = XDONT (NORD)
            iXWRK2 = iX2 (NORD)
         End If
      End Do
!
!
End Subroutine R_inspar
