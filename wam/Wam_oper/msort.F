!--------------------------------------------------------------------

      SUBROUTINE MSORT (IU06, NUM, KICE, IG, IJ)

!-----------------------------------------------------------------

!**** *MSORT* - SORT ICE POINTS

!     R. PORTZ         MPI HAMBURG       OCTOBER 1992


!     PURPOSE.
!     --------
!             SORT ICE-POINTS ACCORDING TO BLOCK AND INDEX
!             IN A BLOCK AND WRITE INFO FILE FOR WAMODEL (TAPE 8)


!**   INTERFACE.
!     ---------
!             *CALL MSORT*

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   PRINTER OUTPUT UNIT.
!      *NUM*       INTEGER   NUMBER OF ICE POINTS IN
!                            MODEL GRID.
!      *NICE*      INTEGER   NUMBER OF ICE POINTS
!      *IG*        INTEGER   BLOCK NUMBER (ICE POINT)
!      *IJ*        INTEGER   NUMBER OF INDEX IN A BLOCK (ICE POINT)

!     METHOD.
!     -------
!             NONE

!    EXTERNAL.
!    ---------

!---------------------------------------------------------------------

!     THIS IS NOW A DUMMY ROUTINE
      RETURN
      END SUBROUTINE MSORT




