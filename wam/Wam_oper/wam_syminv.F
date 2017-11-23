      SUBROUTINE WAM_SYMINV(A,NDIM,N,COND)

!**** *WAM_SYMINV*  - SYMMETRIC POSITIVE MATRIX INVERSION.

!     J. BIDLOT    ECMWF: COPY SYMINV FROM ECLIB INTO WAM CODE.   


!     PURPOSE.                                                          
!     --------

!     INVERTS IN PLACE THE LOWER TRIANGLE OF A (I.E. A(I,J) I.GE.J)
!     A IS A SYMMETRIC POSITIVE DEFINITE MATRIX
!     THE UPPER TRIANGLE OF A (IE. A(I,J) I.LT.J) IS NOT USED OR ALTERED


!**   INTERFACE.                                                        
!     ----------                                                        

!     *CALL* WAM_SYMINV(A ,NDIM ,N ,COND)
     
!     A  REAL  - ARRAY CONTAINING THE LOWER TRIANGLE OF MATRIX INCLUDING
!                THE DIAGONAL.
!                ON OUTPUT, A IS OVERWRITTEN BY THE INVERSE LOWER TRIANGLE
!                MATRIX.
     
!     NDIM INTEGER - FIRST DIMENSION OF A
     
!     N INTEGER - DIMENSION OF MATRIX A
     
!     COND - ESTIMATED CONDITION NUMBER
!            IF A IS NOT POSITIVE-DEFINITE -
!                 COND=-FLOAT(I)
!            WHERE I IS THE RANK OF THE FIRST PRINCIPAL MINOR OF A THAT IS
!            NOT POSITIVE-DEFINITE.
!            if A is positive definite     -
!                 COND=1/(N*N*MAX(A(I,J))*MAX(AINVERTED(I,J)))
!            I.E. IN THE INFINITY NORM, AN UPPER LIMIT ON THE INVERSE OF
!                 NORM(A)*NORM(AINVERTED)
!            COND  IS SMALL FOR ILL-CONDITIONED A.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: NDIM, N, NDIA, I, J, K, NM1, IMX 
      REAL :: COND, EPS, ZMX, X, S
      REAL :: A(NDIM,NDIM)
      REAL :: V(NDIM)

      IF(N.LT.1) RETURN

      EPS=0.
      NDIA=NDIM+1
      IMX=ISAMAX_WAM(N,A(1,1),NDIA)
      ZMX=A(IMX,IMX)*FLOAT(N*N)

!**** 1. ROOT FREE CHOLESKY DECOMPOSITION  A =L D L(TRANSPOSE)
      J=1
      IF(A(1,1).LE.EPS) GOTO 91
      V(1)=1./A(1,1)
      IF(N.EQ.1) GOTO 20
      X=A(2,1)*V(1)
      A(2,2)=A(2,2)-X*A(2,1)
      A(2,1)=X
      IF(A(2,2).LE.EPS) GOTO 91
      V(2)=1./A(2,2)
      IF(N.EQ.2) GOTO 20

      DO I=3,N
        DO J=3,I
          S=A(I,J-1)
          DO K=3,J
            S=S-A(I,K-2)*A(J-1,K-2)
          ENDDO
          A(I,J-1)=S
        ENDDO
        S=A(I,I)
        DO J=2,I
          X=A(I,J-1)*V(J-1)
          S=S-X*A(I,J-1)
          A(I,J-1)=X
        ENDDO
        A(I,I)=S
!        CHECK FOR POSITIVE-DEFINITENESS AND INVERT DIAGONAL MATRIX D.
        IF(A(I,I).LE.EPS) GOTO 91
        V(I)=1./A(I,I)
      ENDDO

!**** 2.  COPY INVERSE OF D WHICH HAS ALREADY BEEN CALCULATED.
20    CONTINUE
      DO J=1,N
        A(J,J)=V(J)
      ENDDO
      IF(N.EQ.1) GOTO 50
!
!**** 3.  INVERSION OF L
      A(2,1)=-A(2,1)
      NM1=N-1
      IF(N.EQ.2) GOTO 40
      DO I=2,NM1
        DO J=2,I
          S=A(I+1,J-1)
          DO K=J,I
            S=S+A(I+1,K)*A(K,J-1)
          ENDDO
          A(I+1,J-1)=-S
        ENDDO
        A(I+1,I)=-A(I+1,I)
      ENDDO

!**** 4.  INV A = INV L(TRANSPOSE) * INV D * INV L
40    CONTINUE
      DO J=2,N
        S=A(J-1,J-1)
        DO I=J,N
          X=A(I,I)*A(I,J-1)
          S=S+A(I,J-1)*X
          A(I,J-1)=X
        ENDDO
        A(J-1,J-1)=S
        IF(J.EQ.N) GOTO 50
        DO I=J,NM1
          S=A(I,J-1)
          DO K=I,NM1
            S=S+A(K+1,I)*A(K+1,J-1)
          ENDDO
          A(I,J-1)=S
        ENDDO
      ENDDO

50    CONTINUE
      IMX=ISAMAX_WAM(N,A(1,1),NDIA)
      COND=1./ABS(A(IMX,IMX)*ZMX)
      RETURN

91    CONTINUE
      COND=-FLOAT(J)
      RETURN

      CONTAINS

! ----------------------------------------------------------------------
      INTEGER FUNCTION ISAMAX_WAM(N,A,M)

!     FIND THE LARGEST ABSOLUTE ELEMENT OF A , SPACED M WORDS APART

      IMPLICIT NONE


      INTEGER :: LARGE, N, M, I, INDEX
      REAL :: A(*)

      LARGE=1
      IF(N.LE.1) GOTO 9
      INDEX=1+M
      DO I=2,N
        IF(ABS(A(INDEX)).GE.ABS(A(LARGE))) LARGE=I
        INDEX=INDEX+M
      ENDDO
9     CONTINUE
      ISAMAX_WAM=LARGE
      END FUNCTION ISAMAX_WAM
      END SUBROUTINE WAM_SYMINV
