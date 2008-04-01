      SUBROUTINE PROPAGS2 (F1, F3, MIJS, MIJL, IG)

! ----------------------------------------------------------------------

!**** *PROPAGS2* - COMPUTATION OF A PROPAGATION TIME STEP.

! ----------------------------------------------------------------------

!**** *PROPAGS2* -  ADVECTION USING THE CONER TRANSPORT SCHEME IN SPACE


!*    PURPOSE.
!     --------

!       COMPUTATION OF A PROPAGATION TIME STEP.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPAGS2(F1, F3, MIJS, MIJL, IG)*
!          *F1*   - SPECTRUM AT TIME T.
!          *F3*   - SPECTRUM AT TIME T+DELT.
!          *MIJS* - INDEX OF FIRST POINT
!          *MIJL* - INDEX OF LAST POINT
!          *IG*   - BLOCK NUMBER.


!     METHOD.
!     -------

!!  need text !!!!!!!!


!     EXTERNALS.
!     ----------


!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : ICASE    ,ISHALLO  ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLATN    ,KLONN    ,SUMWN    ,WLATN    ,
     &            WLONN    ,WKPMN

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F1,F3
      INTEGER :: MIJS,MIJL,IG
      INTEGER :: K,M,IJ
      INTEGER :: KP1,KM1
      REAL :: SUMW 

! ----------------------------------------------------------------------


!*    SPHERICAL OR CARTESIAN GRID?
!     ----------------------------
      IF (ICASE.EQ.1) THEN

!*      SPHERICAL GRID.
!       ---------------

        IF (IREFRA.NE.2 .AND. IREFRA.NE.3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------

          DO K=1,NANG
            KP1 = K+1
            IF (KP1.GT.NANG) KP1 = 1
            KM1 = K-1
            IF (KM1.LT.1) KM1 = NANG

            DO M=1,NFRE

              DO IJ=MIJS,MIJL
                F3(IJ,K,M) = 
     &              (1.-SUMWN(IJ,K,M))* F1(IJ            ,K  ,M)
     &         +      WLATN(IJ,K,M,0) * F1(KLATN(IJ,K,0) ,K  ,M)
     &         +      WLATN(IJ,K,M,1) * F1(KLATN(IJ,K,1) ,K  ,M)
     &         +      WLATN(IJ,K,M,-1)* F1(KLATN(IJ,K,-1),K  ,M)
     &         +      WLONN(IJ,K,M)   * F1(KLONN(IJ,K)   ,K  ,M)
     &         +      WKPMN(IJ,K,M,1) * F1(IJ            ,KP1,M)
     &         +      WKPMN(IJ,K,M,-1)* F1(IJ            ,KM1,M)

              ENDDO

            ENDDO

          ENDDO

        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       -----------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* PROPGAS2:                              *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR DEPTH OR/AND CURRENT REFRACTION !  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1

        ENDIF

      ELSE
!*    CARTESIAN GRID.
!     ---------------
        IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* PROPAGS2:                              *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       ----------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* PROPAGS2:                              *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '* FOR DEPTH OR/AND CURRENT REFRACTION !  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ENDIF

      ENDIF

! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE PROPAGS2 
