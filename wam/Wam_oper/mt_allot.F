      SUBROUTINE MT_ALLOT (F,ETOI, USOI, USMO, THMO,IJS,IJL,NSIZE,NWORK)

! ----------------------------------------------------------------------

!**** *MT_ALLOT* - SET UP TASKS TO RUN UPDATE.

!     B.       HANSEN  MPIfM/E.C.M.W.F.    OCTOBER   1994
!     (  FROM ROUTINE PHYS BY
!        P.A.E.M. JANSSEN K.N.M.I./E.C.M.W.F. SEPTEMBER 1994 )

!*    PURPOSE.
!     --------

!       BY SPLITTING UP INNER DO LOOP IN A NUMBER OF PARTS WE DISTRIBUTE
!       THE WORK OVER A NUMBER OF PROCESSORS.

!**   INTERFACE.
!     ----------

!       *CALL* *MT_ALLOT (F,ETOI, USOI, USMO, THMO, IJS,IJL,IG,NSIZE,
!                         NWORK)*
!          *F*      -  REAL      SPECTRUM ( FIRST GUESS IN INPUT ;
!                                ANALYSIS IN OUTPUT ).
!          *ETOI*   -  REAL      WAVE ENERGY (FROM O.I.).
!          *USOI*   -  REAL      USTAR ( FROM O.I.).
!          *USMO*   -  REAL      USTAR ( FIRST GUESS IN INPUT ;
!                                ANALYSIS IN OUTPUT ).
!          *THMO*   -  REAL      WIND DIRECTION.
!          *IJS*    -  INDEX OF FIRST GRIDPOINT
!          *IJL*    -  INDEX OF LAST GRIDPOINT
!          *IG*     -  BLOCK NUMBER
!          *NSIZE*  -  SIZE OF SUB TASK
!          *NWORK*  -  TOTAL AMOUNT OF SUB TASKS

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ---------

!       *UPDATE*    - ANALYSES THE WAM MODEL SPECTRA.

!     REFERENCE.
!     ----------

!       NONE.                                                     

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.


      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWMLTSK , ONLY : NCHUNK
      USE YOWMPP   , ONLY : NINF     ,NSUP



! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F
      REAL,DIMENSION(NINF:NSUP) :: ETOI,USOI,USMO,THMO

! ----------------------------------------------------------------------


!*    1. INITIALISATION.
!        ---------------

 1000 CONTINUE  
      NN = INCINT(NCHUNK,1)
      IF (NN.LE.NWORK) THEN
        MIJS = IJS + (NN-1)*NSIZE
        MIJL = MIJS + NSIZE-1
        MIJS = MIN(MIJS,IJL)
        MIJL = MIN(MIJL,IJL)

        CALL UPDATE (F, ETOI, USOI, USMO, THMO, MIJS, MIJL)
        GO TO 1000
      ENDIF

! ----------------------------------------------------------------------

      RETURN
      END SUBROUTINE MT_ALLOT
