      SUBROUTINE READBOU (IU09, IU10, IU06)

! ----------------------------------------------------------------------

!**** *READBOU* - READ COARSE AND FINE GRID BOUNDARY INFORMATION FILE.

!     R. PORTZ     MPI          JANUARY 1991

!*    PURPOSE.
!     --------

!       READ COMMON CBOUND AND FBOUND AS WRITTEN BY PREPROC.

!**   INTERFACE.
!     ----------

!       *CALL* *READBOU (IU09, IU10, IU06)*
!          *IU09*    - INPUT  UNIT OF COMMON CBOUND.
!          *IU10*    - INPUT  UNIT OF COMMON FBOUND.
!          *IU06*    - PRINTER OUTPUT UNIT.

!     METHOD.
!     -------

!       SEQUENCIAL UNFORMATED WRITE TO UNIT.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCBOU  , ONLY : IBOUNC   ,NBOUNC   ,IJARC    ,IGARC
      USE YOWFBOU  , ONLY : IBOUNF   ,NBOUNF   ,IJARF    ,IGARF    ,
     &            IBFL     ,IBFR     ,BFW

! ----------------------------------------------------------------------

      IF (IBOUNC.EQ.1) THEN

!     READ COMMON BOUNC

        READ (IU09,ERR=2000) NBOUNC

        ALLOCATE(IGARC(NBOUNC))
        ALLOCATE(IJARC(NBOUNC))

        READ (IU09,ERR=2000) (IGARC(I),I=1,NBOUNC),
     &   (IJARC(I),I=1,NBOUNC)

      ENDIF

      IF (IBOUNF.EQ.1) THEN

!     READ COMMON BOUNF

        READ (IU10,ERR=2000) NBOUNF

        ALLOCATE(IGARF(NBOUNF))
        ALLOCATE(IJARF(NBOUNF))
        ALLOCATE(IBFL(NBOUNF))
        ALLOCATE(IBFR(NBOUNF))
        ALLOCATE(BFW(NBOUNF))

        READ (IU10,ERR=2000) (IGARF(I),I=1,NBOUNF),
     &   (IJARF(I),I=1,NBOUNF),
     &   (IBFL(I),I=1,NBOUNF),
     &   (IBFR(I),I=1,NBOUNF),
     &   (BFW(I),I=1,NBOUNF)

      ENDIF
      RETURN
 2000 CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*    FATAL ERROR IN SUB. READBOU       *'
      WRITE(IU06,*) '*    ===========================       *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*        READ ERROR                    *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*        PROGRAM ABORTS                *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1
      END SUBROUTINE READBOU
