      SUBROUTINE READSARSPEC(LNODATA)

!------------------------------------------------------------------

!     JEAN BIDLOT         ECMWF     MAY 1999.

!     PURPOSE:
!     --------

!     READS SPECTRA IN THE FORMAT OF THE SAR INVERSION PROGRAM.
!     INPUT DATA SHOULD BE IN FILE SFLYYYYMMDDHHmm, where 
!     SFLYYYYMMDDHHmm IS GIVEN BY CDTPRO.
!     FOR EFFICIENCY IT ONLY READS WHAT IT NEEDS.

!     INTERFACE:
!     ----------

!       *CALL* *READSARSPEC*

! ----------------------------------------------------------------------
!     MODULES:
!     --------
      USE YOWICE   , ONLY : FLMIN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : RAD
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWSARAS , ONLY : NSPEC    ,SPEC     ,LONG    ,LAT
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: INUNIT, IFL300, INOSPEC
      INTEGER :: ISPEC, IANG, IFRE, IC, IDOUBLE
      INTEGER :: IANGP1, IANGM1, IFREP1, IFREM1
      INTEGER :: INTDUMMY,JQUAL, IFAIL

      REAL :: TRACK, COST
      REAL :: DUMMY

      CHARACTER* 3 :: YOFID
      CHARACTER*12 :: CDT

      LOGICAL :: ZERO2, LLOPENED, LNODATA, LLTEST

      DATA YOFID /'INV'/
 
!------------------------------------------------------------------
 
      NSPEC = 500

!     note that in parallel run spec, long and lat will be allocated
!     for the other PE's in mpbcastsarin.
      ALLOCATE(SPEC(0:NSPEC,NANG,NFRE,2))
      ALLOCATE(LONG(NSPEC,2))
      ALLOCATE(LAT(NSPEC,2))

      ISPEC = 1
      INOSPEC = 0
      IDOUBLE = 0
      SPEC=FLMIN
      LONG=0.
      LAT=0.
      LNODATA=.FALSE.

!     OPEN INPUT FILE UNIT

      UNIT: DO INUNIT=34,99,1
        INQUIRE ( UNIT=INUNIT, OPENED=LLOPENED)
        IF ( .NOT. LLOPENED) THEN
          WRITE(IU06,*)
     &     ' READSARSPEC : TAKES UNIT ',INUNIT,' AS SAR DATA INPUT UNIT'
          EXIT UNIT 
        ENDIF
      ENDDO UNIT
      IF(INUNIT.GT.99) THEN
        WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) ' +                                        +'
        WRITE(IU06,*) ' + ERROR IN --READSARSPEC--               +'
        WRITE(IU06,*) ' + =========================              +'
        WRITE(IU06,*) ' + RUN OUT OF UNIT TO OPEN SAR DATA INPUT +' 
        WRITE(IU06,*) ' + PROGRAM WILL ABORT                     +' 
        WRITE(IU06,*) ' +                                        +'
        WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        CALL ABORT1
      ENDIF

      LLTEST=.TRUE.
      CALL CONFILE (IU06, INUNIT, CDTPRO, YOFID, IFAIL, LLTEST)
      IF(IFAIL.NE.0) THEN
        WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*) ' +                                     +'
        WRITE(IU06,*) ' + WARNING ERROR IN --READSARSPEC--    +'
        WRITE(IU06,*) ' + ===============================     +'
        WRITE(IU06,*) ' +                                     +'
        WRITE(IU06,*) ' +      NO SAR MEASUREMENTS            +'
        WRITE(IU06,*) ' + MODEL TIME IS  CDTPRO = ', CDTPRO
        WRITE(IU06,*) ' + IFAIL = ', IFAIL
        WRITE(IU06,*) ' + PROCESSING WILL BE CONTINUED        +'
        WRITE(IU06,*) ' +                                     +'
        WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++'
        LNODATA=.TRUE.
        CALL FLUSH(IU06)
        GOTO 5000
      ENDIF


!     READ INPUT FROM SAR INVERSION PROGRAM UNTIL END OF FILE IS 
!     REACHED.

      READ(INUNIT,END=1) CDT
      DO WHILE( .TRUE. )
        READ(INUNIT,END=1) DUMMY 
        READ(INUNIT) DUMMY
        READ(INUNIT) LONG(ISPEC,2),LAT(ISPEC,2)
        READ(INUNIT) DUMMY
        READ(INUNIT) DUMMY
        READ(INUNIT) DUMMY
        READ(INUNIT) DUMMY
        READ(INUNIT) DUMMY
        READ(INUNIT) INTDUMMY, INTDUMMY,JQUAL,INTDUMMY,IFL300
        READ(INUNIT) DUMMY
        READ(INUNIT)
     &      ((SPEC(ISPEC,IANG,IFRE,2),IANG=1,NANG),IFRE=1,NFRE)
        READ(INUNIT) DUMMY 

!       DO NOT KEEP SPECTRA THAT WERE FLAGGED DURING THE INVERSION.
!       -----------------------------------------------------------
        IF(IFL300.LT.10) THEN
          LONG(ISPEC,2) = MOD(LONG(ISPEC,2),360.)
          IF( LONG(ISPEC,2).LT.0 ) LONG(ISPEC,2) = LONG(ISPEC,2) + 360.

!         CHECK FOR DOUBLE OBSERVATIONS
          DO IC=1,ISPEC-1
            IF(LONG(IC,2).EQ.LONG(ISPEC,2).AND.
     &         LAT(IC,2).EQ.LAT(ISPEC,2)) THEN
               IDOUBLE=IDOUBLE+1
               ISPEC = ISPEC-1
               EXIT
            ENDIF
          ENDDO
          ISPEC = ISPEC+1
        ELSE
          INOSPEC = INOSPEC + 1
        ENDIF

        IF( ISPEC.GT.NSPEC ) CALL RESIZE_SARINPUT(ISPEC)

      END DO

1     CALL RESIZE_SARINPUT(ISPEC-1)

      WRITE(IU06,*) 'NUMBER OF VALID SPECTRA READ IN , NSPEC = ',NSPEC
      WRITE(IU06,*) 'NUMBER OF DISCARDED FLAGGED SPECTRA  = ',INOSPEC
      IF(IDOUBLE.GT.0) 
     & WRITE(IU06,*) 'NUMBER OF DOUBLE SPECTRA  = ',IDOUBLE

!     SET MINIMUM SAR SPECTRAL VALUE TO THE SAME VALUE AS IN WAM
!     (ie FLMIN)

      DO IFRE = 1,NFRE
        DO IANG = 1,NANG
          DO ISPEC = 1,NSPEC
            SPEC(ISPEC,IANG,IFRE,2) = MAX(SPEC(ISPEC,IANG,IFRE,2),FLMIN) 
          END DO
        END DO
      END DO

!     SET CUT OFF VALUES CAUSED BY CODING AND DECODING TO ZERO.

      DO IANG = 1,NANG
        IF( IANG.LT.NANG )THEN
          IANGP1 = IANG + 1
        ELSE
          IANGP1 = 1
        END IF
        IF( IANG.GT.1 )THEN
          IANGM1 = IANG - 1
        ELSE
          IANGM1 = NANG
        END IF
        DO IFRE = 1,NFRE
          IFREP1 = MIN(IFRE+1,NFRE)         
          IFREM1 = MAX(IFRE-1,1)
          DO ISPEC = 1,NSPEC
            IF( SPEC(ISPEC,IANGP1,IFRE,2).LE.0.1E-2.AND.
     &          SPEC(ISPEC,IANGM1,IFRE,2).LE.0.1E-2.AND.
     &          SPEC(ISPEC,IANG,IFREP1,2).LE.0.1E-2.AND.
     &          SPEC(ISPEC,IANG,IFREM1,2).LE.0.1E-2
     &          ) SPEC(ISPEC,IANG,IFRE,2) = FLMIN 
          END DO
        END DO
      END DO


!     SMOOTH SAR SPECTRA AROUND 100M WAVE LENGTH CUT OFF.

      DO ISPEC = 1,NSPEC
        ZERO2 = .TRUE.
        IANG = 1
        DO WHILE( IANG.LE.NANG.AND.ZERO2 )
          IFRE = 1
          DO WHILE( IFRE.LE.NFRE.AND.ZERO2 )
            IF( SPEC(ISPEC,IANG,IFRE,2).NE.FLMIN ) ZERO2 = .FALSE.
            IFRE = IFRE + 1
          END DO
          IANG = IANG + 1
        END DO
        IF( ZERO2 ) WRITE(IU06,*) 'WARNING!',
     &    ' INVERTED SAR SPECTRUM NUMBER ',ISPEC,' IS ZERO!'
      END DO

      CLOSE(INUNIT)

5000  CONTINUE

      RETURN

      END SUBROUTINE READSARSPEC
