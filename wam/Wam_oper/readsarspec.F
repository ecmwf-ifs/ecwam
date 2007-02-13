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
      USE YOWICE   , ONLY : FLMINFR
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPCONS , ONLY : ZPI     ,RAD
      USE YOWSARAS , ONLY : NSPEC   ,SPEC     ,LONG     ,LAT      ,
     &            XKLOW    ,XKUP    ,XKCUTOFF ,NSARINVI ,NSARINVR ,
     &            NSAR     ,ISARINV ,RSARINV
      USE YOWSTAT  , ONLY : CDTPRO  ,LSARINV
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: INUNIT, IFL300, INOSPEC
      INTEGER :: ISPEC, IANG, IFRE, IC, IDOUBLE, ISAR, KCOUNT, ISS
      INTEGER :: INTDUMMY,JQUAL, IFAIL

      REAL :: TRACK, ZCOST, USTAR, THW, ZCORWA, ZCORFG, ZCORBE, SCLWSI
      REAL :: XLLOW, XLUP, XHSW, XHSS
      REAL :: DUMMY

      CHARACTER* 3 :: YOFID
      CHARACTER*12 :: CDT

      LOGICAL :: LLOPENED, LNODATA, LLTEST

      DATA YOFID /'INV'/
 
!------------------------------------------------------------------
 
      NSPEC = 1000

!     note that in parallel run spec, long and lat will be allocated
!     for the other PE's in mpbcastsarin.
      ALLOCATE(SPEC(0:NSPEC,NANG,NFRE,2))
      ALLOCATE(LONG(NSPEC,2))
      ALLOCATE(LAT(NSPEC,2))
      ALLOCATE(XKLOW(NSPEC))
      ALLOCATE(XKUP(NSPEC))
      ALLOCATE(XKCUTOFF(NSPEC))

      ISPEC = 1
      INOSPEC = 0
      IDOUBLE = 0
      DO IC=1,2
        DO IANG=1,NANG
          DO IFRE=1,NFRE
            DO ISS=0,NSPEC
            SPEC(ISS,IANG,IFRE,IC) = 0.25*FLMINFR(1,IFRE) 
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      LONG=0.
      LAT=0.
      XKLOW=0.
      XKUP=0.
      XKCUTOFF=0.
      LNODATA=.FALSE.

      IF(.NOT.LSARINV) THEN

!       OPEN INPUT FILE UNIT (if data are not inverted within WAM)

        UNIT: DO INUNIT=34,99,1
          INQUIRE ( UNIT=INUNIT, OPENED=LLOPENED)
          IF ( .NOT. LLOPENED) THEN
            WRITE(IU06,*)
     &     '  READSARSPEC: TAKES UNIT ',INUNIT,' AS SAR DATA INPUT UNIT'
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
      ENDIF


!     READ INPUT FROM SAR INVERSION PROGRAM UNTIL END OF FILE IS 
!     REACHED.

      IF(.NOT.LSARINV) THEN
        READ(INUNIT,END=1) CDT 
        NSAR=99999999
      ENDIF

      DO ISAR = 1, NSAR 
        IF(LSARINV) THEN
          JQUAL = ISARINV(1,ISAR) 
          IFL300 = ISARINV(2,ISAR)

          LONG(ISPEC,2) = RSARINV(1,ISAR)
          LAT(ISPEC,2) =  RSARINV(2,ISAR)
          ZCOST = RSARINV(3,ISAR)
          ZCORWA = RSARINV(4,ISAR)
          ZCORFG = RSARINV(5,ISAR)
          ZCORBE = RSARINV(6,ISAR)
          SCLWSI = RSARINV(7,ISAR)
          XLLOW  = RSARINV(8,ISAR)
          XLUP   = RSARINV(9,ISAR)
          XHSW   = RSARINV(10,ISAR)
          XHSS   = RSARINV(11,ISAR)

          KCOUNT=0
          DO IANG=1,NANG
            DO IFRE=1,NFRE
              KCOUNT=KCOUNT+1
              SPEC(ISPEC,IANG,IFRE,2) = RSARINV(NSARINVR+KCOUNT,ISAR)
            ENDDO
          ENDDO

        ELSE
          READ(INUNIT,END=1) DUMMY 
          READ(INUNIT) XHSW
          READ(INUNIT) LONG(ISPEC,2),LAT(ISPEC,2)
          READ(INUNIT) XHSS
          READ(INUNIT) DUMMY
          READ(INUNIT) DUMMY
          READ(INUNIT) DUMMY
          READ(INUNIT) TRACK,ZCOST,USTAR,THW,ZCORWA,ZCORFG,ZCORBE
          READ(INUNIT) INTDUMMY, INTDUMMY,JQUAL,INTDUMMY,IFL300
          READ(INUNIT) DUMMY
          READ(INUNIT)
     &        ((SPEC(ISPEC,IANG,IFRE,2),IANG=1,NANG),IFRE=1,NFRE)
          READ(INUNIT) DUMMY 
!!!the SAR boundary wave numbers are not yet part of the file input !!! 
          XLLOW  = 1./658.
          XLUP   = 1./65.8
!!! the cut off wave length is not yet part of the file input !!! 
          SCLWSI = 100. 
        ENDIF

!       DO NOT KEEP SPECTRA THAT WERE FLAGGED DURING THE INVERSION
!       OR WHICH MIGHT BE DOUBTFUL
!       OR WHICH DEVIATE TOO MUCH FROM FIRST GUESS (BASIC BACKGROUND
!          CHECK)
!       OR FOR WHICH THE CUT OFF WAVE LENGTH COULD NOT BE DETERMINED.
!       -------------------------------------------------------------
        IF(IFL300.LT.10 .AND. ZCOST.LE.0.5 .AND.
     &     (JQUAL.NE.3 .OR. (JQUAL.EQ.3 .AND. ZCOST.LT.0.1)) .AND.
     &     ZCORFG.GE.0.5 .AND. ZCORBE.GE.0.8 .AND. ZCORWA.LE.0.99 .AND.
     &     ABS(XHSW-XHSS).LT.1.0 .AND.
     &     JQUAL.NE. 7 ) THEN
          LONG(ISPEC,2) = MOD(LONG(ISPEC,2),360.)
          IF( LONG(ISPEC,2).LT.0 ) LONG(ISPEC,2) = LONG(ISPEC,2) + 360.

          XKLOW(ISPEC) = XLLOW
          XKUP(ISPEC) = XLUP
!!!
!!!          WRITE(IU06,'(A21,1x,2(1X,F8.3),1X,I4)') 
!!!     &               'CUT_OFF_WAVE_LENGTH= ',
!!!     &               LONG(ISPEC,2),LAT(ISPEC,2),
!!!     &               NINT(SCLWSI)
!!!          CALL FLUSH(IU06)
!!!
          XKCUTOFF(ISPEC) = ZPI/SCLWSI

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

1     CONTINUE 
      IF( ISPEC.GT.1 ) THEN
        CALL RESIZE_SARINPUT(ISPEC-1)

        WRITE(IU06,*) '  NUMBER OF VALID SPECTRA READ IN , NSPEC = ',
     &  NSPEC
        WRITE(IU06,*) '  NUMBER OF DISCARDED FLAGGED SPECTRA  = ',
     &  INOSPEC
        IF(IDOUBLE.GT.0) 
     &    WRITE(IU06,*) '  NUMBER OF DOUBLE SPECTRA  = ',IDOUBLE

!       SET MINIMUM SAR SPECTRAL VALUE TO THE SAME VALUE AS IN WAM

        DO IFRE = 1,NFRE
          DO IANG = 1,NANG
            DO ISPEC = 1,NSPEC
             SPEC(ISPEC,IANG,IFRE,2)=
     &            MAX(SPEC(ISPEC,IANG,IFRE,2),0.25*FLMINFR(1,IFRE)) 
            END DO
          END DO
        END DO
      ELSE
        LNODATA=.TRUE.
        WRITE(IU06,*) '  NO VALID DATA FOR SAR ASSIMILATION !!!'
        WRITE(IU06,*) '  NOTHING WILL BE DONE. '
      ENDIF

      IF(ALLOCATED(ISARINV)) DEALLOCATE(ISARINV)
      IF(ALLOCATED(RSARINV)) DEALLOCATE(RSARINV)

      IF(.NOT.LSARINV) CLOSE(INUNIT)

5000  CONTINUE

      RETURN

      END SUBROUTINE READSARSPEC
