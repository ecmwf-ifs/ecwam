      SUBROUTINE OUTGRID (FL3, IJS, IJL, IG,U10NEW,THWNEW,USNEW,TAUW,
     &                    NSTART,NEND)
! ----------------------------------------------------------------------

!**** *OUTGRID* - OUTPUT OF WAVE ENERGY,MEAN ANGLE, MEAN FREQUENCY
!****             FRICTION VELOCITY AND WIND DIRECTION IN GRIDDED
!****             FORMAT.

!     P.JANSSEN      KNMI


!     P.LIONELLO  FEB. 87
!                 OUTPUT OF SWELL ENERGY ,MEAN SWELL DIRECTION ,
!                 MEAN SWELL FREQUENCY AND MEAN WIND-SEA WAVE
!                 DIRECTION AT ALL ACTIVE GRID POINTS

!     H.GUNTHER   GKSS/ECMWF     NOVEMBER 1989
!     J. BIDLOT   ECMWF          APRIL 1996
!                 MESSAGE PASSING

!                 OUTPUT FORMAT IS NOW GRIDDED FIELDS.

!*    PURPOSE.
!     --------

!       OUTPUT OF WAVE AND WIND FIELDS INTO COMMON BLOCK.


!**   INTERFACE.
!     ----------

!       *CALL* *OUTGRID (FL3, IJS, IJL, IG,U10NEW,THWNEW,USNEW,TAUW)*
!         *FL3*    - BLOCK OF SPECTRA.
!         *IJS*    - INDEX OF FIRST GRIDPOINT.
!         *IJL*    - INDEX OF LAST GRIDPOINT.
!         *IG*     - BLOCK NUMBER.
!         *U10NEW*    NEW WIND SPEED IN M/S.
!         *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                     NOTATION (POINTING ANGLE OF WIND VECTOR,
!                     CLOCKWISE FROM NORTH).
!         *USNEW*     NEW FRICTION VELOCITY IN M/S.
!         *TAUW*      WAVE STRESS IN (M/S)**2
!         *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!         *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN


!     EXTERNALS.
!     ----------

!       *MAKEGRID*  - GENERATES GRIDDED FIELDS.
!       *PEAKFR*    - COMPUTATION OF PEAK FREQUENCY.
!       *MPGATHERSCFLD* - GATHER SCALAR FIELD DISTRIBUTED ON DIFFERENT
!                         PROCESSES

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : IPFGTBL
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G
      USE YOWICE   , ONLY : LICERUN  ,ICEMASK
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,SMEAN
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : CDTPRO   ,CDTINTT  ,CDTINTS
      USE YOWSWEL  , ONLY : ESWELL   ,FSWELL   ,THSWELL  ,ESEA     ,
     &            FSEA     ,THWISEA

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL3
      REAL,DIMENSION(NINF:NSUP,NBLO) :: TAUW
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW

! ----------------------------------------------------------------------

      REAL,ALLOCATABLE,DIMENSION(:) :: PEAKF,CD,WVSTR,THW,US,U10 

! ----------------------------------------------------------------------

!*    1. INTEGRATED PARAMETERS OF TOTAL SEA.
!        -----------------------------------

      IF (CDTINTT.EQ.CDTPRO) THEN

!*    1.1 COMPUTE PEAK FREQUENCY.
!         -----------------------

        ALLOCATE (PEAKF(NIBLO))
        CALL PEAKFR (FL3, IJS, IJL, PEAKF)

!*    1.2 COMPUTE NORMALISED WAVESTRESS AND DRAG. 
!         ON ICEPOINTS GIVE MISSING VALUE INDICATOR.
!         ---------------------------------------

        ALLOCATE (WVSTR(NIBLO))
        ALLOCATE (CD(NIBLO))
        DO IJ = IJS,IJL
          TAU = USNEW(IJ)**2 +0.0001
          WVSTR(IJ) = TAUW(IJ,IG)/TAU
!!       if the numerical computation of TAU and CD change, a similar
!!       modification has to be put in buildstress where the friction
!!       velocity is determined from U10 and CD.
          CD(IJ) = TAU/(U10NEW(IJ)**2+0.01)
        ENDDO

        IF (LICERUN) THEN
          DO IJ = IJS,IJL
            IF (ICEMASK(IJ,IG).EQ.0) THEN
              EMEAN(IJ)  = ZMISS
              THQ(IJ)    = ZMISS
              PEAKF(IJ)  = ZMISS
              FMEAN(IJ)  = ZMISS
              SMEAN(IJ)  = ZMISS
              THWNEW(IJ) = ZMISS
              USNEW(IJ)  = ZMISS
              CD(IJ)     = ZMISS
              WVSTR(IJ)  = ZMISS   
            ENDIF   
          ENDDO
        ENDIF

!*    1.3 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
!         -------------------------------------------------

        IF(LMESSPASS) THEN
          CALL MPGATHERSCFLD(IPFGTBL(1),1,NSTART,NEND,EMEAN)
          IF(IPFGTBL(1).EQ.IRANK) 
     &     CALL MAKEGRID (EMEAN , WHGTTG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(2),2,NSTART,NEND,THQ)
          IF(IPFGTBL(2).EQ.IRANK) 
     &     CALL MAKEGRID (THQ , WDIRTG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(3),3,NSTART,NEND,FMEAN)
          IF(IPFGTBL(3).EQ.IRANK) 
     &     CALL MAKEGRID (FMEAN , WMNFTG, 1, NEND(NPROC), IG)

          IF(IPFGTBL(4).NE.0) THEN
            ALLOCATE(US(NIBLO))
            US = 0.
            DO IJ = IJS,IJL
              US(IJ) = USNEW(IJ)
            ENDDO
            CALL MPGATHERSCFLD(IPFGTBL(4),4,NSTART,NEND,US)
            IF(IPFGTBL(4).EQ.IRANK) 
     &        CALL MAKEGRID (US, USTARG, 1, NEND(NPROC), IG)
            DEALLOCATE(US)
          ENDIF

          IF(IPFGTBL(5).NE.0) THEN 
            ALLOCATE(THW(NIBLO))
            THW = 0.
            DO IJ = IJS,IJL
              THW(IJ) = THWNEW(IJ) 
            ENDDO
            CALL MPGATHERSCFLD(IPFGTBL(5),5,NSTART,NEND,THW)
            IF(IPFGTBL(5).EQ.IRANK) 
     &        CALL MAKEGRID (THW, UDIRG, 1, NEND(NPROC), IG)
            DEALLOCATE(THW)
          ENDIF

          CALL MPGATHERSCFLD(IPFGTBL(6),6,NSTART,NEND,PEAKF)
          IF(IPFGTBL(6).EQ.IRANK)
     &     CALL MAKEGRID (PEAKF , WPKFTG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(7),7,NSTART,NEND,CD)
          IF(IPFGTBL(7).EQ.IRANK)
     &     CALL MAKEGRID (CD   ,  CDG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(8),8,NSTART,NEND,WVSTR)
          IF(IPFGTBL(8).EQ.IRANK)
     &     CALL MAKEGRID (WVSTR,  TAUWG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(9),9,NSTART,NEND,SMEAN)
          IF(IPFGTBL(9).EQ.IRANK)
     &     CALL MAKEGRID (SMEAN , SMEANG, 1, NEND(NPROC), IG)

          IF(IPFGTBL(10).NE.0) THEN
            ALLOCATE(U10(NIBLO))
            U10 = 0.
            DO IJ = IJS,IJL
              U10(IJ) = U10NEW(IJ)
              IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) U10(IJ)=ZMISS
            ENDDO
            CALL MPGATHERSCFLD(IPFGTBL(10),10,NSTART,NEND,U10)
            IF(IPFGTBL(10).EQ.IRANK) 
     &        CALL MAKEGRID (U10, U10G, 1, NEND(NPROC), IG)
            DEALLOCATE(U10)
          ENDIF

        ELSE
          CALL MAKEGRID (EMEAN , WHGTTG, IJS, IJL, IG)
          CALL MAKEGRID (THQ   , WDIRTG, IJS, IJL, IG)
          CALL MAKEGRID (PEAKF , WPKFTG, IJS, IJL, IG)
          CALL MAKEGRID (FMEAN , WMNFTG, IJS, IJL, IG)
          CALL MAKEGRID (SMEAN , SMEANG, IJS, IJL, IG)
          CALL MAKEGRID (THWNEW, UDIRG , IJS, IJL, IG)
          CALL MAKEGRID (USNEW,  USTARG, IJS, IJL, IG)
          CALL MAKEGRID (CD   ,  CDG   , IJS, IJL, IG)
          CALL MAKEGRID (WVSTR,  TAUWG , IJS, IJL, IG)
          CALL MAKEGRID (U10NEW, U10G  , IJS, IJL, IG)
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    2. INTEGRATED PARAMETERS OF SWELL.
!        -------------------------------

      IF (CDTINTS.EQ.CDTPRO) THEN

!*    2.1 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
!         -------------------------------------------------

        IF (LICERUN) THEN
          DO IJ = IJS,IJL
            IF (ICEMASK(IJ,IG).EQ.0) THEN
              ESEA(IJ)    = ZMISS
              ESWELL(IJ)  = ZMISS
              THWISEA(IJ) = ZMISS
              THSWELL(IJ) = ZMISS
              FSEA(IJ)    = ZMISS
              FSWELL(IJ)  = ZMISS
            ENDIF
          ENDDO
        ENDIF

        IF(LMESSPASS) THEN
          CALL MPGATHERSCFLD(IPFGTBL(11),11,NSTART,NEND,ESEA)
          IF(IPFGTBL(11).EQ.IRANK) 
     &     CALL MAKEGRID (ESEA,   WHGTWG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(12),12,NSTART,NEND,ESWELL)
          IF(IPFGTBL(12).EQ.IRANK) 
     &     CALL MAKEGRID (ESWELL, WHGTSG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(13),13,NSTART,NEND,THWISEA)
          IF(IPFGTBL(13).EQ.IRANK) 
     &     CALL MAKEGRID (THWISEA,WDIRWG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(14),14,NSTART,NEND,THSWELL)
          IF(IPFGTBL(14).EQ.IRANK) 
     &     CALL MAKEGRID (THSWELL,WDIRSG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(15),15,NSTART,NEND,FSEA)
          IF(IPFGTBL(15).EQ.IRANK) 
     &     CALL MAKEGRID (FSEA,   WMNFWG, 1, NEND(NPROC), IG)

          CALL MPGATHERSCFLD(IPFGTBL(16),16,NSTART,NEND,FSWELL)
          IF(IPFGTBL(16).EQ.IRANK) 
     &     CALL MAKEGRID (FSWELL, WMNFSG, 1, NEND(NPROC), IG)

        ELSE
          CALL MAKEGRID(ESWELL, WHGTSG, IJS, IJL, IG)
          CALL MAKEGRID(THSWELL,WDIRSG, IJS, IJL, IG)
          CALL MAKEGRID(FSWELL, WMNFSG, IJS, IJL, IG)
          CALL MAKEGRID(ESEA,   WHGTWG, IJS, IJL, IG)
          CALL MAKEGRID(THWISEA,WDIRWG, IJS, IJL, IG)
          CALL MAKEGRID(FSEA,   WMNFWG, IJS, IJL, IG)
        ENDIF
      ENDIF


! 3. DEALLOCATE LOCAL ARRAYS
! --------------------------
      IF(ALLOCATED(PEAKF)) DEALLOCATE (PEAKF)
      IF(ALLOCATED(CD)) DEALLOCATE(CD)
      IF(ALLOCATED(WVSTR)) DEALLOCATE(WVSTR)

      IF(LMESSPASS) THEN
        CALL MPE_BARRIER(IERR)
        IF(IERR.LT.0)CALL MPEI_ABORT('MPE_BARRIER ERROR in OUTGRID')
      ENDIF

      RETURN
      END SUBROUTINE OUTGRID
