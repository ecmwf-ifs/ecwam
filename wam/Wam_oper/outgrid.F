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

!       *LOC2GLO*   - TRANSFERS A LOCAL BLOCK ARRAY TO A GLOBAL BLOCK 
!                     ARRAY AND IMPOSES ICE MASK.
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
     &            U10G     ,MWP1G    ,MWP2G    ,WSPRDG
      USE YOWICE   , ONLY : LICERUN  ,ICEMASK
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG   ,MWP1SG   ,MWP2SG   ,WSPRDSG  ,
     &            MWP1WG   ,MWP2WG   ,WSPRDWG
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,SMEAN    ,
     &            P1MEAN   ,P2MEAN   ,SPRDMEAN
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NIBLO
      USE YOWPCONS , ONLY : ZMISS    ,EPSUS    ,EPSU10
      USE YOWSTAT  , ONLY : CDTPRO   ,CDTINTS
      USE YOWSWEL  , ONLY : ESWELL   ,FSWELL   ,THSWELL  ,ESEA     ,
     &            FSEA     ,THWISEA  ,P1SWELL  ,P2SWELL  ,SPRDSWELL,
     &            P1SEA    ,P2SEA    ,SPRDSEA
      USE MPL_MODULE

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL3
      REAL,DIMENSION(NINF:NSUP,NBLO) :: TAUW
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW

! ----------------------------------------------------------------------

      REAL,ALLOCATABLE,DIMENSION(:) :: PEAKF,CD,WVSTR
      REAL :: TEMP(NIBLO)

      LOGICAL :: LSQRT

! ----------------------------------------------------------------------

!*    1. INTEGRATED PARAMETERS OF TOTAL SEA.
!        -----------------------------------

        ALLOCATE (CD(IJS:IJL))
        IF(IPFGTBL(8).NE.0) ALLOCATE(WVSTR(IJS:IJL))
        DO IJ = IJS,IJL
          TAU = MAX(USNEW(IJ)**2,EPSUS)
          IF(IPFGTBL(8).NE.0) WVSTR(IJ) = TAUW(IJ,IG)/TAU
!!       if the numerical computation of TAU and CD change, a similar
!!       modification has to be put in buildstress where the friction
!!       velocity is determined from U10 and CD.
          CD(IJ) = TAU/MAX(U10NEW(IJ)**2,EPSU10)
        ENDDO

!*    1.3 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
!         -------------------------------------------------

        IF(LMESSPASS) THEN
          IF(IPFGTBL(1).NE.0) THEN 
            CALL LOC2GLO(EMEAN(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MPGATHERSCFLD(IPFGTBL(1),1,NSTART,NEND,TEMP)
            IF(IPFGTBL(1).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP , WHGTTG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(2).NE.0) THEN 
            CALL LOC2GLO(THQ(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(2),2,NSTART,NEND,TEMP)
            IF(IPFGTBL(2).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP , WDIRTG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(3).NE.0) THEN 
            CALL LOC2GLO(FMEAN(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(3),3,NSTART,NEND,TEMP)
            IF(IPFGTBL(3).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP , WMNFTG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(4).NE.0) THEN
            DO IJ = IJS,IJL
              TEMP(IJ) = USNEW(IJ)
!             cannot set USNEW to ZMISS since it can still be used
!             in the next time step
              IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) THEN
                TEMP(IJ)=ZMISS
              ELSE
                TEMP(IJ) = USNEW(IJ)
              ENDIF
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(4),4,NSTART,NEND,TEMP)
            IF(IPFGTBL(4).EQ.IRANK) 
     &        CALL MAKEGRID (TEMP, USTARG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(5).NE.0) THEN 
            DO IJ = IJS,IJL
!             cannot set THWNEW to ZMISS since it can still be used
!             in the next time step
              IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) THEN
                TEMP(IJ)=ZMISS
              ELSE
                TEMP(IJ) = THWNEW(IJ) 
              ENDIF
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(5),5,NSTART,NEND,TEMP)
            IF(IPFGTBL(5).EQ.IRANK) 
     &        CALL MAKEGRID (TEMP, UDIRG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(6).NE.0) THEN
            ALLOCATE (PEAKF(IJS:IJL))
            CALL PEAKFR (FL3, IJS, IJL, PEAKF)
            CALL LOC2GLO(PEAKF,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(6),6,NSTART,NEND,TEMP)
            IF(IPFGTBL(6).EQ.IRANK)
     &      CALL MAKEGRID (TEMP , WPKFTG, 1, NEND(NPROC), IG, LSQRT)
            DEALLOCATE(PEAKF)
          ENDIF

          IF(IPFGTBL(7).NE.0) THEN 
            CALL LOC2GLO(CD,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(7),7,NSTART,NEND,TEMP)
            IF(IPFGTBL(7).EQ.IRANK)
     &      CALL MAKEGRID (TEMP ,  CDG, 1, NEND(NPROC), IG, LSQRT)
            IF(ALLOCATED(CD)) DEALLOCATE(CD)
          ENDIF

          IF(IPFGTBL(8).NE.0) THEN
            CALL LOC2GLO(WVSTR,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(8),8,NSTART,NEND,TEMP)
            IF(IPFGTBL(8).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,  TAUWG, 1, NEND(NPROC), IG, LSQRT)
            DEALLOCATE(WVSTR)
          ENDIF

          IF(IPFGTBL(9).NE.0) THEN
            CALL LOC2GLO(SMEAN(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(9),9,NSTART,NEND,TEMP)
            IF(IPFGTBL(9).EQ.IRANK)
     &      CALL MAKEGRID (TEMP , SMEANG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(10).NE.0) THEN
            DO IJ = IJS,IJL
!             cannot set U10NEW to ZMISS since it can still be used
!             in the next time step
              IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) THEN
                TEMP(IJ)=ZMISS
              ELSE
                TEMP(IJ) = U10NEW(IJ)
              ENDIF
            ENDDO
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(10),10,NSTART,NEND,TEMP)
            IF(IPFGTBL(10).EQ.IRANK) 
     &        CALL MAKEGRID (TEMP, U10G, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(25).NE.0) THEN
            CALL LOC2GLO(P1MEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(25),25,NSTART,NEND,TEMP)
            IF(IPFGTBL(25).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,  MWP1G, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(26).NE.0) THEN
            CALL LOC2GLO(P2MEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(26),26,NSTART,NEND,TEMP)
            IF(IPFGTBL(26).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,  MWP2G, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(27).NE.0) THEN
            CALL LOC2GLO(SPRDMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(27),27,NSTART,NEND,TEMP)
            IF(IPFGTBL(27).EQ.IRANK)
     &      CALL MAKEGRID (TEMP,  WSPRDG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

        ELSE

          IF(IPFGTBL(1).NE.0) THEN
            CALL LOC2GLO(EMEAN(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MAKEGRID (TEMP , WHGTTG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(2).NE.0) THEN
            CALL LOC2GLO(THQ(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP , WDIRTG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(3).NE.0) THEN
            CALL LOC2GLO(FMEAN(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP , WMNFTG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(4).NE.0) THEN
            DO IJ = IJS,IJL
              TEMP(IJ) = USNEW(IJ)
!             cannot set USNEW to ZMISS since it can still be used
!             in the next time step
              IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) THEN
                TEMP(IJ)=ZMISS
              ELSE
                TEMP(IJ) = USNEW(IJ)
              ENDIF
            ENDDO
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP,  USTARG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(5).NE.0) THEN
            DO IJ = IJS,IJL
!           cannot set THWNEW to ZMISS since it can still be used
!           in the next time step
              IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) THEN
                TEMP(IJ)=ZMISS
              ELSE
                TEMP(IJ) = THWNEW(IJ) 
              ENDIF
            ENDDO
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP, UDIRG , IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(6).NE.0) THEN
            ALLOCATE (PEAKF(IJS:IJL))
            CALL PEAKFR (FL3, IJS, IJL, PEAKF)
            CALL LOC2GLO(PEAKF,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP , WPKFTG, IJS, IJL, IG, LSQRT)
            DEALLOCATE (PEAKF)
          ENDIF

          IF(IPFGTBL(7).NE.0) THEN
            CALL LOC2GLO(CD,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP ,  CDG   , IJS, IJL, IG, LSQRT)
            IF(ALLOCATED(CD)) DEALLOCATE(CD)
          ENDIF

          IF(IPFGTBL(8).NE.0) THEN
            CALL LOC2GLO(WVSTR,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP,  TAUWG , IJS, IJL, IG, LSQRT)
            DEALLOCATE(WVSTR)
          ENDIF

          IF(IPFGTBL(9).NE.0) THEN
            CALL LOC2GLO(SMEAN(IJS),TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP , SMEANG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(10).NE.0) THEN
            DO IJ = IJS,IJL
!           cannot set U10NEW to ZMISS since it can still be used
!           in the next time step
              IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) THEN
                TEMP(IJ)=ZMISS
              ELSE
                TEMP(IJ) = U10NEW(IJ)
              ENDIF
            ENDDO
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP, U10G  , IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(25).NE.0) THEN
            CALL LOC2GLO(P1MEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP,  MWP1G, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(26).NE.0) THEN
            CALL LOC2GLO(P2MEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP,  MWP2G, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(27).NE.0) THEN
            CALL LOC2GLO(SPRDMEAN,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP,  WSPRDG, IJS, IJL, IG, LSQRT)
          ENDIF

        ENDIF

! ----------------------------------------------------------------------

!*    2. INTEGRATED PARAMETERS OF SWELL AND WIND SEA.
!        --------------------------------------------

      IF (CDTINTS.EQ.CDTPRO) THEN

!*    2.1 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
!         -------------------------------------------------

        IF(LMESSPASS) THEN
          IF(IPFGTBL(11).NE.0) THEN
            CALL LOC2GLO(ESEA,TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MPGATHERSCFLD(IPFGTBL(11),11,NSTART,NEND,TEMP)
            IF(IPFGTBL(11).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,   WHGTWG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(12).NE.0) THEN
            CALL LOC2GLO(ESWELL,TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MPGATHERSCFLD(IPFGTBL(12),12,NSTART,NEND,TEMP)
            IF(IPFGTBL(12).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP, WHGTSG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(13).NE.0) THEN
            CALL LOC2GLO(THWISEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(13),13,NSTART,NEND,TEMP)
            IF(IPFGTBL(13).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WDIRWG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(14).NE.0) THEN
            CALL LOC2GLO(THSWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(14),14,NSTART,NEND,TEMP)
            IF(IPFGTBL(14).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WDIRSG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(15).NE.0) THEN
            CALL LOC2GLO(FSEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(15),15,NSTART,NEND,TEMP)
            IF(IPFGTBL(15).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP,WMNFWG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(16).NE.0) THEN
            CALL LOC2GLO(FSWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(16),16,NSTART,NEND,TEMP)
            IF(IPFGTBL(16).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP, WMNFSG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(28).NE.0) THEN
            CALL LOC2GLO(P1SEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(28),28,NSTART,NEND,TEMP)
            IF(IPFGTBL(28).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP, MWP1WG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(29).NE.0) THEN
            CALL LOC2GLO(P1SWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(29),29,NSTART,NEND,TEMP)
            IF(IPFGTBL(29).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP, MWP1SG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(30).NE.0) THEN
            CALL LOC2GLO(P2SEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(30),30,NSTART,NEND,TEMP)
            IF(IPFGTBL(30).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP, MWP2WG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(31).NE.0) THEN
            CALL LOC2GLO(P2SWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(31),31,NSTART,NEND,TEMP)
            IF(IPFGTBL(31).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP, MWP2SG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(32).NE.0) THEN
            CALL LOC2GLO(SPRDSEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(32),32,NSTART,NEND,TEMP)
            IF(IPFGTBL(32).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP, WSPRDWG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

          IF(IPFGTBL(33).NE.0) THEN
            CALL LOC2GLO(SPRDSWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MPGATHERSCFLD(IPFGTBL(33),33,NSTART,NEND,TEMP)
            IF(IPFGTBL(33).EQ.IRANK) 
     &      CALL MAKEGRID (TEMP, WSPRDSG, 1, NEND(NPROC), IG, LSQRT)
          ENDIF

        ELSE
          IF(IPFGTBL(12).NE.0) THEN
            CALL LOC2GLO(ESWELL,TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MAKEGRID(TEMP, WHGTSG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(14).NE.0) THEN
            CALL LOC2GLO(THSWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID(TEMP,WDIRSG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(16).NE.0) THEN
            CALL LOC2GLO(FSWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID(TEMP, WMNFSG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(11).NE.0) THEN
            CALL LOC2GLO(ESEA,TEMP,IJS,IJL,IG)
            LSQRT=.TRUE.
            CALL MAKEGRID(TEMP, WHGTWG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(13).NE.0) THEN
            CALL LOC2GLO(THWISEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID(TEMP,WDIRWG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(13).NE.0) THEN
            CALL LOC2GLO(FSEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID(TEMP, WMNFWG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(28).NE.0) THEN
            CALL LOC2GLO(P1SEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP, MWP1WG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(29).NE.0) THEN
            CALL LOC2GLO(P1SWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP, MWP1SG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(30).NE.0) THEN
            CALL LOC2GLO(P2SEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP, MWP2WG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(31).NE.0) THEN
            CALL LOC2GLO(P2SWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP, MWP2SG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(32).NE.0) THEN
            CALL LOC2GLO(SPRDSEA,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP, WSPRDWG, IJS, IJL, IG, LSQRT)
          ENDIF

          IF(IPFGTBL(33).NE.0) THEN
            CALL LOC2GLO(SPRDSWELL,TEMP,IJS,IJL,IG)
            LSQRT=.FALSE.
            CALL MAKEGRID (TEMP, WSPRDSG, IJS, IJL, IG, LSQRT)
          ENDIF

        ENDIF
      ENDIF

      IF(LMESSPASS) THEN
        CALL MPL_BARRIER(CDSTRING='OUTGRID:')
      ENDIF

      RETURN
      END SUBROUTINE OUTGRID
