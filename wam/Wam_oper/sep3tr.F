      SUBROUTINE SEP3TR (FL1, IJS, IJL, MIJ, U10NEW, THWNEW ,
     &                   ESWELL, FSWELL, THSWELL, FSEA,
     &                   ETNP,
     &                   EMTRAIN  ,THTRAIN  ,PMTRAIN)

! ----------------------------------------------------------------------

!**** *SEP3TR* -  COMPUTE ENERGY, DIRECTION AND PERIOD FOR THE 3 WAVE
!****             TRAINS: SWELL1, SWELL2, SWELL3

!     D.PETTENUZZO     MAY 2012
!     JEAN BIDLOT      JUNE 2013  SIMPLIFIED BY REMOVING THE WIND SEA PART
!                                 AND ONLY PARTITIONED THE SWELL SPECTRUM
!                                 AND ONLY LOOKING IN THE PROGNOSTIC RANGE


!*    PURPOSE.
!     --------

!     CREATED TO TEST WAVE SWELL SPECTRA PARTITIONING INTO 3 WAVE TRAINS 

!**   INTERFACE.
!     ----------

!       *CALL* *SEP3TR (FL1, IJS, IJL, MIJ, U10NEW, THWNEW,
!                       ESWELL, FSWELL, THSWELL, FSEA,
!                       ETNP,
!                       EMTRAIN  ,THTRAIN  ,PMTRAIN)
!          *FL1*    - BLOCK OF SPECTRA
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *U10NEW* - LATESt WIND SPEED.
!          *THWNEW* - LATEST WIND DIRECTION.
!          *ESWELL* - TOTAL SWELL ENERGY
!          *FSWELL* - TOTAL SWELL MEAN FREQUENCY
!          *THSWELL*- TOTAL SWELL MEAN DIRECTION
!          *FSEA*   - WINDSEA MEAN FREQUENCY
!          *ETNP*   - TOTAL ENERGY NOT HELD BY ALL PARTITIONS.
!                     THIS IS POSSIBLE BECAUSE WE SUPPLY SWELL SPECTRA
!                     IN WHICH THE WIND SEA HAS BEEN REMOVED. SOME OF
!                     THE ENERGY LEAKING OUT FROM THE WINDSEA SPECTRUM
!                     MIGHT NOT BE PARTITIONED.
!          *EMTRAIN* - PARTITIONED ENERGY
!          *THTRAIN* - PARTITIONED MEAN DIRECTION
!          *PMTRAIN* - PARTITIONED MEAN PERIOD

!     METHOD.
!     -------

!       HANSON AND PHILLIPS 2001

!     EXTERNALS.
!     ----------

!       *FNDPRT*  - COMPUTE PARTITION MASKS

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : NTRAIN
      USE YOWFRED  , ONLY : FR       ,DFIM     ,C        ,
     &            TH       ,COSTH    ,SINTH    ,FRIC
      USE YOWICE   , ONLY : FLMIN    ,FLMINFR 
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI      ,G        ,EPSMIN
      USE YOWTABL  , ONLY : JUMAX    ,DELU
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IJS, IJL
      INTEGER, DIMENSION(IJS:IJL), INTENT(IN) :: MIJ
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: U10NEW, THWNEW 
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: ESWELL, FSWELL, THSWELL
      REAL, DIMENSION(IJS:IJL), INTENT(IN) :: FSEA
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL, DIMENSION(IJS:IJL), INTENT(OUT) :: ETNP
      REAL, DIMENSION(IJS:IJL,NTRAIN), INTENT(OUT) :: EMTRAIN
      REAL, DIMENSION(IJS:IJL,NTRAIN), INTENT(OUT) :: THTRAIN  ,PMTRAIN

      INTEGER, PARAMETER :: NPMAX=20
      INTEGER :: NPMAX_LOC
      INTEGER :: IJ, M, K, IP
      INTEGER :: ISORT, I, IPLOC
      INTEGER :: IFL, IFH, ITHL, ITHH
      INTEGER :: KM, KP
      INTEGER, DIMENSION(IJS:IJL) :: NPEAK, NPK
      INTEGER, DIMENSION(IJS:IJL) :: FRINVMIJ
      INTEGER, DIMENSION(IJS:IJL) :: MMIN, MMAX 
      INTEGER, DIMENSION(IJS:IJL) :: IPNOW 
      INTEGER, DIMENSION(IJS:IJL) :: JU 
      INTEGER, DIMENSION(IJS:IJL,NTRAIN) :: IENERGY
      INTEGER, DIMENSION(IJS:IJL,NPMAX) :: NFRP, NTHP

      ! relative value above max swell value that is considered above noise level
      REAL, PARAMETER :: XNOISELEVEL=0.005
      ! minimum Hs (m) value to keep a partition HSMIN=HSMIN_INTER+HSMIN_SLOPE*PERIOD
      REAL, PARAMETER :: HSMIN_INTER=0.05
      REAL, PARAMETER :: HSMIN_SLOPE=-0.002
      REAL :: THRS
      REAL :: HSMIN

      REAL :: COSDIFF
      REAL :: DELDW
      REAL :: XJ, DELUINV 
      REAL :: COSDIR, FRLIMIT 
      REAL :: ZHOOK_HANDLE
      REAL :: FLLOWEST
      REAL, DIMENSION(IJS:IJL) :: ENEX, SUMETRAIN
      REAL, DIMENSION(IJS:IJL) :: ETOT
      REAL, DIMENSION(IJS:IJL) :: ETT
      REAL, DIMENSION(IJS:IJL) :: ENMAX, FLNOISE
      REAL, DIMENSION(IJS:IJL,NANG) :: SPRD
      REAL, DIMENSION(IJS:IJL,0:NPMAX) :: DIR,PER,ENE
      REAL, DIMENSION(IJS:IJL,NTRAIN) :: TEMPDIR,TEMPPER,TEMPENE
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: FL3, FLLOW

      LOGICAL :: LLEPSMIN
      LOGICAL :: LLADDPART
      LOGICAL, DIMENSION(IJS:IJL,NTRAIN) :: LPWSECTOR
      LOGICAL, DIMENSION(IJS:IJL,NANG) :: LLCOSDIFF

! ----------------------------------------------------------------------

!*    LOOP THROUGH THE GRID POINTS
!     ----------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SEP3TR',0,ZHOOK_HANDLE)
#endif

      DO IJ=IJS,IJL
        FRINVMIJ(IJ)=1./FR(MIJ(IJ))
      ENDDO

!     TOTAL ENERGY IN THE SWELL SPECTRA
      LLEPSMIN=.FALSE.
      CALL SEMEAN (FL1, IJS, IJL, ETOT, LLEPSMIN)

      DO IJ=IJS,IJL
        ENMAX(IJ)=0.
      ENDDO
!     SMOOTH INPUT SPECTRA (in direction only)
      DO M=1,NFRE
        DO K=1,NANG
          KM=K-1
          IF (KM.LT.1) KM=NANG
          KP=K+1
          IF (KP.GT.NANG) KP=1
          DO IJ=IJS,IJL
!           RE-IMPOSE THE WINDSEA MASK
            IF(FL1(IJ,K,M).LE.0.) THEN
              FL3(IJ,K,M) = 0.
            ELSE
              FL3(IJ,K,M) = 0.10*(FL1(IJ,KM,M)+FL1(IJ,KP,M)) +
     &                      0.80*FL1(IJ,K,M)
              ENMAX(IJ)=MAX(ENMAX(IJ),FL3(IJ,K,M))
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        FLNOISE(IJ) = XNOISELEVEL*ENMAX(IJ)
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          EMTRAIN(IJ,ISORT) = 0. 
          PMTRAIN(IJ,ISORT) = 0.
          IF(ESWELL(IJ).GT.0.) THEN
            THTRAIN(IJ,ISORT) = THSWELL(IJ)
          ELSE
            THTRAIN(IJ,ISORT) = THWNEW(IJ)
          ENDIF
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        NPEAK(IJ)=0
      ENDDO

!*    1. DETERMINATES MAXIMA POSITIONS AND NUMBER
!        ----------------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          COSDIFF=COS(TH(K)-THWNEW(IJ))
          LLCOSDIFF(IJ,K)=(COSDIFF.LT.0.)
          SPRD(IJ,K)=MAX(0.,COSDIFF)**2
        ENDDO
      ENDDO
      DELUINV=1./DELU
      DO IJ=IJS,IJL
        XJ=U10NEW(IJ)*DELUINV
        JU(IJ)=MIN(JUMAX, MAX(NINT(XJ),1))
      ENDDO

      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FLLOW(IJ,K,M)=MAX(FLMINFR(JU(IJ),M)*SPRD(IJ,K),FLMIN)
          ENDDO
        ENDDO
      ENDDO


      DO IJ=IJS,IJL
        OUT: DO M= 2,MIJ(IJ)-1
          IFL    = MAX ( 1 , M-1 )
          IFH    = MIN ( NFRE , M+1 )
          DO K=1,NANG

            FLLOWEST = MAX(FLLOW(IJ,K,M),FLNOISE(IJ))

            IF(FL3(IJ,K,M) .GT. FLLOWEST) THEN
              ITHL   = 1 + MOD(NANG+K-2,NANG)
              ITHH   = 1 + MOD(K,NANG)
              IF ( FL3(IJ,ITHL,M   ) .GT. 0.0 .AND.
     &             FL3(IJ,ITHH,M   ) .GT. 0.0 .AND.
     &             FL3(IJ,K   ,IFL ) .GT. 0.0 .AND.
     &             FL3(IJ,K   ,IFH ) .GT. 0.0 .AND.
     &             FL3(IJ,ITHL,IFL ) .GT. 0.0 .AND.
     &             FL3(IJ,ITHL,IFH ) .GT. 0.0 .AND.
     &             FL3(IJ,ITHH,IFL ) .GT. 0.0 .AND.
     &             FL3(IJ,ITHH,IFH ) .GT. 0.0 ) THEN


                IF ( FL3(IJ,K,M) .GE. FL3(IJ,K   ,IFL ) .AND.
     &               FL3(IJ,K,M) .GE. FL3(IJ,K   ,IFH ) .AND.
     &               FL3(IJ,K,M) .GE. FL3(IJ,ITHL,IFL ) .AND.
     &               FL3(IJ,K,M) .GE. FL3(IJ,ITHL,M   ) .AND.
     &               FL3(IJ,K,M) .GE. FL3(IJ,ITHL,IFH ) .AND.
     &               FL3(IJ,K,M) .GE. FL3(IJ,ITHH,IFL ) .AND.
     &               FL3(IJ,K,M) .GE. FL3(IJ,ITHH,M   ) .AND.
     &               FL3(IJ,K,M) .GE. FL3(IJ,ITHH,IFH ) ) THEN
                    NPEAK(IJ) = NPEAK(IJ) + 1
                    IF(NPEAK(IJ).GT.NPMAX) EXIT OUT
                    NFRP(IJ,NPEAK(IJ)) = M
                    NTHP(IJ,NPEAK(IJ)) = K
                ENDIF
              ENDIF
            ENDIF

          ENDDO

        ENDDO OUT
      ENDDO

      DO IJ=IJS,IJL
        NPEAK(IJ)=MIN(NPEAK(IJ),NPMAX)
      ENDDO

!*    2. GENERATE MASK FOR EACH PARTITION AND COMPUTE STATISTICS
!        -------------------------------------------------------

      CALL FNDPRT(IJS, IJL, NPMAX,
     &            NPEAK, MIJ, NTHP, NFRP,
     &            FLLOW, LLCOSDIFF, FLNOISE,
     &            FL3,
     &            ENE, DIR, PER)


!     REMOVE PARTITION WITH Hs < HSMIN 
!     MEAN PERIOD > 1/FR(MIJ)
      DO IJ=IJS,IJL
        NPK(IJ)=NPEAK(IJ)
        DO IP=1,NPEAK(IJ)
          HSMIN=HSMIN_INTER+HSMIN_SLOPE*PER(IJ,IP)
          THRS=0.0625*HSMIN**2
          IF(ENE(IJ,IP).LT.THRS .OR. PER(IJ,IP).LE.FRINVMIJ(IJ)) THEN
            ENE(IJ,IP)=0.
            DIR(IJ,IP)=0.
            PER(IJ,IP)=0.
            NPK(IJ)=NPK(IJ)-1
          ENDIF
        ENDDO
      ENDDO

!     IF NO SWELL PARTITION BUT TOTAL SWELL EXIST THEN
!     ASSIGN FIRST PARTITION TO TOTAL SWELL
!     IF IT HAS SWELL CHARACTERSITICS
      DO IJ=IJS,IJL
        IF(NPK(IJ).LE.0 .AND. ESWELL(IJ).GT.0.) THEN
          IF(FSWELL(IJ).LT.FSEA(IJ)) THEN
            NPEAK(IJ)=1
            ENE(IJ,1)=ESWELL(IJ)
            DIR(IJ,1)=THSWELL(IJ)
            PER(IJ,1)=1./FSWELL(IJ)
          ENDIF 
        ENDIF 
      ENDDO

!*    5. SORT PARTITIONS ACCORDING TO ENERGY AND ASSIGN THE FIRST NTRAIN TRAINS
!        WE ONLY NEED THE FIRST NTRAIN TO BE SORTED
!        -----------------------------------------------------------------

      NPMAX_LOC=MAXVAL(NPEAK(:))

      DO IJ=IJS,IJL
        ETT(IJ)=ENE(IJ,1)
      ENDDO
      DO IP=2,NPMAX_LOC
        DO IJ=IJS,IJL
          ETT(IJ)=ETT(IJ)+ENE(IJ,IP)
        ENDDO
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          IPNOW(IJ)=0
          ENMAX(IJ)=0.
        ENDDO

        DO IP=1,NPMAX_LOC
          DO IJ=IJS,IJL
            IF (ENE(IJ,IP).GT.ENMAX(IJ)) THEN
              IPNOW(IJ) = IP
              ENMAX(IJ) = ENE(IJ,IP)
            ENDIF
          ENDDO
        ENDDO

        DO IJ=IJS,IJL
          EMTRAIN(IJ,ISORT)=ENE(IJ,IPNOW(IJ))
          THTRAIN(IJ,ISORT)=DIR(IJ,IPNOW(IJ))
          PMTRAIN(IJ,ISORT)=PER(IJ,IPNOW(IJ))
          ENE(IJ,IPNOW(IJ))=0.
          IENERGY(IJ,ISORT)=MIN(IPNOW(IJ),1)
        ENDDO
      ENDDO


!*    6. PRESERVE TOTAL ENERGY
!        ---------------------

!     6.1 Distribute extra energy proportionally to swell trains

      DO IJ=IJS,IJL
        SUMETRAIN(IJ)=MAX(EMTRAIN(IJ,1),EPSMIN)
      ENDDO
      DO ISORT=2,NTRAIN
        DO IJ=IJS,IJL
          SUMETRAIN(IJ)=SUMETRAIN(IJ)+EMTRAIN(IJ,ISORT)
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        ENEX(IJ)=MAX((ETT(IJ)-SUMETRAIN(IJ)),0.)/SUMETRAIN(IJ)
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          EMTRAIN(IJ,ISORT)=
     &      EMTRAIN(IJ,ISORT)+ENEX(IJ)*EMTRAIN(IJ,ISORT)
        ENDDO
      ENDDO

!     6.2 Distribute the non-partioned energy proportionally to swell trains
!         that are in the wind direction

!     non-partitioned energy
      DO IJ=IJS,IJL
        ETNP(IJ)=MAX(ETOT(IJ)-ETT(IJ),0.)
      ENDDO

      DO IJ=IJS,IJL
        SUMETRAIN(IJ)=EPSMIN
      ENDDO
      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          IF(COS(THWNEW(IJ)-THTRAIN(IJ,ISORT)).GT.0.1 .AND.
     &       EMTRAIN(IJ,ISORT).GT.EPSMIN) THEN
            LPWSECTOR(IJ,ISORT)=.TRUE.
            SUMETRAIN(IJ)=SUMETRAIN(IJ)+EMTRAIN(IJ,ISORT)
          ELSE
            LPWSECTOR(IJ,ISORT)=.FALSE.
          ENDIF
        ENDDO
      ENDDO
      DO IJ=IJS,IJL
        ENEX(IJ)=ETNP(IJ)/SUMETRAIN(IJ)
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          IF(LPWSECTOR(IJ,ISORT)) THEN
            EMTRAIN(IJ,ISORT)=
     &      EMTRAIN(IJ,ISORT)+ENEX(IJ)*EMTRAIN(IJ,ISORT)
            ETNP(IJ)=0.
          ENDIF
        ENDDO
      ENDDO


!     RE-ORDER SWELL SYSTEM ACCORDING TO TOTAL ENERGY

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          TEMPENE(IJ,ISORT)=EMTRAIN(IJ,ISORT)
          TEMPDIR(IJ,ISORT)=THTRAIN(IJ,ISORT)
          TEMPPER(IJ,ISORT)=PMTRAIN(IJ,ISORT)
        ENDDO
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          IPNOW(IJ)=0
          ENMAX(IJ)=0.
        ENDDO

        DO IP=1,NTRAIN
          DO IJ=IJS,IJL
            IF (TEMPENE(IJ,IP).GT.ENMAX(IJ)) THEN
              IPNOW(IJ) = IP
              ENMAX(IJ) = TEMPENE(IJ,IP)
            ENDIF
          ENDDO
        ENDDO

        DO IJ=IJS,IJL
          IPLOC=MAX(IPNOW(IJ),1)
          EMTRAIN(IJ,ISORT)=TEMPENE(IJ,IPLOC)
          THTRAIN(IJ,ISORT)=TEMPDIR(IJ,IPLOC)
          PMTRAIN(IJ,ISORT)=TEMPPER(IJ,IPLOC)
          TEMPENE(IJ,IPLOC)=0.
        ENDDO
      ENDDO

!     PREPARE OUTPUT

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          IF(IENERGY(IJ,ISORT).EQ.0) THEN
            EMTRAIN(IJ,ISORT) = 0.0
            THTRAIN(IJ,ISORT) = THWNEW(IJ)
            PMTRAIN(IJ,ISORT) = 0.0 
          ENDIF
        ENDDO
      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('SEP3TR',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE SEP3TR
