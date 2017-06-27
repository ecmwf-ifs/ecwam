      SUBROUTINE SEP3TR (FL1, IJS, IJL, MIJ, U10NEW, THWNEW , ETNP)

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

!       *CALL* *SEP3TR (FL1, IJS, IJL, MIJ, U10NEW, THWNEW, ETNP)
!          *FL1*    - BLOCK OF SPECTRA
!          *IJS*    - INDEX OF FIRST GRIDPOINT
!          *IJL*    - INDEX OF LAST GRIDPOINT
!          *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!          *U10NEW* - LATESt WIND SPEED.
!          *THWNEW* - LATEST WIND DIRECTION.
!          *ETNP*   - TOTAL ENERGY NOT HELD BY ALL PARTITIONS.
!                     THIS IS POSSIBLE BECAUSE WE SUPPLY SWELL SPECTRA
!                     IN WHICH THE WIND SEA HAS BEEN REMOVED. SOME OF
!                     THE ENERGY LEAKING OUT FROM THE WINDSEA SPECTRUM
!                     MIGHT NOT BE PARTITIONED.

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
      USE YOWTRAINS ,ONLY : EMTRAIN  ,THTRAIN  ,PMTRAIN
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, PARAMETER :: NPMAX=20
      INTEGER :: IJ, IJS, IJL, M, K, IP
      INTEGER :: ISORT, I, IPLOC
      INTEGER :: IFL, IFH, ITHL, ITHH
      INTEGER, DIMENSION(IJS:IJL) :: NPEAK
      INTEGER, DIMENSION(IJS:IJL) :: MIJ, MMIN, MMAX 
      INTEGER, DIMENSION(IJS:IJL) :: IPNOW 
      INTEGER, DIMENSION(IJS:IJL) :: JU 
      INTEGER, DIMENSION(IJS:IJL,NTRAIN) :: IENERGY
      INTEGER, DIMENSION(IJS:IJL,NPMAX) :: NFRP, NTHP

      REAL :: COSDIFF
      REAL :: DELDW
      REAL :: XJ, DELUINV 
      REAL :: COSDIR, FRLIMIT 
      REAL :: ZHOOK_HANDLE
      REAL :: FLLOWEST
      REAL :: ZZ(0:1)
      REAL, DIMENSION(IJS:IJL) :: U10NEW, THWNEW 
      REAL, DIMENSION(IJS:IJL) :: ENEX, SUMETRAIN
      REAL, DIMENSION(IJS:IJL) :: ETOT
      REAL, DIMENSION(IJS:IJL) :: ETT, ETNP
      REAL, DIMENSION(IJS:IJL) :: ENMAX
      REAL, DIMENSION(IJS:IJL,NANG) :: SPRD
      REAL, DIMENSION(IJS:IJL,0:NPMAX) :: DIR,PER,ENE
      REAL, DIMENSION(IJS:IJL,NTRAIN) :: TEMPDIR,TEMPPER,TEMPENE
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: FL1
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: FL3
      REAL, DIMENSION(IJS:IJL,NANG,NFRE) :: FLLOW 

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

!     TOTAL ENERGY IN THE SWELL SPECTRA
      LLEPSMIN=.FALSE.
      CALL SEMEAN (FL1, IJS, IJL, ETOT, LLEPSMIN)

!     SMOOTH INPUT SPECTRA
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL3(IJ,K,M) = FL1(IJ,K,M) 
          ENDDO
        ENDDO
      ENDDO
      CALL SMOOTHSARSPEC(IJS,IJL,NANG,NFRE,FL3)
!     RE-IMPOSE THE WINDSEA MASK
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            IF(FL1(IJ,K,M).LE.0.) FL3(IJ,K,M) = 0.
          ENDDO
        ENDDO
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          EMTRAIN(IJ,ISORT) = 0. 
          THTRAIN(IJ,ISORT) = THWNEW(IJ)
          PMTRAIN(IJ,ISORT) = 0.
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
          LLCOSDIFF(IJ,K)=(COSDIFF.LT.-0.4)
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
            FLLOW(IJ,K,M) = MAX(FLMINFR(JU(IJ),M)*SPRD(IJ,K),FLMIN)
          ENDDO
        ENDDO
      ENDDO


      DO IJ=IJS,IJL
        OUT: DO M= 1,MIJ(IJ)-1
          IFL    = MAX ( 1 , M-1 )
          IFH    = MIN ( NFRE , M+1 )
          DO K=1,NANG

            FLLOWEST = FLLOW(IJ,K,M)

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
     &            FLLOW, LLCOSDIFF,
     &            FL3,
     &            ENE, DIR, PER)

!*    5. SORT PARTITIONS ACCORDING TO ENERGY AND ASSIGN THE FIRST NTRAIN TRAINS
!        WE ONLY NEED THE FIRST NTRAIN TO BE SORTED
!        -----------------------------------------------------------------

      DO IJ=IJS,IJL
        ETT(IJ)=ENE(IJ,1)
      ENDDO
      DO IP=2,NPMAX
        DO IJ=IJS,IJL
          ETT(IJ)=ETT(IJ)+ENE(IJ,IP)
        ENDDO
      ENDDO

      DO ISORT=1,NTRAIN
        DO IJ=IJS,IJL
          IPNOW(IJ)=0
          ENMAX(IJ)=0.
        ENDDO

        DO IP=1,NPMAX
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
          IF(COS(THWNEW(IJ)-THTRAIN(IJ,ISORT)).GT.0. .AND.
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
