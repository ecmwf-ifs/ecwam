      SUBROUTINE SEP3TR (FL3, IJS, IJL, USOLD, THWOLD)

! ----------------------------------------------------------------------

!**** *SEP3TR* -  COMPUTE ENERGY, DIRECTION AND PERIOD FOR THE 3 WAVE
!****             TRAINS: WIND WAVE, SWELL1 AND SWELL2

!     D.PETTENUZZO     MAY 2012


!*    PURPOSE.
!     --------

!     CREATED TO TEST WAVE SPECTRA PARTITIONING INTO MORE THAN
!     2 WAVE TRAINS

!**   INTERFACE.
!     ----------

!       *CALL* *SEP3TR (FL3, IJS, IJL, USOLD, THWOLD)*
!          *FL3* - BLOCK OF SPECTRA
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *USOLD*   INTERMEDIATE STORAGE OF FRICTION VELOCITY
!          *THWOLD*  INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                    WIND DIRECTION.

!     METHOD.
!     -------

!       HANSON AND PHILLIPS 2001

!     EXTERNALS.
!     ----------

!       *FNDPRT*  - COMPUTE PARTITION MASKS
!       *PARMEAN* - COMPUTE MEAN PARAMETERS.

! ----------------------------------------------------------------------

      USE YOWFRED  ,  ONLY : FR       ,DFIM     ,C        ,
     &            TH       ,COSTH    ,SINTH    ,FRIC
      USE YOWICE   ,  ONLY : FLMIN
      USE YOWPARAM ,  ONLY : NANG     ,NFRE
      USE YOWPCONS ,  ONLY : ZPI      ,G        ,ZMISS    ,EPSMIN
      USE YOWTEST  ,  ONLY : IU06
      USE YOWTRAINS , ONLY: ETRAIN1  ,THTRAIN1 ,P1TRAIN1 ,
     &                      ETRAIN2  ,THTRAIN2 ,P1TRAIN2 ,
     &                      ETRAIN3  ,THTRAIN3 ,P1TRAIN3
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, PARAMETER :: NPMAX=20
      INTEGER :: IJ, IJS, IJL, M, K, IP
      INTEGER :: NPEAK, ISORT, I, IPNOW
      INTEGER :: IFL, IFH, ITHL, ITHH
      INTEGER,ALLOCATABLE :: WNDI(:)

      REAL :: DELDW, ENMAX, ENEX
      REAL :: COSDIR, FRLIMIT 
      REAL :: FP_DUM
      REAL :: ZHOOK_HANDLE
      REAL,DIMENSION(IJS:IJL,NANG,NFRE) :: FL1,FL3
      REAL,DIMENSION(IJS:IJL) :: THWOLD,USOLD
      INTEGER,DIMENSION(NPMAX) :: NFRP, NTHP
      REAL, DIMENSION(NPMAX) :: DIR,PER,ENE,FP
      REAL, DIMENSION(NPMAX) :: DIRSORT,PERSORT,ENESORT
      REAL,DIMENSION(NANG,NFRE) :: W1, W2, WWMSK, SW1MSK, SW2MSK
      REAL, ALLOCATABLE :: MSK(:,:,:)

! ----------------------------------------------------------------------

!*    LOOP THROUGH THE GRID POINTS
!     ----------------------------

      IF (LHOOK) CALL DR_HOOK('SEP3TR',0,ZHOOK_HANDLE)

      DO IJ=IJS,IJL


!*    1. DETERMINATES MAXIMA POSITIONS AND NUMBER
!        ----------------------------------------

        NPEAK=0

        W1 = 0.

        DO M= 1,NFRE
          IFL    = MAX ( 1 , M-1 )
          IFH    = MIN ( NFRE , M+1 )
          DO K=1,NANG
            ITHL   = 1 + MOD(NANG+K-2,NANG)
            ITHH   = 1 + MOD(K,NANG)
            IF ( FL3(IJ,K,M) .GE. FL3(IJ,K   ,IFL ) .AND.
     &           FL3(IJ,K,M) .GE. FL3(IJ,K   ,IFH ) .AND.
     &           FL3(IJ,K,M) .GE. FL3(IJ,ITHL,IFL ) .AND.
     &           FL3(IJ,K,M) .GE. FL3(IJ,ITHL,M   ) .AND.
     &           FL3(IJ,K,M) .GE. FL3(IJ,ITHL,IFH ) .AND.
     &           FL3(IJ,K,M) .GE. FL3(IJ,ITHH,IFL ) .AND.
     &           FL3(IJ,K,M) .GE. FL3(IJ,ITHH,M   ) .AND.
     &           FL3(IJ,K,M) .GE. FL3(IJ,ITHH,IFH ) .AND.
     &           FL3(IJ,K,M) .GT. FLMIN           ) THEN
                NPEAK  = NPEAK + 1
                IF(NPEAK.GT.NPMAX) EXIT
                NFRP(NPEAK) = M
                NTHP(NPEAK) = K
            ENDIF
          ENDDO
        ENDDO

        IF ( NPEAK .GT. NPMAX ) THEN
            NPEAK  = NPMAX
        ENDIF

        ALLOCATE(MSK(NANG,NFRE,NPEAK))
        ALLOCATE(WNDI(NPEAK))
        MSK=0.
        WNDI=0
        WWMSK=0.


!*    2. GENERATE MASK FOR EACH PARTITION AND COMPUTE STATISTICS
!        -------------------------------------------------------

        DO IP=1, NPEAK
         CALL FNDPRT(NTHP(IP),NFRP(IP),FL3(IJ:IJ,1:NANG,1:NFRE),W1,W2)
  
          DO M=1,NFRE
            DO K=1,NANG
              MSK(K,M,IP)=W2(K,M)
              FL1(IJ,K,M) = FL3(IJ,K,M)*W2(K,M)
            ENDDO
          ENDDO

          CALL PARMEAN(FL1(IJ:IJ,1:NANG,1:NFRE),
     &                 ENE(IP),DIR(IP),PER(IP),FP(IP))

        ENDDO
       

!*    3. DIVIDE SWELL AND WIND WAVE (THE LATTER IS GROUPED)
!        --------------------------------------------------

        DO IP=1, NPEAK
          COSDIR=COS(DIR(IP)-THWOLD(IJ))
          IF (COSDIR.GT.0. .AND. USOLD(IJ).GT.0.) THEN
!!! if it works make all constants into one coefficient
            FRLIMIT= G/(1.5*ZPI*FRIC*USOLD(IJ)*COSDIR)
            IF(FP(IP).GE.FRLIMIT) THEN
              DO M=1,NFRE
                DO K=1,NANG
                WWMSK(K,M) = WWMSK(K,M) + MSK(K,M,IP)
                ENDDO
              ENDDO
              WNDI(IP) = 1
            ENDIF
          ENDIF
        ENDDO
        

!*    4. COMPUTATION OF WIND WAVE OUTPUT PARAMETERS
!        ------------------------------------------

      DO M=1,NFRE
        DO K=1,NANG
          FL1(IJ,K,M) = FL3(IJ,K,M)*WWMSK(K,M)
        ENDDO
      ENDDO

      IF( MAXVAL(WNDI) .GT. 0)THEN
        CALL PARMEAN(FL1(IJ:IJ,1:NANG,1:NFRE),
     &               ETRAIN1(IJ),THTRAIN1(IJ),P1TRAIN1(IJ),FP_DUM)
      ELSE !if there isn't any windsea energy, set the windsea mean direction
           !to the wind direction.
        THTRAIN1(IJ)=THWOLD(IJ)
        ETRAIN1(IJ) =ZMISS
        P1TRAIN1(IJ)=ZMISS
      ENDIF 


!*    5. SORT SWELL PARTITIONS ACCORDING TO ENERGY AND ASSIGN 2ND AND 3RD TRAIN
!        ----------------------------------------------------------------------

      ENESORT=0.
      DIRSORT=0.
      PERSORT=0.
      IPNOW=0

      DO ISORT=1,NPEAK
        ENMAX   = 0.
        DO IP=1, NPEAK
          IF ( ENE(IP).GT.ENMAX .AND. WNDI(IP) .EQ. 0 ) THEN
              IPNOW  = IP
              ENMAX   = ENE(IP)
            ENDIF
        ENDDO

        IF(IPNOW.GT.0) THEN
          ENESORT(ISORT)=ENE(IPNOW)
          PERSORT(ISORT)=PER(IPNOW)
          DIRSORT(ISORT)=DIR(IPNOW)
          ENE(IPNOW)=0.
        ENDIF
        IPNOW  = 0

      ENDDO

      IF(ENESORT(1).GT.0.) THEN
        ETRAIN2(IJ)  = ENESORT(1)
        THTRAIN2(IJ) = DIRSORT(1)
        P1TRAIN2(IJ) = PERSORT(1)
      ELSE
        ETRAIN2(IJ)  = ZMISS
        THTRAIN2(IJ) = ZMISS
        P1TRAIN2(IJ) = ZMISS
      ENDIF
      IF(ENESORT(2).GT.0.) THEN
        ETRAIN3(IJ)  = ENESORT(2)
        THTRAIN3(IJ) = DIRSORT(2)
        P1TRAIN3(IJ) = PERSORT(2)
      ELSE
        ETRAIN3(IJ)  = ZMISS
        THTRAIN3(IJ) = ZMISS
        P1TRAIN3(IJ) = ZMISS
      ENDIF

!*    6. PRESERVE TOTAL ENERGY
!        ---------------------

!     6.1 Compute extra energy

      ENEX=0.
      DO IP=3, NPEAK
        ENEX=ENEX + ENESORT(IP)
      ENDDO

!     6.2 Distribute extra energy proportionally to swell trains

      ETRAIN2(IJ)=
     &  ETRAIN2(IJ)+ENEX *ETRAIN2(IJ)/(ETRAIN2(IJ)+ETRAIN3(IJ))
      ETRAIN3(IJ)=
     &  ETRAIN3(IJ)+ENEX *ETRAIN3(IJ)/(ETRAIN2(IJ)+ETRAIN3(IJ))


!*    7. DEALLOCATE ARRAYS
!        -----------------

        DEALLOCATE(MSK)
        DEALLOCATE(WNDI)


!*    END LOOP THROUGH THE GRID POINTS
!     ----------------------------

      ENDDO

      IF (LHOOK) CALL DR_HOOK('SEP3TR',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SEP3TR
