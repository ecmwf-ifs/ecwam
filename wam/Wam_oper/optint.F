       SUBROUTINE OPTINT(
!           INPUT
     &        NPARTW,NPASC2W,
     &        XMEANEW,
     &        EW,XKW,YKW,
     &        ED,XKD,YKD,
     &        SPECSA,CORRELTAS,
     &        PARTS,
!            OUTPUT
     &        EA,XANGA,XFREA)

!------------------------------------------------------------------

!     PURPOSE
!     -------
!     OPTIMAL INTERPOLATION SCHEME TO COMPUTE CORRECTION FACTORS FOR
!     ALL WAVE SYSTEMS AT WAM POINTS INFLUENCED BY SAR DATA POINTS.

!     AUTHOR
!     ------
!     S. HASSELMANN, MAX PLANCK INSTITUT FUER METEOROLOGIE,
!                    HAMBURG, GER.

!     EXTERNALS
!     ---------
!     SYMINV (ECLIB ROUTINE) 

!     METHOD.
!     -------
!     THE SAR DATA ERROR CORRELATION MATRIX IS COMPUTED FROM AN INVERSE
!     EXPONENTIAL DISTANCE BETWEEN THE CORRELATED POINTS.
!     THE MODEL ERROR CORRELATION MATRIX IS ASSUMED TO BE THE IDENTITY
!     MATRIX. 

!**********************************************************************

!     MODULES:
!     --------

      USE YOWICE   , ONLY : FLMIN
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,R        ,RAD      ,ZPI
      USE YOWSARAS , ONLY : NSPEC    ,NSPECW   ,JPSIW    ,MPARTSW  ,
     &            MAXMPART ,LONG     ,LAT      ,COSSARLAT,SINSARLAT,
     &            NTOSIW1D ,NSIW1D   ,SPECW    ,DIST     ,SARCORDIA
      USE YOWTEST  , ONLY : IU06

   
!     INTERFACE
!     ---------

      IMPLICIT NONE

!     INTEGERS WAM - SAR INDICES
      INTEGER  
     & NPARTW(NSPECW),
!           NUMBER OF PARTITIONINGS IN WAM FIRST GUESS SPECTRUM.
     & NPASC2W(NSPECW,MPARTSW,JPSIW) 
!           INDEX OF PARTITIONING IN SAR SPECTRUM NSIW1D(MSPECW,JPSIW)
!           CROSSASSIGNED TO PARTITIONING MPARTSW OF WAM SPECTRUM MSPEC.

      REAL
     & SARCORLEN,
!           CORRELATION LENGTHS (in radian) !!!
     & WEIGHT(NSPECW,MPARTSW,JPSIW),
!           INTERPOLATION WEIGHTS,
     & EW(NSPECW,MPARTSW),
     & XKW(NSPECW,MPARTSW),
     & YKW(NSPECW,MPARTSW),
!           CHARACTERISTIC PARAMETERS OF WAM FIRST GUESS WAVE 
!           SYSTEMS.
     & ED(NSPEC,0:MPARTSW),
     & XKD(NSPEC,0:MPARTSW),
     & YKD(NSPEC,0:MPARTSW)
!           ERRORS (difference) BETWEEN CHARACTERISTIC PARAMETERS
!           OF SAR WAVE SYSTEMS AND WAM WAVE SYSTEMS AT SAR DATA POINTS
!           (ENERGY AND MEAN WAVE NUMBER IN LAT AND LON DIRECTION).
!           MEAN FREQUENCY OF SAR WAVE SYSTEMS.

      REAL :: XMEANEW(NSPEC,MPARTSW)
!         MODEL FG  MEAN ENERGIES OF PARTITIONINGS AT SAR LOCATIONS
      REAL
     & EA(NSPECW,MPARTSW),
     & XKA(NSPECW,MPARTSW),
     & YKA(NSPECW,MPARTSW),
     & XANGA(NSPECW,MPARTSW),
     & XFREA(NSPECW,MPARTSW)
!           CORRECTION PARAMETERS FOR WAM WAVE SYSTEMS.
      REAL
     & PHI(NSPECW,MPARTSW,JPSIW)
!            CROSS CORRELATION BETWEEN  MODEL AND DATA ERROR.

!!!! in the long run more arrays will have to be allocatable
      REAL, ALLOCATABLE :: SIGMA(:,:)
!            ERROR COVARIANCE MATRIX
!     
      INTEGER
     & IFRE,IANG,IDSAR,IPARTS
      REAL
     &  SPECSA(0:NSPEC,NANG,NFRE)
      INTEGER
     & CORRELTAS(NSPECW,MPARTSW,JPSIW),PARTS(NSPEC,NANG,NFRE)

      REAL :: WSUME,EPS,COND
!            SUMS OF ERRORS.
      REAL, ALLOCATABLE :: WKSPACE(:)
!            INTERMEDIATE STORAGE SPACE.
      REAL :: SIGMASAR(NSPEC,NSPEC)
!     CORRELATION FUNCTION BETWEEN SAR OBS.

      REAL :: CORDIST,STRETCHCOR,YY

      INTEGER :: ISPEC,JSAR,KSAR,IPARTSAR,JPARTSAR,NJS,K,J,KK,JJ
!            LOOP INDICES.
      INTEGER :: NJSAR
!             NUMBER OF SAR DATA POINTS INFLUENCING WAM POINT.


      INTEGER :: MAXMPARTPE
!     MAXIMIUM OF PARTITIONS ON THAT PE.

!------------------------------------------------------------------

!     1. 

      NJSAR=0
      DO ISPEC=1,NSPECW
        NJSAR=MAX(NJSAR,NTOSIW1D(ISPEC))
      END DO
      IF(NJSAR.GT.JPSIW) THEN
        WRITE(IU06,*)'ERROR: NUMBER OF SAR POINTS INFLUENCING WAM'
        WRITE(IU06,*)'POINTS, NJSAR: ',NJSAR
        WRITE(IU06,*)'LARGER THEN DIMENSION OF INVERSE MATRIX.'
        WRITE(IU06,*)'JPSIW: ',JPSIW
        STOP 10
      END IF

!------------------------------------------------------------------

!      2. COMPUTE CROSS CORRELATION BETWEEN MODEL AND DATA ERROR.
!      -----------------------------------------------------------

      SARCORLEN=SARCORDIA/R

      MAXMPARTPE=0
      DO ISPEC=1,NSPECW
        MAXMPARTPE=MAX(MAXMPARTPE,NPARTW(ISPEC))
      ENDDO

      DO JSAR=1,JPSIW
        DO IPARTS=1,MAXMPARTPE
          DO ISPEC=1,NSPECW
            IF(IPARTS.LE.NPARTW(ISPEC).AND.JSAR.LE.NTOSIW1D(ISPEC))THEN
              STRETCHCOR=1.
              IPARTSAR = NPASC2W(ISPEC,IPARTS,JSAR)
              IF(IPARTSAR.GT.0 ) THEN
                YY=XMEANEW(NSIW1D(ISPEC,JSAR),IPARTSAR)
                IF(YY.GT.0.) THEN
                  STRETCHCOR = ABS(EW(ISPEC,IPARTS)-YY)/YY
                  STRETCHCOR = MAX(EXP(-STRETCHCOR),0.1) 
                ENDIF
              ENDIF
              CORDIST=STRETCHCOR*SARCORLEN
              PHI(ISPEC,IPARTS,JSAR) = EXP(-DIST(ISPEC,JSAR)/CORDIST)
            ENDIF
          END DO
        END DO
      END DO


!-------------------------------------------------------------------

!      3. COMPUTE INTERPOLATION WEIGHTS.
!      ---------------------------------

      ALLOCATE(SIGMA(NJSAR,NJSAR))
      ALLOCATE(WKSPACE(NJSAR))


      DO K = 1,NSPEC
        DO J = 1,K-1
          SIGMASAR(K,J) = COS(RAD*(LONG(K,2)-LONG(J,2))) *
     &                      COSSARLAT(J)*COSSARLAT(K) +
     &                      SINSARLAT(J)*SINSARLAT(K)
          SIGMASAR(K,J) = MAX(MIN(SIGMASAR(K,J),1.),-1.) 
          SIGMASAR(K,J) = EXP(-ACOS(SIGMASAR(K,J))/SARCORLEN)
          SIGMASAR(J,K) = SIGMASAR(K,J)
        END DO
        SIGMASAR(K,K) = 1.
      END DO

      WEIGHT = 0.

      DO ISPEC=1,NSPECW

        DO IPARTS=1,NPARTW(ISPEC)
!         create symmetric matrix SIGMA (SYMINV only requires the lower
!         triangular part.)
          NJS=0
          KSAR=0
          DO KK = 1,NTOSIW1D(ISPEC)
            K=NSIW1D(ISPEC,KK)
            IPARTSAR = NPASC2W(ISPEC,IPARTS,KK)
            IF(IPARTSAR.GT.0) THEN
              NJS=NJS+1
              KSAR=KSAR+1
              JSAR=0
              DO JJ = 1,KK-1
                J=NSIW1D(ISPEC,JJ)
                JPARTSAR = NPASC2W(ISPEC,IPARTS,JJ)
                IF(JPARTSAR.GT.0) THEN
                  JSAR=JSAR+1
                  SIGMA(KSAR,JSAR) = SIGMASAR(K,J) 
                ENDIF
              END DO
              SIGMA(KSAR,KSAR) = 1. + SIGMASAR(K,K) 
            ENDIF
          END DO

!         INVERT MATRIX SIGMA.

          IF(NJS.GT.1) THEN
            COND=0.
            EPS=0.
            WKSPACE=0.
            CALL SYMINV(SIGMA,NJSAR,NJS,COND,WKSPACE,EPS)
          ELSE IF (NJS.EQ.1) THEN
            SIGMA(1,1)=1./SIGMA(1,1)
          ENDIF

          DO KSAR=1,NJS
            DO JSAR=1,KSAR-1
              SIGMA(JSAR,KSAR)=SIGMA(KSAR,JSAR)
            END DO
          END DO

          JSAR=0
          DO JJ = 1,NTOSIW1D(ISPEC)
            JPARTSAR = NPASC2W(ISPEC,IPARTS,JJ)
            IF(JPARTSAR.GT.0) THEN
              JSAR=JSAR+1
              WSUME=0.
              KSAR=0
              DO KK = 1,NTOSIW1D(ISPEC)
                K=NSIW1D(ISPEC,KK)
                IPARTSAR = NPASC2W(ISPEC,IPARTS,KK)
                IF(IPARTSAR.GT.0) THEN
                  KSAR=KSAR+1
                  WSUME = WSUME + SIGMA(JSAR,KSAR)*PHI(ISPEC,IPARTS,KK) 
                ENDIF
              ENDDO
              WEIGHT(ISPEC,IPARTS,JJ) = WSUME
            ENDIF
          END DO
        ENDDO
      END DO

      DEALLOCATE(SIGMA)
      DEALLOCATE(WKSPACE)

!------------------------------------------------------------------

!     4.  COMPUTE VALUES TO CORRECT WAM FIRST GUESS WAVE SYSTEMS
!         INFLUENCING WAM POINT.
!     -----------------------------------------------------------

      EA = 0.
      XKA = 0. 
      YKA = 0.

      DO JSAR=1,JPSIW
        DO IPARTS=1,MAXMPARTPE
          DO ISPEC=1,NSPECW
            IF(IPARTS.LE.NPARTW(ISPEC).AND.JSAR.LE.NTOSIW1D(ISPEC))THEN
              IPARTSAR = NPASC2W(ISPEC,IPARTS,JSAR)

              IF(IPARTSAR.GT.0) THEN
                EA(ISPEC,IPARTS) = EA(ISPEC,IPARTS)+
     &                             WEIGHT(ISPEC,IPARTS,JSAR)*
     &                             ED(NSIW1D(ISPEC,JSAR),IPARTSAR)
                XKA(ISPEC,IPARTS) = XKA(ISPEC,IPARTS)+
     &                              WEIGHT(ISPEC,IPARTS,JSAR)*
     &                              XKD(NSIW1D(ISPEC,JSAR),IPARTSAR)
                YKA(ISPEC,IPARTS) = YKA(ISPEC,IPARTS)+
     &                              WEIGHT(ISPEC,IPARTS,JSAR)*
     &                              YKD(NSIW1D(ISPEC,JSAR),IPARTSAR)
              ENDIF
            ENDIF
          END DO
        END DO
      END DO

!!    insure that the increment does no produce negative energy
      DO IPARTS=1,MAXMPART
        DO ISPEC=1,NSPECW
          EA (ISPEC,IPARTS) = MAX(EW(ISPEC,IPARTS)+EA(ISPEC,IPARTS),0.)
!!!! test
!!!! always produce an analysis for wave number even when the
!!!! system is reduced to zero energy. It still needs to be shifted and
!!!! rotated !!!!!
!!!!          IF(EA(ISPEC,IPARTS).GT.0.) THEN
            XKA(ISPEC,IPARTS) = XKW(ISPEC,IPARTS) + XKA(ISPEC,IPARTS)
            YKA(ISPEC,IPARTS) = YKW(ISPEC,IPARTS) + YKA(ISPEC,IPARTS)
!!!!          ELSE
!!!!            XKA(ISPEC,IPARTS) = XKW(ISPEC,IPARTS)
!!!!            YKA(ISPEC,IPARTS) = YKW(ISPEC,IPARTS)
!!!!          ENDIF
        END DO
      END DO

!     COMPUTE MEAN ANALYSED DIRECTIONS AND FREQUENCIES FROM 
!     WAVE NUMBERS.

      DO IPARTS=1,MAXMPARTPE
        DO ISPEC=1,NSPECW
          IF(IPARTS.LE.NPARTW(ISPEC)) THEN
            XANGA(ISPEC,IPARTS) = 
     &               ATAN2(YKA(ISPEC,IPARTS),XKA(ISPEC,IPARTS))
            XFREA(ISPEC,IPARTS) = SQRT(SQRT(XKA(ISPEC,IPARTS)**2+
     &               YKA(ISPEC,IPARTS)**2)*G)/ZPI
          ENDIF
        END DO
      END DO

!------------------------------------------------------------------

!     5. SUPERIMPOSE EXTRA WAVE SYSTEMS FROM SAR RETRIEVAL ON 
!        FIRST GUESS SPECTRUM.
!     ----------------------------------------------------------

!!!   this is no longer used !!!!

      RETURN
      END SUBROUTINE OPTINT
