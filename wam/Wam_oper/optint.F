       SUBROUTINE OPTINT(
!           INPUT
     &        NPARTW,NPASC2W,
     &        EW,XKW,YKW,
     &        ED,XKD,YKD,SARCORLEN,
     &        SPECWA,SPECSA,CORRELTAS,
     &        PARTS,NPARTS,
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
      USE YOWPCONS , ONLY : G        ,RAD      ,ZPI      ,EPSMIN
      USE YOWSARAS , ONLY : NSPEC    ,NSPECW   ,JPSIW    ,MPARTSW  ,
     &            MAXMPART ,LONG     ,LAT      ,COSSARLAT,SINSARLAT,
     &            NTOSIW1D ,NSIW1D   ,SPECW    ,DIST
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
!           CORRELATION LENGTH (in radian) !!!
     & WEIGHT(NSPECW,JPSIW),
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

      REAL
     & EA(NSPECW,MPARTSW),
     & XKA(NSPECW,MPARTSW),
     & YKA(NSPECW,MPARTSW),
     & XANGA(NSPECW,MPARTSW),
     & XFREA(NSPECW,MPARTSW)
!           CORRECTION PARAMETERS FOR WAM WAVE SYSTEMS.
      REAL
     & SIGMAINV,
!            ERROR VARIANCE.
     & PHI(NSPECW,JPSIW)
!            CROSS CORRELATION BETWEEN  MODEL AND DATA ERROR.

!!!! in the long run more arrays will have to be allocatable
      REAL, ALLOCATABLE :: SIGMA(:,:)
!            ERROR COVARIANCE MATRIX
!     
      INTEGER
     & IFRE,IANG,NPARTS(NSPEC),IDSAR,IPARTS
      REAL
     &  SPECWA(0:NSPEC,NANG,NFRE),
     &  SPECSA(0:NSPEC,NANG,NFRE)
      REAL, ALLOCATABLE :: SPECDA(:,:,:)
!     TEMPOARARY ARRAY WHICH WILL CONTAIN SIGMAINV*(SPECWA - SPECSA)
!     FOR 1 to IFRECUT.
      INTEGER
     & CORRELTAS(NSPECW,MPARTSW,JPSIW),PARTS(NSPEC,NANG,NFRE)

      REAL :: WSUME,EPS,COND
!            SUMS OF ERRORS.
      REAL, ALLOCATABLE :: WKSPACE(:)
!            INTERMEDIATE STORAGE SPACE.
      REAL :: SIGMASAR(NSPEC,NSPEC)
!     CORRELATION FUNCTION BETWEEN SAR OBS.
      INTEGER :: ISPEC,JSAR,KSAR,IPARTSAR,NJS,K,J
!            LOOP INDICES.
      INTEGER :: NJSAR
!             NUMBER OF SAR DATA POINTS INFLUENCING WAM POINT.
      INTEGER :: IFRECUT
!     CUF OFF FREQUENCY INDEX WHICH SEPARATES LOW AND HIGH FREQUENCIES

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

!      2. COMPUTE CROSS CORRELATION BETWEEN  MODEL AND DATA ERROR.
!      -----------------------------------------------------------

      DO JSAR=1,NJSAR
        DO ISPEC=1,NSPECW
          PHI(ISPEC,JSAR) = EXP(-DIST(ISPEC,JSAR)/SARCORLEN)
        END DO
      END DO

!-------------------------------------------------------------------

!      3. COMPUTE INTERPOLATION WEIGHTS.
!      ---------------------------------

      ALLOCATE(SIGMA(NJSAR,NJSAR))
      ALLOCATE(WKSPACE(NJSAR))

      DO KSAR=1,NJSAR
        DO ISPEC=1,NSPECW
          WEIGHT(ISPEC,KSAR) = 0.
        END DO
      END DO

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

      SIGMAINV = 0.5

      DO ISPEC=1,NSPECW

        NJS=NTOSIW1D(ISPEC)

!       create symmetric matrix SIGMA (SYMINV only requires the lower
!       triangular part.)
        DO KSAR = 1,NJS
          K=NSIW1D(ISPEC,KSAR)
          DO JSAR = 1,KSAR-1
            SIGMA(KSAR,JSAR) = SIGMASAR(K,NSIW1D(ISPEC,JSAR))
          END DO
          SIGMA(KSAR,KSAR) = 1. + SIGMASAR(K,K) 
        END DO

!       INVERT MATRIX SIGMA.

        IF(NJS.GT.1) THEN
          COND=0.
          EPS=0.
          WKSPACE=0.
          CALL SYMINV(SIGMA,NJSAR,NJS,COND,WKSPACE,EPS)
        ELSE
          SIGMA(1,1)=1./SIGMA(1,1)
        ENDIF

        DO KSAR=1,NJS
          DO JSAR=1,KSAR-1
            SIGMA(JSAR,KSAR)=SIGMA(KSAR,JSAR)
          END DO
        END DO

        DO JSAR=1,NJS
          WSUME=0.
          DO KSAR=1,NJS
            WSUME = WSUME + 
     &      SIGMA(JSAR,KSAR)*PHI(ISPEC,KSAR)
          END DO
          WEIGHT(ISPEC,JSAR) = WSUME
        END DO
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
        DO IPARTS=1,MAXMPART
          DO ISPEC=1,NSPECW
            IF(IPARTS.LE.NPARTW(ISPEC).AND.JSAR.LE.NTOSIW1D(ISPEC))THEN
              IPARTSAR = NPASC2W(ISPEC,IPARTS,JSAR)
              EA(ISPEC,IPARTS) = EA(ISPEC,IPARTS)+WEIGHT(ISPEC,JSAR)*
     &                          ED(NSIW1D(ISPEC,JSAR),IPARTSAR)
              XKA(ISPEC,IPARTS) = XKA(ISPEC,IPARTS)+WEIGHT(ISPEC,JSAR)*
     &                           XKD(NSIW1D(ISPEC,JSAR),IPARTSAR)
              YKA(ISPEC,IPARTS) = YKA(ISPEC,IPARTS)+WEIGHT(ISPEC,JSAR)*
     &                           YKD(NSIW1D(ISPEC,JSAR),IPARTSAR)
            ENDIF
          END DO
        END DO
      END DO

!!    insure that the increment does no produce negative energy
      DO IPARTS=1,MAXMPART
        DO ISPEC=1,NSPECW
          EA (ISPEC,IPARTS) = MAX(EW(ISPEC,IPARTS)+EA(ISPEC,IPARTS),0.)
          IF(EA(ISPEC,IPARTS).GT.0.) THEN
            XKA(ISPEC,IPARTS) = XKW(ISPEC,IPARTS) + XKA(ISPEC,IPARTS)
            YKA(ISPEC,IPARTS) = YKW(ISPEC,IPARTS) + YKA(ISPEC,IPARTS)
          ELSE
            XKA(ISPEC,IPARTS) = XKW(ISPEC,IPARTS)
            YKA(ISPEC,IPARTS) = YKW(ISPEC,IPARTS)
          ENDIF
        END DO
      END DO

!     COMPUTE MEAN ANALYSED DIRECTIONS AND FREQUENCIES FROM 
!     WAVE NUMBERS.

      DO IPARTS=1,MAXMPART
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

!!!   this is not an optimal interpolation  !!!

!!! hardcoded to 12  !!!!!
      IFRECUT=12
 
      DO IFRE=1,IFRECUT
        DO IANG=1,NANG
          DO ISPEC=1,NSPECW
            IF(SPECW(ISPEC,IANG,IFRE).EQ.FLMIN)SPECW(ISPEC,IANG,IFRE)=0.
          END DO
        END DO
      END DO

      ALLOCATE(SPECDA(0:NSPEC,NANG,IFRECUT))

      DO IFRE=1,IFRECUT
        DO IANG=1,NANG
          DO IDSAR=1,NSPEC
            SPECDA(IDSAR,IANG,IFRE) =
     &      SIGMAINV*( SPECSA(IDSAR,IANG,IFRE)-SPECWA(IDSAR,IANG,IFRE) )
          END DO
        END DO
      END DO

      DO JSAR=1,JPSIW
        DO IPARTS=1,MAXMPART
          DO IFRE=1,IFRECUT
            DO IANG=1,NANG
              DO ISPEC=1,NSPECW
                IDSAR=NSIW1D(ISPEC,JSAR)
                IF(CORRELTAS(ISPEC,IPARTS,JSAR).LT.0 .AND.
     &             IDSAR.GT.0 .AND.
     &             PARTS(IDSAR,IANG,IFRE).EQ.IPARTS)
     &          THEN
                  SPECW(ISPEC,IANG,IFRE) = SPECW(ISPEC,IANG,IFRE) 
     &            + WEIGHT(ISPEC,JSAR)* SPECDA(IDSAR,IANG,IFRE)
                ENDIF

               END DO
             END DO
          END DO
        END DO
      END DO

      DO IFRE=1,IFRECUT
        DO IANG=1,NANG
          DO ISPEC=1,NSPECW
            SPECW(ISPEC,IANG,IFRE)=MAX(SPECW(ISPEC,IANG,IFRE),FLMIN)
          END DO
        END DO
      END DO

      DEALLOCATE(SPECDA)

      RETURN
      END SUBROUTINE OPTINT
