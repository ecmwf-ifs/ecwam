      SUBROUTINE MPCRTBL 

! ----------------------------------------------------------------------

!**** *MPCRTBL* -

!     J. BIDLOT     ECMWF    MARCH 1996    MESSAGE PASSING


!*    PURPOSE.
!     --------
!     CREATE A TABLE THAT ASSOCIATES EACH INTEGRATED OUTPUT PARAMETER
!     WITH A SPECIFIC PROCESS. THIS TABLE IS USED TO SPREAD THE
!     RETRIEVAL OF THOSE PARAMETER FIELD IN GRID FORMAT ON AS MANY 
!     PROCESSOR AS POSSIBLE

!**   INTERFACE.
!     ----------


!     METHOD.
!     -------
!     SYSTEMATICALY ASSIGN A PROCESS TO A PARAMETER AND START AGAIN WHEN
!     RUN OUT OF PROCESS.

!     EXTERNALS.
!     ----------
!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : JPPFLAG  ,FFLAG    ,PFLAG    ,GFLAG    ,
     &            NFLAG    ,IPFGTBL  ,NWRTOUTWAM, COUTNAME
      USE YOWMPP   , ONLY : NPROC
      USE YOWSTAT  , ONLY : IASSI


! ----------------------------------------------------------------------
      IMPLICIT NONE
#if defined MODEL_COUPLING_ATM_WAV || defined NETCDF_OUTPUT_WAM
      LOGICAL :: COUPLFLAG(JPPFLAG+1)
#endif

      INTEGER :: IR, IFLAG, IC, ITG

!     1. CREATE TABLE MAPPING OUTPUT INTEGRATED PARAMETER WITH PE RANK
!        -------------------------------------------------------------

#if defined MODEL_COUPLING_ATM_WAV || defined NETCDF_OUTPUT_WAM
      COUPLFLAG=.FALSE.
      COUPLFLAG(1)=.TRUE. ! HS
      COUPLFLAG(2)=.TRUE. ! WDIRTG (Mean wave direction)
      COUPLFLAG(3)=.TRUE.
      COUPLFLAG(6)=.TRUE. ! WPKFTG (wave peak frequency)
      COUPLFLAG(37)=.TRUE.
#endif
#if defined MODEL_COUPLING_OCN_WAV
      COUPLFLAG(27)=.TRUE.   ! DSPR
      COUPLFLAG(25)=.TRUE.   ! T01
      COUPLFLAG(26)=.TRUE.   ! T02
          ! below are for the wave supported stress
      COUPLFLAG(46)=.TRUE.   ! TAUOC
      COUPLFLAG(5)=.TRUE.    ! TWINDDIR
#endif

!     OUTPUT VARIABLE NAMES
      COUTNAME(:) =' VARIABLE NOT DEFINED, see MPCRTBL'
      COUTNAME( 1)=' SIGNIFICANT WAVE HEIGHT .........'
      COUTNAME( 2)=' MEAN WAVE DIRECTION .............'
      COUTNAME( 3)=' WAVE MEAN PERIOD (-1) ...........'
      COUTNAME( 4)=' FRICTION VELOCITY ...............'
      COUTNAME( 5)=' WAVE MODEL WIND DIRECTION........'
      COUTNAME( 6)=' WAVE PEAK FREQUENCY .............'
      COUTNAME( 7)=' DRAG COEFFICIENT ................'
      COUTNAME( 8)=' NORMALISED WAVE STRESS ..........'
      COUTNAME( 9)=' MEAN SQUARE SLOPE................'
      COUTNAME(10)=' WAVE MODEL WIND SPEED............'
      COUTNAME(11)=' WIND SEA WAVE HEIGHT.............'
      COUTNAME(12)=' TOTAL SWELL WAVE HEIGHT..........'
      COUTNAME(13)=' WIND SEA DIRECTION...............'
      COUTNAME(14)=' TOTAL SWELL WAVE DIRECTION.......'
      COUTNAME(15)=' WIND SEA MEAN PERIOD (-1) .......'
      COUTNAME(16)=' TOTAL SWELL MEAN PERIOD (-1) ....'
      COUTNAME(22)=' ALTIMETER WAVE HEIGHT ...........'
      COUTNAME(23)=' CORRECTED ALT WAVE HEIGHT .......'
      COUTNAME(24)=' ALTIMETER RANGE CORRECTION ......'
      COUTNAME(25)=' WAVE MEAN PERIOD (1) ............'
      COUTNAME(26)=' WAVE MEAN PERIOD (2) ............'
      COUTNAME(27)=' WAVE DIRECTIONAL SPREAD .........'
      COUTNAME(28)=' WIND SEA MEAN PERIOD (1) ........'
      COUTNAME(29)=' SWELL MEAN PERIOD (1) ...........'
      COUTNAME(30)=' WIND SEA MEAN PERIOD (2) ........'
      COUTNAME(31)=' SWELL MEAN PERIOD (2) ...........'
      COUTNAME(32)=' WIND SEA DIRECTIONAL SPREAD......'
      COUTNAME(33)=' SWELL DIRECTIONAL SPREAD ........'
      COUTNAME(34)=' WAVE SPECTRAL KURTOSIS ..........'
      COUTNAME(35)=' BENJAMIN-FEIR INDEX .............'
      COUTNAME(36)=' WAVE SPECTRAL PEAKEDNESS ........'
      COUTNAME(37)=' BATHYMETRY ......................'
      COUTNAME(38)=' MAXIMUM WAVE HEIGHT .............'
      COUTNAME(39)=' MAXIMUM WAVE PERIOD .............'
      COUTNAME(40)=' U-COMP SURFACE STOKES DRIFT .....'
      COUTNAME(41)=' V-COMP SURFACE STOKES DRIFT .....'
      COUTNAME(42)=' U-COMP SURFACE CURRENT ..........'
      COUTNAME(43)=' V-COMP SURFACE CURRENT ..........'
      COUTNAME(44)=' NORMALISED ENERGY FLUX TO OCEAN .'
      COUTNAME(45)=' NORMALISED ENERGY FLUX TO WAVES .'
      COUTNAME(46)=' NORMALISED MOMENTUM FLUX TO OCEAN'
      COUTNAME(47)=' SWELL PARTITION 1 WAVE HEIGHT ...'
      COUTNAME(48)=' SWELL PARTITION 1 DIRECTION .....'
      COUTNAME(49)=' SWELL PARTITION 1 MEAN PERIOD (1)'
      COUTNAME(50)=' SWELL PARTITION 2 WAVE HEIGHT ...'
      COUTNAME(51)=' SWELL PARTITION 2 DIRECTION .....'
      COUTNAME(52)=' SWELL PARTITION 2 MEAN PERIOD (1)'
      COUTNAME(53)=' SWELL PARTITION 3 WAVE HEIGHT ...'
      COUTNAME(54)=' SWELL PARTITION 3 DIRECTION .....'
      COUTNAME(55)=' SWELL PARTITION 3 MEAN PERIOD (1)'
      COUTNAME(56)=' MEAN SQUARE STRAIN IN ICE .......'
      COUTNAME(57)=' WAVE HEIGHT WITH PERIOD > 10s ...'
      COUTNAME(58)=' SURFACE AIR DENSITY .............'
      COUTNAME(59)=' CONVECTIVE VELOCITY SCALE........'
      COUTNAME(60)=' SEA ICE COVER ...................'
      COUTNAME(61)=' SEA ICE THICKNESS ...............'
      COUTNAME(62)=' SPECTRAL SKWENESS ...............'
      COUTNAME(63)=' NEMO SST ........................'
      COUTNAME(64)=' NEMO SEA ICE COVER ..............'
      COUTNAME(65)=' NEMO SEA ICE THICKNESS ..........'
      COUTNAME(66)=' NEMO ZONAL CURRENT ..............'
      COUTNAME(67)=' NEMO MERIDIONAL CURRENT .........'
      COUTNAME(68)=' WAVE ENERGY FLUX MAGNITUDE ......'
      COUTNAME(69)=' WAVE ENERGY FLUX DIRECTION ......'
      COUTNAME(70)=' SIG. WAVE HEIGHT 10<=T<=12 ......'
      COUTNAME(71)=' SIG. WAVE HEIGHT 12<=T<=14 ......'
      COUTNAME(72)=' SIG. WAVE HEIGHT 14<=T<=17 ......'
      COUTNAME(73)=' SIG. WAVE HEIGHT 17<=T<=21 ......'
      COUTNAME(74)=' SIG. WAVE HEIGHT 21<=T<=25 ......'
      COUTNAME(75)=' SIG. WAVE HEIGHT 25<=T<=30 ......'

      DO IC=1,5
        ITG=JPPFLAG-5+IC
        COUTNAME(ITG)=' EXTRA FIELD .....................'
      ENDDO

      IR=1
      DO IFLAG=1,JPPFLAG
        IF((IASSI.EQ.1).AND.
     &     ((IFLAG.EQ.1).OR.(IFLAG.EQ.7))) THEN
          IPFGTBL(IFLAG)=IR
          IR=IR+NWRTOUTWAM
          IF(IR.GT.NPROC) IR=1
#if defined MODEL_COUPLING_ATM_WAV || defined NETCDF_OUTPUT_WAM || defined MODEL_COUPLING_OCN_WAV
        ELSEIF(PFLAG(IFLAG).OR.FFLAG(IFLAG).OR.GFLAG(IFLAG).OR.
     &         NFLAG(IFLAG).OR.COUPLFLAG(IFLAG) ) THEN
#else
        ELSEIF(PFLAG(IFLAG).OR.FFLAG(IFLAG).OR.GFLAG(IFLAG).OR.
     &         NFLAG(IFLAG) ) THEN
#endif
          IPFGTBL(IFLAG)=IR
          IR=IR+NWRTOUTWAM
          IF(IR.GT.NPROC) IR=1
        ELSE
          IPFGTBL(IFLAG)=0
        ENDIF
      ENDDO


! FIND P.E. FOR OUTPUT OF RESTART FILES

      IFLAG=JPPFLAG
!c
!cc set to 1 in order to simplify ouput to FDB
      IPFGTBL(JPPFLAG+1)=1

!!!     DO WHILE ((IPFGTBL(IFLAG).EQ.0).AND.(IFLAG.GT.1))
!!!         IFLAG=IFLAG-1
!!!      END DO
!!!      IPFGTBL(JPPFLAG+1)=IPFGTBL(IFLAG)+1
!!!      IF(IPFGTBL(JPPFLAG+1).GT.NPROC)IPFGTBL(JPPFLAG+1)=1

      END SUBROUTINE MPCRTBL
