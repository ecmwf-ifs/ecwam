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
     &            NFLAG    ,IPFGTBL  ,NWRTOUTWAM, NIPRMOUT, ITOBOUT,
     &            NTRAIN   ,LLPARTITION, IFRSTPARTI
      USE YOWMPP   , ONLY : NPROC
      USE YOWSTAT  , ONLY : IASSI


! ----------------------------------------------------------------------
      IMPLICIT NONE
#if defined MODEL_COUPLING_ATM_WAV || defined NETCDF_OUTPUT_WAM
      LOGICAL :: COUPLFLAG(JPPFLAG+1)
#endif

      INTEGER :: IR, IFLAG, ITT

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
      NIPRMOUT=0
      IR=1
      DO IFLAG=1,JPPFLAG
        IF((IASSI.EQ.1).AND.
     &     ((IFLAG.EQ.1).OR.(IFLAG.EQ.7))) THEN
          IPFGTBL(IFLAG)=IR
          NIPRMOUT=NIPRMOUT+1
          ITOBOUT(IFLAG)=NIPRMOUT
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
          NIPRMOUT=NIPRMOUT+1
          ITOBOUT(IFLAG)=NIPRMOUT
          IR=IR+NWRTOUTWAM
          IF(IR.GT.NPROC) IR=1
        ELSE
          IPFGTBL(IFLAG)=0
          ITOBOUT(IFLAG)=0
        ENDIF
      ENDDO


!     WILL THERE BE OUTPUT OF PARTITONED PARAMETERS
      LLPARTITION=.FALSE.
      IFLAG=IFRSTPARTI
      DO ITT=1,3*NTRAIN
        IF(IPFGTBL(IFLAG).NE.0) THEN
          LLPARTITION=.TRUE.
          EXIT
        ENDIF
        IFLAG=IFLAG+1
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
