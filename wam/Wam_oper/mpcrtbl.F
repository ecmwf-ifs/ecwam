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
     &            NFLAG    ,IPFGTBL
      USE YOWMPP   , ONLY : NPROC
      USE YOWSTAT  , ONLY : IASSI

! ----------------------------------------------------------------------

!     1. CREATE TABLE MAPPING OUTPUT INTEGRATED PARAMETER WITH PE RANK
!        -------------------------------------------------------------

      IR=1
      DO IFLAG=1,JPPFLAG
        IF((IASSI.EQ.1).AND.
     &     ((IFLAG.EQ.1).OR.(IFLAG.EQ.7))) THEN
          IPFGTBL(IFLAG)=IR
          IR=IR+1
          IF(IR.GT.NPROC) IR=1
        ELSEIF(PFLAG(IFLAG).OR.FFLAG(IFLAG).OR.GFLAG(IFLAG).OR.
     &         NFLAG(IFLAG) ) THEN
          IPFGTBL(IFLAG)=IR
          IR=IR+1
          IF(IR.GT.NPROC) IR=1
        ELSE
          IPFGTBL(IFLAG)=0
        ENDIF
      ENDDO

! FIND P.E. FOR THE GATHERING OF THE WIND STRESS FIELD
! WHENEVER ASSIMILATION IS DONE BUT NOT OUTPUT OF IT  WAS REQUESTED.

      IF((IASSI.EQ.1).AND.(IPFGTBL(4).EQ.0)) THEN 
        IFLAG=JPPFLAG
        DO WHILE ((IPFGTBL(IFLAG).EQ.0).AND.(IFLAG.GT.1))
          IFLAG=IFLAG-1
        ENDDO
        IPFGTBL(4)=IPFGTBL(IFLAG)+1
        IF(IPFGTBL(4).GT.NPROC)IPFGTBL(4)=1
      ENDIF

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
