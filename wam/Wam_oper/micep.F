!-------------------------------------------------------------------

      SUBROUTINE MICEP (FIELD, IPARAM, ICMSK, KINF,KSUP,IJS,IJL,IG)

!-------------------------------------------------------------------

!**** *MICEP* - DEFINES ICEMASK OF ALL SEA POINTS.

!     R. PORTZ     MPI HAMBURG   OCTOBER 1992
!     J. BIDLOT    ECMWF         JUNE 1996    MESSAGE PASSING
!     P. JANSSEN   ECMWF         OCTOBER 1997 SIMPLIFY BY INTRODUCING
!                                             ICEMASK
!     J. BIDLOT    ECMWF         JANUARY 1998 : CORRECT THRESHOLD FOR
!                                               ICE TO BE LT 271.50
!     J. BIDLOT    ECMWF         FEBRUARY 2000 : USE OF SEA ICE FRACTION
!     J. BIDLOT    ECMWF         AUG 2006 :  


!     PURPOSE.
!     --------
!            SELECTS POINTS (TEMP < 271.5 KELVIN) IN THE TEMPERATURE
!            GRID AND DEFINE AS SEA ICE POINTS.
!            DOES THE SAME BUT WITH THE SEA ICE FRACTION WHERE
!            ANY POINT WITH SEA ICE FRACTION > 0.3 IS CONSIDERED ICE.

!**   INTERFACE
!     ---------
!             *CALL MICEP* *(FIELD, IPARAM, ICMSK, KINF,KSUP,IJS,IJL,IG)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *FIELD*     REAL      GLOBAL FIELD USED TO DETERMINE THE ICE MASK
!      *IPARAM*    INTEGER   GRIB PARAMETER OF FIELD
!      *ICMSK*     INTEGER   ICE MASK 
!      *KINF*      LOWER BOUND ON FIRST DIMENSION OF ICMSK 
!      *KSUP*      UPPER BOUND ON FIRST DIMENSION OF ICMSK 
!      *IJS*       INDEX OF FIRST GRIDPOINT
!      *IJL*       INDEX OF LAST GRIDPOINT
!      *IG*        BLOCK NUMBER


!     METHOD.
!     -------
!             NONE

!    EXTERNAL.
!    ---------
!             *NONE* 

!---------------------------------------------------------------------

      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY     ,NBLO
      USE YOWPCONS , ONLY : ZMISS
      USE YOWTEST  , ONLY : IU06     ,ITEST 

! ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NGX,NGY) :: FIELD
      INTEGER, DIMENSION (KINF:KSUP,NBLO) :: ICMSK 

! ----------------------------------------------------------------------

!     1. INITIALIZATION
!     ------------------

!-----------------------------------------------------------------

!*    2. WE DEFINE ICE FOR THOSE POINTS WHERE SST < 271.5 KELVIN). 
!        OR WITH SEA ICE FRACTION > 0.3.
!        AT THOSE POINTS ICE MASK = 0
!        -------------------------------------------------------

      NICE = 0

      IF (IPARAM.EQ.31) THEN
!       sea ice fraction is only defined between 0 and 1
        DO IJ = IJS, IJL 
          IX = IXLG(IJ,IG)
          IY = NGY- KXLT(IJ,IG) +1
          IF (FIELD(IX,IY) .GT. 0.3 .AND. FIELD(IX,IY) .LE. 1.01  .AND.
     &        FIELD(IX,IY) .NE. ZMISS ) THEN
            NICE = NICE + 1
            ICMSK(IJ,IG) = 0
          ELSE
            ICMSK(IJ,IG) = 1
          ENDIF
        ENDDO
      ELSE IF(IPARAM.EQ.139) THEN
        DO IJ = IJS, IJL 
          IX = IXLG(IJ,IG)
          IY = NGY- KXLT(IJ,IG) +1
            IF (FIELD(IX,IY) .LT. 271.5) THEN
              NICE = NICE + 1
              ICMSK(IJ,IG) = 0
            ELSE
              ICMSK(IJ,IG) = 1
            ENDIF
        ENDDO
      ENDIF 

      IF(ITEST.GE.1) THEN
        IF(IPARAM.EQ.139) THEN
          WRITE(IU06,*)' SUB. MICEP: POINTS WITH TEMP < 271.5 K :',
     &     NICE,' ON PE ',IRANK, ' OUT OF A TOTAL OF ',NPROC, 'PEs'
        ELSEIF (IPARAM.EQ.31) THEN
          WRITE(IU06,*)' SUB. MICEP: ',
     &           'POINTS WITH SEA ICE FRACTION > 0.3: ',
     &     NICE,' ON PE ',IRANK, ' OUT OF A TOTAL OF ',NPROC, 'PEs'
        ELSE
          WRITE(IU06,*)' ********************************************'
          WRITE(IU06,*)' SUB. MICEP: PARAMETER FOR ICE MASK NOT KNOWN'
          WRITE(IU06,*)' SUB. MICEP: NO ICE MASK WAS SET UP'
          WRITE(IU06,*)' ********************************************'
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE MICEP
