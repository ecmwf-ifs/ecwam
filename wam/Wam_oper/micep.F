!-------------------------------------------------------------------

      SUBROUTINE MICEP (IU06, NSTART, NEND, FIELD, IPARAM)

!-------------------------------------------------------------------

!**** *MICEP* - DEFINES ICEMASK OF ALL SEA POINTS.

!     R. PORTZ     MPI HAMBURG   OCTOBER 1992
!     J. BIDLOT    ECMWF         JUNE 1996    MESSAGE PASSING
!     P. JANSSEN   ECMWF         OCTOBER 1997 SIMPLIFY BY INTRODUCING
!                                             ICEMASK
!     J. BIDLOT    ECMWF         JANUARY 1998 : CORRECT THRESHOLD FOR
!                                               ICE TO BE LT 271.50
!     J. BIDLOT    ECMWF         FEBRUARY 2000 : USE OF SEA ICE FRACTION


!     PURPOSE.
!     --------
!            SELECTS POINTS (TEMP < 271.5 KELVIN) IN THE TEMPERATURE
!            GRID AND DEFINE AS SEA ICE POINTS.
!            DOES THE SAME BUT WITH THE SEA ICE FRACTION WHERE
!            ANY POINT WITH SEA ICE FRACTION > 0.3 IS CONSIDERED ICE.

!**   INTERFACE
!     ---------
!             *CALL MICEP* *( IU06, NSTART, NEND, FIELD, IPARAM)*

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   PRINTER OUTPUT UNIT.
!      *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!      *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!      *FIELD*     REAL      GLOBAL FIELD USED TO DETERMINE THE ICE MASK
!      *IPARAM*    INTEGER   GRIB PARAMETER OF FIELD

!     METHOD.
!     -------
!             NONE

!    EXTERNAL.
!    ---------
!             *NONE* 

!---------------------------------------------------------------------

      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : ICEMASK
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWSTAT  , ONLY : IASSI

! ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NGX,NGY) :: FIELD

! ----------------------------------------------------------------------

!     1. INITIALIZATION
!     ------------------

!-----------------------------------------------------------------

!*    2. WE DEFINE ICE FOR THOSE POINTS WHERE SST < 271.5 KELVIN). 
!        OR WITH SEA ICE FRACTION > 0.3.
!        AT THOSE POINTS ICEMASK = 0
!        -------------------------------------------------------

      NICE = 0
      ICEMASK = 1

      DO IG=1,IGL
        IF (LMESSPASS) THEN
!           ICEMASK FOR ALL SEAPOINTS IS NEEDED IN GRDATA
          IF(IASSI.EQ.1) THEN
            IJSTART=NSTART(1)
            IJEND=NEND(NPROC)
          ELSE
            IJSTART=NSTART(IRANK)
            IJEND=NEND(IRANK)
          ENDIF
        ELSE
          IJSTART=IJS(IG)
          IJEND=IJL(IG)
        ENDIF

        DO IJ = IJSTART, IJEND 
          IX = IXLG(IJ,IG)
          IY = NGY- KXLT(IJ,IG) +1
          IF(IPARAM.EQ.139) THEN
            IF (FIELD(IX,IY) .LT. 271.5) THEN
              NICE = NICE + 1
              ICEMASK(IJ,IG) = 0
            ENDIF
          ELSEIF (IPARAM.EQ.31) THEN
!           sea ice fraction is only defined between 0 and 1
            IF (FIELD(IX,IY) .GT. 0.3 .AND. FIELD(IX,IY) .LE. 1.01) THEN
              NICE = NICE + 1
              ICEMASK(IJ,IG) = 0
            ENDIF
          ENDIF 
        ENDDO
      ENDDO

      IF(IASSI.EQ.1) THEN
        IF(IPARAM.EQ.139) THEN
          WRITE(IU06,*)' SUB. MICEP: POINTS WITH TEMP < 271.5 K : ',
     &             NICE,' ON ALL ',NPROC, 'PEs'
        ELSEIF (IPARAM.EQ.31) THEN
          WRITE(IU06,*)' SUB. MICEP: ',
     &             'POINTS WITH SEA ICE FRACTION > 0.3: ',
     &             NICE,' ON ALL ',NPROC, 'PEs'
        ELSE
          WRITE(IU06,*)' ********************************************'
          WRITE(IU06,*)' SUB. MICEP: PARAMETER FOR ICE MASK NOT KNOWN'
          WRITE(IU06,*)' SUB. MICEP: NO ICEMASK WAS SET UP'
          WRITE(IU06,*)' ********************************************'
        ENDIF
      ELSE
        IF(IPARAM.EQ.139) THEN
          WRITE(IU06,*)' SUB. MICEP: POINTS WITH TEMP < 271.5 K :',
     &       NICE,' ON PE ',IRANK, ' OUT OF A TOTAL OF ',NPROC, 'PEs'
        ELSEIF (IPARAM.EQ.31) THEN
          WRITE(IU06,*)' SUB. MICEP: ',
     &             'POINTS WITH SEA ICE FRACTION > 0.3: ',
     &       NICE,' ON PE ',IRANK, ' OUT OF A TOTAL OF ',NPROC, 'PEs'
        ELSE
          WRITE(IU06,*)' ********************************************'
          WRITE(IU06,*)' SUB. MICEP: PARAMETER FOR ICE MASK NOT KNOWN'
          WRITE(IU06,*)' SUB. MICEP: NO ICEMASK WAS SET UP'
          WRITE(IU06,*)' ********************************************'
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE MICEP
