!-------------------------------------------------------------------

      SUBROUTINE MICEP (IU06, NSTART, NEND, TEMPG)

!-------------------------------------------------------------------

!**** *MICEP* - DEFINE BLOCK AND INDEX OF EACH ICE POINT.

!     R. PORTZ     MPI HAMBURG   OCTOBER 1992
!     J. BIDLOT    ECMWF         JUNE 1996    MESSAGE PASSING
!     P. JANSSEN   ECMWF         OCTOBER 1997 SIMPLIFY BY INTRODUCING
!                                             ICEMASK
!     J. BIDLOT    ECMWF           JANURAY 1998 : CORRECT THRESHOLD FOR
!                                                 ICE TO BE LT 271.50


!     PURPOSE.
!     --------
!            SELECT POINTS (TEMP < 271.5 KELVIN) IN THE TEMPERATURE
!            GRID AND DEFINE AS SEA ICE POINTS.



!**   INTERFACE
!     ---------
!             *CALL MICEP* *( IU06, NICE, NSTART, NEND)*

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   PRINTER OUTPUT UNIT.
!      *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!      *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!      *TEMPG*     REAL      GLOBAL TEMPRERATURES (KELVIN).


!     METHOD.
!     -------
!             NONE


!    EXTERNAL.
!    ---------
!             *NONE* 

!---------------------------------------------------------------------

      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : ICEMASK
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWSTAT  , ONLY : IASSI
      USE YOWTEMP  , ONLY : NXT      ,NYT      ,TMOWEP   ,TMONOP   ,
     &            DELLA    ,DELLO

! ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NXT,NYT) :: TEMPG

! ----------------------------------------------------------------------

!     1. INITIALIZATION
!     ------------------

!-----------------------------------------------------------------

!*    2. WE DEFINE ICE FOR THOSE POINTS WHERE SST < 271.5 KELVIN). 
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
          DO IJ = IJSTART, IJEND 
            IX = IXLG(IJ,IG)
            KX = KXLT(IJ,IG)
            XI = AMOWEP + REAL(IX-1)*ZDELLO(KX)-TMOWEP
            XI = MOD(XI+720.,360.)
            XI = XI/DELLO+1.00001
            XJ = TMONOP - AMOSOP - REAL(KX-1)*XDELLA
            XJ = XJ/DELLA+1.00001

            I  = NINT(XI)
            J  = NINT(XJ)
            IF (TEMPG(I,J) .LT. 271.5) THEN
              NICE = NICE + 1
              ICEMASK(IJ,IG) = 0
            ENDIF
          ENDDO
        ELSE
          DO IJ = IJS(IG),IJL(IG)
            IX = IXLG(IJ,IG)
            KX = KXLT(IJ,IG)
            XI = AMOWEP + REAL(IX-1)*ZDELLO(KX)-TMOWEP
            XI = MOD(XI+720.,360.)
            XI = XI/DELLO+1.00001
            XJ = TMONOP - AMOSOP - REAL(KX-1)*XDELLA
            XJ = XJ/DELLA+1.00001

            I  = NINT(XI)
            J  = NINT(XJ)
            IF (TEMPG(I,J) .LT. 271.5) THEN
              NICE = NICE + 1
              ICEMASK(IJ,IG) = 0
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      IF(IASSI.EQ.1) THEN
        WRITE(IU06,*)'   POINTS WITH TEMP < 271.5 K:',NICE,
     &             ' ON ALL ',NPROC, 'PEs'
      ELSE
        WRITE(IU06,*)'   POINTS WITH TEMP < 271.5 K:',NICE,
     &             ' ON PE ',IRANK, ' OUT OF A TOTAL OF ',NPROC, 'PEs'
      ENDIF

      RETURN
      END SUBROUTINE MICEP
