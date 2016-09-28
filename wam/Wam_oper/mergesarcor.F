       SUBROUTINE MERGESARCOR

!---------------------------------------------------------------------

!     PURPOSE
!     -------

!     MERGE ANALYSED SPECTRA INTO BLOCKED GRID SPECTRA.
!     MERGE ANALYSED WINDS INTO WIND BLOCKED ARRAYS.


!     AUTHOR: RENATE BROKOPF, MPM HAMBURG, GERMANY. MAY 1994.
!             JEAN BIDLOT ECMWF MAY 1999.
!     

!     EXTERNALS
!     ---------
       

!***********************************************************************

!     MODULES:
!     --------

      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL      ,IJL2
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO    ,NBLO
      USE YOWSARAS , ONLY : NW1D     ,SPECW    ,UW10     ,THWW     ,
     &              NSPECW
      USE YOWSPEC,   ONLY : FL1
      USE YOWTEST  , ONLY : IU06
      USE YOWUNIT  , ONLY : IU11     ,IU12     ,IU13     ,IU14
      USE YOWWIND  , ONLY : WSPMIN

!-----------------------------------------------------------------------

!     INTERFACE
!     ---------

      IMPLICIT NONE

!     VARIABLES
!     ---------

      INTEGER :: JBL,JP,JFRE,JANG,IFIL

      REAL :: UST(2), USTM
     
!     
!-----------------------------------------------------------------------

!      1. UPDATING SPECTRA 
!         ----------------

      IF(NSPECW.GT.0) THEN

        IF(IGL.GT.1) THEN
          REWIND IU11
          REWIND IU12
          REWIND IU13
          REWIND IU14
          IFIL=IU11
          IU11=IU12
          IU12=IFIL
          IFIL=IU13
          IU13=IU14
          IU14=IFIL
          READ(IU11) FL1
        ENDIF

        JBL=1

            DO JFRE=1,NFRE
              DO JANG=1,NANG
                DO JP = IJS(JBL),IJL(JBL)
                  IF(NW1D(JP,JBL).NE.0)
     &            FL1(JP,JANG,JFRE)=SPECW(NW1D(JP,JBL),JANG,JFRE)       
                ENDDO
              ENDDO
            ENDDO


        IF(IGL.GT.1) THEN
          CLOSE(IU11)
          CLOSE(IU12)
        ENDIF
      ENDIF


!     2. UPDATING WINDS
!        --------------

!     no longer used.

      RETURN
      END SUBROUTINE MERGESARCOR
