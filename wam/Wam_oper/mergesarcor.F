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
       
!     *SPLITBL*

!***********************************************************************

!     MODULES:
!     --------

      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL      ,IJL2
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO    ,NBLO
      USE YOWSARAS , ONLY : NW1D     ,SPECW    ,UW10     ,THWW
      USE YOWSPEC,   ONLY : U10OLD   ,THWOLD   ,FL1      ,FL2
      USE YOWTEST  , ONLY : IU06
      USE YOWUNIT  , ONLY : IU11     ,IU12     ,IU13     ,IU14

!-----------------------------------------------------------------------

!     INTERFACE
!     ---------

      IMPLICIT NONE

!     VARIABLES
!     ---------

      INTEGER JBL,JP,JFRE,JANG,IFIL
!     
!-----------------------------------------------------------------------

!      1. UPDATING SPECTRA 
!         ----------------

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

      DO JBL=1,IGL

!       ODD BLOCKS
        IF(MOD(JBL,2).EQ.1) THEN

          IF(JBL.NE.IGL) READ(IU11) FL2

          DO JFRE=1,NFRE
            DO JANG=1,NANG
              DO JP = IJS(JBL),IJL(JBL)
                IF(NW1D(JP,JBL).NE.0)
     &          FL1(JP,JANG,JFRE)=SPECW(NW1D(JP,JBL),JANG,JFRE)       
              ENDDO
            ENDDO
          ENDDO

          IF(JBL.NE.1) THEN
             CALL SPLITBL (FL1, IJS(JBL), IJL2(JBL), IU14,
     &                     NIBLO,NANG,NFRE)
             WRITE(IU12) FL1
          ENDIF

!       EVEN BLOCKS
        ELSE

          IF(JBL.NE.IGL) READ(IU11) FL1

          DO JFRE=1,NFRE
            DO JANG=1,NANG
              DO JP = IJS(JBL),IJL(JBL)
                IF(NW1D(JP,JBL).NE.0)
     &          FL2(JP,JANG,JFRE)=SPECW(NW1D(JP,JBL),JANG,JFRE)       
              ENDDO
            ENDDO
          ENDDO
          IF(JBL.NE.1) THEN
            CALL SPLITBL (FL2, IJS(JBL), IJL2(JBL), IU14,
     &                    NIBLO,NANG,NFRE)
             
          ENDIF
        ENDIF
      ENDDO

      IF(IGL.GT.1) THEN
        CLOSE(IU11)
        CLOSE(IU12)
      ENDIF


!     2. UPDATING WINDS
!        --------------

      DO JBL=1,IGL
         DO JP = IJS(JBL),IJL(JBL)
           IF(NW1D(JP,JBL).GT.0) THEN
             U10OLD(JP,JBL) = UW10(NW1D(JP,JBL))
             THWOLD(JP,JBL) = THWW(NW1D(JP,JBL))
           ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE MERGESARCOR
