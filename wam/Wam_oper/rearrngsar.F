       SUBROUTINE REARRNGSAR(TH0) 
 
!---------------------------------------------------------------------

!     PURPOSE
!     -------
!     TO ADJUST THE SAR SPECTRA TO THE FIRST DIRECTION OF THE MODEL
!     (TH0) AND TO FETCH THE MODEL SPECTRA AND MODEL WINDS THE CLOSEST
!     TO SAR POINT LOCATIONS. 

!     AUTHOR
!     ------
!     J. BIDLOT  ECMWF JUNE 1999.

!     EXTERNALS
!     ---------
!     *FINDB*   

!**********************************************************************

!     MODULES:
!     --------
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : CICOVER  ,CITHRSH
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NPROC    ,IRANK    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI
      USE YOWSPEC,   ONLY : NSTART   ,NEND     ,U10OLD   ,USOLD    ,
     &            THWOLD   ,FL1      ,FL2
      USE YOWSARAS , ONLY : NSPEC    ,IGSAR    ,IJSAR    ,SPEC     ,
     &            LONG     ,LAT      ,U10      ,USSAR    ,THW
      USE YOWSTAT  , ONLY : LALTAS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU11

!----------------------------------------------------------------------

      IMPLICIT NONE

      REAL :: TH0 ! FIRST DIRECTION IN WAM SPECTRUM.

!     VARIABLES
!     ---------

      REAL, ALLOCATABLE :: FL(:,:)
      REAL :: FTH, ADIF, BDIF
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISPECPE 
      INTEGER :: ISPEC, K, M, KC, KC1, INC
      INTEGER :: JBL,JP,JFRE,JANG,NSPECPE

!     ROTATE INVERTED SAR SPECTRA TO THE FIRST DIRECTION OF WAM (TH0)
!     --------------------------------------------------------------- 
!     !!!! NOTE IT IS ASSUMED THAT INPUT INVERTED SPECTRA ARE DEFINED
!          WITH 0 AS THEIR FIRST DIRECTION.

      FTH = MOD(-TH0+ZPI,ZPI)
      FTH = FTH * REAL(NANG) / ZPI
      INC = INT(FTH)
      ADIF = FTH - INC
      BDIF = 1. - ADIF
      ALLOCATE(FL(NANG,NFRE))

      DO ISPEC = 1,NSPEC
        DO K=1,NANG
          KC  = K  - INC
          IF (KC .LT. 1) KC  = KC + NANG
          KC1 = KC - 1
          IF (KC1.LT. 1) KC1 = KC1+ NANG
          DO M=1,NFRE
            FL(K,M) = BDIF*SPEC(ISPEC,KC,M,2)+ADIF*SPEC(ISPEC,KC1,M,2)
          ENDDO
        ENDDO
        DO M=1,NFRE
          DO K=1,NANG
            SPEC(ISPEC,K,M,2) = FL(K,M)
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(FL)

!     GET WAM SPECTRA AT OBSERVATION LOCATION
!     ---------------------------------------

      ALLOCATE(IGSAR(NSPEC))
      ALLOCATE(IJSAR(NSPEC))
      ALLOCATE(ISPECPE(NSPEC))
      ISPECPE=0

!     FIND BLOCK POSITION OF SAR OBSERVATION POINTS

      CALL FINDB(NSPEC,NSPEC,LAT(1,2),LONG(1,2),IGSAR,IJSAR,
     &           NPROC,NSTART,NEND,IRANK)

!     REPLACE FIRST GUESS SPECTRA USED IN THE SAR INVERSION BY WAM
!     SPECTRA FOUND AT THE CLOSEST GRID POINT

      NSPECPE=0

!     READ FIRST BLOCK
      IF(IGL.GT.1) THEN
        REWIND IU11
        READ(IU11) FL1
      ENDIF
      DO JBL=1,IGL
!       ODD BLOCKS
        IF(MOD(JBL,2).EQ.1) THEN
          IF(JBL.NE.IGL) READ(IU11) FL2
          DO ISPEC = 1, NSPEC 
            JP=IJSAR(ISPEC)
            IF(JP.NE.0) THEN
              IF(IGSAR(ISPEC).EQ.JBL) THEN
                IF(CICOVER(JP,JBL).LE.CITHRSH) THEN
                  DO JFRE=1,NFRE
                    DO JANG=1,NANG
                      SPEC(ISPEC,JANG,JFRE,1)=FL1(JP,JANG,JFRE)
                    ENDDO
                  ENDDO
                  LAT(ISPEC,1)=AMOSOP+REAL(KXLT(JP,JBL)-1)*XDELLA
                  LONG(ISPEC,1)=AMOWEP+
     &                        REAL(IXLG(JP,JBL)-1)*ZDELLO(KXLT(JP,JBL))
                ELSE
                  DO JFRE=1,NFRE
                    DO JANG=1,NANG
                      SPEC(ISPEC,JANG,JFRE,1)=SPEC(ISPEC,JANG,JFRE,2)
                    ENDDO
                  ENDDO
                  LAT(ISPEC,1)=LAT(ISPEC,2)
                  LONG(ISPEC,1)=LONG(ISPEC,2)
                ENDIF
                NSPECPE=NSPECPE+1
                ISPECPE(NSPECPE)=ISPEC
              ENDIF
            ELSE
              DO JFRE=1,NFRE
                DO JANG=1,NANG
                  SPEC(ISPEC,JANG,JFRE,1)=SPEC(ISPEC,JANG,JFRE,2)
                ENDDO
              ENDDO
              LAT(ISPEC,1)=LAT(ISPEC,2)
              LONG(ISPEC,1)=LONG(ISPEC,2)
            ENDIF

          ENDDO
!       EVEN BLOCKS
        ELSE
          IF(JBL.NE.IGL) READ(IU11) FL1
          DO ISPEC = 1, NSPEC 
            JP=IJSAR(ISPEC)
            IF(CICOVER(MAX(JP,1),JBL).LE.CITHRSH .AND. JP.NE.0) THEN
              IF(IGSAR(ISPEC).EQ.JBL) THEN
                DO JFRE=1,NFRE
                  DO JANG=1,NANG
                    SPEC(ISPEC,JANG,JFRE,1)=FL2(JP,JANG,JFRE)
                  ENDDO
                ENDDO
                LAT(ISPEC,1)=AMOSOP+REAL(KXLT(JP,JBL)-1)*XDELLA
                LONG(ISPEC,1)=AMOWEP+
     &                        REAL(IXLG(JP,JBL)-1)*ZDELLO(KXLT(JP,JBL))
                NSPECPE=NSPECPE+1
                ISPECPE(NSPECPE)=ISPEC
              ENDIF
            ELSE
              DO JFRE=1,NFRE
                DO JANG=1,NANG
                  SPEC(ISPEC,JANG,JFRE,1)=SPEC(ISPEC,JANG,JFRE,2)
                ENDDO
              ENDDO
              LAT(ISPEC,1)=LAT(ISPEC,2)
              LONG(ISPEC,1)=LONG(ISPEC,2)
            ENDIF
          ENDDO
        ENDIF
      ENDDO

! REPLACE WINDS FROM SAR INPUT BY MODEL WINDS

      ALLOCATE(U10(NSPEC))
      ALLOCATE(USSAR(NSPEC))
      ALLOCATE(THW(NSPEC))

      DO JBL=1,IGL
        DO ISPEC = 1, NSPEC 
          JP=IJSAR(ISPEC)
          IF(IGSAR(ISPEC).EQ.JBL) THEN
            U10(ISPEC) = U10OLD(JP,JBL)
            USSAR(ISPEC) = USOLD(JP,JBL)
            THW(ISPEC) = THWOLD(JP,JBL)
          ELSE
            U10(ISPEC) = 0.
            USSAR(ISPEC) = 0.
            THW(ISPEC) = 0.
          ENDIF
        ENDDO
      ENDDO

!     EXCHANGE THE DIFFERENT CONTRIBUTIONS OF THE MODEL VALUES AT SAR
!     LOCATION WITH ALL THE OTHER PE'S.

      IF(LMESSPASS) THEN
        IF(ITEST.GT.0)
     &  WRITE(IU06,*)'-----------MPEXCHNGSARIN BEGINS------------'
        CALL MPEXCHNGSARIN(NSPECPE,ISPECPE,KTAG)
        IF(ITEST.GT.0)
     &  WRITE(IU06,*)'-----------MPEXCHNGSARIN ENDS------------'
      ENDIF

      DEALLOCATE(ISPECPE)

      RETURN
      END SUBROUTINE REARRNGSAR 
