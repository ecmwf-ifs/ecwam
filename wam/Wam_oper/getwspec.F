       SUBROUTINE GETWSPEC(SARCORLEN)
 
!---------------------------------------------------------------------

!     PURPOSE
!     -------
!     TO COLLECT THE MODEL SPECTRA INFLUENCED BY SAR DATA POINTS 
!     AND ORGANIZE THEM IN A ONE DIMENSIONAL ARRAY.

!     AUTHOR
!     ------
!     S.HASSELMANN, MPM, HAMBURG,GERMANY
!     R.BROKOPF DTO.
!     P.LIONELLO (UNIV. OF PADUA) , HAMBURG, 11/93
!     J. BIDLOT  ECMWF MAY 1999.

!     EXTERNALS
!     ---------
!       *RESIZE_GETWSPEC*

!**********************************************************************

!     MODULES:
!     --------
      USE YOWGRID  , ONLY : SINPH    ,COSPH    ,IGL      ,IJS      ,
     &            IJL
      USE YOWICE   , ONLY : ICEMASK
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,ZDELLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NBLO     ,NGY
      USE YOWPCONS , ONLY : RAD      ,DEG      ,R
      USE YOWSPEC,   ONLY : U10OLD   ,THWOLD   ,FL1      ,FL2
      USE YOWSARAS , ONLY : PPSARCF  ,SARCORDIS,NSPEC    ,NSPECW   ,
     &            JPSIW    ,NW1D     ,NTOSIW1D ,IGSAR    ,IJSAR    ,     
     &            LONG     ,LAT      ,COSSARLAT,SINSARLAT,
     &            NSIW1D   ,SPECW    ,DIST     ,
     &            UW10     ,THWW
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06
      USE YOWUNIT  , ONLY : IU11

!----------------------------------------------------------------------
!     INTERFACE
!     ---------
      IMPLICIT NONE

!           
!     VARIABLES
!     ---------
      INTEGER :: JBL,JP,JSAR,JFRE,JANG,ICOUNT,MAXJPSIW,LLAT,NSIW1D0
      INTEGER :: ICOUNTOLD,IC,ICC
      INTEGER,DIMENSION(NSPEC) :: IS,KS 
      INTEGER,ALLOCATABLE :: INDEX(:)

      REAL :: DELLON
      REAL :: SARCORLEN, XMAX, DEFAULT_DIST
      REAL,ALLOCATABLE,DIMENSION(:) :: ADS
     
      LOGICAL,ALLOCATABLE,DIMENSION(:) :: LOSAR 
!--------------------------------------------------------------------

!     1.0 ALLOCATE ARRAYS FROM MODULE
!     --------------------------------

      ALLOCATE(NW1D(NINF:NSUP,NBLO))

      NSPECW=INT((1.5/XDELLA)**2*4400./NPROC)
      JPSIW=10
      NSIW1D0=0

      ALLOCATE(NTOSIW1D(NSPECW))
      ALLOCATE(NSIW1D(NSPECW,JPSIW))
      NSIW1D=NSIW1D0
      ALLOCATE(DIST(NSPECW,JPSIW))

             
!     2.0 LABEL WAM SPECTRA THAT ARE INFLUENCED BY THE SAR 
!         OBSERVATIONS AND ASSOCIATE THEM TO SAR WITHIN A 
!         SARCORLEN DISTANCE (in radian).

      SARCORLEN=SARCORDIS/R
      XMAX=PPSARCF*SARCORLEN

!     MAXIMUN EXTEND OF XMAX IN NUMBER OF GRID POINTS 
      LLAT = NINT(DEG*XMAX/XDELLA)+1

      DEFAULT_DIST=10*XMAX

      DIST=DEFAULT_DIST
      NTOSIW1D=0
      MAXJPSIW=0
      NW1D=0
      ICOUNT=0

      DO JSAR=1,NSPEC
        IF(IJSAR(JSAR).NE.0) THEN
          IS(JSAR)=IXLG(IJSAR(JSAR),IGSAR(JSAR))
          KS(JSAR)=KXLT(IJSAR(JSAR),IGSAR(JSAR))
        ELSE
          IS(JSAR)=0
          KS(JSAR)=0
        ENDIF
      ENDDO


!     FIND POINTS WHICH FALL WITHIN XMAX FROM SAR OBSERVATIONS

      ALLOCATE(INDEX((2*LLAT+3)**2))

      DO JBL=1,IGL

        ALLOCATE(ADS(IJS(JBL):IJL(JBL)))
        ALLOCATE(LOSAR(IJS(JBL):IJL(JBL)))
        DO JP=IJS(JBL),IJL(JBL)
          LOSAR(JP)=.FALSE.
        ENDDO

        DO JSAR=1,NSPEC
!         exclude sar point that are obviously outside the sub region
          IF(KS(JSAR).GE.KXLT(IJS(JBL),JBL)-LLAT .AND.
     &       KS(JSAR).LE.KXLT(IJL(JBL),JBL)+LLAT .AND.
     &       IS(JSAR).GT.0 .AND. KS(JSAR).GT.0 ) THEN

            DO JP=IJS(JBL),IJL(JBL)
              IF(ICEMASK(JP,JBL).NE.0 .AND.
     &           ABS(KS(JSAR)-KXLT(JP,JBL)).LE.LLAT) THEN

                DELLON=RAD*(LONG(JSAR,2)
     &                     -REAL(IXLG(JP,JBL)-1)*ZDELLO(KXLT(JP,JBL)))
                ADS(JP)=COS(DELLON)*COSSARLAT(JSAR)*COSPH(KXLT(JP,JBL))+
     &                              SINSARLAT(JSAR)*SINPH(KXLT(JP,JBL))
                ADS(JP)=ACOS(MAX(MIN(ADS(JP),1.),-1.))
              ELSE
                ADS(JP)=DEFAULT_DIST
              ENDIF
            ENDDO

            IC=0
            DO JP=IJS(JBL),IJL(JBL)
              IF(ADS(JP).LE.XMAX.AND.LOSAR(JP))THEN
                IC=IC+1
                INDEX(IC)=NW1D(JP,JBL)
              ENDIF
            ENDDO
            DO ICC=1,IC
              NTOSIW1D(INDEX(ICC)) = NTOSIW1D(INDEX(ICC)) + 1
            ENDDO

            DO JP=IJS(JBL),IJL(JBL)
              IF(ADS(JP).LE.XMAX.AND.LOSAR(JP))THEN
                MAXJPSIW=MAX(MAXJPSIW,NTOSIW1D(NW1D(JP,JBL)))
              ENDIF
            ENDDO

            IC=0
            DO JP=IJS(JBL),IJL(JBL)
              IF(ADS(JP).LE.XMAX.AND..NOT.LOSAR(JP))THEN
                IC=IC+1
                INDEX(IC)=JP
              ENDIF
            ENDDO

            ICOUNTOLD=ICOUNT
            ICOUNT=ICOUNTOLD+IC
            IF(ICOUNT.GT.NSPECW)
     &         CALL RESIZE_GETWSPEC(ICOUNT, JPSIW, DEFAULT_DIST,
     &                              NSIW1D0)
            DO ICC=1,IC
              JP=INDEX(ICC)
              LOSAR(JP)=.TRUE.
              NW1D(JP,JBL)=ICOUNTOLD+ICC
              NTOSIW1D(ICOUNTOLD+ICC)=1
            ENDDO

            IF(MAXJPSIW.GT.JPSIW)
     &         CALL RESIZE_GETWSPEC(NSPECW, MAXJPSIW, DEFAULT_DIST,
     &                              NSIW1D0)

            DO JP=IJS(JBL),IJL(JBL)
              IF(ADS(JP).LE.XMAX)THEN
                NSIW1D(NW1D(JP,JBL),NTOSIW1D(NW1D(JP,JBL)))=JSAR
                DIST(NW1D(JP,JBL),NTOSIW1D(NW1D(JP,JBL)))=ADS(JP)
              ENDIF
            ENDDO

          ENDIF
        ENDDO
        DEALLOCATE(ADS)
        DEALLOCATE(LOSAR)
      ENDDO

      DEALLOCATE(INDEX)

      CALL RESIZE_GETWSPEC(ICOUNT, MAXJPSIW, DEFAULT_DIST, NSIW1D0)
      WRITE(IU06,*)' '
      WRITE(IU06,*)'  IN GETWSPEC: '
      WRITE(IU06,*)'  NUMBER OF MODEL PTS INFLUENCED BY SAR PTS ',NSPECW
      WRITE(IU06,*)'  MAX NUMBER OF SAR PTS INFLUENCING A MODEL PTS ',
     &             MAXJPSIW
      WRITE(IU06,*)' '


!      3. COLLECTING THE WAM SPECTRA INFLUENCED BY SAR DATA POINTS.
!      ------------------------------------------------------------

      ALLOCATE(SPECW(0:NSPECW,NANG,NFRE))

!     READ FIRST BLOCK
      IF(IGL.GT.1) THEN
        REWIND IU11
        READ(IU11) FL1
      ENDIF
      DO JBL=1,IGL
!       ODD BLOCKS
        IF(MOD(JBL,2).EQ.1) THEN
          IF(JBL.NE.IGL) READ(IU11) FL2
            DO JFRE=1,NFRE
              DO JANG=1,NANG
                DO JP = IJS(JBL),IJL(JBL)
                  IF(NW1D(JP,JBL).GT.0) THEN
                    SPECW(NW1D(JP,JBL),JANG,JFRE)=FL1(JP,JANG,JFRE)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
!       EVEN BLOCKS
        ELSE
          IF(JBL.NE.IGL) READ(IU11) FL1
            DO JFRE=1,NFRE
              DO JANG=1,NANG
                DO JP = IJS(JBL),IJL(JBL)
                  IF(NW1D(JP,JBL).GT.0) THEN
                    SPECW(NW1D(JP,JBL),JANG,JFRE)=FL2(JP,JANG,JFRE)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
        ENDIF
      ENDDO

!     4. READ U10 WINDS AND DEFINE U10 AND THW FROM WAM SPECTRA
!     --------------------------------------------------------

      ALLOCATE(UW10(NSPECW))
      ALLOCATE(THWW(NSPECW))

      DO JBL=1,IGL
         DO JP = IJS(JBL),IJL(JBL)
           IF(NW1D(JP,JBL).GT.0) THEN
             UW10(NW1D(JP,JBL)) = U10OLD(JP,JBL)       
             THWW(NW1D(JP,JBL)) = THWOLD(JP,JBL)       
           ENDIF
         ENDDO
      ENDDO


      RETURN
      END SUBROUTINE GETWSPEC
