SUBROUTINE GETWSPEC(IJSG, IJLG, BLK2GLO, FF_NOW, FL1, NW1D)
 
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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, FORCING_FIELDS

      USE YOWICE   , ONLY : CITHRSH
      USE YOWGRID  , ONLY : SINPH, COSPH, NPROMA_WAM, NCHNK,            &
     &                      ICHNKFROMIJ, IPRMFROMIJ
      USE YOWMAP   , ONLY : KXLTMIN  ,KXLTMAX  ,IRGG     ,AMOWEP   ,    &
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,ZDELLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : RAD      ,DEG      ,R
      USE YOWSARAS , ONLY : PPSARCF  ,SARCORDIA,NSPEC    ,              &
     &            NSPECW   ,JPSIW    ,NTOSIW1D ,                        &
     &            IJSAR    ,LONG     ,LAT      ,COSSARLAT,SINSARLAT,    &
     &            NSIW1D   ,SPECW    ,DIST     ,                        &
     &            UW10     ,USW      ,THWW
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------
!     INTERFACE
!     ---------

      IMPLICIT NONE

#include "resize_getwspec.intfb.h"


      INTEGER(KIND=JWIM), INTENT(IN) :: IJSG, IJLG
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
      INTEGER(KIND=JWIM), INTENT(OUT) :: NW1D(IJSG:IJLG)

!           
!     VARIABLES
!     ---------
      INTEGER(KIND=JWIM) :: JP, JSAR, JFRE, JANG, ICOUNT,          &
     &                      MAXJPSIW, LLAT, NSIW1D0
      INTEGER(KIND=JWIM) :: ICOUNTOLD, IC, ICC
      INTEGER(KIND=JWIM) :: IP, IK
      INTEGER(KIND=JWIM), DIMENSION(NSPEC) :: IS, KS 
      INTEGER(KIND=JWIM), ALLOCATABLE :: INDEX(:)

      REAL(KIND=JWRB) :: DELLON, XI
      REAL(KIND=JWRB) :: XMAX, DEFAULT_DIST
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:) :: ADS
     
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LOSAR 
!--------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GETWSPEC',0,ZHOOK_HANDLE)


!     1.0 ALLOCATE ARRAYS FROM MODULE
!     --------------------------------

      IF (XDELLA < 0.2_JWRB) THEN
        NSPECW=400000/NPROC
      ELSEIF (XDELLA < 0.3_JWRB) THEN
        NSPECW=20000/NPROC
      ELSE
        NSPECW=INT((1.5_JWRB/XDELLA)**2*9000.0_JWRB/NPROC)
      ENDIF

      JPSIW=12
      NSIW1D0=0

      ALLOCATE(NTOSIW1D(NSPECW))
      ALLOCATE(NSIW1D(NSPECW,JPSIW))
      NSIW1D=NSIW1D0
      ALLOCATE(DIST(NSPECW,JPSIW))

             
!     2.0 LABEL WAM SPECTRA THAT ARE INFLUENCED BY THE SAR 
!         OBSERVATIONS AND ASSOCIATE THEM TO SAR WITHIN A 
!         XMAX DISTANCE (in radian).

!     MAXIMUN EXTEND OF XMAX IN NUMBER OF GRID POINTS 
      XMAX=PPSARCF*SARCORDIA/R

      LLAT = NINT(DEG*XMAX/XDELLA)+1

      IF (XMAX > 0.0_JWRB) THEN
        DEFAULT_DIST=10.0_JWRB*XMAX
      ELSE
        DEFAULT_DIST=10.0_JWRB
      ENDIF

      DIST=DEFAULT_DIST
      NTOSIW1D=0
      MAXJPSIW=0
      NW1D(:)=0
      ICOUNT=0

      DO JSAR=1,NSPEC
        IF (IJSAR(JSAR) /= 0) THEN
          IS(JSAR)=BLK2GLO%IXLG(IJSAR(JSAR))
          KS(JSAR)=BLK2GLO%KXLT(IJSAR(JSAR))
        ELSE
          IS(JSAR)=0
          KS(JSAR)=0
        ENDIF
      ENDDO


!     FIND POINTS WHICH FALL WITHIN XMAX FROM SAR OBSERVATIONS

      ALLOCATE(INDEX(1+IJLG-IJSG))

        ALLOCATE(ADS(IJSG:IJLG))
        ALLOCATE(LOSAR(IJSG:IJLG))
        DO JP = IJSG, IJLG
          LOSAR(JP)=.FALSE.
        ENDDO

        DO JSAR=1,NSPEC
!         exclude sar point that are obviously outside the sub region
          IF (KS(JSAR) >= KXLTMIN(IRANK)-LLAT .AND.                     &
     &        KS(JSAR) <= KXLTMAX(IRANK)+LLAT .AND.                     &
     &        IS(JSAR) > 0 .AND. KS(JSAR) > 0) THEN

            DO JP = IJSG, IJLG
              IK = ICHNKFROMIJ(JP)
              IP = IPRMFROMIJ(JP)
              IF (FF_NOW%CICOVER(IP, IK) <= CITHRSH .AND. ABS(KS(JSAR)-BLK2GLO%KXLT(JP)) <= LLAT) THEN

                XI = AMOWEP + REAL(BLK2GLO%IXLG(JP)-1)*ZDELLO(BLK2GLO%KXLT(JP))
                DELLON=RAD*(LONG(JSAR,2) - XI)
                ADS(JP)=COS(DELLON)*COSSARLAT(JSAR)*COSPH(BLK2GLO%KXLT(JP))+ SINSARLAT(JSAR)*SINPH(BLK2GLO%KXLT(JP))
                ADS(JP)=ACOS(MAX(MIN(ADS(JP),1.0_JWRB),-1.0_JWRB))
              ELSE
                ADS(JP)=DEFAULT_DIST
              ENDIF
            ENDDO

!           the point is already influences by an observation,
!           increment its number by 1.
            IC=0
            DO JP = IJSG, IJLG
              IF (ADS(JP) <= XMAX .AND. LOSAR(JP)) THEN
                IC=IC+1
                INDEX(IC)=NW1D(JP)
              ENDIF
            ENDDO
            DO ICC=1,IC
              NTOSIW1D(INDEX(ICC)) = NTOSIW1D(INDEX(ICC)) + 1
            ENDDO

            DO JP = IJSG, IJLG
              IF (ADS(JP) <= XMAX .AND. LOSAR(JP)) THEN
                MAXJPSIW=MAX(MAXJPSIW,NTOSIW1D(NW1D(JP)))
              ENDIF
            ENDDO

!           the point is not yet influenced by an observation
!           add it to the list
            IC=0
            DO JP = IJSG, IJLG
              IF (ADS(JP) <= XMAX .AND. .NOT.LOSAR(JP)) THEN
                IC=IC+1
                INDEX(IC)=JP
              ENDIF
            ENDDO

            ICOUNTOLD=ICOUNT
            ICOUNT=ICOUNTOLD+IC
            IF (ICOUNT > NSPECW) CALL RESIZE_GETWSPEC(ICOUNT, JPSIW, DEFAULT_DIST, NSIW1D0) 
            DO ICC=1,IC
              JP=INDEX(ICC)
              LOSAR(JP)=.TRUE.
              NW1D(JP)=ICOUNTOLD+ICC
              NTOSIW1D(ICOUNTOLD+ICC)=1
            ENDDO

            IF (MAXJPSIW > JPSIW) CALL RESIZE_GETWSPEC(NSPECW, MAXJPSIW, DEFAULT_DIST, NSIW1D0)

            DO JP = IJSG, IJLG
              IF (ADS(JP) <= XMAX) THEN
                NSIW1D(NW1D(JP),NTOSIW1D(NW1D(JP)))=JSAR
                DIST(NW1D(JP),NTOSIW1D(NW1D(JP)))=ADS(JP)
              ENDIF
            ENDDO

          ENDIF
        ENDDO

        DEALLOCATE(ADS)
        DEALLOCATE(LOSAR)

        DO JP = IJSG, IJLG
           IF (NW1D(JP) > 0) THEN
             MAXJPSIW=MAX(MAXJPSIW,NTOSIW1D(NW1D(JP)))
           ENDIF
        ENDDO

      DEALLOCATE(INDEX)

      IF (ICOUNT <= 0 .OR. MAXJPSIW <= 0) THEN
        NSPECW = 0
      ELSE 
        CALL RESIZE_GETWSPEC(ICOUNT, MAXJPSIW, DEFAULT_DIST, NSIW1D0)
      ENDIF
      WRITE(IU06,*)' '
      WRITE(IU06,*)'  IN GETWSPEC: '
      WRITE(IU06,*)'  NUMBER OF MODEL PTS INFLUENCED BY SAR PTS ',NSPECW
      WRITE(IU06,*)'  MAX NUMBER OF SAR PTS INFLUENCING A MODEL POINT ', &
     &             MAXJPSIW
      WRITE(IU06,*)' '


!      3. COLLECTING THE WAM SPECTRA INFLUENCED BY SAR DATA POINTS.
!      ------------------------------------------------------------

      IF (NSPECW > 0 ) THEN
        ALLOCATE(SPECW(0:NSPECW,NANG,NFRE))
        SPECW(:,:,:)=0.0_JWRB

        DO JFRE=1,NFRE
          DO JANG=1,NANG
            DO JP = IJSG, IJLG
              IF (NW1D(JP) > 0) THEN
                IK = ICHNKFROMIJ(JP)
                IP = IPRMFROMIJ(JP)
                SPECW(NW1D(JP),JANG,JFRE) = FL1(IP, JANG, JFRE, IK)
              ENDIF
            ENDDO
          ENDDO
        ENDDO

!       4. READ U10 WINDS AND DEFINE U10 AND THW FROM WAM SPECTRA
!       --------------------------------------------------------

        ALLOCATE(UW10(NSPECW))
        ALLOCATE(USW(NSPECW))
        ALLOCATE(THWW(NSPECW))

        DO JP = IJSG, IJLG
          IF (NW1D(JP) > 0) THEN
            IK = ICHNKFROMIJ(JP)
            IP = IPRMFROMIJ(JP)
            UW10(NW1D(JP)) = FF_NOW%WSWAVE(IP, IK)
            USW(NW1D(JP)) = FF_NOW%UFRIC(IP, IK)
            THWW(NW1D(JP)) = FF_NOW%WDWAVE(IP, IK)
          ENDIF
        ENDDO

      ENDIF

IF (LHOOK) CALL DR_HOOK('GETWSPEC',1,ZHOOK_HANDLE)

END SUBROUTINE GETWSPEC
