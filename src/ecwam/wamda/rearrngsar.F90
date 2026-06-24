SUBROUTINE REARRNGSAR(TH0, BLK2GLO, FF_NOW, FL1)
 
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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, FORCING_FIELDS

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, ICHNKFROMIJ, IPRMFROMIJ 
      USE YOWICE   , ONLY : CITHRSH
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,XDELLA   ,ZDELLO
      USE YOWMPP   , ONLY : NPROC    ,IRANK    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI
      USE YOWSPEC,   ONLY : NSTART   ,NEND
      USE YOWSARAS , ONLY : NSPEC    ,IJSAR    ,SPEC     ,              &
     &                      LONG     ,LAT      ,U10      ,USSAR    ,THW
      USE YOWSTAT  , ONLY : LALTAS

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK


!----------------------------------------------------------------------

      IMPLICIT NONE

#include "findb.intfb.h"
#include "mpexchngsarin.intfb.h"

      REAL(KIND=JWRB) :: TH0 ! FIRST DIRECTION IN WAM SPECTRUM.
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1


!     VARIABLES
!     ---------

      REAL(KIND=JWRB), ALLOCATABLE :: FL(:,:)
      REAL(KIND=JWRB) :: FTH, ADIF, BDIF
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: ISPECPE 
      INTEGER(KIND=JWIM) :: ISPEC, K, M, KC, KC1, INC
      INTEGER(KIND=JWIM) :: IP, IK
      INTEGER(KIND=JWIM) :: JP, JFRE, JANG, NSPECPE

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('REARRNGSAR',0,ZHOOK_HANDLE)


!     ROTATE INVERTED SAR SPECTRA TO THE FIRST DIRECTION OF WAM (TH0)
!     --------------------------------------------------------------- 
!     !!!! NOTE IT IS ASSUMED THAT INPUT INVERTED SPECTRA ARE DEFINED
!          WITH 0 AS THEIR FIRST DIRECTION.

      FTH = MOD(-TH0+ZPI,ZPI)
      FTH = FTH * REAL(NANG) / ZPI
      INC = INT(FTH)
      ADIF = FTH - INC
      BDIF = 1.0_JWRB - ADIF
      ALLOCATE(FL(NANG,NFRE))

      DO ISPEC = 1,NSPEC
        DO K=1,NANG
          KC  = K  - INC
          IF (KC < 1) KC  = KC + NANG
          KC1 = KC - 1
          IF (KC1 < 1) KC1 = KC1+ NANG
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

      ALLOCATE(IJSAR(NSPEC))
      ALLOCATE(ISPECPE(NSPEC))
      ISPECPE=0

!     FIND BLOCK POSITION OF SAR OBSERVATION POINTS

      CALL FINDB(NSPEC,NSPEC,LAT(1,2),LONG(1,2),IJSAR,NPROC,NSTART,NEND,IRANK,IRANK)

!     REPLACE FIRST GUESS SPECTRA USED IN THE SAR INVERSION BY WAM
!     SPECTRA FOUND AT THE CLOSEST GRID POINT

      NSPECPE=0

!     READ FIRST BLOCK
          DO ISPEC = 1, NSPEC 
            JP=IJSAR(ISPEC)
            IF (JP /= 0) THEN

                IK = ICHNKFROMIJ(JP)
                IP = IPRMFROMIJ(JP)

                IF (FF_NOW%CICOVER(IP,IK) <= CITHRSH) THEN
                  DO JFRE=1,NFRE
                    DO JANG=1,NANG
                      SPEC(ISPEC,JANG,JFRE,1)=FL1(IP,JANG,JFRE,IK)
                    ENDDO
                  ENDDO
                  LAT(ISPEC,1)=AMOSOP+REAL(BLK2GLO%KXLT(JP)-1)*XDELLA
                  LONG(ISPEC,1)=AMOWEP+                                 &
     &                        REAL(BLK2GLO%IXLG(JP)-1)*ZDELLO(BLK2GLO%KXLT(JP))
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

! REPLACE WINDS FROM SAR INPUT BY MODEL WINDS

      ALLOCATE(U10(NSPEC))
      ALLOCATE(USSAR(NSPEC))
      ALLOCATE(THW(NSPEC))

        DO ISPEC = 1, NSPEC 
          JP=IJSAR(ISPEC)
          IF (JP /= 0) THEN
            IK = ICHNKFROMIJ(JP)
            IP = IPRMFROMIJ(JP) 
            U10(ISPEC) = FF_NOW%WSWAVE(IP, IK)
            USSAR(ISPEC) = FF_NOW%UFRIC(IP, IK)
            THW(ISPEC) = FF_NOW%WDWAVE(IP, IK)
          ELSE
            U10(ISPEC) = 0.0_JWRB
            USSAR(ISPEC) = 0.0_JWRB
            THW(ISPEC) = 0.0_JWRB
          ENDIF
        ENDDO

!     EXCHANGE THE DIFFERENT CONTRIBUTIONS OF THE MODEL VALUES AT SAR
!     LOCATION WITH ALL THE OTHER PE'S.

      CALL MPEXCHNGSARIN(NSPECPE,ISPECPE,KTAG)

      DEALLOCATE(ISPECPE)

IF (LHOOK) CALL DR_HOOK('REARRNGSAR',1,ZHOOK_HANDLE)

END SUBROUTINE REARRNGSAR 
