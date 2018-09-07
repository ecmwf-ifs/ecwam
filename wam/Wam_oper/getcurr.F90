      SUBROUTINE GETCURR(LWCUR, IREAD)

!****  *GETCURR* - READS SURFACE CURRENTS FROM FILE (IF UNCOUPLED)
!                  OR EXTRACT THEM FROM FORCING_FIELDS DATA STRUCTURE (IF COUPLED)

!     PURPOSE.
!     --------

!       GET SURFACE CURRENTS AND COMPUTE THE ASSOCIATED
!       REFRACTION TERMS.

!*    INTERFACE.
!     ----------

!     CALL *GETCURR*(LWCUR, IREAD)

!     *LWCUR*   -  LOGICAL CONTROLLING THE PRESENCE OF MEANINGFUL
!                  SURFACE CURRENTS IN ARRAY FORCING_FIELDS DATA STRUCTURE   
!     *IREAD*      INTEGER   PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                            (IF NEEDED).


!     METHOD.
!     -------

!     IF UNCOUPLED
!     READS CURRENTS FROM FILE "currents" WHEN IT IS NEEDED
!     IF COUPLED
!     GET CURRENTS AS PROCESSED BY *IFSTOWAM*
!     THEN 
!     CALL PROPDOT FOR THE CALCULATION OF THE REFRACTION TERMS.
!     ALSO SET LLCHKCFLA=.TRUE. TO CHECK ON THE CFL CRITERIA.

!     EXTERNALS.
!     ----------

!       *INCDATE*
!       *ABORT1*
!       *CURRENT2WAM*
!       *PROPDOT*

!     REFERENCES.
!     ----------- 

!       NONE

! ------------------------------------------------------------------- 

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOU, LWNEMOCOUCUR
      USE YOWNEMOFLDS,ONLY: NEMOUCUR, NEMOVCUR
      USE YOWCURR  , ONLY : U        ,V        ,CDTCUR   ,IDELCUR  ,    &
     &             LLCHKCFL,LLCHKCFLA, CURRENT_MAX
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NIBLO    ,NBLO
      USE YOWSTAT  , ONLY : CDTPRO   ,IREFRA   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE YOWWIND  , ONLY : FIELDG   ,LLNEWCURR 
      USE UNSTRUCT_CURR, ONLY : SET_CURTXY, SET_CURTXY_SINGLEFILE
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! --------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "current2wam.intfb.h"
#include "incdate.intfb.h"
#include "mpexchng.intfb.h"
#include "propdot.intfb.h"
#include "wamcur.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      LOGICAL, INTENT(IN) :: LWCUR

      INTEGER(KIND=JWIM) :: IG
      INTEGER(KIND=JWIM) :: LIU
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, IJ, IX, IY

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      CHARACTER(LEN=14) :: CDATEIN, CDTNEWCUR
      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: LLCURRENT

! --------------------------------------------------------------------- 

      IF (LHOOK) CALL DR_HOOK('GETCURR',0,ZHOOK_HANDLE)

!     1.0  GET NEW CURRENTS IF IT IS REQUIRED.
!          ----------------------------------

      CALL GSTATS(1984,0)

      IF (LLNEWCURR) THEN
        IF ( (LWCOU .AND. LWCUR ) .OR.                                  &
     &        IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN

          IF(.NOT.ALLOCATED(U)) ALLOCATE(U(NINF-1:NSUP,NBLO))
          IF(.NOT.ALLOCATED(V)) ALLOCATE(V(NINF-1:NSUP,NBLO))

          CDTNEWCUR=CDTCUR

          IF(.NOT.LWCOU) CALL INCDATE(CDTNEWCUR,IDELCUR/2)

          IF(CDTPRO.GE.CDTNEWCUR) THEN

            LLCURRENT=.FALSE.

            CALL INCDATE(CDTCUR,IDELCUR)

            IF(LWCOU) THEN
!             CURRENTS FROM COUPLING INTERFACE 
!             --------------------------------
              IF(LWCUR.AND..NOT.LWNEMOCOUCUR) THEN

                IG=1
!!!               from NINF to NSUP as the halo has to be included
!!!               for the calculation of the gradients !!!
! Mod for OPENMP
                  CALL GSTATS(1444,0)
                  NPROMA=NPROMA_WAM
!$OMP           PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
                  DO JKGLO=NINF,NSUP,NPROMA
                    KIJS=JKGLO
                    KIJL=MIN(KIJS+NPROMA-1,NSUP)
                    CALL WAMCUR (U(KIJS,1), V(KIJS,1), KIJS, KIJL)
                  ENDDO
!$OMP           END PARALLEL DO
                  CALL GSTATS(1444,1)
                  U(NINF-1,IG)=0.0_JWRB
                  V(NINF-1,IG)=0.0_JWRB

                IF (ITEST.GE.1) THEN
                   WRITE (IU06,*) ' '
                   WRITE (IU06,*) '  SUB. GETCURR:',                    &
     &             ' INPUT CURRENTS FIELDS CONVERTED TO BLOCKS'
                   CALL FLUSH(IU06)
                ENDIF
!             CURRENTS FROM NEMO
!             ------------------
              ELSEIF(LWNEMOCOUCUR) THEN
                WRITE(IU06,*)'NEMO CURRENTS'!
                IG=1
                DO IJ = IJS(IG),IJL(IG)
                  IX = IFROMIJ(IJ,IG)
                  IY = JFROMIJ(IJ,IG)
                  IF (FIELDG(IX,IY)%LKFR .LE. 0.0_JWRB ) THEN
!                   if lake cover = 0, we assume open ocean point, then get currents directly from NEMO 
                    U(IJ,IG) = SIGN(MIN(ABS(NEMOUCUR(IJ)),CURRENT_MAX),NEMOUCUR(IJ))
                    V(IJ,IG) = SIGN(MIN(ABS(NEMOVCUR(IJ)),CURRENT_MAX),NEMOVCUR(IJ))
                  ELSE
!                   no currents over lakes and land
                    U(IJ,IG) = 0.0_JWRB
                    V(IJ,IG) = 0.0_JWRB
                  ENDIF
                ENDDO
                U(NINF-1,IG)=0.0_JWRB
                V(NINF-1,IG)=0.0_JWRB
                CALL MPEXCHNG(U,1,1)
                CALL MPEXCHNG(V,1,1)
              ELSE
                U=0.0_JWRB
                V=0.0_JWRB
                WRITE(IU06,*)' '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)'     IREFRA = ',IREFRA
                WRITE(IU06,*)'     LWCUR IS FALSE !!!!' 
                WRITE(IU06,*)'     CURRENTS ARE SET TO 0. '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)' '
                CALL FLUSH(IU06)
              ENDIF
            ELSE
!             CURRENTS FROM INPUT FILE
!             ------------------------
              FILNM = 'currents'
              LIU   = LEN_TRIM(FILNM)
              FILNM=FILNM(1:LIU)
              INQUIRE(FILE=FILNM,EXIST=LLCURRENT)

              IF (LLCURRENT) THEN
                WRITE(IU06,*) ' '
                WRITE(IU06,*) ' SUB. GETCURR: GETTING OCEAN CURRENTS',  &
     &                        ' FOR DATE ',CDTCUR
                CALL FLUSH(IU06)

                IF (LLUNSTR) THEN
                  CALL SET_CURTXY_SINGLEFILE
                ELSE
                  CALL CURRENT2WAM (FILNM,IREAD,CDATEIN)
                END IF
                

                IF(CDATEIN.NE.CDTCUR) THEN
                WRITE (IU06,*) ' **************************************'
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' * PROBLEM IN GETCURR :               *'
                WRITE (IU06,*) ' * THE REQUESTED DATE FOR THE CURRENTS*'
                WRITE (IU06,*) ' * DOES NOT CORRESPOND TO THE DECODED *'
                WRITE (IU06,*) ' * DATE !!!!                          *'
                WRITE (IU06,*) ' * CDTCUR =',CDTCUR 
                WRITE (IU06,*) ' * CDATEIN=',CDATEIN
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' **************************************'
                CALL ABORT1
                ENDIF
              ELSE
                U=0.0_JWRB
                V=0.0_JWRB
                WRITE(IU06,*)' '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)'     FILE ',FILNM,' NOT FOUND '
                WRITE(IU06,*)'     CURRENTS ARE SET TO 0. '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)' '
                CALL FLUSH(IU06)
              ENDIF

            ENDIF


!           COMPUTE REFRACTION TERMS
            IF (IREFRA .NE. 0) THEN
              CALL PROPDOT
              IF (ITEST.GE.2) THEN
                WRITE(IU06,*) ' SUB. GETCURR: REFRACTION ',             &
     &         ' TERMS INITIALIZED'
                CALL FLUSH(IU06)
              END IF
            END IF

            LLCHKCFLA=.TRUE.

!           SET LOGICAL TO RECOMPUTE THE WEIGHTS IN CTUW.
            LUPDTWGHT=.TRUE.

          ELSE
            LLCHKCFLA=.FALSE.
          ENDIF
          IF (LLUNSTR) THEN
            CALL SET_CURTXY
          END IF
        ENDIF
      ELSE
        LLCHKCFLA=.FALSE.
      ENDIF

      CALL GSTATS(1984,1)

      IF (LHOOK) CALL DR_HOOK('GETCURR',1,ZHOOK_HANDLE)

      END SUBROUTINE GETCURR 
