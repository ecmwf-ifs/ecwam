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
     &             LLCHKCFLA, CURRENT_MAX
      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE_RED ,NIBLO
      USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
      USE YOWSTAT  , ONLY : CDTPRO   ,IREFRA   ,NPROMA_WAM
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE YOWWIND  , ONLY : FIELDG   ,LLNEWCURR 

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

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

      INTEGER(KIND=JWIM) :: LIU
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, IJ, IX, IY

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1) :: OLDU, OLDV

      CHARACTER(LEN=14) :: CDATEIN, CDTNEWCUR
      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: LLCURRENT
      LOGICAL :: LLUPDATE

! --------------------------------------------------------------------- 

      IF (LHOOK) CALL DR_HOOK('GETCURR',0,ZHOOK_HANDLE)

!     1.0  GET NEW CURRENTS IF IT IS REQUIRED.
!          ----------------------------------

      CALL GSTATS(1984,0)


      IF (LLNEWCURR) THEN
        IF ( (LWCOU .AND. LWCUR ) .OR. LWNEMOCOUCUR .OR. IREFRA == 2 .OR. IREFRA == 3) THEN

          IF (.NOT.ALLOCATED(U)) THEN 
            ALLOCATE(U(NINF:NSUP+1))
            U(:)=0.0_JWRB
          ENDIF

          IF (.NOT.ALLOCATED(V)) THEN
            ALLOCATE(V(NINF:NSUP+1))
            V(:)=0.0_JWRB
          ENDIF

          CDTNEWCUR=CDTCUR

          IF (.NOT.LWCOU) CALL INCDATE(CDTNEWCUR,IDELCUR/2)

          IF (CDTPRO >= CDTNEWCUR) THEN

            OLDU(:)=U(:)
            OLDV(:)=V(:)

            LLCURRENT=.FALSE.

            CALL INCDATE(CDTCUR,IDELCUR)

            IF (LWCOU) THEN
!             CURRENTS FROM COUPLING INTERFACE 
!             --------------------------------
              IF (LWCUR .AND. .NOT.LWNEMOCOUCUR) THEN

! Mod for OPENMP
                  CALL GSTATS(1444,0)
                  NPROMA=NPROMA_WAM
!$OMP           PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
                  DO JKGLO = IJS, IJL, NPROMA
                    KIJS=JKGLO
                    KIJL=MIN(KIJS+NPROMA-1,IJL)
                    CALL WAMCUR (U(KIJS), V(KIJS), KIJS, KIJL)
                  ENDDO
!$OMP           END PARALLEL DO
                  CALL GSTATS(1444,1)
                  U(NSUP+1)=0.0_JWRB
                  V(NSUP+1)=0.0_JWRB

!!!               The halo has to be included for the calculation of the gradients !!!
                  CALL MPEXCHNG(U,1,1)
                  CALL MPEXCHNG(V,1,1)

!             CURRENTS FROM NEMO
!             ------------------
              ELSEIF (LWNEMOCOUCUR) THEN
                WRITE(IU06,*)' NEMO CURRENTS OBTAINED'!
                CALL GSTATS(1444,0)
                NPROMA=NPROMA_WAM
!$OMP           PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL,IX,IY)
                DO JKGLO = IJS, IJL, NPROMA
                  KIJS=JKGLO
                  KIJL=MIN(KIJS+NPROMA-1,IJL)
                  DO IJ=KIJS,KIJL
                    IX = IFROMIJ(IJ)
                    IY = JFROMIJ(IJ)
                    IF (FIELDG(IX,IY)%LKFR <=  0.0_JWRB ) THEN
!                     if lake cover = 0, we assume open ocean point, then get currents directly from NEMO 
                      U(IJ) = SIGN(MIN(ABS(NEMOUCUR(IJ)),CURRENT_MAX),NEMOUCUR(IJ))
                      V(IJ) = SIGN(MIN(ABS(NEMOVCUR(IJ)),CURRENT_MAX),NEMOVCUR(IJ))
                    ELSE
!                     no currents over lakes and land
                      U(IJ) = 0.0_JWRB
                      V(IJ) = 0.0_JWRB
                    ENDIF
                  ENDDO
                ENDDO
!$OMP           END PARALLEL DO
                CALL GSTATS(1444,1)

                U(NSUP+1)=0.0_JWRB
                V(NSUP+1)=0.0_JWRB

!!!             The halo has to be included for the calculation of the gradients !!!
                CALL MPEXCHNG(U,1,1)
                CALL MPEXCHNG(V,1,1)
              ELSE
                U(:)=0.0_JWRB
                V(:)=0.0_JWRB
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

                CALL CURRENT2WAM (FILNM,IREAD,CDATEIN)
                

                IF (CDATEIN /= CDTCUR) THEN
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
                U(:)=0.0_JWRB
                V(:)=0.0_JWRB
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
!           ------------------------
!           CHECK IF UPDATE IS NEEDED
            LLUPDATE=.FALSE.
            DO IJ = NINF, NSUP
              IF ( U(IJ) /= OLDU(IJ) .OR. V(IJ) /= OLDV(IJ) ) THEN
                LLUPDATE=.TRUE.
                EXIT
              ENDIF 
            ENDDO

            IF (LLUPDATE) THEN
              IF (IREFRA /= 0) THEN
                IF (.NOT.ALLOCATED(THDC)) ALLOCATE(THDC(IJS:IJL,NANG))
                IF (.NOT.ALLOCATED(THDD)) ALLOCATE(THDD(IJS:IJL,NANG))
                IF (.NOT.ALLOCATED(SDOT)) ALLOCATE(SDOT(IJS:IJL,NANG,NFRE_RED))

                CALL GSTATS(1444,0)
                NPROMA=NPROMA_WAM
!$OMP           PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
                DO JKGLO = IJS, IJL, NPROMA
                  KIJS=JKGLO
                  KIJL=MIN(KIJS+NPROMA-1,IJL)
                  CALL PROPDOT(KIJS, KIJL, THDC(KIJS:KIJL,:), THDD(KIJS:KIJL,:), SDOT(KIJS:KIJL,:,:))
                ENDDO
!$OMP           END PARALLEL DO
                CALL GSTATS(1444,1)

                WRITE(IU06,*) ' SUB. GETCURR: REFRACTION TERMS INITIALIZED'
                CALL FLUSH(IU06)
              END IF

              LLCHKCFLA=.TRUE.

!             SET LOGICAL TO RECOMPUTE THE WEIGHTS IN CTUW.
              LUPDTWGHT=.TRUE.

            ELSE
              LLCHKCFLA=.FALSE.
            ENDIF

          ELSE
            LLCHKCFLA=.FALSE.
          ENDIF
        ENDIF
      ELSE
        LLCHKCFLA=.FALSE.
      ENDIF

      CALL GSTATS(1984,1)

      IF (LHOOK) CALL DR_HOOK('GETCURR',1,ZHOOK_HANDLE)

END SUBROUTINE GETCURR 
