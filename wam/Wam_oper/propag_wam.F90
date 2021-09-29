SUBROUTINE PROPAG_WAM (IJS, IJL, WVENVI, WVPRPT, FL1)

! ----------------------------------------------------------------------

!**** *PROPAG_WAM* - WAVE PROPGATION

!*    PURPOSE.
!     --------

!     PROPAGATION

!**   INTERFACE.
!     ----------

!     *CALL* *PROPAG_WAM (IJS, IJL, WVENVI, WVPRPT, FL1)
!          *IJS*       - INDEX OF FIRST GRIDPOINT.
!          *IJL*       - INDEX OF LAST GRIDPOINT.
!          *WVENVI*    - WAVE ENVIRONMENT FIELDS
!          *WVPRPT*    - WAVE PROPERTIES FIELDS
!          *FL1*       - SPECTRUM


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : ENVIRONMENT, FREQUENCY

      USE YOWCURR  , ONLY : LLCHKCFL ,LLCHKCFLA
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,NIBLO
      USE YOWREFD  , ONLY : LLUPDTTD ,THDD     ,THDC     ,SDOT
      USE YOWSHAL  , ONLY : BATHYMAX
      USE YOWSTAT  , ONLY : IPROPAGS ,NPROMA_WAM
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE UNWAM    , ONLY : PROPAG_UNWAM
      USE YOWUNPOOL ,ONLY : LLUNSTR

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "ctuwupdt.intfb.h"
#include "mpexchng.intfb.h"
#include "proenvhalo.intfb.h"
#include "propags.intfb.h"
#include "propags1.intfb.h"
#include "propags2.intfb.h"
#include "propdot.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      TYPE(ENVIRONMENT), DIMENSION(IJS:IJL), INTENT(IN) :: WVENVI
      TYPE(FREQUENCY), DIMENSION(IJS:IJL,NFRE), INTENT(IN) :: WVPRPT
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M, J
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, MTHREADS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

!     Spectra extended with the halo exchange for the propagation
!     But limited to NFRE_RED frequencies
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NANG, NFRE_RED) :: FL1_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED) :: WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1) :: DEPTH_EXT, UCUR_EXT, VCUR_EXT

      LOGICAL :: L1STCALL

      DATA L1STCALL / .TRUE. /

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',0,ZHOOK_HANDLE)

ASSOCIATE(DEPTH => WVENVI%DEPTH, &
 &        UCUR => WVENVI%UCUR, &
 &        VCUR => WVENVI%VCUR, &
 &        WAVNUM => WVPRPT%WAVNUM, &
 &        CGROUP => WVPRPT%CGROUP, &
 &        OMOSNH2KD => WVPRPT%OMOSNH2KD)


      IF (NIBLO > 1) THEN

        NPROMA=NPROMA_WAM
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,KIJS,KIJL,K,M,IJ) 
        DO JKGLO=IJS,IJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,IJL)
          DO M=1,NFRE_RED
            DO K=1,NANG
              DO IJ=KIJS,KIJL
                FL1_EXT(IJ,K,M) = FL1(IJ,K,M)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

!       SET THE DUMMY LAND POINT TO 0.
        FL1_EXT(NSUP+1,:,:) = 0.0_JWRB 


        IF (LLUNSTR) THEN 

           CALL PROPAG_UNWAM(FL1_EXT, FL1)

        ELSE


!          OBTAIN INFORMATION AT NEIGHBORING GRID POINTS (HALO)
!          ----------------------------------------------------
           CALL MPEXCHNG(FL1_EXT,NANG,NFRE_RED)


           CALL GSTATS(1430,0)

!          UPDATE DOT THETA TERM
!          ---------------------

           IF (LLUPDTTD) THEN
             IF (.NOT.ALLOCATED(THDC)) ALLOCATE(THDC(IJS:IJL,NANG))
             IF (.NOT.ALLOCATED(THDD)) ALLOCATE(THDD(IJS:IJL,NANG))
             IF (.NOT.ALLOCATED(SDOT)) ALLOCATE(SDOT(IJS:IJL,NANG,NFRE_RED))

!            NEED HALO VALUES
             CALL  PROENVHALO (IJS, IJL, NINF, NSUP,                  &
&                              WAVNUM, CGROUP, OMOSNH2KD,             &
&                              DEPTH, U, V,                           &
&                              WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                              DEPTH_EXT, UCUR_EXT, VCUR_EXT )


!            DOT THETA TERM:

             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO = IJS, IJL, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPDOT(KIJS, KIJL, NINF, NSUP,                                     &
     &                      WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT,                      &
     &                      DEPTH_EXT, UCUR_EXT, VCUR_EXT,                                    &
     &                      THDC(KIJS:KIJL,:), THDD(KIJS:KIJL,:), SDOT(KIJS:KIJL,:,:))
             ENDDO
!$OMP        END PARALLEL DO


             LLUPDTTD = .FALSE.
           ENDIF


!          IPROPAGS = 2 is the default option (the other options are kept but usually not used)
!          -----------------------------------------------------------------------------------
           SELECT CASE (IPROPAGS)

           CASE(2)

             IF (LUPDTWGHT) THEN
!              NEED HALO VALUES
               CALL  PROENVHALO (IJS, IJL, NINF, NSUP,                  &
&                                WAVNUM, CGROUP, OMOSNH2KD,             &
&                                DEPTH, UCUR, VCUR,                     &
&                                WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                                DEPTH_EXT, UCUR_EXT, VCUR_EXT )

!              COMPUTES ADVECTION WEIGTHS AND CHECK CFL CRITERIA
               CALL CTUWUPDT(IJS, IJL, NINF, NSUP,       &
&                            CGROUP_EXT, OMOSNH2KD_EXT,  &
&                            DEPTH_EXT, UCUR_EXT, VCUR_EXT )

               LUPDTWGHT=.FALSE.
             ENDIF

             MTHREADS=1
!$           MTHREADS=OMP_GET_MAX_THREADS()
             NPROMA=(IJL-IJS+1)/MTHREADS + 1
!$OMP        PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS,IJL,NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPAGS2(FL1_EXT, FL1, NINF, NSUP, IJS, IJL, KIJS, KIJL)
             ENDDO
!$OMP        END PARALLEL DO


           CASE(1)
             IF (L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.

!            NEED HALO VALUES
             CALL  PROENVHALO (IJS, IJL, NINF, NSUP,                  &
&                              WAVNUM, CGROUP, OMOSNH2KD,             &
&                              DEPTH, U, V,                           &
&                              WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                              DEPTH_EXT, UCUR_EXT, VCUR_EXT )

             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS,IJL,NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPAGS1(FL1_EXT, FL1, NINF, NSUP, IJS, IJL, KIJS, KIJL, &
&                            DEPTH,                                          &
&                            CGROUP_EXT, OMOSNH2KD_EXT,                      &
&                            UCUR_EXT, VCUR_EXT,                             &
&                            L1STCALL)
             ENDDO
!$OMP       END PARALLEL DO

           CASE(0)
             IF (L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.

!            NEED HALO VALUES
             CALL  PROENVHALO (IJS, IJL, NINF, NSUP,                  &
&                              WAVNUM, CGROUP, OMOSNH2KD,             &
&                              DEPTH, U, V,                           &
&                              WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                              DEPTH_EXT, UCUR_EXT, VCUR_EXT )

             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS,IJL,NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPAGS(FL1_EXT, FL1, NINF, NSUP, IJS, IJL, KIJS, KIJL, &
&                           DEPTH,                                          &
&                           CGROUP_EXT, OMOSNH2KD_EXT,                      &
&                           UCUR_EXT, VCUR_EXT,                             &
&                           L1STCALL)
             ENDDO
!$OMP        END PARALLEL DO
           END SELECT 


           CALL GSTATS(1430,1)

        ENDIF  ! end propagation

      ENDIF ! more than one grid point

      L1STCALL=.FALSE.
      LLCHKCFL=.FALSE.

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',1,ZHOOK_HANDLE)

END SUBROUTINE PROPAG_WAM
