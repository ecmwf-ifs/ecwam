SUBROUTINE PROPAG_WAM (IJS, IJL, WAVNUM, CGROUP, OMOSNH2KD, &
&                      DEPTH, U, V,                         &
&                      FL1)

! ----------------------------------------------------------------------

!**** *PROPAG_WAM* - WAVE PROPGATION

!*    PURPOSE.
!     --------

!     PROPAGATION

!**   INTERFACE.
!     ----------

!     *CALL* *PROPAG_WAM (IJS, IJL, WAVNUM, CGROUP, OMOSNH2KD,
!                         DEPTH, U, V,
!                         FL1)*
!          *IJS*       - INDEX OF FIRST GRIDPOINT.
!          *IJL*       - INDEX OF LAST GRIDPOINT.
!          *WAVNUM*    - WAVE NUMBER.
!          *CGROUP*    - GROUP SPPED.
!          *OMOSNH2KD* - OMEGA / SINH(2KD)
!          *DEPTH*     - WATER DEPTH
!          *U*         - U-COMPONENT OF SURFACE CURRENT
!          *V*         - V-COMPONENT OF SURFACE CURRENT
!          *FL1*       - SPECTRUM


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCURR  , ONLY : LLCHKCFL ,LLCHKCFLA
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,NIBLO
      USE YOWREFD  , ONLY : LLUPDTTD ,THDD     ,THDC     ,SDOT
      USE YOWSHAL  , ONLY : BATHYMAX ,WAVNUM_LAND, CGROUP_LAND, OMOSNH2KD_LAND
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
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NFRE), INTENT(IN) :: WAVNUM, CGROUP, OMOSNH2KD
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: DEPTH, U, V
      REAL(KIND=JWRB), DIMENSION(IJS:IJL, NANG, NFRE), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M, J
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA, MTHREADS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

!     Spectra extended with the halo exchange for the propagation
!     But limited to NFRE_RED frequencies
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NANG, NFRE_RED) :: FL1_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NFRE_RED) :: WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1) :: DEPTH_EXT, U_EXT, V_EXT

      LOGICAL :: L1STCALL

      DATA L1STCALL / .TRUE. /

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',0,ZHOOK_HANDLE)


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
&                              DEPTH_EXT, U_EXT, V_EXT )


!            DOT THETA TERM:

             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO = IJS, IJL, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPDOT(KIJS, KIJL, NINF, NSUP,                                     &
     &                      WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT,                      &
     &                      DEPTH_EXT, U_EXT, V_EXT,                                    &
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
!!!
write(*,*) 'debile, before call to PROENVHALO ',IJS, IJL, NINF, NSUP
write(*,*) 'WAVNUM ',WAVNUM 
write(*,*) 'CGROUP ',CGROUP
write(*,*) 'OMOSNH2KD ',OMOSNH2KD
write(*,*) 'DEPTH ',DEPTH
write(*,*) 'U ',U
write(*,*) 'V ',V

               CALL  PROENVHALO (IJS, IJL, NINF, NSUP,                  &
&                                WAVNUM, CGROUP, OMOSNH2KD,             &
&                                DEPTH, U, V,                           &
&                                WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                                DEPTH_EXT, U_EXT, V_EXT )

!              COMPUTES ADVECTION WEIGTHS AND CHECK CFL CRITERIA
               CALL CTUWUPDT(IJS, IJL, NINF, NSUP,       &
&                            CGROUP_EXT, OMOSNH2KD_EXT,  &
&                            DEPTH_EXT, U_EXT, V_EXT )

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
             CALL  PROENVHALO (IJS, IJL, NINF, NSUP,      &
&                  WAVNUM, CGROUP, OMOSNH2KD,             &
&                  DEPTH, U, V,                           &
&                  WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                  DEPTH_EXT, U_EXT, V_EXT )

             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS,IJL,NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPAGS1(FL1_EXT, FL1, NINF, NSUP, IJS, IJL, KIJS, KIJL, &
&                            CGROUP_EXT, OMOSNH2KD_EXT,                      &
&                            U_EXT, V_EXT,                                   &
&                            L1STCALL)
             ENDDO
!$OMP       END PARALLEL DO

           CASE(0)
             IF (L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.

!            NEED HALO VALUES
             CALL  PROENVHALO (IJS, IJL, NINF, NSUP,      &
&                  WAVNUM, CGROUP, OMOSNH2KD,             &
&                  DEPTH, U, V,                           &
&                  WAVNUM_EXT, CGROUP_EXT, OMOSNH2KD_EXT, &
&                  DEPTH_EXT, U_EXT, V_EXT )

             NPROMA=NPROMA_WAM
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
             DO JKGLO=IJS,IJL,NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1,IJL)
               CALL PROPAGS(FL1_EXT, FL1, NINF, NSUP, IJS, IJL, KIJS, KIJL, &
&                           CGROUP_EXT, OMOSNH2KD_EXT,                      &
&                           U_EXT, V_EXT,                                   &
&                           L1STCALL)
             ENDDO
!$OMP        END PARALLEL DO
           END SELECT 


           CALL GSTATS(1430,1)

        ENDIF  ! end propagation

      ENDIF ! more than one grid point

      L1STCALL=.FALSE.
      LLCHKCFL=.FALSE.

      IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',1,ZHOOK_HANDLE)

END SUBROUTINE PROPAG_WAM
