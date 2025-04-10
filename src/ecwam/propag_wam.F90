! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PROPAG_WAM (BLK2GLO, WAVNUM, CGROUP, OMOSNH2KD, FL1, &
&  DEPTH, DELLAM1, COSPHM1, UCUR, VCUR)

! ----------------------------------------------------------------------

!**** *PROPAG_WAM* - WAVE PROPGATION

!*    PURPOSE.
!     --------

!     PROPAGATION

!**   INTERFACE.
!     ----------

!     *CALL* *PROPAG_WAM (BLK2GLO, WVENVI, WVPRPT, FL1)
!          *BLK2GLO*   - BLOCK TO GRID TRANSFORMATION
!          *WVENVI*    - WAVE ENVIRONMENT FIELDS
!          *WVPRPT*    - WAVE PROPERTIES FIELDS
!          *FL1*       - SPECTRUM


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY

      USE YOWCURR  , ONLY : LLCHKCFL ,LLCHKCFLA
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK
      USE YOWMAP   , ONLY : NIBLO
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,LLUNSTR
      USE YOWREFD  , ONLY : LLUPDTTD ,THDD     ,THDC     ,SDOT
      USE YOWSTAT  , ONLY : IPROPAGS ,IFRELFMAX, DELPRO_LF, IDELPRO
      USE YOWUBUF  , ONLY : LUPDTWGHT
#ifdef WAM_HAVE_UNWAM
      USE UNWAM    , ONLY : PROPAG_UNWAM
#endif

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR
      USE YOWABORT , ONLY : WAM_ABORT
      USE OML_MOD  , ONLY : OML_GET_MAX_THREADS

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "ctuwupdt.intfb.h"
#include "mpexchng.intfb.h"
#include "proenvhalo.intfb.h"
#include "propags.intfb.h"
#include "propags1.intfb.h"
#include "propags2.intfb.h"
#include "propdot.intfb.h"

      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NFRE, NCHNK), INTENT(IN) :: WAVNUM, CGROUP, OMOSNH2KD
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: DEPTH, DELLAM1, COSPHM1, UCUR, VCUR


      INTEGER(KIND=JWIM) :: IJ, K, M, J, II
      INTEGER(KIND=JWIM) :: JKGLO, NPROMA, MTHREADS
      INTEGER(KIND=JWIM) :: NSTEP_LF, ISUBST
      INTEGER(KIND=JWIM) :: IJSG, IJLG, ICHNK, KIJS, KIJL, IJSB, IJLB
      INTEGER(KIND=JWIM) :: ND3SF1, ND3EF1, ND3S, ND3E

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_MPI

!     Spectra extended with the halo exchange for the propagation
!     But limited to NFRE_RED frequencies
!!! the advection schemes are still written in block structure
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NANG, NFRE_RED) :: FL1_EXT, FL3_EXT
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, 3*NFRE_RED + 5) :: BUFFER_EXT

      LOGICAL :: L1STCALL

      DATA L1STCALL / .TRUE. /

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',0,ZHOOK_HANDLE)


!$acc data present(FL1, WAVNUM, CGROUP, OMOSNH2KD, DEPTH, DELLAM1,COSPHM1,UCUR,VCUR,BLK2GLO) CREATE(FL1_EXT,FL3_EXT) &
!$acc & create(BUFFER_EXT)
      IF (NIBLO > 1) THEN

        IJSG = IJFROMCHNK(1,1)
        IJLG = IJSG + SUM(KIJL4CHNK) - 1

        MTHREADS=1
#ifndef _OPENACC
        MTHREADS=OML_GET_MAX_THREADS()
#endif
        NPROMA=(IJLG-IJSG+1)/MTHREADS + 1


!!! the advection schemes are still written in block structure
!!! mapping chuncks to block ONLY for actual grid points !!!!
#ifdef _OPENACC
        !$acc parallel loop private(KIJS, IJSB, KIJL, IJLB)
#else
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB, M, K)
#endif /*_OPENACC*/
        DO ICHNK = 1, NCHNK
          KIJS = 1
          IJSB = IJFROMCHNK(KIJS, ICHNK)
          KIJL = KIJL4CHNK(ICHNK)
          IJLB = IJFROMCHNK(KIJL, ICHNK)
          !$acc loop independent collapse(2)
          DO M = 1, NFRE_RED
            DO K = 1, NANG
              FL1_EXT(IJSB:IJLB, K, M) = FL1(KIJS:KIJL, K, M, ICHNK)
            ENDDO
          ENDDO
        ENDDO
#ifdef _OPENACC
        !$acc end parallel loop
#else
!$OMP   END PARALLEL DO
#endif /*_OPENACC*/

!       SET THE DUMMY LAND POINT TO 0.
        !$acc kernels
        FL1_EXT(NSUP+1,:,:) = 0.0_JWRB 
        !$acc end kernels


        IF (LLUNSTR) THEN 

        WRITE(NULERR,*) '!!! ********************************* !!'
        WRITE(NULERR,*) '!!! in PROPAG_WAM Not yet ready !!!' 
        WRITE(NULERR,*) '!!! PROPAG_UNWAM will need to be adapted to new FL1 structure !!!' 
        WRITE(NULERR,*) '!!! ********************************* !!'
        CALL ABORT1

!!!!!           CALL PROPAG_UNWAM(FL1_EXT, FL1)

        ELSE


!          OBTAIN INFORMATION AT NEIGHBORING GRID POINTS (HALO)
!          ----------------------------------------------------
           IF (LHOOK) CALL DR_HOOK('MPI_TIME',0,ZHOOK_HANDLE_MPI)
           CALL MPEXCHNG(FL1_EXT, NANG, 1, NFRE_RED)
           IF (LHOOK) CALL DR_HOOK('MPI_TIME',1,ZHOOK_HANDLE_MPI)


           CALL GSTATS(1430,0)

!          UPDATE DOT THETA TERM
!          ---------------------

           IF (LLUPDTTD) THEN
             IF (.NOT.ALLOCATED(THDC)) ALLOCATE(THDC(IJSG:IJLG, NANG))
             IF (.NOT.ALLOCATED(THDD)) ALLOCATE(THDD(IJSG:IJLG, NANG))
             IF (.NOT.ALLOCATED(SDOT)) ALLOCATE(SDOT(IJSG:IJLG, NANG, NFRE_RED))

!            NEED HALO VALUES
             CALL  PROENVHALO (NINF, NSUP,                            &
&                              WAVNUM, CGROUP, OMOSNH2KD,            &
&                              DEPTH, DELLAM1, COSPHM1, UCUR, VCUR,   &
&                              BUFFER_EXT)


!            DOT THETA TERM:

#ifdef _OPENACC
             !$acc data create(THDC,THDD,SDOT) present(BUFFER_EXT,BLK2GLO)
#else
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
#endif
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPDOT(KIJS, KIJL, NINF, NSUP,                                &
     &                      BLK2GLO,                                               &
     &                      BUFFER_EXT(:,1:NFRE_RED), BUFFER_EXT(:,NFRE_RED+1:2*NFRE_RED), &
     &                      BUFFER_EXT(:,2*NFRE_RED+1:3*NFRE_RED),      &
     &                      BUFFER_EXT(:,3*NFRE_RED+2), BUFFER_EXT(:,3*NFRE_RED+3),    &
     &                      BUFFER_EXT(:,3*NFRE_RED+4), BUFFER_EXT(:,3*NFRE_RED+5),    &
     &                      THDC(KIJS:KIJL,:), THDD(KIJS:KIJL,:), SDOT(KIJS:KIJL,:,:))
             ENDDO
#ifdef _OPENACC
             !$acc end data
#else
!$OMP        END PARALLEL DO
#endif

             LLUPDTTD = .FALSE.
           ENDIF


!          IPROPAGS = 2 is the default option (the other options are kept but usually not used)
!          -----------------------------------------------------------------------------------
           SELECT CASE (IPROPAGS)

           CASE(2)

             IF (LUPDTWGHT) THEN
!              NEED HALO VALUES
               CALL  PROENVHALO (NINF, NSUP,                            &
&                                WAVNUM, CGROUP, OMOSNH2KD,            &
&                                DEPTH, DELLAM1, COSPHM1, UCUR, VCUR,   &
&                                BUFFER_EXT )

!              COMPUTES ADVECTION WEIGTHS AND CHECK CFL CRITERIA
               CALL CTUWUPDT(IJSG, IJLG, NINF, NSUP,                      &
&                            BLK2GLO,                                     &
&                            BUFFER_EXT(:,NFRE_RED+1:2*NFRE_RED), BUFFER_EXT(:,2*NFRE_RED+1:3*NFRE_RED),        &
&                            BUFFER_EXT(:,3*NFRE_RED+2), BUFFER_EXT(:,3*NFRE_RED+3), &
&                            BUFFER_EXT(:,3*NFRE_RED+4), BUFFER_EXT(:,3*NFRE_RED+5) )

               LUPDTWGHT=.FALSE.
             ENDIF


             ND3SF1=1
             ND3EF1=NFRE_RED
             ND3S=1
             ND3E=NFRE_RED

#ifndef _OPENACC
!$OMP        PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO, KIJS, KIJL)
#endif /*_OPENACC*/
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS2(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL, NANG, ND3SF1, ND3EF1, ND3S, ND3E)
             ENDDO
#ifndef _OPENACC             
!$OMP        END PARALLEL DO
#endif /*_OPENACC*/

!            SUB TIME STEPPING FOR FAST WAVES (only if IFRELFMAX > 0)
             IF (IFRELFMAX > 0 .AND. IFRELFMAX < NFRE_RED) THEN
               NSTEP_LF = NINT(REAL(IDELPRO, JWRB)/DELPRO_LF)
               ISUBST = 2  ! The first step was done as part of the previous call to PROPAGS2

               ND3SF1=1
               ND3EF1=IFRELFMAX+1
               ND3S=1
               ND3E=IFRELFMAX

               DO WHILE (ISUBST <= NSTEP_LF)

#ifdef _OPENACC
!$acc kernels loop private(KIJS, KIJL)
#else
!$OMP            PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO, KIJS, KIJL, M, K, IJ)
#endif /*_OPENACC*/
                 DO JKGLO = IJSG, IJLG, NPROMA
                   KIJS=JKGLO
                   KIJL=MIN(KIJS+NPROMA-1, IJLG)
                   !$acc loop independent collapse(3)
                   DO M = ND3S, ND3E
                     DO K = 1, NANG
                       DO IJ = KIJS, KIJL
                         FL1_EXT(IJ, K, M) = FL3_EXT(IJ, K, M)
                       ENDDO
                     ENDDO
                   ENDDO
                 ENDDO
#ifdef _OPENACC
!$acc end kernels
#else
!$OMP            END PARALLEL DO
#endif /*_OPENACC*/

!                OBTAIN INFORMATION AT NEIGHBORING GRID POINTS (HALO)
                 IF (LHOOK) CALL DR_HOOK('MPI_TIME',0,ZHOOK_HANDLE_MPI)
                 CALL MPEXCHNG(FL1_EXT(:,:,ND3S:ND3E), NANG, ND3S, ND3E)
                 IF (LHOOK) CALL DR_HOOK('MPI_TIME',1,ZHOOK_HANDLE_MPI)

#ifndef _OPENACC
!$OMP            PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO, KIJS, KIJL)
#endif /*_OPENACC*/
                 DO JKGLO = IJSG, IJLG, NPROMA
                   KIJS=JKGLO
                   KIJL=MIN(KIJS+NPROMA-1, IJLG)
                   CALL PROPAGS2(FL1_EXT(:,:,ND3SF1:ND3EF1), FL3_EXT(:,:,ND3S:ND3E), &
                  &              NINF, NSUP, KIJS, KIJL, NANG, ND3SF1, ND3EF1, ND3S, ND3E)
                 ENDDO
#ifndef _OPENACC
!$OMP            END PARALLEL DO
#endif /*_OPENACC*/

                 ISUBST = ISUBST + 1

               ENDDO
             
ENDIF  ! end sub time steps (if needed)

           CASE(1)
#ifdef _OPENACC
           CALL WAM_ABORT("PROPAG_WAM: BRANCH NOT YET PORTED FOR GPU EXECUTION")
#endif
             IF (L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.

!            NEED HALO VALUES
             CALL  PROENVHALO (NINF, NSUP,                            &
&                              WAVNUM, CGROUP, OMOSNH2KD,            &
&                              DEPTH, DELLAM1, COSPHM1, UCUR, VCUR,   &
&                              BUFFER_EXT)

!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS1(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL,       &
&                            BLK2GLO,                                        &
&                            BUFFER_EXT(:,3*NFRE_RED+3),                       &
&                            BUFFER_EXT(:,NFRE_RED+1:2*NFRE_RED), BUFFER_EXT(:,2*NFRE_RED+1:3*NFRE_RED),           &
&                            BUFFER_EXT(:,3*NFRE_RED+1), BUFFER_EXT(:,3*NFRE_RED+2), &
&                            BUFFER_EXT(:,3*NFRE_RED+4), BUFFER_EXT(:,3*NFRE_RED+5), &
&                            L1STCALL)
             ENDDO
!$OMP       END PARALLEL DO

           CASE(0)
             IF (L1STCALL .OR. LLCHKCFLA) LLCHKCFL=.TRUE.

!            NEED HALO VALUES
             CALL  PROENVHALO (NINF, NSUP,                            &
&                              WAVNUM, CGROUP, OMOSNH2KD,            &
&                              DEPTH, DELLAM1, COSPHM1, UCUR, VCUR,   &
&                              BUFFER_EXT)

!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL,       &
&                           BLK2GLO,                                        &
&                           BUFFER_EXT(:,3*NFRE_RED+3),                       &
&                           BUFFER_EXT(:,NFRE_RED+1:2*NFRE_RED), BUFFER_EXT(:,2*NFRE_RED+1:3*NFRE_RED), &
&                           BUFFER_EXT(:,3*NFRE_RED+1), BUFFER_EXT(:,3*NFRE_RED+2), &
&                           BUFFER_EXT(:,3*NFRE_RED+4), BUFFER_EXT(:,3*NFRE_RED+5), &
&                           L1STCALL)
             ENDDO
!$OMP        END PARALLEL DO
           END SELECT 


!!! the advection schemes are still written in block structure
!!!  So need to convert back to the nproma_wam chuncks
#ifdef _OPENACC
        !$acc kernels loop independent private(KIJS, IJSB, KIJL, IJLB)
#else
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB, M, K, II, J)
#endif /*_OPENACC*/
          DO ICHNK = 1, NCHNK
            KIJS = 1
            IJSB = IJFROMCHNK(KIJS, ICHNK)
            KIJL = KIJL4CHNK(ICHNK)
            IJLB = IJFROMCHNK(KIJL, ICHNK)
            !$acc loop independent collapse(3)
            DO M = 1, NFRE_RED
              DO K = 1, NANG
                 DO J = KIJS, KIJL
                   II = IJSB + J - KIJS
                   FL1(J, K, M, ICHNK) = FL3_EXT(II, K, M)
                 ENDDO
              ENDDO
            ENDDO

            IF (KIJL < NPROMA_WAM) THEN
              !!! make sure fictious points keep values of the first point in the chunk
              !$acc loop independent collapse(3)
              DO M = 1, NFRE_RED
                DO K = 1, NANG
                  DO J = KIJL+1,NPROMA_WAM
                    FL1(J, K, M, ICHNK) = FL1(1, K, M, ICHNK)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

          ENDDO
#ifdef _OPENACC
        !$acc end kernels
#else
!$OMP     END PARALLEL DO
#endif /*_OPENACC*/

           CALL GSTATS(1430,1)

        ENDIF  ! end propagation

      ENDIF ! more than one grid point
!$acc end data

      L1STCALL=.FALSE.
      LLCHKCFL=.FALSE.

IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',1,ZHOOK_HANDLE)

END SUBROUTINE PROPAG_WAM
