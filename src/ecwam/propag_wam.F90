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


      INTEGER(KIND=JWIM) :: IJ, K, M, J, IO
      INTEGER(KIND=JWIM) :: JKGLO, NPROMA, MTHREADS
      INTEGER(KIND=JWIM) :: NSTEP_LF, ISUBST
      INTEGER(KIND=JWIM) :: IJSG, IJLG, ICHNK, KIJS, KIJL, IJSB
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


#ifdef OMPGPU
!$omp target data map(to:FL1,WAVNUM,CGROUP,OMOSNH2KD,DEPTH,DELLAM1,COSPHM1,UCUR,VCUR,BLK2GLO) &
!$omp & map(alloc:FL1_EXT,FL3_EXT,BUFFER_EXT)
#else
!$acc data present(FL1, WAVNUM, CGROUP, OMOSNH2KD, DEPTH, DELLAM1,COSPHM1,UCUR,VCUR,BLK2GLO) &
!$acc & create(FL1_EXT,FL3_EXT,BUFFER_EXT)
#endif
      IF (NIBLO > 1) THEN

        IJSG = IJFROMCHNK(1,1)
        IJLG = IJSG + SUM(KIJL4CHNK) - 1

        MTHREADS=1
#ifndef WAM_GPU
        MTHREADS=OML_GET_MAX_THREADS()
#endif
        NPROMA=(IJLG-IJSG+1)/MTHREADS + 1


!!! the advection schemes are still written in block structure
!!! mapping chuncks to block ONLY for actual grid points !!!!
#ifdef WAM_GPU
#ifdef OMPGPU
        !$omp target teams distribute parallel do collapse(4)
#else
        !$acc parallel loop gang vector collapse(4)
#endif
#else
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, M, K, IJ)
#endif
        DO ICHNK = 1, NCHNK
          DO M = 1, NFRE_RED
            DO K = 1, NANG
              DO IJ = 1, NPROMA_WAM
                IF (IJ <= KIJL4CHNK(ICHNK)) THEN
                  FL1_EXT(IJ-1+IJFROMCHNK(1, ICHNK), K, M) = FL1(IJ, K, M, ICHNK)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
#ifdef WAM_GPU
#ifdef OMPGPU
        !$omp end target teams distribute parallel do
#else
        !$acc end parallel loop
#endif
#else
!$OMP   END PARALLEL DO
#endif

!       SET THE DUMMY LAND POINT TO 0.
#ifdef OMPGPU
        !$omp target teams distribute parallel do collapse(2)
#else
        !$acc kernels
#endif
        DO M = 1, NFRE_RED
          DO K = 1, NANG
            FL1_EXT(NSUP+1,K,M) = 0.0_JWRB 
          ENDDO
        ENDDO
#ifdef OMPGPU
        !$omp end target teams distribute parallel do
#else
        !$acc end kernels
#endif
        
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
             IF (.NOT.ALLOCATED(THDC)) THEN
               ALLOCATE(THDC(IJSG:IJLG, NANG))
#ifdef OMPGPU
               !$omp target enter data map(alloc: THDC)
#else
               !$acc enter data create(THDC)
#endif
             ENDIF
             IF (.NOT.ALLOCATED(THDD)) THEN
               ALLOCATE(THDD(IJSG:IJLG, NANG))
#ifdef OMPGPU
               !$omp target enter data map(alloc: THDD)
#else
               !$acc enter data create(THDD)
#endif
             ENDIF
             IF (.NOT.ALLOCATED(SDOT)) THEN
               ALLOCATE(SDOT(IJSG:IJLG, NANG, NFRE_RED))
#ifdef OMPGPU
               !$omp target enter data map(alloc: SDOT)
#else
               !$acc enter data create(SDOT)
#endif
             ENDIF

!            NEED HALO VALUES
             CALL  PROENVHALO (NINF, NSUP,                            &
&                              WAVNUM, CGROUP, OMOSNH2KD,            &
&                              DEPTH, DELLAM1, COSPHM1, UCUR, VCUR,   &
&                              BUFFER_EXT)


!            DOT THETA TERM:

#ifdef WAM_GPU
#ifdef OMPGPU
             !$omp target data map(to:BUFFER_EXT,BLK2GLO)
#else
             !$acc data present(BUFFER_EXT,BLK2GLO)
#endif
#else
!$OMP        PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
#endif
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPDOT(KIJS, KIJL, NINF, NSUP,                                &
     &                      BLK2GLO, IJSG, IJLG,                                           &
     &                      BUFFER_EXT(:,1:NFRE_RED), BUFFER_EXT(:,NFRE_RED+1:2*NFRE_RED), &
     &                      BUFFER_EXT(:,2*NFRE_RED+1:3*NFRE_RED),      &
     &                      BUFFER_EXT(:,3*NFRE_RED+2), BUFFER_EXT(:,3*NFRE_RED+3),    &
     &                      BUFFER_EXT(:,3*NFRE_RED+4), BUFFER_EXT(:,3*NFRE_RED+5),    &
     &                      THDC, THDD, SDOT)
             ENDDO
#ifdef WAM_GPU
#ifdef OMPGPU
             !$omp end target data
#else
             !$acc end data
#endif
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

#ifndef WAM_GPU
!$OMP        PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO, KIJS, KIJL)
#endif
             DO JKGLO = IJSG, IJLG, NPROMA
               KIJS=JKGLO
               KIJL=MIN(KIJS+NPROMA-1, IJLG)
               CALL PROPAGS2(FL1_EXT, FL3_EXT, NINF, NSUP, KIJS, KIJL, NANG, ND3SF1, ND3EF1, ND3S, ND3E)
             ENDDO
#ifndef WAM_GPU
!$OMP        END PARALLEL DO
#endif

!            SUB TIME STEPPING FOR FAST WAVES (only if IFRELFMAX > 0)
             IF (IFRELFMAX > 0 .AND. IFRELFMAX < NFRE_RED) THEN
               NSTEP_LF = NINT(REAL(IDELPRO, JWRB)/DELPRO_LF)
               ISUBST = 2  ! The first step was done as part of the previous call to PROPAGS2

               ND3SF1=1
               ND3EF1=IFRELFMAX+1
               ND3S=1
               ND3E=IFRELFMAX

               DO WHILE (ISUBST <= NSTEP_LF)

#ifdef WAM_GPU
#ifdef OMPGPU
!$omp target teams distribute parallel do collapse(4) private(IJ)
#else
!$acc parallel loop gang vector collapse(4) private(IJ)
#endif
#else
!$OMP            PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO, M, K, IO, IJ)
#endif
                 DO JKGLO = IJSG, IJLG, NPROMA
                   DO M = ND3S, ND3E
                     DO K = 1, NANG
                       DO IO = 0, NPROMA-1
                         IJ = JKGLO + IO
                         IF (IJ <= IJLG) FL1_EXT(IJ, K, M) = FL3_EXT(IJ, K, M)
                       ENDDO
                     ENDDO
                   ENDDO
                 ENDDO
#ifdef WAM_GPU
#ifdef OMPGPU
!$omp end target teams distribute parallel do
#else
!$acc end parallel loop
#endif
#else
!$OMP            END PARALLEL DO
#endif

!                OBTAIN INFORMATION AT NEIGHBORING GRID POINTS (HALO)
                 IF (LHOOK) CALL DR_HOOK('MPI_TIME',0,ZHOOK_HANDLE_MPI)
                 CALL MPEXCHNG(FL1_EXT(:,:,ND3S:ND3E), NANG, ND3S, ND3E)
                 IF (LHOOK) CALL DR_HOOK('MPI_TIME',1,ZHOOK_HANDLE_MPI)

#ifndef WAM_GPU
!$OMP            PARALLEL DO SCHEDULE(STATIC,1) PRIVATE(JKGLO, KIJS, KIJL)
#endif
                 DO JKGLO = IJSG, IJLG, NPROMA
                   KIJS=JKGLO
                   KIJL=MIN(KIJS+NPROMA-1, IJLG)
                   CALL PROPAGS2(FL1_EXT(:,:,ND3SF1:ND3EF1), FL3_EXT(:,:,ND3S:ND3E), &
                  &              NINF, NSUP, KIJS, KIJL, NANG, ND3SF1, ND3EF1, ND3S, ND3E)
                 ENDDO
#ifndef WAM_GPU
!$OMP            END PARALLEL DO
#endif

                 ISUBST = ISUBST + 1

               ENDDO
             
ENDIF  ! end sub time steps (if needed)

           CASE(1)
#ifdef WAM_GPU
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
!!! fictitious points (J > KIJL) keep the values of the first point in the chunk
#ifdef WAM_GPU
#ifdef OMPGPU
        !$omp target teams distribute parallel do collapse(4) private(IJSB)
#else
        !$acc parallel loop gang vector collapse(4) private(IJSB)
#endif
#else
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, M, K, J, IJSB)
#endif
          DO ICHNK = 1, NCHNK
            DO M = 1, NFRE_RED
              DO K = 1, NANG
                DO J = 1, NPROMA_WAM
                  IJSB = IJFROMCHNK(1, ICHNK)
                  IF (J <= KIJL4CHNK(ICHNK)) THEN
                    FL1(J, K, M, ICHNK) = FL3_EXT(IJSB + J - 1, K, M)
                  ELSE
                    FL1(J, K, M, ICHNK) = FL3_EXT(IJSB, K, M)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
#ifdef WAM_GPU
#ifdef OMPGPU
        !$omp end target teams distribute parallel do
#else
        !$acc end parallel loop
#endif
#else
!$OMP     END PARALLEL DO
#endif

           CALL GSTATS(1430,1)

        ENDIF  ! end propagation

      ENDIF ! more than one grid point
#ifdef OMPGPU
!$omp end target data
#else
!$acc end data
#endif

      L1STCALL=.FALSE.
      LLCHKCFL=.FALSE.

IF (LHOOK) CALL DR_HOOK('PROPAG_WAM',1,ZHOOK_HANDLE)

END SUBROUTINE PROPAG_WAM
