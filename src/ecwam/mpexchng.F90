! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPEXCHNG(FLD, NDIM2, ND3S, ND3E)


!****  *MPEXCHNG* - EXHANGES MESSAGE BETWEEN THE PROCESS IRANK
!****               AND THE ALL PE's THAT CONTAINS INFORMATION NEEDED
!                   FOR THE ADVECTION OR THE CALCULATION OF THE
!                   CURRENT VELOCITY GRADIENT. 

!     J. BIDLOT    ECMWF   MARCH 1996  MESSAGE PASSING
!                  ECMWF   JANUARY 2004 MADE THE INTERFACE MORE GENERAL

!     PURPOSE.
!     --------
!     EXHANGE MESSAGE BETWEEN ONE PROCESS AND ITS NEIGHBOURS (HALO). 

!*    INTERFACE.
!     ----------
!     CALL *MPEXCHNG*(FLD, NDIM2, ND2S, ND3E)

!      *FLD*       ARRAY TO EXCHANGE 
!                  THE FIRST DIMENSION IS ALWAYS GIVEN AS NINF:NSUP+1
!                  BECAUSE OF THE HALO CONFIGURATION (SEE MPDECOMP)
!      *NDIM2*     SECOND DIMENSION OF FLD.
!      *ND3S*      THIRD DIMENSION: ND3S:ND3E
!      *ND3E*      

!     METHOD.
!     -------

!     NOW UPDATED TO USE NON-BLOCKING SENDS AND RECEIVES


!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_SEND
!         MPL_RECV
!         MPL_WAIT

!     REFERENCES.
!     -----------
!     CHAPTER 4 OF
!     USING MPI, PORTABLE PARALLEL PROGRAMMING WITH THE MESSAGE PASSING
!     INTERFACE. W.CROPP, E LUSK, A SKJELLUM. MIT PRESS 1995

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,    &
     &          KTAG
      USE YOWSPEC, ONLY   : NTOPE    ,NTOPEMAX ,IJTOPE   ,NGBTOPE  ,    &
     &          NTOPELST   ,NFROMPE  ,NFROMPEMAX,NIJSTART,NGBFROMPE,    &
     &          NFROMPELST

      USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_RECV, MPL_SEND, MPL_WAIT, &
                           & JP_NON_BLOCKING_STANDARD, MPL_COMM_OML
#ifdef WITH_GPU_AWARE_MPI
      USE OML_MOD   , ONLY : OML_MY_THREAD
      USE MPI_F08   , ONLY : MPI_ISEND, MPI_IRECV, MPI_COMM, MPI_REQUEST
#ifdef WAM_HAVE_SINGLE_PRECISION
      USE MPI_F08   , ONLY : ECWAM_MPI_DATATYPE => MPI_REAL4
#else
      USE MPI_F08   , ONLY : ECWAM_MPI_DATATYPE => MPI_REAL8
#endif
#endif

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: NDIM2, ND3S, ND3E
      REAL(KIND=JWRB), DIMENSION(NINF:NSUP+1, NDIM2, ND3S:ND3E), INTENT(INOUT) :: FLD

      INTEGER(KIND=JWIM) :: INGB, M, K, IH, IJ, KCOUNT, IPROC, IR
      INTEGER(KIND=JWIM) :: NBUFMAX
      INTEGER(KIND=JWIM) :: NDIM3
      INTEGER(KIND=JWIM), DIMENSION(NGBTOPE+NGBFROMPE) :: IREQ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: ZDUM(2)
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFS(:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFR(:,:)
#ifdef WITH_GPU_AWARE_MPI
      TYPE(MPI_COMM) :: ICOMM
      TYPE(MPI_REQUEST) :: IREQUEST_LOCAL
#endif

      LOGICAL :: LLOK
      INTEGER(KIND=JWIM) :: IERROR

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MPEXCHNG',0,ZHOOK_HANDLE)

      IF (NPROC <= 1) THEN
        IF (LHOOK) CALL DR_HOOK('MPEXCHNG',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

      CALL GSTATS_BARRIER(736)

      NDIM3 = ND3E-ND3S+1

      NBUFMAX=MAX(NTOPEMAX,NFROMPEMAX)*NDIM2*NDIM3
      ALLOCATE(ZCOMBUFS(NBUFMAX,NGBTOPE))
      ALLOCATE(ZCOMBUFR(NBUFMAX,NGBFROMPE))
!$acc enter data create(ZCOMBUFS,ZCOMBUFR)

!     PACK SEND BUFFERS FOR NGBTOPE NEIGHBOURING PE's
!     -------------------------------------------------
      CALL GSTATS(1892,0)
#ifdef _OPENACC
!$acc kernels loop independent private(IPROC) present(ZCOMBUFS,FLD) &
!$acc copyin(NTOPELST,NTOPE,IJTOPE)
      DO INGB=1,NGBTOPE !Total number of PE's to which information will be sent
        IPROC=NTOPELST(INGB)  !To which PE to send informations
          !$acc loop independent collapse(3) private(IJ,KCOUNT,M,K,IH)
          DO M = ND3S, ND3E
            DO K = 1, NDIM2
              DO IH = 1, NTOPE(IPROC) !How many halo points to be sent
                IJ=IJTOPE(IH,IPROC) !The index of which points to send
                KCOUNT = (M - ND3S) * (NDIM2 * NTOPE(IPROC)) + (K - 1) * NTOPE(IPROC) + IH
                ZCOMBUFS(KCOUNT,INGB)=FLD(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$acc end kernels
#else
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(INGB,IPROC,KCOUNT,M,K,IH,IJ)
       DO INGB=1,NGBTOPE
         IPROC=NTOPELST(INGB)
         KCOUNT=0
         DO M = ND3S, ND3E
           DO K = 1, NDIM2
             DO IH = 1, NTOPE(IPROC)
               IJ=IJTOPE(IH,IPROC)
               KCOUNT=KCOUNT+1
               ZCOMBUFS(KCOUNT,INGB)=FLD(IJ,K,M)
             ENDDO
           ENDDO
         ENDDO
       ENDDO
!$OMP END PARALLEL DO
#endif /*_OPENACC*/

      CALL GSTATS(1892,1)

!     DO NON BLOCKING SENDS AND RECVS

      IR=0
      CALL GSTATS(676,0)

      DO INGB=1,NGBFROMPE
        IR=IR+1
        IPROC=NFROMPELST(INGB)
        KCOUNT=NDIM3*NDIM2*NFROMPE(IPROC)
#ifdef WITH_GPU_AWARE_MPI
        ICOMM%MPI_VAL=MPL_COMM_OML(OML_MY_THREAD())
!$acc host_data use_device(ZCOMBUFR)
        CALL MPI_IRECV(ZCOMBUFR(1:KCOUNT,INGB),KCOUNT,                 &
     &     ECWAM_MPI_DATATYPE,IPROC-1, KTAG,                           &
     &     ICOMM,IREQUEST_LOCAL, IERROR)
!$acc end host_data
        IREQ(IR) = IREQUEST_LOCAL%MPI_VAL
#else
        CALL MPL_RECV(ZCOMBUFR(1:KCOUNT,INGB),KSOURCE=IPROC,KTAG=KTAG,  &
     &     KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),         &
     &     CDSTRING='MPEXCHNG:')
#endif
      ENDDO

      DO INGB=1,NGBTOPE
        IR=IR+1
        IPROC=NTOPELST(INGB)
        KCOUNT=NDIM3*NDIM2*NTOPE(IPROC)

#ifdef WITH_GPU_AWARE_MPI
        ICOMM%MPI_VAL=MPL_COMM_OML(OML_MY_THREAD())
!$acc host_data use_device(ZCOMBUFS)
        CALL MPI_ISEND(ZCOMBUFS(1:KCOUNT,INGB),KCOUNT,                 &
     &     ECWAM_MPI_DATATYPE,IPROC-1, KTAG,                           &
     &     ICOMM,IREQUEST_LOCAL, IERROR)
!$acc end host_data
        IREQ(IR) = IREQUEST_LOCAL%MPI_VAL
#else              
!$acc update self(ZCOMBUFS)
        CALL MPL_SEND(ZCOMBUFS(1:KCOUNT,INGB),KDEST=IPROC,KTAG=KTAG,    &
     &     KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),         &
     &     CDSTRING='MPEXCHNG:')
#endif
      ENDDO

!     NOW WAIT FOR ALL TO COMPLETE

      CALL MPL_WAIT(KREQUEST=IREQ(1:IR),CDSTRING='MPEXCHNG:')
#ifndef WITH_GPU_AWARE_MPI
!$acc update device(ZCOMBUFR)
#endif

      CALL GSTATS(676,1)

!     DECODE THE RECEIVED BUFFERS

      CALL GSTATS(1893,0)
#ifdef _OPENACC
      !$acc kernels loop independent private(IPROC) present(ZCOMBUFR,FLD) &
      !$acc copyin(NFROMPELST,NFROMPE,NIJSTART)
      DO INGB=1,NGBFROMPE
        IPROC=NFROMPELST(INGB)
        !$acc loop vector independent collapse(3) private(IJ,KCOUNT,M,K,IH)
        DO M = ND3S, ND3E
          DO K = 1, NDIM2
            DO IH = 1, NFROMPE(IPROC)
              IJ=NIJSTART(IPROC)+IH-1
              KCOUNT = (M - ND3S) * (NDIM2 * NFROMPE(IPROC)) + (K - 1) * NFROMPE(IPROC) + IH
              FLD(IJ,K,M)=ZCOMBUFR(KCOUNT,INGB)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !$acc end kernels
#else
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(INGB,IPROC,KCOUNT,M,K,IH,IJ)
      DO INGB=1,NGBFROMPE
        IPROC=NFROMPELST(INGB)
        KCOUNT=0
        DO M = ND3S, ND3E
          DO K = 1, NDIM2
            DO IH = 1, NFROMPE(IPROC)
              IJ=NIJSTART(IPROC)+IH-1
              KCOUNT=KCOUNT+1
              FLD(IJ,K,M)=ZCOMBUFR(KCOUNT,INGB)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
#endif /*_OPENACC*/
      CALL GSTATS(1893,1)

      KTAG=KTAG+1

!$acc exit data delete(ZCOMBUFS,ZCOMBUFR)
      DEALLOCATE(ZCOMBUFS)
      DEALLOCATE(ZCOMBUFR)

      IF (LHOOK) CALL DR_HOOK('MPEXCHNG',1,ZHOOK_HANDLE)

      END SUBROUTINE MPEXCHNG
