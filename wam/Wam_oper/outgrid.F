      SUBROUTINE OUTGRID
! ----------------------------------------------------------------------

!**** *OUTGRID* - OUTPUT OF WAVE ENERGY,MEAN ANGLE, MEAN FREQUENCY
!****             FRICTION VELOCITY AND WIND DIRECTION IN GRIDDED
!****             FORMAT.

!*    PURPOSE.
!     --------

!       OUTPUT OF WAVE AND WIND FIELDS INTO COMMON BLOCK.


!**   INTERFACE.
!     ----------

!        *CALL* *OUTGRID

!     EXTERNALS.
!     ----------

!       *MAKEGRID*  - GENERATES GRIDDED FIELDS.

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : JPPFLAG  ,IPFGTBL ,NIPRMOUT ,ITOBOUT  ,BOUT
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC
      USE YOWINTP  , ONLY : GOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,MPMAXLENGTH ,KTAG
      USE YOWPARAM , ONLY : NIBLO    ,NGX      ,NGY
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSPEC  , ONLY : NBLKS    ,NBLKE
      USE YOWTEST  , ONLY : IU06
      USE OUTPUT_STRUCT, ONLY : ARR_OUT_RECV, LocalPosICT
      USE OUTPUT_STRUCT, ONLY : INITIAL_OUTPUT_INITS_NEXTGEN
      USE OUTPUT_STRUCT, ONLY : SET_UP_ARR_OUT_RECV
      USE YOWUNPOOL, ONLY : NIBLO_FD, NIBLO_OUT, OUT_METHOD, LLUNSTR

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "makegrid.intfb.h"

      INTEGER :: IG
      INTEGER :: IJ, ITT, ICT, ITG, IFLD, IR, IPR, ICOUNT
      INTEGER :: NFLDTOT, NFLDPPEMAX
      INTEGER, DIMENSION(NPROC) :: ICNT, NFLDPPE
      INTEGER, DIMENSION(NPROC+JPPFLAG) :: IREQ
      INTEGER, DIMENSION(NPROC) :: ISENDCOUNTS,IRECVCOUNTS

      REAL :: ZHOOK_HANDLE
      REAL :: ZDUM(2)
      REAL, ALLOCATABLE, DIMENSION(:) :: GTEMP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ZSENDBUF, ZRECVBUF 

      LOGICAL, SAVE :: HaveMPI_arrays = .FALSE.
!----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTGRID',0,ZHOOK_HANDLE)

      IG=1

!     FIND THE NUMBER OF FIELDS EACH PE HAS TO DEAL WITH
      NFLDPPE(:)=0
      NFLDTOT=0
      NFLDPPEMAX=1
      DO ICT=1,JPPFLAG
        IF(IPFGTBL(ICT).NE.0) THEN 
          NFLDPPE(IPFGTBL(ICT))=NFLDPPE(IPFGTBL(ICT))+1
          NFLDTOT=NFLDTOT+1
          NFLDPPEMAX=MAX(NFLDPPEMAX,NFLDPPE(IPFGTBL(ICT)))
        ENDIF
      ENDDO
      IF (HaveMPI_arrays .eqv. .FALSE.) THEN
        IF (LLUNSTR .and. (OUT_METHOD .eq. 2)) THEN
           CALL INITIAL_OUTPUT_INITS_NEXTGEN
        END IF
      END IF
      HaveMPI_arrays=.TRUE.
      IF(NFLDTOT.EQ.0) THEN
#ifdef ECMWF
        IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)
#endif
        RETURN
      ENDIF
 

!     SENDING TO RELEVANT PE'S
!     ------------------------

      DO IPR=1,NPROC
        ISENDCOUNTS(IPR)=NFLDPPE(IPR)*(NBLKE(IRANK)-NBLKS(IRANK)+1)
      ENDDO

      DO IPR=1,NPROC
        IRECVCOUNTS(IPR)=NFLDPPE(IRANK)*(NBLKE(IPR)-NBLKS(IPR)+1)
      ENDDO

      IF ((.NOT.LLUNSTR).OR.(OUT_METHOD.eq.1)) THEN
      
!     LOADING THE COMMUNICATION BUFFER
        ALLOCATE(ZSENDBUF(NFLDPPEMAX * MPMAXLENGTH,NPROC))

        ICNT(:)=0
        DO ICT=1,JPPFLAG
          IPR=IPFGTBL(ICT)
          IF(IPR.GT.0) THEN 
            DO IJ=IJSLOC,IJLLOC
              ICNT(IPR) = ICNT(IPR) + 1
              ZSENDBUF(ICNT(IPR),IPR) = BOUT(IJ,ITOBOUT(ICT))
            ENDDO
          ENDIF
        ENDDO
      ELSE
!!!!! this will need to be adapted to use BOUT
        WRITE(0,*) '!!! ********************************* !!'
        WRITE(0,*) '!!! in outgrid. Not yet ready !!!' 
        WRITE(0,*) '!!! ********************************* !!'
        CALL ABORT1

        CALL SET_UP_ARR_OUT_RECV(IJSLOC, IJLLOC, BOUT, NFLDPPE)
      END IF


!     GLOBAL EXCHANGE

      IF ((.NOT.LLUNSTR).OR.(OUT_METHOD.eq.1)) THEN
        CALL GSTATS(693,0)
        IR=0
        IF(NFLDPPE(IRANK).GT.0) THEN
          ALLOCATE(ZRECVBUF(NFLDPPE(IRANK)*MPMAXLENGTH,NPROC))
          DO IPR=1,NPROC
            IR=IR+1
            CALL MPL_RECV(ZRECVBUF(1:IRECVCOUNTS(IPR),IPR),KSOURCE=IPR,
     &        KTAG=KTAG,
     &        KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),
     &        CDSTRING='OUTGRID:')
          ENDDO
        ENDIF
        DO IPR=1,NPROC
          IF(NFLDPPE(IPR).GT.0) THEN
            IR=IR+1
            CALL MPL_SEND(ZSENDBUF(1:ICNT(IPR),IPR),KDEST=IPR,
     &        KTAG=KTAG,
     &        KMP_TYPE=JP_NON_BLOCKING_STANDARD,KREQUEST=IREQ(IR),
     &        CDSTRING='OUTGRID:')
          ENDIF
        ENDDO
!       NOW WAIT FOR ALL TO COMPLETE
        CALL MPL_WAIT(KREQUEST=IREQ(1:IR),CDSTRING='OUTGRID:')
        CALL GSTATS(693,1)
      END IF


!     RETRIEVE THE INFORMATION
!     ------------------------

      IF (LLUNSTR .and. (OUT_METHOD .eq. 2)) THEN
        NIBLO_OUT = NIBLO_FD
      ELSE
        NIBLO_OUT = NIBLO
      END IF
      ALLOCATE(GTEMP(NIBLO_OUT))


!     RECEIVING AND TRANSFERING FROM BLOCK TO GRID
!     --------------------------------------------

      IF(NFLDPPE(IRANK).GT.0) THEN
        IF(ALLOCATED(GOUT)) DEALLOCATE(GOUT)
        ALLOCATE(GOUT(NFLDPPE(IRANK),NGX,NGY))
      ENDIF

      ICNT(:)=0
      ICT=1
      IFLD=1
      DO WHILE ( IFLD.LE.NFLDPPE(IRANK) .AND. ICT .LE. JPPFLAG )
        IF(IPFGTBL(ICT).EQ.IRANK ) THEN
          IF (.NOT.(LLUNSTR) .OR. (OUT_METHOD .eq. 1)) THEN
            DO IPR=1,NPROC
              ICOUNT=ICNT(IPR)
              DO IJ=NBLKS(IPR),NBLKE(IPR)
                ICOUNT=ICOUNT+1
                GTEMP(IJ)=ZRECVBUF(ICOUNT,IPR)
              ENDDO
              ICNT(IPR)=ICOUNT
            ENDDO
          ELSE
            GTEMP=ARR_OUT_RECV(LocalPosICT(ICT),:)
          END IF

          CALL MAKEGRID (GTEMP,GOUT(IFLD,:,:),IG,ZMISS)
          IFLD=IFLD+1

        ENDIF ! (IPFGTBL) 

        ICT=ICT+1
      ENDDO  ! WHILE

      DEALLOCATE(GTEMP)
      IF(ALLOCATED(ZRECVBUF)) DEALLOCATE(ZRECVBUF)
      IF(ALLOCATED(ZSENDBUF)) DEALLOCATE(ZSENDBUF)

      IF (LHOOK) CALL DR_HOOK('OUTGRID',1,ZHOOK_HANDLE)

      END SUBROUTINE OUTGRID
