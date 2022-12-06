SUBROUTINE MPFLDTOIFS(IJS, IJL, BLK2GLO, NWVFIELDS, BLOCK,       &
 &                    GRID, DEFVAL, MASK_OUT, LLGLOBAL)

!****  *MPFLDTOIFS* - TRANSFORMS BLOCK DATA TO GRID DATA FOR
!****                 FIELDS THAT WILL BE RETURNED TO IFS.                

!     J. BIDLOT    ECMWF AUGUST 2008

!     PURPOSE.
!     --------

!     TRANSFORMS BLOCK DATA TO GRID DATA.

!*    INTERFACE.
!     ----------

!     CALL *MPFLDTOIFS(IJS, IJL, BLK2GLO, NWVFIELDS, BLOCK,
!    &                 GRID, DEFVAL, MASK_OUT )*
!         *IJS*        - BLOCK INDEX OF FIRST GRIDPOINT.
!         *IJL*        - BLOCK INDEX OF LAST GRIDPOINT.
!         *BLK2GLO*    - BLOCK TO GRID TRANSFORMATION
!         *NWVFIELDS*  - TOTAl NUMBER OF FIELDS RETURNED TO IFS
!         *BLOCK*      - FIELDS IN BLOCK FORM (INPUT) !     CALL *MPFLDTOIFS(IJS, IJL, BLK2GLO, NWVFIELDS, BLOCK,
!    &                 GRID, DEFVAL, MASK_OUT )*
!         *IJS*        - BLOCK INDEX OF FIRST GRIDPOINT.
!         *IJL*        - BLOCK INDEX OF LAST GRIDPOINT.
!
!                        ONLY DEFINED LOCALLY.
!         *GRID*       - FIELDS IN GRID FORM (OUTPUT)
!                        DEFINED  GLOBALLY OR JUST ON THE
!                        SUBAREA NEEDED (SEE MASK_OUT) DEPENDING
!                        ON LWCOUNORMS OR LMASK_OUT_NOT_SET
!         *DEFVAL*     - DEFAULT VALUE TO ASSIGN EACH FIELD WHEN
!                        NOT DEFINED BY WAVE MODEL.
!         *MASK_OUT*   - MASK POINTING TO VALUES OF FIELDS THAT ARE
!                        NEEDED BY IFS ON CURRENT TASK.
!         *LLGLOBAL*   - TRUE IF GRID HAS BEEN DEFINED GLOBALLY (OUTPUT).


!     METHOD.
!     -------
!     MPL ALLGATHERV OF ARRAY FIELD TO ALL PE's FOR FIRST CALL
!     FOR SUBSEQUENT CALLS COMMUNICATION TAKES PLACE WITH ONLY THE
!     TASKS THAT WE NEED TO SEND DATA TO.

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWCOUP  , ONLY : LWCOUNORMS,LMASK_OUT_NOT_SET,               &
     &            LMASK_TASK_STR,                                       &
     &            LFROMTASK,IJFROMTASK,NFROMTASKS,ISTFROMTASK,          &
     &            LTOTASK  ,IJTOTASK  ,NTOTASKS  ,ISTTOTASK
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_ALLGATHERV, MPL_SEND, MPL_RECV, MPL_WAIT, &
                           & JP_NON_BLOCKING_STANDARD, JP_BLOCKING_STANDARD
!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      TYPE(WVGRIDGLO), DIMENSION(NIBLO), INTENT(IN) :: BLK2GLO
      INTEGER(KIND=JWIM), INTENT(IN) :: NWVFIELDS
      REAL(KIND=JWRB), DIMENSION(IJS:IJL,NWVFIELDS), INTENT(IN) :: BLOCK
      REAL(KIND=JWRB), DIMENSION(:,:,:), INTENT(INOUT) :: GRID
      REAL(KIND=JWRB), DIMENSION(NWVFIELDS), INTENT(IN) :: DEFVAL 
      INTEGER(KIND=JWIM), INTENT(IN) :: MASK_OUT(NGX,NGY)
      LOGICAL, INTENT(OUT) :: LLGLOBAL


      INTEGER(KIND=JWIM) :: IFLD, IP, IP1, IJ, I, J, IX, IY
      INTEGER(KIND=JWIM) :: ITAG, IREQ, ICOUNT, IC, IST, IEND
      INTEGER(KIND=JWIM), DIMENSION(NPROC) :: KRECVCOUNTS
      INTEGER(KIND=JWIM) :: ISENDREQ(NPROC)

      REAL(KIND=JWRB), ALLOCATABLE :: ZBUFS(:) ,ZBUFR(:)
      REAL(KIND=JWRB), ALLOCATABLE :: ZSENDBUF(:), ZRECVBUF(:) 
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE1

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MPFLDTOIFS',0,ZHOOK_HANDLE)

ASSOCIATE(IXLG => BLK2GLO%IXLG, &
 &        KXLT => BLK2GLO%KXLT)

      CALL GSTATS_BARRIER(734)

      LLGLOBAL=(LWCOUNORMS .OR. LMASK_OUT_NOT_SET)

      CALL GSTATS(1503,0)
!     DEFAULT VALUES
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IFLD,J,I) COLLAPSE(2)
      DO IFLD=1,NWVFIELDS
        DO J = 1,NGY
          DO I = 1,NGX
            GRID(I,J,IFLD) = DEFVAL(IFLD) 
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1503,1)


      CALL GSTATS(686,0)

!     LOOP OVER ALL INPUT FIELDS

      DO IFLD=1,NWVFIELDS

        IF (NPROC > 1) THEN

!         GLOBAL EXCHANGE

          IF (LLGLOBAL) THEN

            ALLOCATE(ZBUFS(NIBLO))
            ALLOCATE(ZBUFR(NIBLO))

            DO IJ=IJS,IJL
              ZBUFS(IJ)=BLOCK(IJ,IFLD)
            ENDDO

            DO IP=1,NPROC
              KRECVCOUNTS(IP)=NEND(IP) - NSTART(IP)+1
            ENDDO

            CALL MPL_ALLGATHERV(ZBUFS(NSTART(IRANK):NEND(IRANK)),       &
     &                          ZBUFR(1:NIBLO),KRECVCOUNTS,             &
     &                          CDSTRING='MPFLDTOIFS:')

!           TRANSFORM FROM BLOCK TO GRID
            DO IP=1,NPROC
              DO IJ = NSTART(IP), NEND(IP)
               IX = IXLG(IJ)
               IY = NGY- KXLT(IJ) +1
               GRID(IX,IY,IFLD) = ZBUFR(IJ)
              ENDDO
            ENDDO

            DEALLOCATE(ZBUFS)
            DEALLOCATE(ZBUFR)

            LMASK_OUT_NOT_SET=.FALSE.

          ELSE

!           INITIALISE THE TASK STRUCTURE FOR OPTIMAL DATA EXCHANGE
!           FOR FIELDS RETURNED TO IFS.

            IF (LMASK_TASK_STR) THEN

              LMASK_TASK_STR=.FALSE.

              IF (ALLOCATED(LFROMTASK)) DEALLOCATE(LFROMTASK)
              ALLOCATE(LFROMTASK(NPROC))

              DO IP=1,NPROC
                LFROMTASK(IP)=0
                DO IJ = NSTART(IP), NEND(IP)
                  IX = IXLG(IJ)
                  IY = NGY- KXLT(IJ) +1
                  IF (MASK_OUT(IX,IY) == 1) LFROMTASK(IP)=LFROMTASK(IP)+1
                ENDDO
              ENDDO

              NFROMTASKS=SUM(LFROMTASK)
              IF (ALLOCATED(IJFROMTASK)) DEALLOCATE(IJFROMTASK)
              ALLOCATE(IJFROMTASK(NFROMTASKS))

              IF (ALLOCATED(ISTFROMTASK)) DEALLOCATE(ISTFROMTASK)
              ALLOCATE(ISTFROMTASK(NPROC))
              DO IP=1,NPROC
                ISTFROMTASK(IP)=0
              ENDDO

              ICOUNT=0
              DO IP=1,NPROC
                IF (LFROMTASK(IP) > 0) THEN
                  ISTFROMTASK(IP)=ICOUNT+1
                  DO IJ = NSTART(IP), NEND(IP)
                    IX = IXLG(IJ)
                    IY = NGY- KXLT(IJ) +1
                    IF (MASK_OUT(IX,IY) == 1) THEN
                      ICOUNT=ICOUNT+1
                      IJFROMTASK(ICOUNT)=IJ
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO

              IF (ALLOCATED(LTOTASK)) DEALLOCATE(LTOTASK)
              ALLOCATE(LTOTASK(NPROC))
              DO IP=1,NPROC
                LTOTASK(IP)=0
              ENDDO

!             SEND TO ALL OTHER TASKS TO LET THEM KNOW IF CONTRIBUTIONS
!             FROM THEM ARE NEEDED.
              ITAG=1
              DO IP=1,NPROC
                CALL MPL_SEND(LFROMTASK(IP),KDEST=IP,KTAG=ITAG,         &
     &                        KMP_TYPE=JP_NON_BLOCKING_STANDARD,        &
     &                        KREQUEST=ISENDREQ(IP),                    &
     &                        CDSTRING='MPFLDTOIFS: SEND COUNT ' )
              ENDDO
!             RECEIVE INFORMATION ON WHICH TASKS WILL NEED TO BE SENT
!             SOME CONTRIBUTIONS.
              DO IP=1,NPROC
                CALL MPL_RECV(LTOTASK(IP),KSOURCE=IP,KTAG=ITAG,         &
     &                        KMP_TYPE=JP_BLOCKING_STANDARD,            &
     &                        CDSTRING='MPFLDTOIFS: RECV COUNT ' )
              ENDDO

              CALL MPL_WAIT(KREQUEST=ISENDREQ(1:NPROC),                 &
     &                      CDSTRING='MPFLDTOIFS: WAIT SEND COUNT')


              NTOTASKS=SUM(LTOTASK)

              IF (ALLOCATED(IJTOTASK)) DEALLOCATE(IJTOTASK)
              ALLOCATE(IJTOTASK(NTOTASKS))

              IF (ALLOCATED(ISTTOTASK)) DEALLOCATE(ISTTOTASK)
              ALLOCATE(ISTTOTASK(NPROC))
              DO IP=1,NPROC
                ISTTOTASK(IP)=0
              ENDDO

!             SEND TO ALL OTHER TASKS TO LET THEM KNOW WHICH CONTRIBUTIONS
!             FROM THEM ARE NEEDED.
              ITAG=2
              IREQ=0

              DO IP=1,NPROC
                IF (LFROMTASK(IP) > 0) THEN
                  IREQ=IREQ+1

                  IST=ISTFROMTASK(IP)
                  IEND=IST+LFROMTASK(IP)-1
                  CALL MPL_SEND(IJFROMTASK(IST:IEND),                   &
     &                          KDEST=IP,KTAG=ITAG,                     &
     &                          KMP_TYPE=JP_NON_BLOCKING_STANDARD,      &
     &                          KREQUEST=ISENDREQ(IREQ),                &
     &                          CDSTRING='MPFLDTOIFS: SEND IJ ' )
                ENDIF
              ENDDO
!             RECEIVE INFORMATION ON WHAT CONTRIBUTIONS WILL NEED TO BE 
!             SENT.
              IST=1
              DO IP=1,NPROC
                IF (LTOTASK(IP) > 0) THEN
                  IEND=IST+LTOTASK(IP)-1
                  ISTTOTASK(IP)=IST
                  CALL MPL_RECV(IJTOTASK(IST:IEND),                     &
     &                          KSOURCE=IP,KTAG=ITAG,                   &
     &                          KMP_TYPE=JP_BLOCKING_STANDARD,          &
     &                          CDSTRING='MPFLDTOIFS: RECV IJ ' )
                  IST=IST+LTOTASK(IP)
                ENDIF
              ENDDO

              IF ( IREQ > 0 )THEN
                CALL MPL_WAIT(KREQUEST=ISENDREQ(1:IREQ),                &
     &                        CDSTRING='MPFLDTOIFS: WAIT SEND IJ')
              ENDIF

!             DO IP1=1,NPROC
!               IF (IRANK == IP1)THEN
!                 DO IP=1,NPROC
!                   IF ( LTOTASK(IP)>0 .OR. LFROMTASK(IP)>0 )THEN
!                     WRITE(0,                                    &
!    &                  '("MPFLDTOIFS: IP=",I6," LTOTASK=",I10,   &
!    &                  " LFROMTASK=",I10)') IP,LTOTASK(IP),      &
!    &                  LFROMTASK(IP)
!                   ENDIF
!                 ENDDO
!               ENDIF
!               CALL MPL_BARRIER(CDSTRING='MPFLDTOIFS')
!             ENDDO

            ENDIF


!           EXCHANGE INFORMATION

!           SEND LOCAL TASK CONTRIBUTION TO TASKS THAT NEED IT

            ALLOCATE(ZSENDBUF(MAX(1,NTOTASKS)))
            DO IP=1,NPROC
              IF (LTOTASK(IP) > 0) THEN
                IST=ISTTOTASK(IP)
                IEND=IST+LTOTASK(IP)-1
                DO IC=IST,IEND
                  ZSENDBUF(IC)=BLOCK(IJTOTASK(IC),IFLD)
                ENDDO
              ENDIF
            ENDDO

            ITAG=3
            IREQ=0

            DO IP=1,NPROC
              IF (LTOTASK(IP) > 0) THEN
                IREQ=IREQ+1
                IST=ISTTOTASK(IP)
                IEND=IST+LTOTASK(IP)-1
                CALL MPL_SEND(ZSENDBUF(IST:IEND),                       &
     &                        KDEST=IP,KTAG=ITAG,                       &
     &                        KMP_TYPE=JP_NON_BLOCKING_STANDARD,        &
     &                        KREQUEST=ISENDREQ(IREQ),                  &
     &                        CDSTRING='MPFLDTOIFS: SEND DATA')
              ENDIF
            ENDDO

            ALLOCATE(ZRECVBUF(NFROMTASKS))

!           RECEIVE INFORMATION FROM OTHER TASKS (IF NEEDED)
            DO IP=1,NPROC

              IF (LFROMTASK(IP) > 0) THEN
                IST=ISTFROMTASK(IP)
                IEND=IST+LFROMTASK(IP)-1
                CALL MPL_RECV(ZRECVBUF(IST:IEND),                       &
     &                        KSOURCE=IP,KTAG=ITAG,                     &
     &                        KMP_TYPE=JP_BLOCKING_STANDARD,            &
     &                        CDSTRING='MPFLDTOIFS: RECV DATA ')

              ENDIF
            ENDDO

!           ENSURE ALL SENDS ARE FINISHED.

            IF (IREQ > 0) THEN
              CALL MPL_WAIT(KREQUEST=ISENDREQ(1:IREQ),                  &
     &                      CDSTRING='MPFLDTOIFS: WAIT SEND DATA ')
            ENDIF

            DO IP=1,NPROC
              IF (LFROMTASK(IP) > 0) THEN
                IST=ISTFROMTASK(IP)
                IEND=IST+LFROMTASK(IP)-1
!               TRANSFORM FROM BLOCK TO GRID
                DO IC = IST,IEND
                  IJ=IJFROMTASK(IC)
                  IX = IXLG(IJ)
                  IY = NGY- KXLT(IJ) +1
                  GRID(IX,IY,IFLD) = ZRECVBUF(IC)
                ENDDO
              ENDIF
            ENDDO

            DEALLOCATE(ZSENDBUF)
            DEALLOCATE(ZRECVBUF)

          ENDIF

        ELSE ! ONE PE

!         TRANSFORM FROM BLOCK TO GRID
          DO IP=1,NPROC
            DO IJ = NSTART(IP), NEND(IP)
             IX = IXLG(IJ)
             IY = NGY- KXLT(IJ) +1
             GRID(IX,IY,IFLD) = BLOCK(IJ,IFLD)
            ENDDO
          ENDDO

        ENDIF

      ENDDO

      CALL GSTATS(686,1)
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MPFLDTOIFS',1,ZHOOK_HANDLE)

END SUBROUTINE MPFLDTOIFS