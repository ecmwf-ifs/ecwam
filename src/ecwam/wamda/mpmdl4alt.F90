      SUBROUTINE MPMDL4ALT(IJS, IJL, HSMOD, CICOVER, &
     &                     MINIJS, MAXIJL, WHMOD, CICVR)
! ----------------------------------------------------------------------

!****  *MPMDL4ALT* - GATHER MODEL VALUES NEEDED FOR ALTIMETER ASSIMILATION 
!                    ON GIVEN PE.


!     PURPOSE.                                                          
!     --------                                                          
!     GATHERS MODEL VALUES NEEDED FOR ALTIMETER ASSIMILATION ON PE IRANK.

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *MPMDL4ALT(IJS, IJL, HSMOD, CICOVER, &
!    &                  MINIJS, MAXIJL, WHMOD, CICVR)

!                       HSMOD, CICOVER HAVE LOCAL POINT INDEX 
!                       WHMOD, CICVR HAVE GLOBAL POINT INDEX
!     METHOD.                                                           
!     -------                                                           

!     EXTERNALS.                                                        
!     ----------                                                        

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWABORT , ONLY : WAM_ABORT
      USE YOWALTAS , ONLY : INTLMAX
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG    ,MPMAXLENGTH
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWPARAM , ONLY : LLUNSTR
#ifdef WAM_HAVE_UNWAM
      USE YOWPD    , ONLY : RANK
#endif
      USE MPL_MODULE,ONLY : MPL_SEND ,MPL_RECV ,MPL_WAIT, JP_NON_BLOCKING_STANDARD, &
                          & MPL_ABORT
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL
      INTEGER(KIND=JWIM), INTENT(IN) :: MINIJS, MAXIJL

      REAL(KIND=JWRB),DIMENSION(IJS:IJL), INTENT(IN) :: HSMOD, CICOVER
      REAL(KIND=JWRB),DIMENSION(MINIJS:MAXIJL), INTENT(OUT) :: WHMOD, CICVR

      INTEGER(KIND=JWIM) :: IJ, IR, IPR, IP
      INTEGER(KIND=JWIM) :: NDISPE, IDISPE
      INTEGER(KIND=JWIM) :: MBUFLEN 
      INTEGER(KIND=JWIM) :: ICOUNT, KRCOUNT, KRTAG
      INTEGER(KIND=JWIM),DIMENSION(NPROC) :: ISENDLEN
      INTEGER(KIND=JWIM),ALLOCATABLE :: ISENDREQ(:)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),ALLOCATABLE :: ZCOMBUFR(:)
      REAL(KIND=JWRB),ALLOCATABLE :: ZCOMBUFS(:,:)

!     -------------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('MPMDL4ALT',0,ZHOOK_HANDLE)

      WHMOD(MINIJS:MAXIJL)=0.0_JWRB
      CICVR(MINIJS:MAXIJL)=0.0_JWRB


!     DISTRIBUTE MODEL VALUES TO WHERE IT IS NEEDED

      MBUFLEN=2*MPMAXLENGTH

      IF (NPROC > 1) THEN
        NDISPE=0
        DO IR = 1, NPROC
          IF (INTLMAX(IR) == 1) THEN
            NDISPE=NDISPE+1
            ISENDLEN(IR)=MBUFLEN
          ELSE
            ISENDLEN(IR)=0
          ENDIF
        ENDDO
      ENDIF


      IF (NPROC > 1 .AND. NDISPE > 0) THEN
        KTAG=KTAG+1
!       PACKING THE CONTRIBUTION WHICH MIGHT BE COMMON
        ALLOCATE(ZCOMBUFS(MBUFLEN,NDISPE))
        ALLOCATE(ISENDREQ(NDISPE))
        IDISPE=0
        DO IR = 1, NPROC
          IF (INTLMAX(IR) == 1) THEN
            IDISPE=IDISPE+1
            ICOUNT = 0
            DO IJ = IJS, IJL
              ICOUNT = ICOUNT + 1
              ZCOMBUFS(ICOUNT,IDISPE) = HSMOD(IJ)
            ENDDO
            DO IJ = IJS, IJL
              ICOUNT = ICOUNT + 1
              ZCOMBUFS(ICOUNT,IDISPE) = CICOVER(IJ)
            ENDDO
          ENDIF
        ENDDO

!       SENDING THE CONTRIBUTION WHICH MIGHT BE COMMON
!       SEND NON BLOCKING THE BUFFERS

        IDISPE=0
        CALL GSTATS(618,0)
        DO IR = 1,NPROC
          IF (INTLMAX(IR) == 1) THEN
            IDISPE=IDISPE+1
            ICOUNT = ISENDLEN(IR) 
            CALL MPL_SEND(ZCOMBUFS(1:ICOUNT,IDISPE),                     &
     &                    KDEST=IR,KTAG=KTAG,                            &
     &                    KMP_TYPE=JP_NON_BLOCKING_STANDARD,             &
     &                    KREQUEST=ISENDREQ(IDISPE),                     &
     &                    CDSTRING='MPMDL4ALT:')
          ENDIF
        ENDDO
        CALL GSTATS(618,1)

!       RECEIVING THE CONTRIBUTIONS WHICH MIGHT BE COMMON
!       (in whatever order they arrive.)
        ALLOCATE(ZCOMBUFR(MBUFLEN))
        DO IPR =1,NDISPE
          CALL GSTATS(618,0)
          CALL MPL_RECV(ZCOMBUFR(1:MBUFLEN),KFROM=IR,KTAG=KTAG,         &
     &           KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='MPMDL4ALT:')
          CALL GSTATS(618,1)
          IF (KRTAG /= KTAG) CALL MPL_ABORT                             &
     &    ('MPL_RECV ERROR in MPMDL4ALT: MISMATCHED TAGS' )

          ICOUNT = 0
          IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
            DO IP=1,RANK(IR)%NP
              ICOUNT = ICOUNT + 1
              WHMOD(RANK(IR)%IPLG(IP))=ZCOMBUFR(ICOUNT)
            ENDDO
            DO IP=1,RANK(IR)%NP
              ICOUNT = ICOUNT + 1
              CICVR(RANK(IR)%IPLG(IP))=ZCOMBUFR(ICOUNT)
            ENDDO
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
          ELSE
            DO IJ = NSTART(IR), NEND(IR)
              ICOUNT = ICOUNT + 1
              WHMOD(IJ)=ZCOMBUFR(ICOUNT)
            ENDDO
            DO IJ = NSTART(IR), NEND(IR)
              ICOUNT = ICOUNT + 1
              CICVR(IJ)=ZCOMBUFR(ICOUNT)
            ENDDO
          ENDIF
        ENDDO
        DEALLOCATE(ZCOMBUFR)

!       WAIT ANY OUTSTANDING SENDS TO COMPLETE
        CALL GSTATS(618,0)
        IF (MBUFLEN > 0) THEN
          CALL MPL_WAIT(KREQUEST=ISENDREQ,                              &
     &                  CDSTRING='MPMDL4ALT')
        ENDIF
        CALL GSTATS(618,1)
        IF (ALLOCATED(ISENDREQ)) DEALLOCATE(ISENDREQ)
        KTAG = KTAG + 1
        IF (ALLOCATED(ZCOMBUFS)) DEALLOCATE(ZCOMBUFS)

      ELSE
        DO IJ = IJSLOC, IJLLOC
          WHMOD(IJ+IJGLOBAL_OFFSET)=HSMOD(IJ)
        ENDDO
        DO IJ = IJSLOC, IJLLOC
          CICVR(IJ+IJGLOBAL_OFFSET)=CICOVER(IJ)
        ENDDO
      ENDIF

      IF (LHOOK) CALL DR_HOOK('MPMDL4ALT',1,ZHOOK_HANDLE)

      END SUBROUTINE MPMDL4ALT
