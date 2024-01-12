! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GETBOBSTRCT(IREAD, NPR, MAXLEN, NEWIJ2IJ)

!****  *GETBOBSTRCT* - DETERMINES OBSTRUCTION COEFFICIENTS FROM BINARY INPUT 

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWMAP   , ONLY : NIBLO
      USE YOWMPP   , ONLY : IRANK
      USE YOWPARAM , ONLY : NFRE_RED ,LL1D
      USE YOWSTAT  , ONLY : IPROPAGS ,LSUBGRID
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : OBSLAT   ,OBSLON   ,OBSCOR   ,OBSRLAT  ,OBSRLON 
      USE YOWUNIT  , ONLY : IU08

      USE MPL_MODULE, ONLY : MPL_SCATTERV
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD         !! PE READING THE GRIB FILE
      INTEGER(KIND=JWIM), INTENT(IN) :: NPR     !! NUMBER OF SUBDOMAINS (USUALLY THE NUMBER OF PE'S )
      INTEGER(KIND=JWIM), INTENT(IN) :: MAXLEN  !! MAXIMUM NUMBER OF POINTS IN ANY SUB DOMAIN
      INTEGER(KIND=JWIM), DIMENSION(0:NIBLO), INTENT(IN) :: NEWIJ2IJ 

      INTEGER(KIND=JWIM) :: IJ, M, K, IP, IC, ICC, NIJ 
      INTEGER(KIND=JWIM) :: ITAG
      INTEGER(KIND=JWIM) :: MPLENGTH, KCOUNT
      INTEGER(KIND=JWIM), DIMENSION(NPR) :: ICOUNTS
      INTEGER(KIND=JWIM), ALLOCATABLE :: ICOMBUF_S(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: ICOMBUF_R(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IDUM(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: KDUM(:)
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: KOBSLON, KOBSLAT
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: KOBSCOR
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: KOBSRLON, KOBSRLAT

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GETBOBSTRCT',0,ZHOOK_HANDLE)

        ITAG = 0

      IF (LSUBGRID) THEN
!       READ IU08
!       =========

        IF (.NOT.LL1D .AND. NPR > 1) ALLOCATE(KDUM(NIBLO))

        MPLENGTH=MAXLEN
        ALLOCATE(ICOMBUF_S(MPLENGTH*NPR))
        ALLOCATE(ICOMBUF_R(MPLENGTH))

        ALLOCATE(IDUM(NIBLO))

        ALLOCATE(KOBSLON(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
        ALLOCATE(KOBSLAT(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
        IF (IPROPAGS == 1) THEN
          ALLOCATE(KOBSRLAT(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
          ALLOCATE(KOBSRLON(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
        ENDIF
        IF (IPROPAGS == 2) THEN
          ALLOCATE(KOBSCOR(NSTART(IRANK):NEND(IRANK),NFRE_RED,4))
        ENDIF


        DO M=1,NFRE_RED ! loop over frequencies

!         READING KOBSLAT
!         ---------------
          DO IC=1,2
            ITAG=ITAG+1
            IF (IRANK == IREAD) THEN

              CALL GSTATS(1771,0)
              READ (IU08(IPROPAGS)) IDUM 
              CALL GSTATS(1771,1)

!             RELABELLING OF THE ARRAY
              IF (.NOT.LL1D .AND. NPR > 1) THEN
                DO NIJ=NSTART(1),NEND(NPR)
                  KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                ENDDO
                DO NIJ=NSTART(1),NEND(NPR)
                  IDUM(NIJ)=KDUM(NIJ)
                ENDDO
              ENDIF

!             FILL THE SEND BUFFER
              DO IP=1,NPR
                KCOUNT=(IP-1)*MPLENGTH
                DO IJ=NSTART(IP),NEND(IP)
                  KCOUNT=KCOUNT+1
                  ICOMBUF_S(KCOUNT)=IDUM(IJ)
                ENDDO
              ENDDO
            ENDIF

            IF (NPR > 1) THEN
              CALL GSTATS(694,0)
              ICOUNTS(:)=MPLENGTH
              CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD, KSENDBUF=ICOMBUF_S,     &
     &                          KSENDCOUNTS=ICOUNTS,CDSTRING='GETBOBSTRCT 1:')
              CALL GSTATS(694,1)
            ENDIF

!           KEEP THE RELEVANT PART OF KOBSLAT
            IF (IRANK == IREAD) THEN
              KCOUNT=(IRANK-1)*MPLENGTH
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                KOBSLAT(IJ,M,IC)=ICOMBUF_S(KCOUNT)
              ENDDO
            ELSE
              KCOUNT=0
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                KOBSLAT(IJ,M,IC)=ICOMBUF_R(KCOUNT)
              ENDDO
            ENDIF

          ENDDO

!         READING KOBSLON
!         ---------------

          DO IC=1,2

            ITAG=ITAG+1
            IF (IRANK == IREAD) THEN
              CALL GSTATS(1771,0)
              READ (IU08(IPROPAGS)) IDUM 
              CALL GSTATS(1771,1)

!             RELABELLING OF THE ARRAY
              IF (.NOT.LL1D .AND. NPR > 1) THEN
                DO NIJ=NSTART(1),NEND(NPR)
                   KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                ENDDO
                DO NIJ=NSTART(1),NEND(NPR)
                  IDUM(NIJ)=KDUM(NIJ)
                ENDDO
              ENDIF
!             FILL THE SEND BUFFER
              DO IP=1,NPR
                KCOUNT=(IP-1)*MPLENGTH
                DO IJ=NSTART(IP),NEND(IP)
                  KCOUNT=KCOUNT+1
                  ICOMBUF_S(KCOUNT)=IDUM(IJ)
                ENDDO
              ENDDO
            ENDIF

            IF (NPR > 1) THEN
              CALL GSTATS(694,0)
              ICOUNTS(:)=MPLENGTH
              CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD, KSENDBUF=ICOMBUF_S,     & 
     &                          KSENDCOUNTS=ICOUNTS,CDSTRING='GETBOBSTRCT 2:')
              CALL GSTATS(694,1)
            ENDIF

!           KEEP THE RELEVANT PART OF KOBSLON
            IF (IRANK == IREAD) THEN
              KCOUNT=(IRANK-1)*MPLENGTH
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                KOBSLON(IJ,M,IC)=ICOMBUF_S(KCOUNT)
              ENDDO
            ELSE
              KCOUNT=0
              DO IJ=NSTART(IRANK),NEND(IRANK)
                KCOUNT=KCOUNT+1
                KOBSLON(IJ,M,IC)=ICOMBUF_R(KCOUNT)
              ENDDO
            ENDIF

          ENDDO


!         READING KOBSRLAT
!         ----------------
          IF (IPROPAGS == 1) THEN

            DO IC=1,2
              ITAG=ITAG+1

              IF (IRANK == IREAD) THEN
                CALL GSTATS(1771,0)
                READ (IU08(IPROPAGS)) IDUM 
                CALL GSTATS(1771,1)

!               RELABELLING OF THE ARRAY
                IF (.NOT.LL1D .AND. NPR > 1) THEN
                  DO NIJ=NSTART(1),NEND(NPR)
                    KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                  ENDDO
                  DO NIJ=NSTART(1),NEND(NPR)
                    IDUM(NIJ)=KDUM(NIJ)
                  ENDDO
                ENDIF

!               FILL THE SEND BUFFER
                DO IP=1,NPR
                  KCOUNT=(IP-1)*MPLENGTH
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    ICOMBUF_S(KCOUNT)=IDUM(IJ)
                  ENDDO
                ENDDO
              ENDIF

              IF (NPR > 1) THEN
                CALL GSTATS(694,0)
                ICOUNTS(:)=MPLENGTH
                CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD, KSENDBUF=ICOMBUF_S,     &
     &                            KSENDCOUNTS=ICOUNTS,CDSTRING='GETBOBSTRCT 3:')
                CALL GSTATS(694,1)
              ENDIF

!             KEEP THE RELEVANT PART OF KOBSRLAT
              IF (IRANK == IREAD) THEN
                KCOUNT=(IRANK-1)*MPLENGTH
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSRLAT(IJ,M,IC)=ICOMBUF_S(KCOUNT)
                ENDDO
              ELSE
                KCOUNT=0
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSRLAT(IJ,M,IC)=ICOMBUF_R(KCOUNT)
                ENDDO
              ENDIF

            ENDDO
          ENDIF

!         READING KOBSRLON
!         ----------------
          IF (IPROPAGS == 1) THEN

            DO IC=1,2
              ITAG=ITAG+1

              IF (IRANK == IREAD) THEN
                CALL GSTATS(1771,0)
                READ (IU08(IPROPAGS)) IDUM 
                CALL GSTATS(1771,1)

!               RELABELLING OF THE ARRAY
                IF (.NOT.LL1D .AND. NPR > 1) THEN
                  DO NIJ=NSTART(1),NEND(NPR)
                     KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                  ENDDO
                  DO NIJ=NSTART(1),NEND(NPR)
                    IDUM(NIJ)=KDUM(NIJ)
                  ENDDO
                ENDIF

!               FILL THE SEND BUFFER
                DO IP=1,NPR
                  KCOUNT=(IP-1)*MPLENGTH
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    ICOMBUF_S(KCOUNT)=IDUM(IJ)
                  ENDDO
                ENDDO
              ENDIF

              IF (NPR > 1) THEN
                CALL GSTATS(694,0)
                ICOUNTS(:)=MPLENGTH
                CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD, KSENDBUF=ICOMBUF_S,     &
     &                            KSENDCOUNTS=ICOUNTS,CDSTRING='GETBOBSTRCT 4:')
                CALL GSTATS(694,1)
              ENDIF

!             KEEP THE RELEVANT PART OF KOBSRLON
              IF (IRANK == IREAD) THEN
                KCOUNT=(IRANK-1)*MPLENGTH
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSRLON(IJ,M,IC)=ICOMBUF_S(KCOUNT)
                ENDDO
              ELSE
                KCOUNT=0
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSRLON(IJ,M,IC)=ICOMBUF_R(KCOUNT)
                ENDDO
              ENDIF

            ENDDO
          ENDIF


!         READING KOBSCOR
!         ---------------
          IF (IPROPAGS == 2) THEN
            DO IC=1,4
              ITAG=ITAG+1

              IF (IRANK == IREAD) THEN
                CALL GSTATS(1771,0)
                READ (IU08(IPROPAGS)) IDUM 
                CALL GSTATS(1771,1)

!               RELABELLING OF THE ARRAY
                IF (.NOT.LL1D .AND. NPR > 1 ) THEN
                  DO NIJ=NSTART(1),NEND(NPR)
                   KDUM(NIJ)=IDUM(NEWIJ2IJ(NIJ))
                  ENDDO
                  DO NIJ=NSTART(1),NEND(NPR)
                    IDUM(NIJ)=KDUM(NIJ)
                  ENDDO
                ENDIF

!               FILL THE SEND BUFFER
                DO IP=1,NPR
                  KCOUNT=(IP-1)*MPLENGTH
                  DO IJ=NSTART(IP),NEND(IP)
                    KCOUNT=KCOUNT+1
                    ICOMBUF_S(KCOUNT)=IDUM(IJ)
                  ENDDO
                ENDDO
              ENDIF

              IF (NPR > 1) THEN
                CALL GSTATS(694,0)
                ICOUNTS(:)=MPLENGTH
                CALL MPL_SCATTERV(ICOMBUF_R,KROOT=IREAD, KSENDBUF=ICOMBUF_S,     &
     &                            KSENDCOUNTS=ICOUNTS,CDSTRING='GETBOBSTRCT 5:')
                CALL GSTATS(694,1)
              ENDIF           

!             KEEP THE RELEVANT PART OF KOBSCOR
              IF (IRANK == IREAD) THEN
                KCOUNT=(IRANK-1)*MPLENGTH
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSCOR(IJ,M,IC)=ICOMBUF_S(KCOUNT)
                ENDDO
              ELSE
                KCOUNT=0
                DO IJ=NSTART(IRANK),NEND(IRANK)
                  KCOUNT=KCOUNT+1
                  KOBSCOR(IJ,M,IC)=ICOMBUF_R(KCOUNT)
                ENDDO
              ENDIF
            ENDDO

          ENDIF

        ENDDO ! end loop on frequencies


        IF (ALLOCATED(ICOMBUF_S)) DEALLOCATE(ICOMBUF_S)
        IF (ALLOCATED(ICOMBUF_R)) DEALLOCATE(ICOMBUF_R)

        DEALLOCATE(IDUM)

        IF (.NOT.LL1D .AND. NPR > 1 ) DEALLOCATE(KDUM)

        WRITE(IU06,*) ' WAVE MODEL PREPROC UBUF INFORMATION READ IN  (second part)'
        CALL FLUSH (IU06)

      ENDIF  !! LSUBGRID (reading of subgrid information)


!     OBSTRUCTION COEFFICIENTS

!     NOTE: THE VALUE OF OBSLON WILL BE RESET IN THE FIRST
!     CALL TO PROPAGS TO CONTAIN THE OBSTRUCTION TIME THE GROUP
!     VELOCITY At THE INTERFACE !!!!!!!!
      IF (ALLOCATED(OBSLON)) DEALLOCATE(OBSLON)
      ALLOCATE(OBSLON(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))

!     NOTE: THE VALUE OF OBSLAT WILL BE RESET IN THE FIRST
!     CALL TO PROPAGS TO CONTAIN THE OBSTRUCTION TIME THE GROUP
!     VELOCITY At THE INTERFACE !!!!!!!!
      IF (ALLOCATED(OBSLAT)) DEALLOCATE(OBSLAT)
      ALLOCATE(OBSLAT(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))

      CALL GSTATS(1497,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICC,IC,M,IJ)
      DO ICC=1,2
        IC=ICC
        DO M=1,NFRE_RED
          DO IJ=NSTART(IRANK),NEND(IRANK)
            IF (.NOT. LSUBGRID) THEN
              OBSLON(IJ,M,IC)=1.0_JWRB
            ELSEIF (KOBSLON(IJ,M,IC) == 0) THEN
              OBSLON(IJ,M,IC)=0.0_JWRB
            ELSEIF (MOD(KOBSLON(IJ,M,IC),1000) == 0) THEN
              OBSLON(IJ,M,IC)=1.0_JWRB
            ELSE
              OBSLON(IJ,M,IC)=REAL(KOBSLON(IJ,M,IC),JWRB)*0.001_JWRB
            ENDIF
          ENDDO
        ENDDO

        DO M=1,NFRE_RED
          DO IJ=NSTART(IRANK),NEND(IRANK)
            IF (.NOT. LSUBGRID) THEN
              OBSLAT(IJ,M,IC)=1.0_JWRB
            ELSEIF (KOBSLAT(IJ,M,IC) == 0) THEN
              OBSLAT(IJ,M,IC)=0.0_JWRB
            ELSEIF (MOD(KOBSLAT(IJ,M,IC),1000) == 0) THEN
              OBSLAT(IJ,M,IC)=1.0_JWRB
            ELSE
              OBSLAT(IJ,M,IC)=REAL(KOBSLAT(IJ,M,IC),JWRB)*0.001_JWRB
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL GSTATS(1497,1)

      IF ( ALLOCATED(KOBSLON)) DEALLOCATE(KOBSLON)
      IF ( ALLOCATED(KOBSLAT)) DEALLOCATE(KOBSLAT)

      IF (IPROPAGS == 1) THEN
!       NOTE: THE VALUE OF OBSRLON WILL NOT BE RESET IN THE FIRST
        IF (ALLOCATED(OBSRLON)) DEALLOCATE(OBSRLON)
        ALLOCATE(OBSRLON(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
!       NOTE: THE VALUE OF OBSRLAT WILL NOT BE RESET IN THE FIRST
        IF (ALLOCATED(OBSRLAT)) DEALLOCATE(OBSRLAT)
        ALLOCATE(OBSRLAT(NSTART(IRANK):NEND(IRANK),NFRE_RED,2))
        CALL GSTATS(1497,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICC,IC,M,IJ)
        DO ICC=1,2
          IC=ICC
          DO M=1,NFRE_RED
            DO IJ=NSTART(IRANK),NEND(IRANK)
              IF (.NOT. LSUBGRID) THEN
                OBSRLON(IJ,M,IC)=1.0_JWRB
              ELSEIF (KOBSRLON(IJ,M,IC) == 0) THEN
                OBSRLON(IJ,M,IC)=0.0_JWRB
              ELSEIF (MOD(KOBSRLON(IJ,M,IC),1000) == 0) THEN
                OBSRLON(IJ,M,IC)=1.0_JWRB
              ELSE
                OBSRLON(IJ,M,IC)=REAL(KOBSRLON(IJ,M,IC),JWRB)*0.001_JWRB
              ENDIF
            ENDDO
          ENDDO

          DO M=1,NFRE_RED
            DO IJ=NSTART(IRANK),NEND(IRANK)
              IF (.NOT. LSUBGRID) THEN
                OBSRLAT(IJ,M,IC)=1.0_JWRB
              ELSEIF (KOBSRLAT(IJ,M,IC) == 0) THEN
                OBSRLAT(IJ,M,IC)=0.0_JWRB
               ELSEIF (MOD(KOBSRLAT(IJ,M,IC),1000) == 0) THEN
              OBSRLAT(IJ,M,IC)=1.0_JWRB
              ELSE
                OBSRLAT(IJ,M,IC)=REAL(KOBSRLAT(IJ,M,IC),JWRB)*0.001_JWRB
              ENDIF
            ENDDO
          ENDDO

        ENDDO
!$OMP END PARALLEL DO
        CALL GSTATS(1497,1)
        IF (ALLOCATED(KOBSRLON)) DEALLOCATE(KOBSRLON)
        IF (ALLOCATED(KOBSRLAT)) DEALLOCATE(KOBSRLAT)
      ENDIF

      IF (IPROPAGS == 2) THEN
        IF (ALLOCATED(OBSCOR)) DEALLOCATE(OBSCOR)
        ALLOCATE(OBSCOR(NSTART(IRANK):NEND(IRANK),NFRE_RED,4))
        CALL GSTATS(1497,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICC,IC,M,IJ)
        DO ICC=1,4
          IC=ICC
          DO M=1,NFRE_RED
            DO IJ=NSTART(IRANK),NEND(IRANK)
              IF (.NOT. LSUBGRID) THEN
                OBSCOR(IJ,M,IC)=1.0_JWRB
              ELSEIF (KOBSCOR(IJ,M,IC) == 0) THEN
                OBSCOR(IJ,M,IC)=0.0_JWRB
              ELSEIF (MOD(KOBSCOR(IJ,M,IC),1000) == 0) THEN
                OBSCOR(IJ,M,IC)=1.0_JWRB
              ELSE
                OBSCOR(IJ,M,IC)=REAL(KOBSCOR(IJ,M,IC),JWRB)*0.001_JWRB
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
        CALL GSTATS(1497,1)
        IF (ALLOCATED(KOBSCOR)) DEALLOCATE(KOBSCOR)
      ENDIF

      IF (IRANK == IREAD) CLOSE (UNIT=IU08(IPROPAGS))


WRITE(IU06,*) ''
WRITE(IU06,*) ' OBSTRUCTION COEFFICIENTS FROM BINARY INPUT COMPUTED'
CALL FLUSH(IU06)

IF (LHOOK) CALL DR_HOOK('GETBOBSTRCT',1,ZHOOK_HANDLE)

END SUBROUTINE GETBOBSTRCT
