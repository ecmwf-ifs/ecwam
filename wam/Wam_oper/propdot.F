      SUBROUTINE PROPDOT(SL)

! ----------------------------------------------------------------------

!**** *PROPDOT* - PROPAGATION DOT TERMS FROM DEPTH AND CURRENT GRADIENT.

!     H. GUNTHER   GKSS/ECMWF   17/02/91
!     J. BIDLOT    ECMWF        FEBRUARY 97 MESSAGE PASSING

!*    PURPOSE.
!     --------

!       COMPUTATION OF COMMON REFDOT FOR PROPAGATION.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPDOT*

!     METHOD.
!     -------

!       IN A LOOP OVER THE BLOCKS THE COMMON UBUF IS READ,
!       THE DEPTH AND CURRENT GRADIENTS ARE COMPUTED,
!       COMMON REFDOT (DEPTH AND CURRENT REFRACTION FOR THETA DOT)
!       IS COMPUTED AND WRITTEN TO MASS STORAGE (IU16).
!       IN CASE OF CURRENT REFRACTION THE COMPLETE SIGMA DOT TERM
!       IS COMPUTED AND WRITTEN TO IU16 ADDITIONALLY.
!       WRITE OPERATIONS ARE NOT DONE FOR COMMON UBUF AND REFDOT
!       IF THIS IS A ONE BLOCK MODEL.

!     EXTERNALS.
!     ----------

!       *GRADI*     - COMPUTES DEPTH AND CURRENT GRADIENTS.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : COSTH    ,SINTH
      USE YOWGRID  , ONLY : COSPH    ,IGL      ,IJS      ,IJL      ,IJLT
      USE YOWMAP   , ONLY : KXLT
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLD    ,NIBLC
      USE YOWREFD  , ONLY : THDD     ,THDC
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,
     &            TCGOND   ,TFAK     ,TSIHKD   ,INDEP
      USE YOWSTAT  , ONLY : ICASE    ,ISHALLO  ,IREFRA
      USE YOWUBUF  , ONLY : KLAT     ,KLON
      USE YOWUNIT  , ONLY : IU08     ,IU16

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: SL

! ----------------------------------------------------------------------

      REAL, ALLOCATABLE, DIMENSION(:) :: DDPHI, DDLAM, DUPHI, DULAM
      REAL, ALLOCATABLE, DIMENSION(:) :: DVPHI, DVLAM, DCO, OMDD

      CHARACTER*30 FILENAME

! ----------------------------------------------------------------------

      ALLOCATE(DDPHI(NIBLD), DDLAM(NIBLD), DUPHI(NIBLC), DULAM(NIBLC))
      ALLOCATE(DVPHI(NIBLC), DVLAM(NIBLC), DCO(NIBLD), OMDD(NIBLC))

      IF(.NOT.ALLOCATED(THDC)) ALLOCATE(THDC(NIBLC,NANG))
      IF(.NOT.ALLOCATED(THDD)) ALLOCATE(THDD(NIBLD,NANG))

!*    1. IF CARTESIAN PROPAGATION SET COSINE OF LAT TO 1.
!         -----------------------------------------------

      IF (ICASE.NE.1) THEN
        DO IJ = 1,NIBLD
          DCO(IJ) = 1.
        ENDDO
      ENDIF

!*    2. LOOP OVER BLOCKS.
!        -----------------

      IF (IGL.NE.1) REWIND IU08
      DO IG = 1,IGL

!*    2.1 IF MULTI BLOCK VERSION.
!         -----------------------

        IF (IGL.NE.1) THEN

!*    2.1.1 READ COMMON BLOCK UBUF.
!           -----------------------

          READ (IU08) KLAT, KLON
!SHALLOW

!*    2.1.2 COMPUTE SHALLOW WATER TABLE INDICES.
!           ------------------------------------

          IF (ISHALLO.NE.1) THEN
            DO IJ=NINF,IJLT(IG)
              XD = LOG(DEPTH(IJ,IG)/DEPTHA)/LOG(DEPTHD)+1.
              ID = NINT(XD)
              ID = MAX(ID,1)
              INDEP(IJ) = MIN(ID,NDEPTH)
            ENDDO
          ENDIF
!SHALLOW
        ENDIF

!*    2.2 DEPTH AND CURRENT GRADIENTS.
!         ----------------------------

        CALL GRADI (IG, IREFRA, DDPHI, DDLAM, DUPHI,
     &   DULAM, DVPHI, DVLAM)

!*    2.3 COSINE OF LATITUDES IF SPHERICAL PROPAGATION.
!         ---------------------------------------------

        IF (ICASE.EQ.1) THEN
          DO IJ = IJS(IG),IJL(IG)
            JH = KXLT(IJ,IG)
            DCO(IJ) = 1./COSPH(JH)
          ENDDO
        ENDIF

!*    2.4 DEPTH GRADIENT PART OF SIGMA DOT.
!         ---------------------------------

        IF (ISHALLO.NE.1 .AND. IREFRA.EQ.2) THEN
          DO IJ = IJS(IG),IJL(IG)
            OMDD(IJ) = V(IJ,IG)*DDPHI(IJ) +
     &       U(IJ,IG)*DDLAM(IJ)*DCO(IJ)
          ENDDO
        ENDIF

!*    2.5. LOOP OVER DIRECTIONS.
!          ---------------------

        DO K=1,NANG
          SD = SINTH(K)
          CD = COSTH(K)

!*    2.5.1. DEPTH GRADIENT OF THETA DOT.
!            ----------------------------

          IF (ISHALLO.NE.1) THEN
            DO IJ = IJS(IG),IJL(IG)
              THDD(IJ,K) = SD*DDPHI(IJ) - CD*DDLAM(IJ)*DCO(IJ)
            ENDDO
          ENDIF

!*    2.5.2 SIGMA DOT AND THETA DOT PART FROM CURRENT GRADIENT.
!           ---------------------------------------------------

          IF (IREFRA.EQ.2) THEN
            SS  = SD**2
            SC  = SD*CD
            CC  = CD**2
            DO IJ = IJS(IG),IJL(IG)
              SL(IJ,K,NFRE) = -SC*DUPHI(IJ) - CC*DVPHI(IJ)
     &         - (SS*DULAM(IJ) + SC*DVLAM(IJ))*DCO(IJ)
              THDC(IJ,K) =  SS*DUPHI(IJ) + SC*DVPHI(IJ)
     &         - (SC*DULAM(IJ) + CC*DVLAM(IJ))*DCO(IJ)
            ENDDO

!*    2.5.3 LOOP OVER FREQUENCIES.
!           ----------------------

            IF (ISHALLO.NE.1) THEN
              DO M=1,NFRE
                DO IJ=IJS(IG),IJL(IG)
                  SL(IJ,K,M) = (SL(IJ,K,NFRE)*TCGOND(INDEP(IJ),M)
     &             + OMDD(IJ)*TSIHKD(INDEP(IJ),M))
     &             * TFAK(INDEP(IJ),M)
                ENDDO

!*    BRANCH BACK TO 2.5.3 FOR NEXT FREQUENCY.

              ENDDO
            ENDIF
          ENDIF

!*    BRANCH BACK TO 2.5 FOR NEXT DIRECTION.

        ENDDO

!*    2.6. WRITE COMPUTED DOT TERMS.
!           ------------------------

        IF (IGL.NE.1) THEN
          IF (ISHALLO.NE.1) THEN
            IF (IREFRA.EQ.1) THEN
              WRITE(IU16) THDD
            ELSE
              WRITE(IU16) THDD, THDC, SL
            ENDIF
          ELSE
            WRITE(IU16) THDC, ((SL(IJ,K,NFRE),IJ=1,NIBLC),K=1,NANG)
          ENDIF
        ELSE
          IF(LMESSPASS) THEN
            FILENAME='fort.16.%p'
            CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAME,1)
            OPEN(IU16,FILE=FILENAME,STATUS='UNKNOWN')

            WRITE(IU16) SL
            CLOSE(IU16)
          ELSE
            IF (ISHALLO.NE.1) THEN
              IF (IREFRA.EQ.2) THEN
                WRITE(IU16) SL
              ENDIF
            ELSE
              WRITE(IU16) ((SL(IJ,K,NFRE),IJ=1,NIBLC),K=1,NANG)
            ENDIF
          ENDIF
        ENDIF

!*    BRANCH BACK TO 2. FOR NEXT BLOCK.

      ENDDO

      DEALLOCATE(DDPHI, DDLAM, DUPHI, DULAM)
      DEALLOCATE(DVPHI, DVLAM, DCO, OMDD)
      IF (IGL.NE.1) THEN
        DEALLOCATE(THDC)
        DEALLOCATE(THDD)
      ENDIF

      RETURN
      END SUBROUTINE PROPDOT
