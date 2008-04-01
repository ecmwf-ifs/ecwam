      SUBROUTINE BOUINPT (FL3, IG, IU02)

! ----------------------------------------------------------------------

!**** *BOUINPT* - BOUNDARY VALUE INPUT INTO A BLOCK.

!     H. GUNTHER    GKSS/ECMWF         JANUARY 1991

!*    PURPOSE.
!     --------

!       INPUT AND SPACE INTERPOLATION OF BOUNDARY SPECTRA.

!**   INTERFACE.
!     ----------

!       *CALL* *BOUINPT (FL3, IG, IU02)*
!          *FL3*     REAL         BLOCK OF SPECTRA.
!          *IG*      INTEGER      BLOCK NUMBER.
!          *IU02*    INTEGER      UNIT FOR INPUT OF BOUNDARY VALUES.

!     EXTERNALS.
!     ----------

!       *ABORT*     - TERMINATES PROCESSING.
!       *GSFILE*    - GETS OR SAVES A FILE.
!       *INTSPEC*   - INTERPOLATE A SPECTRUM.

!     METHOD.
!     -------

!       IN THE FIRST CALL OF THE SUB. THE FILE HEADER IS READ
!       AND THE CONSISTENCY IS CHECKED. THE SUB. READS A COMPLETE
!       SET OF BOUNDARY VALUES EACH PROPAGATION TIMESTEP WHEN
!       IT IS CALLED FOR THE FIRST BLOCK. THE SPECTRA REQUIRED FOR
!       THE ACTUAL BLOCK ARE INTERPOLATED AND STORED IN THE BLOCK.
!       INDICES AND WEIGHTS NECESSARY FOR THE INTERPOLATION AND
!       STORAGE ARE PRECOMPUTED IN PROG. PREPROC.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWFBOU  , ONLY : NBOUNF   ,IJARF    ,IGARF    ,IBFL     ,
     &            IBFR     ,BFW
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,CDTRES   ,IDELPRO
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      DIMENSION FL3(NINF-1:NSUP,NANG,NFRE)

! ----------------------------------------------------------------------
      REAL, ALLOCATABLE :: F1(:,:,:)
      REAL, ALLOCATABLE :: FL(:,:)
      REAL, ALLOCATABLE :: FMEAN1(:), EMEAN1(:), THQ1(:) 

      SAVE F1, FMEAN1, EMEAN1, THQ1, NBOINP

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *F1*        REAL      SPECTRA FROM COARSE GRID.
!      *FMEAN1*    REAL      MEAN FREQUENCIES FROM COARSE GRID.
!      *EMEAN1*    REAL      TOTAL ENERGIES FROM COARSE GRID.
!      *THQ1*      REAL      MEAN DIRECTIONS FROM COARSE GRID (RAD).
!      *NBOINP*    REAL      NUMBER OF INPUT SPECTRA.
!      *FL*        REAL      INTERPOLATED SPECTRUM.

      CHARACTER*14 CDTLAST, CDATE1, ZERO
      DATA IJAF, CDTLAST /0, '00000000000000'/  

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IJFA*      REAL      NUMBER OF FIRST INPUT POINT IN BLOCK.
!      *CDTLAST*   CHAR*14   DATE OF LAST FETCH FILE.


! ----------------------------------------------------------------------

!*    1. INITIALIZATION.
!        ---------------


!*    1.1 FIRST CALL OF SUB?
!         ------------------

      ZERO = ' '
      IF (CDTLAST.EQ.ZERO) THEN

!*    1.1.1 INITIALISE FOR LANDPOINTS.
!           --------------------------

        DO M=1,NFRE
          DO K=1,NANG
            F1(K,M,0) = 0.
          ENDDO
        ENDDO
        FMEAN1(0) = 0.
        EMEAN1(0) = 0.
        THQ1(0) = 0.
      ENDIF

!*    1.2 HAS A NEW BOUNDARY VALUE INPUT FILE TO BE FETCH?
!         ------------------------------------------------

      IF (CDTLAST.LT.CDTRES) THEN

!*    1.2.1 FETCH INPUT FILE.
!           -----------------

        CDTLAST = CDTRES
        CALL GSFILE (IU06, IU02, 0, CDTRES, CDATEF, 'FBI', 'G')

!*    1.2.2 READ BOUNDARY FILE HEADER.
!           --------------------------

        READ (IU02, ERR=5000, END=5000)
     &   XANG, XFRE, TH0, FR1, CO, XNBO, XDELIN
        KL = NINT(XANG)
        ML = NINT(XFRE)
        NBOINP  = NINT(XNBO)
        IDELINP = NINT(XDELIN)

        IF(.NOT.ALLOCATED(F1)) ALLOCATE(F1(NANG,NFRE,0:NBOINP))
        IF(.NOT.ALLOCATED(FL)) ALLOCATE(FL(NANG,NFRE))
        IF(.NOT.ALLOCATED(FMEAN1)) ALLOCATE(FMEAN1(0:NBOINP))
        IF(.NOT.ALLOCATED(EMEAN1)) ALLOCATE(EMEAN1(0:NBOINP))
        IF(.NOT.ALLOCATED(THQ1)) ALLOCATE(THQ1(0:NBOINP))

        IF (ITEST.GT.3) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' BOUNDARY VALUE INPUT FILE HEADER IS:'
          WRITE (IU06,*) ' NO. OF DIRECTIONS IS    KL     = ', KL
          WRITE (IU06,*) ' NO. OF FREQUENCIES IS   ML     = ', ML
          WRITE (IU06,*) ' FIRST DIRECTION IS      TH0    = ', TH0
          WRITE (IU06,*) ' FIRST FREQUENCY IS      FR1    = ', FR1
          WRITE (IU06,*) ' FREQUENCY RATIO IS      CO     = ', CO
          WRITE (IU06,*) ' NO. OF BOUNDRAY POINTS  NBOINP = ', NBOINP
          WRITE (IU06,*) ' TIME STEP OF DATA IS    IDELINP= ', IDELINP
        ENDIF

!*    1.2.3 CHECK CONSISTENCY.
!           ------------------

        IF (KL.NE.NANG .OR. ML.NE.NFRE .OR.
     &   FR1.NE.FR(1) .OR. TH0.NE.TH(1) .OR. IDELINP.NE.IDELPRO) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*    FATAL ERROR SUB. BOUINP           *'
          WRITE (IU06,*) '*    =======================           *'
          WRITE (IU06,*) '* VALUES IN BOUNDARY FILE HEADER ARE   *'
          WRITE (IU06,*) '* INCONSISTENT WITH MODEL SET-UP.      *'
          WRITE (IU06,*) '* MODEL VALUES ARE:                    *'
          WRITE (IU06,*) '* NO. OF DIRECTIONS     NANG   = ', NANG
          WRITE (IU06,*) '* NO. OF FREQUENCIES    NFRE   = ', NFRE
          WRITE (IU06,*) '* FIRST DIRECTION       TH0    = ', TH0
          WRITE (IU06,*) '* FIRST FREQUENCY       FR(1)  = ', FR(1)
          WRITE (IU06,*) '* PROPAGATION TIMESTEP  IDELPRO= ', IDELPRO
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.    *'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          CALL ABORT1
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    2. CALL OF SUB FOR FIRST BLOCK?
!        ----------------------------

      IF (IG.EQ.1) THEN
        IJAF = 1

!*    2.1 READ BOUNDARY SPECTRA.
!         ----------------------

 2100   CONTINUE
        DO IJ=1,NBOINP
          READ (IU02, ERR=5001, END=5001)
     &     DUM1, DUM2, CDATE1, EMEAN1(IJ),
     &     THQ1(IJ), FMEAN1(IJ), DUM3, DUM4
          READ (IU02, ERR=5002, END=5002)
     &     ((F1(K,M,IJ),K=1,NANG),M=1,NFRE)
        ENDDO
        IF (ITEST.GT.3) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' BOUNDARY VALUES READ CDATE1 = ', CDATE1
        ENDIF

!*    2.2 CHECK DATES.
!         ------------

        IF (CDATE1.LT.CDTPRO) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*    WARNING ERROR SUB. BOUINP         *'
          WRITE (IU06,*) '*    =========================         *'
          WRITE (IU06,*) '* BOUNDARY VALUE INPUT DATE IS BEFORE  *'
          WRITE (IU06,*) '* MODEL DATE.                          *'
          WRITE (IU06,*) '* MODEL DATE IS          CDTPRO = ', CDTPRO
          WRITE (IU06,*) '* BOUNDARY VALUE DATE IS CDATE1 = ', CDATE1
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '* PROGRAM WILL TRY NEXT INPUT.         *'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          GOTO 2100
        ELSEIF (CDATE1.GT.CDTPRO) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*    FATAL ERROR SUB. BOUINP           *'
          WRITE (IU06,*) '*    =======================           *'
          WRITE (IU06,*) '* DATES DO NOT MATCH.                  *'
          WRITE (IU06,*) '* MODEL DATE IS          CDTPRO = ', CDTPRO
          WRITE (IU06,*) '* BOUNDARY VALUE DATE IS CDATE1 = ', CDATE1
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.    *'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          CALL ABORT1
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    3. LOOP OVER BOUNDARY POINTS.
!        --------------------------

      DEL12 = 1.
      DO IJ=IJAF,NBOUNF

!*    3.1 GET INDICES.
!         ------------

        IGF = IGARF(IJ)
        IF (IGF.EQ.IG) THEN
          IJF = IJARF(IJ)
          IBCL = IBFL(IJ)
          IBCR = IBFR(IJ)
          DEL1L = BFW(IJ)

!*    3.2 CHECK INTERPOLATION WEIGHT.
!         ---------------------------

          IF (DEL1L.GT.0.) THEN

!*    3.2.1. WEIGHT IS GT ZERO INTERPOLATE.
!            ------------------------------

            CALL INTSPEC (NFRE, NANG, NFRE, NANG, FR, DEL12, DEL1L,
     &       F1(1,1,IBCL), FMEAN1(IBCL), EMEAN1(IBCL), THQ1(IBCL),
     &       F1(1,1,IBCR), FMEAN1(IBCR), EMEAN1(IBCR), THQ1(IBCR),
     &       FL, FMEAN, EMEAN, THQ)
            DO M=1,NFRE
              DO K=1,NANG
                FL3(IJF,K,M) = FL(K,M)
              ENDDO
            ENDDO
          ELSE

!*    3.2.2. WEIGHT IS ZERO COPY LEFT POINT.
!            -------------------------------

            DO M=1,NFRE
              DO K=1,NANG
                FL3(IJF,K,M) = F1(K,M,IBCL)
              ENDDO
            ENDDO
          ENDIF
        ELSEIF (IGF.GT.IG) THEN
          IJAF = IJ
          GOTO 4000
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    4. RETURN TO CALLING PROG.
!        -----------------------

 4000 CONTINUE

      RETURN

! ----------------------------------------------------------------------

!*    5. ERROR MESSAGES.
!        ---------------

 5000 CONTINUE
      WRITE (IU06,*) '*******************************************'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*      FATAL ERROR SUB. BOUNINP.          *'
      WRITE (IU06,*) '*      =========================          *'
      WRITE (IU06,*) '* END OF FILE OR READ ERROR.              *'
      WRITE (IU06,*) '* PROGRAM TRIES TO READ                   *'
      WRITE (IU06,*) '* HEADER OF BOUNDARY VALUES               *'
      WRITE (IU06,*) '* UNIT IS IU02 = ', IU02
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*******************************************'
      CALL ABORT1
 5001 CONTINUE
      WRITE (IU06,*) '*******************************************'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*      FATAL ERROR SUB. BOUNINP.          *'
      WRITE (IU06,*) '*      =========================          *'
      WRITE (IU06,*) '* END OF FILE OR READ ERROR.              *'
      WRITE (IU06,*) '* PROGRAM TRIES TO READ                   *'
      WRITE (IU06,*) '* A SPECTRUM HEADER                       *'
      WRITE (IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
      WRITE (IU06,*) '* UNIT IS          IU02 = ', IU02
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*******************************************'
      CALL ABORT1
 5002 CONTINUE
      WRITE (IU06,*) '*******************************************'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*      FATAL ERROR SUB. BOUNINP.          *'
      WRITE (IU06,*) '*      =========================          *'
      WRITE (IU06,*) '* END OF FILE OR READ ERROR.              *'
      WRITE (IU06,*) '* PROGRAM TRIES TO READ A SPECTRUM        *'
      WRITE (IU06,*) '* DATE IS        CDATE1 = ', CDATE1
      WRITE (IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
      WRITE (IU06,*) '* UNIT IS          IU02 = ', IU02
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
      WRITE (IU06,*) '*                                         *'
      WRITE (IU06,*) '*******************************************'
      CALL ABORT1

      END SUBROUTINE BOUINPT
