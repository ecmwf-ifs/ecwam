      SUBROUTINE MINTF

! ----------------------------------------------------------------------

!**** *MINTF* - MAKE INTERPOLATION TABLES FOR BOUNDARY INPUT.

!     R. PORTZ     MPI         15/01/1991
!     H. GUNTHER   GKSS/ECMWF  15/01/1991

!*    PURPOSE.
!     -------

!       GENERATE SPACE INTERPOLATION TABLES USED FOR BOUNDARY
!       VALUE INPUT INTO A FINE GRID MODEL.

!**   INTERFACE.
!     ----------

!       *CALL* *MINT*

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCPBO  , ONLY : NBOUNC   ,DLAMAC   ,DPHIAC   ,BLATC    ,
     &            BLNGC
      USE YOWFPBO  , ONLY : IBFL     ,IBFR     ,BFW
      USE YOWMAP   , ONLY : NX       ,NY       ,AMOWEP   ,AMOSOP   ,
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO

! ----------------------------------------------------------------------

!*    1. RATIOS OF GRID INCREMENTS.
!        --------------------------

      IDELLA = NINT(DPHIAC/XDELLA)
      IDELLO = NINT(DLAMAC/XDELLO)

!*    2. SOUTHERN MOST LATITUDE OF FINE GRID.
!        ------------------------------------

      NI    = IDELLO - 1
      PHI   = AMOSOP

!*    2.1 LOOP OVER COARSE GRID POINTS.
!         -----------------------------

      DO I = 1, NX, IDELLO

!*    2.2 INTERPOLATION WEIGHT FOR INTERMEDIATE POINTS.
!         ---------------------------------------------

        IF (I.NE.NX) THEN
          DO N = 1, NI
            BFW (I+N) = REAL(N) / REAL(IDELLO)
          ENDDO
        ENDIF

!*    2.3 INDICES OF COARSE GRID OUTPUT POINTS.
!         -------------------------------------

        XLAMDA = AMOWEP + REAL(I-1) * XDELLO
        DO M=1,NBOUNC
          IF (BLATC(M).EQ.PHI .AND. BLNGC(M).EQ.XLAMDA) THEN
            IBFL(I) = M
            IBFR(I) = M
            DO N = 1, NI
              IF (I.NE.1) THEN
                IBFR(I-N) = M
              ENDIF
              IF (I.NE.NX) THEN
                IBFL(I+N) = M
              ENDIF
            ENDDO
            GOTO 2101
          ENDIF
        ENDDO
 2101   continue
      ENDDO

!*    3. NORTHERN MOST LATITUDE OF FINE GRID.
!        ------------------------------------

      PHI = AMONOP
      IS  = NX + 2*(NY-2) + 1
      IE  = 2*(NX+NY-2)

!*    3.1 LOOP OVER COARSE GRID POINTS.
!         -----------------------------

      DO I = IS, IE, IDELLO

!*    3.2 INTERPOLATION WEIGHT FOR INTERMEDIATE POINTS.
!         ---------------------------------------------

        IF (I.NE.IE) THEN
          DO N = 1, NI
            BFW (I+N) = REAL(N) / REAL(IDELLO)
          ENDDO
        ENDIF

!*    3.3 INDICES OF COARSE GRID OUTPUT POINTS.
!         -------------------------------------

        XLAMDA = AMOWEP + (I-IS) * XDELLO
        DO M=1,NBOUNC
          IF (BLATC(M).EQ.PHI .AND. BLNGC(M).EQ.XLAMDA) THEN
            IBFL(I) = M
            IBFR(I) = M
            DO N = 1, NI
              IF (I.NE.IS) THEN
                IBFR(I-N) = M
              ENDIF
              IF (I.NE.IE) THEN
                IBFL(I+N) = M
              ENDIF
            ENDDO
            GOTO 3101
          ENDIF
        ENDDO
 3101   continue
      ENDDO

!*    4. WESTERN MOST LONGITUDE OF FINE GRID.
!        ------------------------------------

      XLAMDA = AMOWEP
      NI = IDELLA - 1
      NSTEP = 2 * IDELLA
      IE = NX + 2*(NY-2) - 1
      IS = NX + 2*NI + 1
      K = 1

!*    4.1 WEIGHTS AND LEFT INDICES FOR FIRST COARSE GRID SECTION.
!        --------------------------------------------------------

      DO N = 1, NI
        IBFL(NX-1+2*N) = IBFL(1)
        BFW (NX-1+2*N) = REAL(N) / REAL(IDELLA)
      ENDDO

!*    4.2 LOOP OVER COARSE GRID POINTS.
!         -----------------------------

      DO I = IS, IE, NSTEP

!*    4.3 INTERPOLATION WEIGHT FOR INTERMEDIATE POINTS.
!         ---------------------------------------------

        DO N = 1, NI
          BFW (I+2*N) = REAL(N) / REAL(IDELLA)
        ENDDO

!*    4.4 INDICES OF COARSE GRID OUTPUT POINTS.
!         -------------------------------------

        K = K + 1
        PHI = AMOSOP + (K-1) * DPHIAC
        DO M=1,NBOUNC
          IF (BLNGC(M).EQ.XLAMDA .AND. BLATC(M).EQ.PHI) THEN
            IBFL(I) = M
            IBFR(I) = M
            DO N = 1, NI
              IBFL(I+2*N) = M
              IBFR(I-2*N) = M
            ENDDO
            GOTO 4201
          ENDIF
        ENDDO
 4201   continue
      ENDDO

!*    4.5 RIGHT INDICES FOR LAST COARSE GRID SECTION.
!         -------------------------------------------

      IF (IBFR(K) .NE. 0) THEN
        K = NX + (2* (NY-2)) + 1
        DO N = 1, NI
          IBFR(K-2*N) = IBFR(K)
        ENDDO
      ENDIF

!*    5. EASTERN MOST LONGITUDE OF FINE GRID.
!        ------------------------------------

      XLAMDA = AMOEAP
      IS = NX + 2*NI + 2
      IE = NX + 2*(NY-2)
      K = 1

!*    5.1 WEIGHTS AND LEFT INDICES FOR FIRST COARSE GRID SECTION.
!        --------------------------------------------------------

      DO N = 1, NI
        IBFL(NX+2*N) = IBFL(NX)
        BFW (NX+2*N) = REAL(N) / REAL(IDELLA)
      ENDDO

!*    5.2 LOOP OVER COARSE GRID POINTS.
!         -----------------------------

      DO I = IS, IE, NSTEP

!*    5.3 INTERPOLATION WEIGHT FOR INTERMEDIATE POINTS.
!         ---------------------------------------------

        DO N = 1, NI
          BFW (I+2*N) = REAL(N) / REAL(IDELLA)
        ENDDO

!*    5.4 INDICES OF COARSE GRID OUTPUT POINTS.
!         -------------------------------------

        K = K + 1
        PHI = AMOSOP + (K-1) * DPHIAC
        DO M=1,NBOUNC
          IF (BLNGC(M).EQ.XLAMDA .AND. BLATC(M).EQ.PHI) THEN
            IBFL(I) = M
            IBFR(I) = M
            DO N = 1, NI
              IBFL(I+2*N) = M
              IBFR(I-2*N) = M
            ENDDO
            GOTO 5201
          ENDIF
        ENDDO
 5201   continue
      ENDDO

!*    5.5 RIGHT INDICES FOR LAST COARSE GRID SECTION.
!         -------------------------------------------

      IF (IBFR(K) .NE. 0) THEN
        K = 2*(NX + NY - 2)
        M = NX + 2*NY - 2
        DO N = 1, NI
          IBFR(M-2*N) = IBFR(K)
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE MINTF
