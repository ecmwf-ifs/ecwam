      SUBROUTINE PROPAGS (F1, F3, MIJS, MIJL, IG, SL, FCONST)

! ----------------------------------------------------------------------

!**** *PROPAGS* - COMPUTATION OF A PROPAGATION TIME STEP.

!     S.D. HASSELMANN.
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER

!     MODIFIED BY   H. GUNTHER   01/06/90    -   LAND POINTS ARE TAKEN
!                             OUT OF BLOCKS AND REFRACTION INTEGRATION
!                             CORRECTED FOR N-S AND S-N PROPAGATION.

!     K.P. HUBBERT                /07/89    -   DEPTH AND CURRENT
!     S. HASSELMANN   MPIFM       /04/90        REFRACTION SHALLOW

!     H. GUNTHER   GKSS/ECMWF   17/01/91    -   MODIFIED FOR CYCLE_4

!     J. BIDLOT    ECMWF   FEBRUARY 1997    - ADD SL AND FCONST
!     J. BIDLOT    ECMWF   APRIL 1997       - MODIFIED ADVECTION SCHEME 

!*    PURPOSE.
!     --------

!       COMPUTATION OF A PROPAGATION TIME STEP.

!**   INTERFACE.
!     ----------

!       *CALL* *PROPAGS(F1, F3, MIJS, MIJL, IG, SL, FCONST)*
!          *F1*   - SPECTRUM AT TIME T.
!          *F3*   - SPECTRUM AT TIME T+DELT.
!          *MIJS* - INDEX OF FIRST POINT
!          *MIJL* - INDEX OF LAST POINT
!          *IG*   - BLOCK NUMBER.
!          *SL*   - SIGMA DOT 
!          *FCONST* - SIGMA/SINH(2KD)

!     METHOD.
!     -------

!       FIRST ORDER FLUX SCHEME.

!     EXTERNALS.
!     ----------

!       *DOTDC*     - READ DOT TERMS FOR REFRACTION AND SCATTER TABLE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : FR       ,GOM      ,DELTH    ,FRATIO   ,
     &            COSTH    ,SINTH
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,
     &            IGL      ,IJLT
      USE YOWMAP   , ONLY : KXLT
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO    ,NIBLD    ,
     &            NIBLC
      USE YOWPCONS , ONLY : PI       ,ZPI      ,R
      USE YOWREFD  , ONLY : THDD     ,THDC
      USE YOWSHAL  , ONLY : NDEPTH   ,TCGOND   ,INDEP
      USE YOWSTAT  , ONLY : IDELPRO  ,ICASE    ,ISHALLO  ,IREFRA
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2    ,WADV     ,
     &            LSAMEDEPTH

! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: F1,F3,SL
      REAL,DIMENSION(NINF:NSUP,NFRE):: FCONST

      REAL,ALLOCATABLE,DIMENSION(:) :: DELLA0,DCO,DP1,DP2,DPN,DPN2,DPS
      REAL,ALLOCATABLE,DIMENSION(:) :: DEA,DLE,DLW,DPH,DLA
      REAL,ALLOCATABLE,DIMENSION(:) :: DTP,DTM,DTC,DRGP,DRGM
      REAL,ALLOCATABLE,DIMENSION(:) :: DOP,DOM,DRCP,DRCM
      REAL,ALLOCATABLE,DIMENSION(:) :: DRDP,DRDM
      REAL,ALLOCATABLE,DIMENSION(:,:) :: CGOND
      REAL,ALLOCATABLE,DIMENSION(:,:) :: WADVM1
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: CGKLON, CGKLAT, CGKLAT2

! ----------------------------------------------------------------------

      ALLOCATE(DELLA0(NINF-1:IJLT(IG)))
      ALLOCATE(DPH(NINF-1:IJLT(IG)))
      ALLOCATE(DLA(NINF-1:IJLT(IG)))
      ALLOCATE(DTP(MIJS:MIJL))
      ALLOCATE(DTM(MIJS:MIJL))
      ALLOCATE(DTC(MIJS:MIJL))


!*    0. SPECTRUM AT LAND TO ZERO.
!        -------------------------

      DO M=1,NFRE
        DO K=1,NANG
          F3(NINF-1,K,M) = 0.
        ENDDO
      ENDDO

      DELPRO = REAL(IDELPRO)   

!*    0.0 GROUP VELOCITIES. (if shallow water on)
!         -----------------
      IF (ISHALLO.NE.1) THEN
        ALLOCATE(CGOND(NINF-1:IJLT(IG),NFRE))
        DO M=1,NFRE
          CGOND(NINF-1,M) = TCGOND(NDEPTH,M)
          DO IJ=NINF,IJLT(IG)
            CGOND(IJ,M) = TCGOND(INDEP(IJ),M)
          ENDDO
        ENDDO

        ALLOCATE(CGKLON(MIJS:MIJL,NFRE,2))
        ALLOCATE(CGKLAT(MIJS:MIJL,NFRE,2))
        ALLOCATE(CGKLAT2(MIJS:MIJL,NFRE,2))

        DO IC=1,2
          DO M=1,NFRE
            DO IJ=MIJS,MIJL
              CGKLON(IJ,M,IC) = CGOND(KLON(IJ,IC),M)
              CGKLAT(IJ,M,IC) = CGOND(KLAT(IJ,IC),M)
              CGKLAT2(IJ,M,IC) = CGOND(KLAT2(IJ,IC),M)
            ENDDO
          ENDDO
        ENDDO

      ENDIF


!*    0.1 READ REFRACTION DOT TERMS.
!         --------------------------

      IF (IREFRA.NE.0) THEN
        IF(.NOT.ALLOCATED(THDC)) ALLOCATE(THDC(NIBLC,NANG))
        IF(.NOT.ALLOCATED(THDD)) ALLOCATE(THDD(NIBLD,NANG))
        CALL DOTDC (MIJS, MIJL, IG, IGL, ISHALLO, IREFRA, SL, FCONST)
      ENDIF

!*    0.2 SPHERICAL OR CARTESIAN GRID?
!         ----------------------------

      IF (ICASE.EQ.1) THEN

!*    0.2.1 SPHERICAL GRID.
!           ---------------
        ALLOCATE(DCO(NINF-1:IJLT(IG)))
        ALLOCATE(DP1(MIJS:MIJL))
        ALLOCATE(DP2(MIJS:MIJL))

!*    0.2.1.1 COSINE OF LATITUDE.
!             -------------------

        DELLA0(NINF-1) = 0.
        DO IJ = NINF,IJLT(IG)
          JH = KXLT(IJ,IG)
          DCO(IJ) = 1./COSPH(JH)
          DELLA0(IJ) = DELPRO/DELLAM(JH)
        ENDDO

!*    0.2.1.2 COMPUTE COS PHI FACTOR FOR ADJOINING GRID POINT.
!             ------------------------------------------------

        DO IJ = MIJS,MIJL
          JH = KLAT(IJ,1)
          IF (JH.EQ.NINF-1) THEN
            DP1(IJ) = 1.
          ELSE
            DP1(IJ) = DCO(IJ)/DCO(JH)
          ENDIF
          JH = KLAT(IJ,2)
          IF (JH.EQ.NINF-1) THEN
            DP2(IJ) = 1.
          ELSE
            DP2(IJ) = DCO(IJ)/DCO(JH)
          ENDIF
        ENDDO

        ALLOCATE(WADVM1(MIJS:MIJL,2))
        DO IC=1,2
          DO IJ = MIJS,MIJL
            WADVM1(IJ,IC) = 1. - WADV(IJ,IC)
          ENDDO
        ENDDO
        IF (IREFRA.NE.2) THEN

!*       BRANCH TO 3. IF WITHOUT REFRACTION OR DEPTH.
!        --------------------------------------------

          GOTO 3000
        ELSE

!*       BRANCH TO 4. IF DEPTH AND CURRENT REFRACTION.
!        ---------------------------------------------

          GOTO 4000
        ENDIF
      ELSE

!*    0.2.2 CARTESIAN GRID.
!           ---------------

        DELLA0(NINF-1) = 0.
        DO IJ = NINF,IJLT(IG)
          JH = KXLT(IJ,IG)
          DELLA0(IJ) = DELPRO/DELLAM(JH)
        ENDDO

!*    0.2.2.1 BRANCH TO 2. IF DEPTH AND CURRENT REFRACTION.
!             ---------------------------------------------

        IF (IREFRA.EQ.2) GOTO 2000
      ENDIF

! ----------------------------------------------------------------------

!*    1. PROPAGATION FOR CARTESIAN GRID
!*       WITHOUT REFRACTION OR DEPTH REFRATION.
!        --------------------------------------

      DELPH0 = DELPRO/DELPHI
      DELTH0 = 0.25*DELPRO/DELTH

      IF (ISHALLO.EQ.1) ALLOCATE(DEA(MIJS:MIJL))

      IF (ISHALLO.NE.1.AND.IREFRA.EQ.1) THEN
        ALLOCATE(DRDP(NIBLD))
        ALLOCATE(DRDM(NIBLD))
      ENDIF

!*    1.1 LOOP OVER DIRECTIONS.
!         ---------------------

      DO K=1,NANG
        SD = SINTH(K)
        CD = COSTH(K)

!*    1.1.1 INDEX FOR ADJOINING POINTS.
!           ---------------------------

        IF (SD.LT.0) THEN
          IJLA = 2
        ELSE
          IJLA = 1
        ENDIF
        IF (CD.LT.0) THEN
          IJPH = 2
        ELSE
          IJPH = 1
        ENDIF

        IF (ISHALLO.EQ.1) THEN

!*    1.1.2 DEEP WATER.
!           -----------

          SD = ABS(SD)
          CD = ABS(CD)
          SDA = ABS(SD)
          CDA = ABS(CD)

!*    1.1.2.1 LOOP OVER FREQUENCIES.
!             ----------------------

          DO M=1,NFRE

!*    1.1.2.1.1 LOOP OVER GRIDPOINTS.
!               ---------------------

            DNO = CD*GOM(M)*DELPH0
            DO IJ = MIJS,MIJL
              DEA(IJ) = SD*GOM(M)*DELLA0(IJ)
            ENDDO
            DO IJ = MIJS,MIJL
              DTT = 1.-(SDA*DELLA0(IJ)+CDA*DELPH0)*GOM(M)
              F3(IJ,K,M) = DTT * F1(IJ,K,M )
     &         + DNO * F1(KLAT(IJ,IJPH),K  ,M)
     &         + DEA(IJ) * F1(KLON(IJ,IJLA),K  ,M)
            ENDDO

!*    BRANCH BACK TO 1.1.2.1 FOR NEXT FREQUENCY.

          ENDDO
        ELSE
!SHALLOW

!*    1.1.3 SHALLOW WATER.
!           --------------

          SD = 0.5*SD
          CD = 0.5*CD

!*    1.1.3.1 DEPTH REFRACTION.
!             -----------------

          IF(IREFRA.EQ.1) THEN
            KP1 = K+1
            IF (KP1.GT.NANG) KP1 = 1
            KM1 = K-1
            IF (KM1.LT.1) KM1 = NANG
            DO IJ = MIJS,MIJL
              DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
              DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
            ENDDO
          ENDIF

!*    1.1.3.2 LOOP OVER FREQUENCIES.
!             ----------------------

          DO M=1,NFRE

!*    1.1.3.2.2 WEIGHTS IN INTEGRATION SCHEME.
!               ------------------------------

            IF (SD.GE.0.) THEN
              DO IJ=MIJS,MIJL
                SDD = SD*DELLA0(IJ)
                DLA(IJ) = SDD*(CGOND(KLON(IJ,1),M) + CGOND(IJ,M))
                DTC(IJ) = SDD*(CGOND(KLON(IJ,2),M) + CGOND(IJ,M))
              ENDDO
            ELSE
              DO IJ=MIJS,MIJL
                SDD = SD*DELLA0(IJ)
                DLA(IJ) =-SDD*(CGOND(KLON(IJ,2),M) + CGOND(IJ,M))
                DTC(IJ) =-SDD*(CGOND(KLON(IJ,1),M) + CGOND(IJ,M))
              ENDDO
            ENDIF

            IF (CD.GE.0.) THEN
              DO IJ=MIJS,MIJL
                CDD = CD*DELPH0
                DPH(IJ) = CDD*(CGOND(KLAT(IJ,1),M) + CGOND(IJ,M))
                DTC(IJ) = DTC(IJ)
     &           + CDD*(CGOND(KLAT(IJ,2),M) + CGOND(IJ,M))
              ENDDO
            ELSE
              DO IJ=MIJS,MIJL
                CDD = CD*DELPH0
                DPH(IJ) =-CDD*(CGOND(KLAT(IJ,2),M) + CGOND(IJ,M))
                DTC(IJ) = DTC(IJ)
     &           -CDD*(CGOND(KLAT(IJ,1),M) + CGOND(IJ,M))
              ENDDO
            ENDIF
            IF (IREFRA.EQ.1) THEN
              DO IJ = MIJS,MIJL
                DTHP = FCONST(IJ,M)*DRDP(IJ)
                DTHM = FCONST(IJ,M)*DRDM(IJ)
                DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
                DTP(IJ) = -DTHP+ABS(DTHP)
                DTM(IJ) =  DTHM+ABS(DTHM)
              ENDDO
            ENDIF

!*    1.1.3.2.3 LOOP OVER GRIDPOINTS.
!               ---------------------

            DO IJ = MIJS,MIJL
              F3(IJ,K,M) = (1.-DTC(IJ))*F1(IJ,K,M )
     &         + DPH(IJ) * F1(KLAT(IJ,IJPH),K  ,M)
     &         + DLA(IJ) * F1(KLON(IJ,IJLA),K  ,M)
            ENDDO
            IF (IREFRA.EQ.1) THEN
              DO IJ = MIJS,MIJL
                F3(IJ,K,M) = F3(IJ,K,M )
     &           + DTP(IJ) * F1(IJ,KP1,M)
     &           + DTM(IJ) * F1(IJ,KM1,M)
              ENDDO
            ENDIF

!*    BRANCH BACK TO 1.1.3.2 FOR NEXT FREQUENCY.

          ENDDO
!SHALLOW
        ENDIF

!*    BRANCH BACK TO 1.1 FOR NEXT DIRECTION.

      ENDDO

!*    1.2 END OF PROPAGATION FOR CARTESIAN GRID
!*        WITHOUT REFRACTION OR DEPTH REFRACTION, RETURN.
!         -----------------------------------------------

      DEALLOCATE(DELLA0)
      DEALLOCATE(DPH,DLA)
      DEALLOCATE(DTP,DTM,DTC)
      IF (ISHALLO.NE.1) DEALLOCATE(CGOND)
      IF(ALLOCATED(DEA))DEALLOCATE(DEA)
      IF (ISHALLO.NE.1.AND.IREFRA.EQ.1) THEN
        DEALLOCATE(DRDP)
        DEALLOCATE(DRDM)
      ENDIF
      IF (IGL.NE.1) THEN
        IF(ALLOCATED(THDC))DEALLOCATE(THDC)
        IF(ALLOCATED(THDD))DEALLOCATE(THDD)
      ENDIF

      RETURN

! ----------------------------------------------------------------------

!*    2. PROPAGATION FOR CARTESIAN GRID
!*       WITH DEPTH AND CURRENT REFRACTION.
!        ----------------------------------

 2000 CONTINUE

      ALLOCATE(DPN(MIJS:MIJL))
      ALLOCATE(DPS(MIJS:MIJL))
      ALLOCATE(DLE(MIJS:MIJL))
      ALLOCATE(DLW(MIJS:MIJL))

      ALLOCATE(DOP(NIBLC))
      ALLOCATE(DOM(NIBLC))
      ALLOCATE(DRCP(NIBLC))
      ALLOCATE(DRCM(NIBLC))
      IF (ISHALLO.NE.1) THEN
        ALLOCATE(DRDP(NIBLD))
        ALLOCATE(DRDM(NIBLD))
      ENDIF

      DELPH0 = 0.25*DELPRO/DELPHI
      DELTH0 = 0.25*DELPRO/DELTH
      DELLA0 = 0.25*DELLA0
      DELFR0 = 0.25*DELPRO/((1.-FRATIO)*ZPI)

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

      DO K=1,NANG
        KP1 = K+1
        IF (KP1.GT.NANG) KP1 = 1
        KM1 = K-1
        IF (KM1.LT.1) KM1 = NANG
        SD = SINTH(K)
        CD = COSTH(K)

!*    2.1.1 DEPTH REFRACTION IF SHALLOW WATER.
!           ----------------------------------

        IF (ISHALLO.NE.1) THEN
          DO IJ = MIJS,MIJL
            DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
            DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
          ENDDO
        ENDIF

!*    2.1.2 CURRENT REFRACTION.
!           -------------------

        DO IJ = MIJS,MIJL
          DRCP(IJ) = (THDC(IJ,K) + THDC(IJ,KP1))*DELTH0
          DRCM(IJ) = (THDC(IJ,K) + THDC(IJ,KM1))*DELTH0
        ENDDO

!*    2.1.3 LOOP OVER FREQUENCIES.
!           ----------------------

        DO M=1,NFRE
          IF (ISHALLO.EQ.1) THEN

!*    2.1.3.1 DEEP WATER.
!             -----------

            MP1 = MIN(NFRE,M+1)
            MM1 = MAX(1,M-1)
            DFP = PI*(1.+FRATIO)*DELFR0

!*    2.1.3.1.1 GROUP VELOCITIES.
!               -----------------

            CGS = GOM(M)*SD
            CGC = GOM(M)*CD

!*    2.1.3.1.2 WEIGHTS IN INTEGRATION SCHEME.
!               ------------------------------

            DLA(NINF-1) = CGS*DELLA0(NINF-1)
            DPH(NINF-1) = CGC*DELPH0
            DO IJ=NINF,IJLT(IG)
              DLA(IJ) = (U(IJ,IG) + CGS)*DELLA0(IJ)
              DPH(IJ) = (V(IJ,IG) + CGC)*DELPH0
            ENDDO
            DO IJ=MIJS,MIJL
              DLWE = DLA(IJ) + DLA(KLON(IJ,1))
              DLEA = DLA(IJ) + DLA(KLON(IJ,2))
              DLE(IJ) = -DLEA+ABS(DLEA)
              DLW(IJ) =  DLWE+ABS(DLWE)
              DTC(IJ) =  DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

              DPSO = DPH(IJ) + DPH(KLAT(IJ,1))
              DPNO = DPH(IJ) + DPH(KLAT(IJ,2))
              DPN(IJ) = -DPNO+ABS(DPNO)
              DPS(IJ) =  DPSO+ABS(DPSO)
              DTC(IJ) =  DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

              DTHP = DRCP(IJ)
              DTHM = DRCM(IJ)
              DTP(IJ) = -DTHP+ABS(DTHP)
              DTM(IJ) =  DTHM+ABS(DTHM)
              DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)

              DTHP    = SL(IJ,K,NFRE) * DFP
              DTC(IJ) = DTC(IJ) + 2.* ABS(DTHP)
              DOP(IJ) = (-DTHP+ABS(DTHP))/FRATIO
              DOM(IJ) = ( DTHP+ABS(DTHP))*FRATIO
            ENDDO
          ELSE
!SHALLOW

!*    2.1.3.2 SHALLOW WATER.
!             --------------

            MP1 = MIN(NFRE,M+1)
            MM1 = MAX(1,M-1)
            DFP = DELFR0/FR(M)
            DFM = DELFR0/FR(MM1)


!*    2.1.3.2.2 WEIGHTS IN INTEGRATION SCHEME.
!               ------------------------------

            DLA(NINF-1) = SD*CGOND(NINF-1,M)*DELLA0(NINF-1)
            DPH(NINF-1) = CD*CGOND(NINF-1,M)*DELPH0
            DO IJ=NINF,IJLT(IG)
              DLA(IJ) = (U(IJ,IG) + SD*CGOND(IJ,M))*DELLA0(IJ)
              DPH(IJ) = (V(IJ,IG) + CD*CGOND(IJ,M))*DELPH0
            ENDDO
            DO IJ=MIJS,MIJL
              DLWE = DLA(IJ) + DLA(KLON(IJ,1))
              DLEA = DLA(IJ) + DLA(KLON(IJ,2))
              DLE(IJ) = -DLEA+ABS(DLEA)
              DLW(IJ) =  DLWE+ABS(DLWE)
              DTC(IJ) = DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

              DPSO = DPH(IJ) + DPH(KLAT(IJ,1))
              DPNO = DPH(IJ) + DPH(KLAT(IJ,2))
              DPN(IJ) = -DPNO+ABS(DPNO)
              DPS(IJ) =  DPSO+ABS(DPSO)
              DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

              DTHP = FCONST(IJ,M)*DRDP(IJ) + DRCP(IJ)
              DTHM = FCONST(IJ,M)*DRDM(IJ) + DRCM(IJ)
              DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
              DTP(IJ) = -DTHP+ABS(DTHP)
              DTM(IJ) =  DTHM+ABS(DTHM)

              DTHP = (SL(IJ,K,M) + SL(IJ,K,MP1))*DFP
              DTHM = (SL(IJ,K,M) + SL(IJ,K,MM1))*DFM
              DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
              DOP(IJ) = (-DTHP+ABS(DTHP))/FRATIO
              DOM(IJ) = ( DTHM+ABS(DTHM))*FRATIO
            ENDDO
!SHALLOW
          ENDIF

!*    2.1.3.3 LOOP OVER GRIDPOINTS.
!             ---------------------

          DO IJ = MIJS,MIJL
            F3(IJ,K,M) = (1.-DTC(IJ))*F1(IJ,K,M )
     &       + DPN(IJ) * F1(KLAT(IJ,2),K  ,M)
     &       + DPS(IJ) * F1(KLAT(IJ,1),K  ,M)
     &       + DLE(IJ) * F1(KLON(IJ,2),K  ,M)
     &       + DLW(IJ) * F1(KLON(IJ,1),K  ,M)
     &       + DTP(IJ) * F1(IJ        ,KP1,M)
     &       + DTM(IJ) * F1(IJ        ,KM1,M)
     &       + DOP(IJ) * F1(IJ        ,K  ,MP1)
     &       + DOM(IJ) * F1(IJ        ,K  ,MM1)
          ENDDO

!*    BRANCH BACK TO 2.1.3 FOR NEXT FREQUENCY.

        ENDDO

!*    BRANCH BACK TO 2.1 FOR NEXT DIRECTION.

      ENDDO

!*    2.2 END OF PROPAGATION FOR CARTESIAN GRID
!*        WITH DEPTH AND CURRENT REFRACTION, RETURN.
!         ------------------------------------------

      DEALLOCATE(DELLA0,DPN,DPS)
      DEALLOCATE(DLE,DLW,DPH,DLA)
      DEALLOCATE(DTP,DTM,DTC)
      IF (ISHALLO.NE.1) DEALLOCATE(CGOND)

      DEALLOCATE(DOP,DOM,DRCP,DRCM)
      IF (ISHALLO.NE.1) THEN
        DEALLOCATE(DRDP)
        DEALLOCATE(DRDM)
      ENDIF

      IF (IGL.NE.1) THEN
        IF(ALLOCATED(THDC))DEALLOCATE(THDC)
        IF(ALLOCATED(THDD))DEALLOCATE(THDD)
      ENDIF

      RETURN

! ----------------------------------------------------------------------

!*    3. PROPAGATION FOR SPHERICAL LATITUDE/LONGITUDE GRID
!*       WITHOUT OR DEPTH REFRACTION.
!        -------------------------------------------------

 3000 CONTINUE

      ALLOCATE(DPN(MIJS:MIJL))
      ALLOCATE(DPN2(MIJS:MIJL))
      ALLOCATE(DLE(MIJS:MIJL))
      ALLOCATE(DRGP(MIJS:MIJL))
      ALLOCATE(DRGM(MIJS:MIJL))

      IF (ISHALLO.NE.1.AND.IREFRA.EQ.1) THEN
        ALLOCATE(DRDP(NIBLD))
        ALLOCATE(DRDM(NIBLD))
      ENDIF

      DELTH0 = 0.25*DELPRO/DELTH
      DELPH0 = 0.5*DELPRO/DELPHI

!*    3.1 LOOP OVER DIRECTIONS.
!         ---------------------

      DO K=1,NANG
        KP1 = K+1
        IF (KP1.GT.NANG) KP1 = 1
        KM1 = K-1
        IF (KM1.LT.1) KM1 = NANG
        SD = SINTH(K)
        CD = COSTH(K)
        SDA = ABS(SD)
        SDA2 = 0.5*SDA
        CDA = ABS(CD)
        DELPH0_CDA = DELPH0*CDA

!*    3.1.1 COMPUTE GRID REFRACTION.
!           ------------------------

        SP  = DELTH0*(SINTH(K)+SINTH(KP1))/R
        SM  = DELTH0*(SINTH(K)+SINTH(KM1))/R
        DO IJ = MIJS,MIJL
          JH = KXLT(IJ,IG)
          TANPH = SINPH(JH)*DCO(IJ)
          DRGP(IJ) = TANPH*SP
          DRGM(IJ) = TANPH*SM
        ENDDO

!*    3.1.2 INDEX FOR ADJOINING POINTS.
!           ---------------------------

        IF (SD.LT.0) THEN
          IJLA = 2
        ELSE
          IJLA = 1
        ENDIF
        IF (CD.LT.0) THEN
          IJPH = 2
        ELSE
          IJPH = 1
        ENDIF

        IF (ISHALLO.EQ.1) THEN

!*    3.1.3 DEEP WATER.
!           -----------

!*    3.1.3.1 LAT / LONG WEIGHTS IN INTEGRATION SCHEME.
!             -----------------------------------------

          DO IJ=MIJS,MIJL
            DLE(IJ) = DCO(IJ)*SDA*DELLA0(IJ)
          ENDDO
          IF (CD.GT.0.) THEN
            DO IJ=MIJS,MIJL
              DTC(IJ) = DLE(IJ) + DELPH0_CDA*(DP2(IJ) + 1.)
              DPN(IJ) = DELPH0*CDA*(DP1(IJ) + 1.)
            ENDDO
          ELSE
            DO IJ=MIJS,MIJL
              DTC(IJ) = DLE(IJ) + DELPH0_CDA*(DP1(IJ) + 1.)
              DPN(IJ) = DELPH0_CDA*(DP2(IJ) + 1.)
            ENDDO
          ENDIF

!*    3.1.3.2 REFRACTION WEIGHTS IN INTEGRATION SCHEME.
!             -----------------------------------------

          DO IJ=MIJS,MIJL
            DTHP = DRGP(IJ)
            DTHM = DRGM(IJ)
            DTC(IJ) = DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
            DTP(IJ) = -DTHP+ABS(DTHP)
            DTM(IJ) =  DTHM+ABS(DTHM)
          ENDDO

!*    3.1.3.3 LOOP OVER FREQUENCIES.
!             ----------------------

          DO M=1,NFRE

!*    3.1.3.3.1 LOOP OVER GRIDPOINTS.
!               ---------------------

            DO IJ = MIJS,MIJL
              DTT = 1. - DTC(IJ)*GOM(M)
              F3(IJ,K,M) = DTT*F1(IJ,K,M ) + GOM(M) *
     &          (DPN(IJ) * WADV(IJ,IJPH) * F1(KLAT(IJ,IJPH),K  ,M)
     &         + DPN(IJ) * WADVM1(IJ,IJPH) * F1(KLAT2(IJ,IJPH),K,M)
     &         + DLE(IJ) * F1(KLON(IJ,IJLA),K  ,M)
     &         + DTP(IJ) * F1(IJ           ,KP1,M)
     &         + DTM(IJ) * F1(IJ           ,KM1,M))
            ENDDO

!*    BRANCH BACK TO 3.1.3.3 FOR NEXT FREQUENCY.

          ENDDO
        ELSE
!SHALLOW

!*    3.1.4 SHALLOW WATER.
!           --------------


!*    3.1.4.1 COMPUTE DEPTH REFRACTION.
!             -------------------------

          IF (IREFRA.EQ.1) THEN
            DO IJ = MIJS,MIJL
              DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
              DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
            ENDDO
          ENDIF

!*    3.1.4.2 LOOP OVER FREQUENCIES.
!             ----------------------

          DO M=1,NFRE


!*    3.1.4.3.2 LAT / LONG WEIGHTS IN INTEGRATION SCHEME.
!               -----------------------------------------

            IF (SD.GT.0.) THEN
              DO IJ=MIJS,MIJL
                IF(LSAMEDEPTH(IJ)) THEN
                  DTC(IJ) = 1. - DCO(IJ)*SDA*DELLA0(IJ)*CGOND(IJ,M)
                  DLE(IJ) = DCO(IJ)*SDA*DELLA0(IJ)*CGOND(IJ,M)
                ELSE
                  DTC(IJ) = 1. - DCO(IJ)*SDA2*DELLA0(IJ)*
     &             (CGKLON(IJ,M,2) + CGOND(IJ,M))
                  DLE(IJ) = DCO(IJ)*SDA2*DELLA0(IJ)*
     &             (CGKLON(IJ,M,1) + CGOND(IJ,M))
                ENDIF
              ENDDO
            ELSE
              DO IJ=MIJS,MIJL
                IF(LSAMEDEPTH(IJ)) THEN
                  DTC(IJ) = 1. - DCO(IJ)*SDA*DELLA0(IJ)*CGOND(IJ,M)
                  DLE(IJ) = DCO(IJ)*SDA*DELLA0(IJ)*CGOND(IJ,M)
                ELSE
                  DTC(IJ) = 1. - DCO(IJ)*SDA2*DELLA0(IJ)*
     &             (CGKLON(IJ,M,1) + CGOND(IJ,M))
                  DLE(IJ) = DCO(IJ)*SDA2*DELLA0(IJ)*
     &             (CGKLON(IJ,M,2) + CGOND(IJ,M))
                ENDIF
              ENDDO
            ENDIF
            IF (CD.GT.0.) THEN
              DO IJ=MIJS,MIJL
                IF(LSAMEDEPTH(IJ)) THEN
                  DTC(IJ) = DTC(IJ)-DELPH0_CDA*CGOND(IJ,M)*(DP2(IJ)+1.)
                  DPN(IJ) = DELPH0_CDA*CGOND(IJ,M)*(DP1(IJ)+1.)
                  DPN2(IJ)= WADVM1(IJ,1)*DPN(IJ)
                  DPN(IJ) = WADV(IJ,1)*DPN(IJ)
                ELSE
                  DTC(IJ) = DTC(IJ) - DELPH0_CDA*( DP2(IJ)*
     &                  (  WADV(IJ,2)*CGKLAT(IJ,M,2)
     &                   + WADVM1(IJ,2)*CGKLAT2(IJ,M,2) )
     &                   + CGOND(IJ,M) )
                  DPN(IJ) = DELPH0_CDA*WADV(IJ,1)*
     &                      (CGKLAT(IJ,M,1)*DP1(IJ) + CGOND(IJ,M))
                  DPN2(IJ) = DELPH0_CDA*WADVM1(IJ,1)*
     &                      (CGKLAT2(IJ,M,1)*DP1(IJ) + CGOND(IJ,M))
                ENDIF
              ENDDO
            ELSE
              DO IJ=MIJS,MIJL
                IF(LSAMEDEPTH(IJ)) THEN
                  DTC(IJ) = DTC(IJ)-DELPH0_CDA*CGOND(IJ,M)*(DP1(IJ)+1.)
                  DPN(IJ) = DELPH0_CDA*CGOND(IJ,M)*(DP2(IJ)+1.)
                  DPN2(IJ)= WADVM1(IJ,2)*DPN(IJ)
                  DPN(IJ) = WADV(IJ,2)*DPN(IJ)
                ELSE
                  DTC(IJ) = DTC(IJ) - DELPH0_CDA*( DP1(IJ)*
     &                  (  WADV(IJ,1)*CGKLAT(IJ,M,1)
     &                   + WADVM1(IJ,1)*CGKLAT2(IJ,M,1))
     &                   + CGOND(IJ,M))
                  DPN(IJ) = DELPH0_CDA*WADV(IJ,2)*
     &                     (CGKLAT(IJ,M,2)*DP2(IJ) + CGOND(IJ,M))
                  DPN2(IJ) = DELPH0_CDA*WADVM1(IJ,2)*
     &                     (CGKLAT2(IJ,M,2)*DP2(IJ) + CGOND(IJ,M))
                ENDIF
              ENDDO
            ENDIF

!*    3.1.4.2.3 REFRACTION WEIGHTS IN INTEGRATION SCHEME.
!               -----------------------------------------

            IF (IREFRA.EQ.0) THEN
              DO IJ=MIJS,MIJL
                DTHP = DRGP(IJ)*CGOND(IJ,M)
                DTHM = DRGM(IJ)*CGOND(IJ,M)
                DTC(IJ) = DTC(IJ) - DTHP-ABS(DTHP)+DTHM-ABS(DTHM)
                DTP(IJ) = -DTHP+ABS(DTHP)
                DTM(IJ) =  DTHM+ABS(DTHM)
              ENDDO
            ELSE
              DO IJ=MIJS,MIJL
                DTHP = DRGP(IJ)*CGOND(IJ,M)+FCONST(IJ,M)*DRDP(IJ)
                DTHM = DRGM(IJ)*CGOND(IJ,M)+FCONST(IJ,M)*DRDM(IJ)
                DTC(IJ) = DTC(IJ) - DTHP-ABS(DTHP)+DTHM-ABS(DTHM)
                DTP(IJ) = -DTHP+ABS(DTHP)
                DTM(IJ) =  DTHM+ABS(DTHM)
              ENDDO
            ENDIF

!*    3.1.4.2.4 LOOP OVER GRIDPOINTS.
!               ---------------------

            DO IJ = MIJS,MIJL
              F3(IJ,K,M) = DTC(IJ)*F1(IJ,K,M )
     &         + DPN(IJ) * F1(KLAT(IJ,IJPH),K  ,M)
     &         + DPN2(IJ)* F1(KLAT2(IJ,IJPH),K  ,M)
     &         + DLE(IJ) * F1(KLON(IJ,IJLA),K  ,M)
     &         + DTP(IJ) * F1(IJ           ,KP1,M)
     &         + DTM(IJ) * F1(IJ           ,KM1,M)
            ENDDO

!*    BRANCH BACK TO 3.1.4.2 FOR NEXT FREQUENCY.

          ENDDO
!SHALLOW
        ENDIF

!*    BRANCH BACK TO 3.1 FOR NEXT DIRECTION.

      ENDDO

!*    3.2 END OF PROPAGATION FOR SPHERICAL GRID
!*        WITHOUT REFRACTION OR DEPTH REFRACTION, RETURN.
!         -----------------------------------------------

      DEALLOCATE(DELLA0,DCO,DP1,DP2,DPN,DPN2)
      DEALLOCATE(DLE,DPH,DLA)
      DEALLOCATE(DTP,DTM,DRGP,DRGM,DTC)
      DEALLOCATE(WADVM1)
      IF (ISHALLO.NE.1) THEN
        DEALLOCATE(CGKLON)
        DEALLOCATE(CGKLAT)
        DEALLOCATE(CGKLAT2)
      ENDIF
      IF (ISHALLO.NE.1) DEALLOCATE(CGOND)
      IF (ISHALLO.NE.1.AND.IREFRA.EQ.1) THEN
        DEALLOCATE(DRDP)
        DEALLOCATE(DRDM)
      ENDIF
      IF (IGL.NE.1) THEN
        IF(ALLOCATED(THDC))DEALLOCATE(THDC)
        IF(ALLOCATED(THDD))DEALLOCATE(THDD)
      ENDIF

      RETURN

! ----------------------------------------------------------------------

!*    4. PROPAGATION FOR SPHERICAL LATITUDE/LONGITUDE GRID
!*       WITH DEPTH AND CURRENT REFRACTION.
!        -------------------------------------------------

 4000 CONTINUE

      ALLOCATE(DPN(MIJS:MIJL))
      ALLOCATE(DPS(MIJS:MIJL))
      ALLOCATE(DLE(MIJS:MIJL))
      ALLOCATE(DLW(MIJS:MIJL))

      ALLOCATE(DOP(NIBLC))
      ALLOCATE(DOM(NIBLC))
      ALLOCATE(DRCP(NIBLC))
      ALLOCATE(DRCM(NIBLC))
      IF (ISHALLO.NE.1) THEN
        ALLOCATE(DRDP(NIBLD))
        ALLOCATE(DRDM(NIBLD))
      ENDIF


      ALLOCATE(DRGP(MIJS:MIJL))
      ALLOCATE(DRGM(MIJS:MIJL))

      DELPH0 = 0.25*DELPRO/DELPHI
      DELTH0 = 0.25*DELPRO/DELTH
      DELLA0 = 0.25*DELLA0
      DELFR0 = 0.25*DELPRO/((1.-FRATIO)*ZPI)

!*    4.1 LOOP OVER DIRECTIONS.
!         ---------------------

      DO K=1,NANG
        KP1 = K+1
        IF (KP1.GT.NANG) KP1 = 1
        KM1 = K-1
        IF (KM1.LT.1) KM1 = NANG
        SD = SINTH(K)
        CD = COSTH(K)

!*    4.1.1 COMPUTE GRID REFRACTION.
!           ------------------------

        SP = DELTH0*(SINTH(K)+SINTH(KP1))/R
        SM = DELTH0*(SINTH(K)+SINTH(KM1))/R
        DO IJ = MIJS,MIJL
          JH = KXLT(IJ,IG)
          TANPH = SINPH(JH)*DCO(IJ)
          DRGP(IJ) = TANPH*SP
          DRGM(IJ) = TANPH*SM
        ENDDO

!*    4.1.2 COMPUTE DEPTH REFRACTION.
!           -------------------------

        IF (ISHALLO.NE.1) THEN
          DO IJ = MIJS,MIJL
            DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
            DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
          ENDDO
        ENDIF

!*    4.1.3 COMPUTE CURRENT REFRACTION.
!           ---------------------------

        DO IJ = MIJS,MIJL
          DRCP(IJ) = (THDC(IJ,K) + THDC(IJ,KP1))*DELTH0
          DRCM(IJ) = (THDC(IJ,K) + THDC(IJ,KM1))*DELTH0
        ENDDO

!*    4.1.4 LOOP OVER FREQUENCIES.
!           ----------------------

        DO M=1,NFRE
          MP1 = MIN(NFRE,M+1)
          MM1 = MAX(1,M-1)
          IF (ISHALLO.EQ.1) THEN

!*    4.1.4.1 DEEP WATER.
!             -----------

!*    4.1.4.1.1 GROUP VELOCITIES.
!               -----------------

            DFP = PI*(1.+FRATIO)*DELFR0
            CGS = GOM(M)*SD
            CGC = GOM(M)*CD

!*    4.1.4.1.2 WEIGHTS IN INTEGRATION SCHEME.
!               ------------------------------

            DLA(NINF-1) = CGS*DELLA0(NINF-1)
            DPH(NINF-1) = CGC*DELPH0
            DO IJ=NINF,IJLT(IG)
              DLA(IJ) = (U(IJ,IG) + CGS)*DELLA0(IJ)*DCO(IJ)
              DPH(IJ) = (V(IJ,IG) + CGC)*DELPH0
            ENDDO
            DO IJ=MIJS,MIJL
              DLWE = DLA(IJ) + DLA(KLON(IJ,1))
              DLEA = DLA(IJ) + DLA(KLON(IJ,2))
              DLE(IJ) = -DLEA+ABS(DLEA)
              DLW(IJ) =  DLWE+ABS(DLWE)
              DTC(IJ) =  DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

              DPSO = DPH(IJ) + DPH(KLAT(IJ,1))*DP1(IJ)
              DPNO = DPH(IJ) + DPH(KLAT(IJ,2))*DP2(IJ)
              DPN(IJ) = -DPNO+ABS(DPNO)
              DPS(IJ) =  DPSO+ABS(DPSO)
              DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

              DTHP = DRGP(IJ)*GOM(M) + DRCP(IJ)
              DTHM = DRGM(IJ)*GOM(M) + DRCM(IJ)
              DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
              DTP(IJ) = -DTHP+ABS(DTHP)
              DTM(IJ) =  DTHM+ABS(DTHM)

              DTHP =  SL(IJ,K,NFRE) * DFP
              DTC(IJ) =  DTC(IJ) + 2. * ABS(DTHP)
              DOP(IJ) = (-DTHP+ABS(DTHP))/FRATIO
              DOM(IJ) = ( DTHP+ABS(DTHP))*FRATIO
            ENDDO
          ELSE
!SHALLOW

!*    4.1.4.2 SHALLOW WATER.
!             --------------

            DFP = DELFR0/FR(M)
            DFM = DELFR0/FR(MM1)

!*    4.1.4.2.2 LON/LAT/DIR WEIGHTS IN INTEGRATION SCHEME.
!               ------------------------------------------

            DLA(0) = SD*CGOND(0,M)*DELLA0(NINF-1)
            DPH(0) = CD*CGOND(0,M)*DELPH0
            DLA(NINF-1) = SD*CGOND(NINF-1,M)*DELLA0(NINF-1)
            DPH(NINF-1) = CD*CGOND(NINF-1,M)*DELPH0
            DO IJ=NINF,IJLT(IG)
              DLA(IJ) =(U(IJ,IG)+SD*CGOND(IJ,M))*DELLA0(IJ)*DCO(IJ)
              DPH(IJ) =(V(IJ,IG)+CD*CGOND(IJ,M))*DELPH0
            ENDDO
            DO IJ=MIJS,MIJL
              DLWE = DLA(IJ) + DLA(KLON(IJ,1))
              DLEA = DLA(IJ) + DLA(KLON(IJ,2))
              DLE(IJ) = -DLEA+ABS(DLEA)
              DLW(IJ) =  DLWE+ABS(DLWE)
              DTC(IJ) =  DLEA+ABS(DLEA)-DLWE+ABS(DLWE)

              DPSO = DPH(IJ) + DPH(KLAT(IJ,1))*DP1(IJ)
              DPNO = DPH(IJ) + DPH(KLAT(IJ,2))*DP2(IJ)
              DPN(IJ) = -DPNO+ABS(DPNO)
              DPS(IJ) =  DPSO+ABS(DPSO)
              DTC(IJ) = DTC(IJ) + DPNO+ABS(DPNO)-DPSO+ABS(DPSO)

              DTHP=DRGP(IJ)*CGOND(IJ,M)
     &         +FCONST(IJ,M)*DRDP(IJ)+DRCP(IJ)
              DTHM=DRGM(IJ)*CGOND(IJ,M)
     &         +FCONST(IJ,M)*DRDM(IJ)+DRCM(IJ)
              DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
              DTP(IJ) = -DTHP+ABS(DTHP)
              DTM(IJ) =  DTHM+ABS(DTHM)

              DTHP = (SL(IJ,K,M) + SL(IJ,K,MP1))*DFP
              DTHM = (SL(IJ,K,M) + SL(IJ,K,MM1))*DFM
              DTC(IJ) =  DTC(IJ) + DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
              DOP(IJ) = (-DTHP+ABS(DTHP))/FRATIO
              DOM(IJ) = ( DTHM+ABS(DTHM))*FRATIO
            ENDDO
!SHALLOW
          ENDIF

!*    4.1.4.3 LOOP OVER GRIDPOINTS.
!             ---------------------

          DO IJ = MIJS,MIJL
            F3(IJ,K,M) = (1.-DTC(IJ))*F1(IJ,K,M )
     &       + DPN(IJ) * F1(KLAT(IJ,2),K  ,M)
     &       + DPS(IJ) * F1(KLAT(IJ,1),K  ,M)
     &       + DLE(IJ) * F1(KLON(IJ,2),K  ,M)
     &       + DLW(IJ) * F1(KLON(IJ,1),K  ,M)
     &       + DTP(IJ) * F1(IJ        ,KP1,M)
     &       + DTM(IJ) * F1(IJ        ,KM1,M)
     &       + DOP(IJ) * F1(IJ        ,K  ,MP1)
     &       + DOM(IJ) * F1(IJ        ,K  ,MM1)
          ENDDO

!*    BRANCH BACK TO 4.1.4 FOR NEXT FREQUENCY.

        ENDDO

!*    BRANCH BACK TO 4.2 FOR NEXT DIRECTION.

      ENDDO

!*    4.4 END OF PROPAGATION FOR SPHERICAL GRID
!*        WITH DEPTH AND CURRENT REFRACTION, RETURN.
!         ------------------------------------------

      DEALLOCATE(DELLA0,DCO,DP1,DP2,DPN,DPS)
      DEALLOCATE(DLE,DLW,DPH,DLA)
      DEALLOCATE(DTP,DTM,DRGP,DRGM,DTC)
      DEALLOCATE(WADVM1)
      IF (ISHALLO.NE.1) DEALLOCATE(CGOND)

      DEALLOCATE(DOP,DOM,DRCP,DRCM)

      IF (ISHALLO.NE.1) THEN
        DEALLOCATE(DRDP)
        DEALLOCATE(DRDM)
      ENDIF

      IF (IGL.NE.1) THEN
        IF(ALLOCATED(THDC))DEALLOCATE(THDC)
        IF(ALLOCATED(THDD))DEALLOCATE(THDD)
      ENDIF

      RETURN
      END SUBROUTINE PROPAGS
