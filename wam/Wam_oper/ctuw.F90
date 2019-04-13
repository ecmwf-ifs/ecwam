      SUBROUTINE CTUW (MIJS, MIJL, IG, LCFLFAIL)
! ----------------------------------------------------------------------

!**** *CTUW* - COMPUTATION OF THE CONER TRANSPORT SCHEME WEIGHTS.


!*    PURPOSE.
!     --------

!       COMPUTATION OF THE CORNER TRANSPORT UPSTREAM WEIGHTS
!       USED IN THE PROPAGATION FOR A GIVEN TIME STEP.

!**   INTERFACE.
!     ----------

!       *CALL* *CTUW(MIJS, MIJL, IG, LCFLFAIL)*
!          *MIJS*     - INDEX OF FIRST POINT.
!          *MIJL*     - INDEX OF LAST POINT.
!          *IG*       - BLOCK NUMBER.
!          *LCFLFAIL* - TRUE IF CFL CRITERION WAS VIOLATED.

!     METHOD.
!     -------


!     EXTERNALS.
!     ----------


!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : FR       ,GOM      ,DELTH    ,FRATIO   ,    &
     &            COSTH    ,SINTH
      USE YOWGRID  , ONLY : SINPH    ,COSPH    ,COSPHM1
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE YOWMAP   , ONLY : IXLG     , KXLT     ,IRGG    ,IPER     ,    &
     &            XDELLA   ,ZDELLO   ,AMOWEP   ,AMOSOP 
      USE YOWMPP   , ONLY : NINF
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY
      USE YOWPCONS , ONLY : PI       ,ZPI      ,R        ,CIRC
      USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
      USE YOWSHAL  , ONLY : NDEPTH   ,TCGOND   ,INDEP    ,DEPTH   ,     &
     &               TSIHKD
      USE YOWSTAT  , ONLY : IDELPRO  ,ICASE    ,ISHALLO  ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,WLAT     ,              &
     &            KCOR     ,WCOR     ,                                  &
     &            SUMWN    ,WLATN    ,WLONN    ,WCORN    ,              &
     &            WKPMN    ,WMPMN    ,OBSLAT   ,OBSLON   ,OBSCOR   ,    &
     &            LLWLATN  ,LLWLONN  ,LLWCORN  ,LLWKPMN  ,LLWMPMN  ,    &
     &            JXO      ,JYO      ,KCR      ,                        &
     &            LSAMEDEPTH

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: MIJS,MIJL,IG
      LOGICAL, DIMENSION(MIJS:MIJL), INTENT(INOUT) :: LCFLFAIL

      INTEGER(KIND=JWIM) :: IP,IJP
      INTEGER(KIND=JWIM) :: K,M,IJ,IC,IX,KY,KK,KKM
      INTEGER(KIND=JWIM) :: KP1,KM1,JH
      INTEGER(KIND=JWIM) :: ICL,ICR,ICC, ICRM, JCR
      INTEGER(KIND=JWIM) :: MP1, MM1
      INTEGER(KIND=JWIM) :: ISAMESIGN
      INTEGER(KIND=JWIM) :: ISSU(2),ISSV(2)

      REAL(KIND=JWRB) :: DELPRO, DTNEW
      REAL(KIND=JWRB) :: DXP, DYP, XM 
      REAL(KIND=JWRB) :: CMTODEG
      REAL(KIND=JWRB) :: GRIDAREA
      REAL(KIND=JWRB) :: CGYP, CGTH
      REAL(KIND=JWRB) :: DELTH0, DELFR0, SP, SM, DTHP, DTHM, DFP, DFM 
      REAL(KIND=JWRB) :: TANPH
      REAL(KIND=JWRB) :: DXX, DYY
      REAL(KIND=JWRB) :: UU, VV, UREL, VREL
      REAL(KIND=JWRB) :: XLAT, XLON
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(2) :: ADXP, ADYP
      REAL(KIND=JWRB), DIMENSION(2) :: DXUP, DXDW, DYUP, DYDW
      REAL(KIND=JWRB), DIMENSION(4) :: WEIGHT
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL) :: DRGP,DRGM
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL) :: DRDP,DRDM
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL) :: DRCP,DRCM
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL,2) :: CGX, CGY
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL,2) :: DP 
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL,2) :: WLATM1 
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL,4) :: WCORM1 
      REAL(KIND=JWRB), DIMENSION(MIJS:MIJL,NFRE) :: CGR
      REAL(KIND=JWRB), ALLOCATABLE :: SIGSI(:,:) 


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CTUW',0,ZHOOK_HANDLE)

      DELPRO = REAL(IDELPRO)   
      CMTODEG = 360.0_JWRB/CIRC

      DO IJ=MIJS,MIJL
        LCFLFAIL(IJ)=.FALSE.
      ENDDO

      IF (ISHALLO.NE.1) THEN
        DO M=1,NFRE
          DO IJ=MIJS,MIJL
            CGR(IJ,M)=TCGOND(INDEP(IJ),M)
          ENDDO
        ENDDO
      ELSE
        DO M=1,NFRE
          DO IJ=MIJS,MIJL
            CGR(IJ,M)=GOM(M)
          ENDDO
        ENDDO
      ENDIF

      DO IC=1,2
        DO IJ = MIJS,MIJL
          IF (KLAT(IJ,IC,1).GT. NINF-1 .AND.                            &
     &        KLAT(IJ,IC,2).GT. NINF-1) THEN
!           BOTH CLOSEST AND SECOND CLOSEST POINTS ARE OVER THE OCEAN
            WLATM1(IJ,IC) = 1. - WLAT(IJ,IC)
          ELSE IF (KLAT(IJ,IC,1).EQ. NINF-1) THEN
!           ADAPT CORNER POINT INTERPOLATION WEIGHT IF LAND IS PRESENT
!           CLOSEST POINT IS OVER LAND
            IF(WLAT(IJ,IC).LE.0.75_JWRB) WLAT(IJ,IC)=0.0_JWRB
            WLATM1(IJ,IC) = 1.0_JWRB - WLAT(IJ,IC)
          ELSE
!           ADAPT CORNER POINT INTERPOLATION WEIGHT IF LAND IS PRESENT
!           SECOND CLOSEST POINT IS OVER LAND
            IF(WLAT(IJ,IC).GE.0.5_JWRB) WLAT(IJ,IC)=1.0_JWRB
            WLATM1(IJ,IC) = 1.0_JWRB - WLAT(IJ,IC)
          ENDIF
        ENDDO
      ENDDO

      DO ICR=1,4
        DO IJ = MIJS,MIJL
          IF(KCOR(IJ,ICR,1).GT.NINF-1 .AND.                             &
     &       KCOR(IJ,ICR,2).GT.NINF-1) THEN
!           BOTH CLOSEST AND SECOND CLOSEST CORNER POINTS ARE OVER THE OCEAN
            WCORM1(IJ,ICR) = 1.0_JWRB - WCOR(IJ,ICR)
          ELSE IF(KCOR(IJ,ICR,1).EQ.NINF-1) THEN
!           ADAPT CORNER POINT INTERPOLATION WEIGHT IF LAND IS PRESENT
!           CLOSEST CORNER POINT IS OVER LAND
            IF(WCOR(IJ,ICR).LE.0.75_JWRB) WCOR(IJ,ICR)=0.0_JWRB
            WCORM1(IJ,ICR) = 1.0_JWRB - WCOR(IJ,ICR)
          ELSE
!           ADAPT CORNER POINT INTERPOLATION WEIGHT IF LAND IS PRESENT
!           SECOND CLOSEST CORNER POINT IS OVER LAND
            IF(WCOR(IJ,ICR).GT.0.5_JWRB) WCOR(IJ,ICR)=1.0_JWRB 
            WCORM1(IJ,ICR) = 1.0_JWRB - WCOR(IJ,ICR)
          ENDIF
        ENDDO
      ENDDO

      DO ICL=1,2
        DO IC=1,2
          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=MIJS,MIJL
                WLATN(IJ,K,M,IC,ICL)=0.0_JWRB
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO IC=1,2
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=MIJS,MIJL
              WLONN(IJ,K,M,IC)=0.0_JWRB
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO ICL=1,2
        DO ICR=1,4
          DO M=1,NFRE
            DO K=1,NANG
              DO IJ=MIJS,MIJL
                WCORN(IJ,K,M,ICR,ICL)=0.0_JWRB
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO

!*    ADVECTION IN PHYSICAL SPACE
!     =========================== 

!*    SPHERICAL OR CARTESIAN GRID?
!     ----------------------------
      IF (ICASE.EQ.1) THEN

!*      SPHERICAL GRID.
!       ---------------


!*        COMPUTE COS PHI FACTOR FOR ADJOINING GRID POINT.
!         (for all grid points)
          DO IC=1,2
            DO IJ = MIJS,MIJL
              KY=KXLT(IJ,IG)
              KK=KY+2*IC-3
              KKM=MAX(1,MIN(KK,NGY))
              DP(IJ,IC) = COSPH(KKM)*COSPHM1(IJ,IG)
            ENDDO
          ENDDO

!         FIND THE RELATIVE WEIGHT IN
!         THE CONER TRANSPORT UPSTREAM SCHEME.

!*        LOOP OVER DIRECTIONS.
!         ---------------------
          DO K=1,NANG

!*          LOOP OVER FREQUENCIES.
!           ----------------------
            DO M=1,NFRE

!             FIND MEAN GROUP VELOCITY COMPONENTS FOR DIRECTION TH(K)+180
!             -----------------------------------------------------------
              IF (ISHALLO.NE.1) THEN
!             SHALLOW WATER
                DO IC=1,2
                  DO IJ=MIJS,MIJL
                    IF(LSAMEDEPTH(IJ)) THEN
                      CGX(IJ,IC)=CGR(IJ,M)*SINTH(K)*COSPHM1(IJ,IG)
                      CGY(IJ,IC)=0.5_JWRB*CGR(IJ,M)*COSTH(K)*(1.+DP(IJ,IC))
                    ELSE
                      CGX(IJ,IC)=                                          &
     &                   0.5_JWRB*(CGR(IJ,M)+TCGOND(INDEP(KLON(IJ,IC)),M)) &
     &                      *SINTH(K)*COSPHM1(IJ,IG)
!                     IRREGULAR GRID
                      IF(IRGG.EQ.1) THEN
                        CGYP=WLAT(IJ,IC)*TCGOND(INDEP(KLAT(IJ,IC,1)),M)+       &
     &                   (1.0_JWRB-WLAT(IJ,IC))*TCGOND(INDEP(KLAT(IJ,IC,2)),M)
                      ELSE
!                     REGULAR GRID
                        CGYP=TCGOND(INDEP(KLAT(IJ,IC,1)),M)
                      ENDIF
                      CGY(IJ,IC)=0.5_JWRB*(CGR(IJ,M)+DP(IJ,IC)*CGYP)*COSTH(K)
                    ENDIF
                  ENDDO
                ENDDO
              ELSE
!             DEEP WATER
                DO IC=1,2
                  DO IJ=MIJS,MIJL
                    CGX(IJ,IC)=GOM(M)*SINTH(K)*COSPHM1(IJ,IG)
                    CGY(IJ,IC)=0.5_JWRB*GOM(M)*COSTH(K)*(1.0_JWRB+DP(IJ,IC))
                  ENDDO
                ENDDO
              ENDIF


!             LOOP OVER GRID POINTS
!             ---------------------
              DO IJ=MIJS,MIJL
                IX=IXLG(IJ,IG)
                KY=KXLT(IJ,IG)

!               FLUX VELOCITUES AT THE GRID BOX INTERFACE 

                DO IC=1,2

                  IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3 ) THEN
                    UU=U(IJ,IG)*COSPHM1(IJ,IG)
                    UREL=CGX(IJ,IC)+UU
                    ISSU(IC)=ISAMESIGN(UREL,CGX(IJ,IC))
                    VV=V(IJ,IG)*0.5_JWRB*(1.0_JWRB+DP(IJ,IC))
                    VREL=CGY(IJ,IC)+VV
                    ISSV(IC)=ISAMESIGN(VREL,CGY(IJ,IC))
                  ELSE
                    UREL=CGX(IJ,IC)
                    ISSU(IC)=1
                    VREL=CGY(IJ,IC)
                    ISSV(IC)=1
                  ENDIF
                  DXP=-DELPRO*UREL*CMTODEG
                  DYP=-DELPRO*VREL*CMTODEG 

                  ADXP(IC)=ABS(DXP)
                  ADYP(IC)=ABS(DYP)

!                 BASIC CFL CHECKS (IN EACH DIRECTION)
!                 ----------------
                  IF(ADXP(IC).GT.ZDELLO(KY))THEN
                    WRITE (IU06,*) '********************************'
                    WRITE (IU06,*) '* CTUW:                        *'
                    WRITE (IU06,*) '* CFL VIOLATED IN X DIRECTION. *'
                    WRITE (IU06,*) '* ADXP SHOULD BE < ZDELLO, BUT *'
                    WRITE (IU06,*) '* ADXP = ',ADXP(IC),IC
                    WRITE (IU06,*) '* ZDELLO = ',ZDELLO(KY)
                    DTNEW=ZDELLO(KY)*DELPRO/ADXP(IC)
                    WRITE (IU06,*) '* TIME STEP SHOULD BE REDUCED TO',  &
     &                              DTNEW
                    WRITE (IU06,*) '*                              *'
                    WRITE (IU06,*) '********************************'
                    LCFLFAIL(IJ)=.TRUE.
                  ENDIF
                  IF(ADYP(IC).GT.XDELLA)THEN
                    WRITE (IU06,*) '********************************'
                    WRITE (IU06,*) '* CTUW:                        *'
                    WRITE (IU06,*) '* CFL VIOLATED IN Y DIRECTION. *'
                    WRITE (IU06,*) '* ADYP SHOULD BE < XDELLA, BUT *'
                    WRITE (IU06,*) '* ADYP = ',ADYP(IC),IC
                    WRITE (IU06,*) '* XDELLA = ',XDELLA
                    DTNEW=XDELLA*DELPRO/ADYP(IC)
                    WRITE (IU06,*) '* TIME STEP SHOULD BEREDUCED TO',   &
     &                              DTNEW
                    WRITE (IU06,*) '*                              *'
                    WRITE (IU06,*) '********************************'
                    LCFLFAIL(IJ)=.TRUE.
                  ENDIF

                  DXUP(IC)=ADXP(IC)*ISSU(IC)
                  DXDW(IC)=ADXP(IC)*(1-ISSU(IC))
                  DYUP(IC)=ADYP(IC)*ISSV(IC)
                  DYDW(IC)=ADYP(IC)*(1-ISSV(IC))
    
                ENDDO

!                 GET ADVECTION WEIGHT FOR ALL NEIGHBOURING GRID POINTS

                  DXX=ZDELLO(KY)-DXUP(JXO(K,2))-DXDW(JXO(K,1))
                  DYY=XDELLA-DYUP(JYO(K,2))-DYDW(JYO(K,1))

                  GRIDAREA =  ZDELLO(KY)*XDELLA

!                 WEIGHTED CONTRIBUTION FROM NORTH-SOUTH DIRECTION (WLATN)

                  WEIGHT(JYO(K,1))=DXX*DYUP(JYO(K,1))/GRIDAREA
                  WEIGHT(JYO(K,2))=DXX*DYDW(JYO(K,2))/GRIDAREA
                  DO IC=1,2
                    WLATN(IJ,K,M,IC,1)=WLAT(IJ,IC)*WEIGHT(IC)
                    WLATN(IJ,K,M,IC,2)=WLATM1(IJ,IC)*WEIGHT(IC)
                  ENDDO

!                 WEIGHTED CONTRIBUTION FROM EAST-WEST DIRECTION (WLONN)

                  WLONN(IJ,K,M,JXO(K,1))=DYY*DXUP(JXO(K,1))/GRIDAREA
                  WLONN(IJ,K,M,JXO(K,2))=DYY*DXDW(JXO(K,2))/GRIDAREA


!                 CONTRIBUTION FROM CORNERS (KCOR)
                  WEIGHT(1)=DXUP(JXO(K,1))*DYUP(JYO(K,1))/GRIDAREA
                  WEIGHT(2)=DXDW(JXO(K,2))*DYUP(JYO(K,1))/GRIDAREA
                  WEIGHT(3)=DXUP(JXO(K,1))*DYDW(JYO(K,2))/GRIDAREA
                  WEIGHT(4)=DXDW(JXO(K,2))*DYDW(JYO(K,2))/GRIDAREA
                  DO ICR=1,4
                    WCORN(IJ,K,M,ICR,1)=WCOR(IJ,KCR(K,ICR))*WEIGHT(ICR)
                    WCORN(IJ,K,M,ICR,2)=WCORM1(IJ,KCR(K,ICR))*WEIGHT(ICR)
                  ENDDO

!                 CONTRIBUTIONS FOR IJ
                  SUMWN(IJ,K,M)=(ZDELLO(KY)*                            &
     &                            (DYDW(JYO(K,1))+DYUP(JYO(K,2))) +     &
     &                           XDELLA*                                &
     &                            (DXUP(JXO(K,2))+DXDW(JXO(K,1))) -     &
     &                           (DXDW(JXO(K,1))+DXUP(JXO(K,2)))*       &
     &                           (DYDW(JYO(K,1))+DYUP(JYO(K,2)))  )     &
     &                           /GRIDAREA

                ENDDO  ! END LOOP OVER GRID POINTS


            ENDDO  ! END LOOP OVER FREQUENCIES

          ENDDO  ! END LOOP OVER DIRECTIONS


      ELSE
!*    CARTESIAN GRID.
!     ---------------
        IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* CTUW:                                  *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '* FOR DEPTH OR/AND CURRENT REFRACTION !  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       ----------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* CTUW:                                  *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR  CARTESIAN GRID !                  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ENDIF

      ENDIF


!*    REFRACTIONS TERMS
!     ================= 

      DELTH0 = 0.25*DELPRO/DELTH

!*    GET SCATTER SIGMA/ SINH (2*K*D) TABLE 
!     -------------------------------------
      IF (IREFRA.NE.0 .AND. ISHALLO.NE.1 ) THEN
        ALLOCATE(SIGSI(MIJS:MIJL,NFRE))
        DO M=1,NFRE
          DO IJ=MIJS,MIJL
            SIGSI(IJ,M) = TSIHKD(INDEP(IJ),M)
          ENDDO
        ENDDO
      ENDIF

!*    LOOP OVER DIRECTIONS.
!     ---------------------

      DO K=1,NANG
        KP1 = K+1
        IF (KP1.GT.NANG) KP1 = 1
        KM1 = K-1
        IF (KM1.LT.1) KM1 = NANG

!*      COMPUTE GRID REFRACTION.
!       ------------------------
        SP  = DELTH0*(SINTH(K)+SINTH(KP1))/R
        SM  = DELTH0*(SINTH(K)+SINTH(KM1))/R

        DO IJ = MIJS,MIJL
          JH=KXLT(IJ,IG)
          TANPH = SINPH(JH)/COSPH(JH)
          DRGP(IJ) = TANPH*SP
          DRGM(IJ) = TANPH*SM
        ENDDO

!*      COMPUTE DEPTH REFRACTION.
!       -------------------------
        IF (IREFRA.EQ.1 .AND. ISHALLO.NE.1) THEN
          DO IJ = MIJS,MIJL
            DRDP(IJ) = (THDD(IJ,K) + THDD(IJ,KP1))*DELTH0
            DRDM(IJ) = (THDD(IJ,K) + THDD(IJ,KM1))*DELTH0
          ENDDO
        ELSE
          DO IJ = MIJS,MIJL
            DRDP(IJ) =  0.
            DRDM(IJ) =  0.
          ENDDO
        ENDIF

!*      COMPUTE CURRENT REFRACTION.
!       ---------------------------

        IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3 ) THEN

!!!!debile
          THDC(1849111,:)=0.0_JWRB
   
          DO IJ = MIJS,MIJL
            DRCP(IJ) = (THDC(IJ,K) + THDC(IJ,KP1))*DELTH0
            DRCM(IJ) = (THDC(IJ,K) + THDC(IJ,KM1))*DELTH0
          ENDDO
        ELSE
          DO IJ = MIJS,MIJL
            DRCP(IJ) = 0. 
            DRCM(IJ) = 0.
          ENDDO
        ENDIF


!*      REFRACTION WEIGHTS IN INTEGRATION SCHEME.
!       -----------------------------------------

!*      DEEP WATER OR
!*      SHALLOW WATER (NO DEPTH REFRACTION).
!       ------------------------------------
        IF (ISHALLO.EQ.1 .OR. IREFRA.EQ.0) THEN
          DO M=1,NFRE
            DO IJ=MIJS,MIJL
              DTHP = DRGP(IJ)*CGR(IJ,M) + DRCP(IJ)
              DTHM = DRGM(IJ)*CGR(IJ,M) + DRCM(IJ)
              WKPMN(IJ,K,M,0)=(DTHP+ABS(DTHP))+(ABS(DTHM)-DTHM)
              WKPMN(IJ,K,M,1)=-DTHP+ABS(DTHP)
              WKPMN(IJ,K,M,-1)=DTHM+ABS(DTHM)
            ENDDO
          ENDDO
        ELSE
!*      SHALLOW WATER AND DEPTH REFRACTION.
!       -----------------------------------
          DO M=1,NFRE
            DO IJ=MIJS,MIJL
              DTHP = DRGP(IJ)*CGR(IJ,M)+SIGSI(IJ,M)*DRDP(IJ)+DRCP(IJ)
              DTHM = DRGM(IJ)*CGR(IJ,M)+SIGSI(IJ,M)*DRDM(IJ)+DRCM(IJ)
              WKPMN(IJ,K,M,0)=(DTHP+ABS(DTHP))+(ABS(DTHM)-DTHM)
              WKPMN(IJ,K,M,1)=-DTHP+ABS(DTHP)
              WKPMN(IJ,K,M,-1)=DTHM+ABS(DTHM)
            ENDDO
          ENDDO
        ENDIF

!*      COMPUTE FREQUENCY SHIFTING DUE TO CURRENTS.
!       -------------------------------------------

        IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3 ) THEN

          DELFR0 = 0.25_JWRB*DELPRO/((FRATIO-1)*ZPI)

!*        DEEP WATER 
!         ----------
          IF (ISHALLO.EQ.1) THEN
            DO M=1,NFRE
              DFP = PI*(1.+FRATIO)*DELFR0
              DO IJ=MIJS,MIJL
                DTHP    = SDOT(IJ,K,NFRE) * DFP
                WMPMN(IJ,K,M,0) =2.0_JWRB* ABS(DTHP)
                WMPMN(IJ,K,M,1) =(-DTHP+ABS(DTHP))/FRATIO
                WMPMN(IJ,K,M,-1)=( DTHP+ABS(DTHP))*FRATIO
              ENDDO
            ENDDO
          ELSE
!*        SHALLOW WATER
!         -------------
            DO M=1,NFRE
              MP1 = MIN(NFRE,M+1)
              MM1 = MAX(1,M-1)
              DFP = DELFR0/FR(M)
              DFM = DELFR0/FR(MM1)

              DO IJ=MIJS,MIJL
                DTHP = (SDOT(IJ,K,M) + SDOT(IJ,K,MP1))*DFP
                DTHM = (SDOT(IJ,K,M) + SDOT(IJ,K,MM1))*DFM
                WMPMN(IJ,K,M,0) =(DTHP+ABS(DTHP))+(ABS(DTHM)-DTHM)
                WMPMN(IJ,K,M,1) =(-DTHP+ABS(DTHP))/FRATIO
                WMPMN(IJ,K,M,-1)=(DTHM+ABS(DTHM))*FRATIO
              ENDDO
            ENDDO
          ENDIF

        ENDIF

      ENDDO  ! END LOOP ON DIRECTIONS

      IF (IREFRA.NE.0 .AND. ISHALLO.NE.1 ) DEALLOCATE(SIGSI)
 

!     CHECK THAT WEIGHTS ARE LESS THAN 1
!     AND COMPUTE THEIR SUM AND CHECK IT IS LESS THAN 1 AS WELL
!!!   THE SUM IS NEEDED LATER ON !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=MIJS,MIJL
            DO IC=1,2
              DO ICL=1,2
              IF(WLATN(IJ,K,M,IC,ICL).GT.1.0 .OR.                       &
     &           WLATN(IJ,K,M,IC,ICL).LT.0.0      ) THEN
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '* CTUW:                           *'
                WRITE (IU06,*) '* CFL VIOLATED IN Y DIRECTION     *'
                WRITE (IU06,*) '* WLATN SHOULD BE < 1 AND > 0, BUT*'
                WRITE (IU06,*) '* WLATN(IJ,K,M,IC,ICL)=',IJ,K,M,IC,ICL, &
     &                            WLATN(IJ,K,M,IC,ICL)
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                LCFLFAIL(IJ)=.TRUE.
              ENDIF
              ENDDO
            ENDDO

            DO IC=1,2
              IF(WLONN(IJ,K,M,IC).GT.1.0 .OR.                           &
     &           WLONN(IJ,K,M,IC).LT.0.0      ) THEN
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '* CTUW:                           *'
                WRITE (IU06,*) '* CFL VIOLATED IN X DIRECTION     *'
                WRITE (IU06,*) '* WLONN SHOULD BE < 1 AND > 0, BUT*'
                WRITE (IU06,*) '* WLONN(IJ,K,M,IC)= ',IJ,K,M,IC,        &
     &                          WLONN(IJ,K,M,IC)
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                LCFLFAIL(IJ)=.TRUE.
              ENDIF
            ENDDO

            DO ICR=1,4
              DO ICL=1,2
              IF(WCORN(IJ,K,M,ICR,ICL).GT.1.0 .OR.                      &
     &           WCORN(IJ,K,M,ICR,ICL).LT.0.0     ) THEN
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '* CTUW:                           *'
                WRITE (IU06,*) '* CFL VIOLATED IN CORNER DIRECTION*'
                WRITE (IU06,*) '* WCORN SHOULD BE < 1 AND > 0, BUT*'
                WRITE (IU06,*) '*WCORN(IJ,K,M,ICR,ICL)=',IJ,K,M,ICR,ICL, &
     &                            WCORN(IJ,K,M,ICR,ICL)
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                LCFLFAIL(IJ)=.TRUE.
              ENDIF
              ENDDO
            ENDDO

            DO IC=-1,1
              IF(WKPMN(IJ,K,M,IC).GT.1.0 .OR.                           &
     &           WKPMN(IJ,K,M,IC).LT.0.0      ) THEN
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '* CTUW:                           *'
                WRITE (IU06,*) '* CFL VIOLATED IN DIRECTION SPACE *'
                WRITE (IU06,*) '* WKPMN SHOULD BE < 1 AND > 0, BUT*'
                WRITE (IU06,*) '* WKPMN(IJ,K,M,IC)= ',IJ,K,M,IC,        &
     &                            WKPMN(IJ,K,M,IC)
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                LCFLFAIL(IJ)=.TRUE.
              ENDIF
            ENDDO

            SUMWN(IJ,K,M)=SUMWN(IJ,K,M)+WKPMN(IJ,K,M,0)


            IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3 ) THEN

              DO IC=-1,1
                IF(WMPMN(IJ,K,M,IC).GT.1.0 .OR.                         &
     &             WMPMN(IJ,K,M,IC).LT.0.0      ) THEN
                  WRITE (IU06,*) '***********************************'
                  WRITE (IU06,*) '* CTUW:                           *'
                  WRITE (IU06,*) '* CFL VIOLATED IN FREQUENCY SPACE *'
                  WRITE (IU06,*) '* WMPMN SHOULD BE < 1 AND > 0, BUT*'
                  WRITE (IU06,*) '* WMPMN(IJ,K,M,IC)= ',IJ,K,M,IC,      &
     &                              WMPMN(IJ,K,M,IC)
                  WRITE (IU06,*) '*                                 *'
                  WRITE (IU06,*) '***********************************'

                  WRITE (0,*) '* CTUW: CFL VIOLATED IN FREQUENCY*',IJ,K,M,IC,WMPMN(IJ,K,M,IC),U(IJ,IG),V(IJ,IG)

                  LCFLFAIL(IJ)=.TRUE.
                ENDIF
              ENDDO

              SUMWN(IJ,K,M)=SUMWN(IJ,K,M)+WMPMN(IJ,K,M,0)
            ENDIF

!           SUM < 1  ?
            IF(SUMWN(IJ,K,M).GT.1.0 .OR. SUMWN(IJ,K,M).LT.0.0) THEN
              IX=IXLG(IJ,IG)
              KY=KXLT(IJ,IG)
              XLON=AMOWEP+(IX-1)*ZDELLO(KY)
              XLAT=AMOSOP+(KY-1)*XDELLA
              WRITE (IU06,*) '***********************************'
              WRITE (IU06,*) '* CTUW:                           *'
              WRITE (IU06,*) '* CFL VIOLATED                    *'
              WRITE (IU06,*) '* SUMW SHOULD BE < 1 AND > 0, BUT*'
              WRITE (IU06,*) '* IJ, SUMWN(IJ) = ',IJ,SUMWN(IJ,K,M)
              WRITE (IU06,*) '* XLAT= ',XLAT,' XLON= ',XLON 
              WRITE (IU06,*) '* DEPTH= ',DEPTH(IJ,IG)
              WRITE (0,*) '* CTUW: SUMW SHOULD BE < 1 AND > 0, BUT*',IJ,SUMWN(IJ,K,M),XLAT,XLON,DEPTH(IJ,IG),U(IJ,IG),V(IJ,IG)
              DO IP=1,2
              DO IC=1,2
                IJP = KLAT(IJ,IC,IP)
                IF(IJP.NE.NINF-1) THEN
                  WRITE (IU06,*) '* DEPTH= ',IJP,IP,IC,DEPTH(IJP,IG)
                  WRITE (IU06,*) '*     U= ',U(IJP,IG)
                  WRITE (IU06,*) '*     V= ',V(IJP,IG)
                ELSE
                  WRITE (IU06,*) '* DEPTH= ',IJP,IP,IC,'LAND'
                ENDIF
              ENDDO
              ENDDO
              DO IC=1,2
                IJP = KLON(IJ,IC)
                IF(IJP.NE.NINF-1) THEN
                  WRITE (IU06,*) '* DEPTH= ',IJP,IC,DEPTH(IJP,IG)
                  WRITE (IU06,*) '*     U= ',U(IJP,IG)
                  WRITE (IU06,*) '*     V= ',V(IJP,IG)
                ELSE
                  WRITE (IU06,*) '* DEPTH= ',IJP,IC,'LAND'
                ENDIF
              ENDDO

              IF (IREFRA.EQ.2 .OR. IREFRA.EQ.3 ) THEN
              WRITE (IU06,*) '* U = ',U(IJ,IG),' V = ',V(IJ,IG)
              ENDIF
              WRITE (IU06,*) '*                                 *'
              WRITE (IU06,*) '***********************************'
              LCFLFAIL(IJ)=.TRUE.
            ENDIF

          ENDDO  ! END LOOP OVER GRID POINTS
        ENDDO  ! END LOOP OVER FREQUENCIES
      ENDDO  ! END LOOP OVER DIRECTIONS

      DO IJ=MIJS,MIJL
        IF(LCFLFAIL(IJ)) THEN
          IF (LHOOK) CALL DR_HOOK('CTUW',1,ZHOOK_HANDLE)
          RETURN
        ENDIF
      ENDDO


!!!!!!INCLUDE THE BLOCKING COEFFCIENTS INTO THE WEIGHTS OF THE
!     SURROUNDING POINTS.

      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=MIJS,MIJL

!           POINTS ON SURROUNDING LATITUDES 
            DO IC=1,2
              DO ICL=1,2
                WLATN(IJ,K,M,IC,ICL)=                                   &
     &                   WLATN(IJ,K,M,IC,ICL)*OBSLAT(IJ,M,IC)
              ENDDO
            ENDDO

!           POINTS ON SURROUNDING LONGITUDE
            DO IC=1,2
              WLONN(IJ,K,M,IC)=WLONN(IJ,K,M,IC)*OBSLON(IJ,M,IC)
            ENDDO

!           SURROUNDING CORNER POINTS
            DO ICR=1,4
              DO ICL=1,2
                WCORN(IJ,K,M,ICR,ICL)=                                  &
     &                   WCORN(IJ,K,M,ICR,ICL)*OBSCOR(IJ,M,KCR(K,ICR))
              ENDDO
            ENDDO

          ENDDO  ! END LOOP OVER GRID POINTS
        ENDDO  ! END LOOP ON FREQUENCIES
      ENDDO  ! END LOOP OVER DIRECTIONS

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CTUW',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE CTUW 

      INTEGER FUNCTION ISAMESIGN(A,B)
!       =1 IF A AND B HAVE THE SAME SIGN, 0 OTHERWISE
        IMPLICIT NONE
        REAL :: A,B
        IF(SIGN(1.,A).EQ.SIGN(1.,B)) THEN
          ISAMESIGN=1
        ELSE
          ISAMESIGN=0
        ENDIF
      END FUNCTION ISAMESIGN
