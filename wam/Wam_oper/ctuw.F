      SUBROUTINE CTUW (MIJS, MIJL, IG, LCFLFAIL)
! ----------------------------------------------------------------------

!**** *CTUW* - COMPUTATION OF THE CONER TRANSPORT SCHEME WEIGHTS.


!*    PURPOSE.
!     --------

!       COMPUTATION OF THE CORNER TRANSPORT UPSTREAM WEIGHT
!       USE IN THE PROPAGATION FOR A GIVEN TIME STEP.

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

      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : FR       ,GOM      ,DELTH    ,FRATIO   ,
     &            COSTH    ,SINTH
      USE YOWGRID  , ONLY : SINPH    ,COSPH    ,COSPHM1
      USE YOWMAP   , ONLY : IXLG     , KXLT     ,IRGG    ,IPER     ,
     &            XDELLA   ,ZDELLO   ,IJFROMIK  ,IFROMIJ ,KFROMIJ
      USE YOWMPP   , ONLY : NINF
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY
      USE YOWPCONS , ONLY : PI       ,ZPI      ,R        ,CIRC
      USE YOWREFD  , ONLY : THDD     ,THDC     ,SDOT
      USE YOWSHAL  , ONLY : NDEPTH   ,TCGOND   ,INDEP    ,DEPTH   ,
     &               TSIHKD
      USE YOWSTAT  , ONLY : IDELPRO  ,ICASE    ,ISHALLO  ,IREFRA
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2    ,WLAT     ,
     &             KCOR    ,
     &             KLATN   ,KLONN    ,SUMWN    ,WLATN    ,WLONN    ,
     &             WKPMN   ,OBSLAT   ,OBSLON   ,OBSCOR   ,LSAMEDEPTH

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: MIJS,MIJL,IG
      INTEGER :: K,M,IJ,IC,IX,KY,KK,KKM,IM2,IM3
      REAL, DIMENSION(MIJS:MIJL) :: IM 
      INTEGER :: JX,JY,JXO,JYO
      INTEGER :: IJP,IXP,KYP,ICR
      INTEGER :: KP1,KM1,JH
      INTEGER :: ICC,ICL 
      INTEGER :: ISIDE, ISDXP
 
      REAL :: DELPRO, DTNEW
      REAL :: DXP, DYP, XM 
      REAL, DIMENSION(MIJS:MIJL) :: XLON, XL
      REAL :: XLONP, XDIFF
      REAL, DIMENSION(2) :: ADXP, ADYP
      REAL :: XLEFT, XMID, XRIGHT, XLP, XRP
      REAL :: CMTODEG
      REAL :: GRIDAREA
      REAL :: CGYP, CGTH
      REAL :: DELTH0, SP, SM, DTHP, DTHM 
      REAL :: TANPH
      REAL :: OFFSET, AOFFSET, AREAM, AREAC, WP1, WM1
      REAL :: DXX, DYY
      REAL, DIMENSION(MIJS:MIJL,NFRE) :: CGR
      REAL, DIMENSION(MIJS:MIJL,2) :: CGX, CGY
      REAL, DIMENSION(MIJS:MIJL) :: DRGP,DRGM
      REAL, DIMENSION(MIJS:MIJL) :: DRDP,DRDM
      REAL :: DP(MIJS:MIJL,2)
      REAL, ALLOCATABLE :: SIGSI(:,:) 

      LOGICAL:: LCFLFAIL
      LOGICAL:: LLNOTOK
      LOGICAL :: LLNOTNINF

! ----------------------------------------------------------------------

      DELPRO = REAL(IDELPRO)   
      CMTODEG = 360./CIRC

      LCFLFAIL=.FALSE.

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

!*    ADVECTION IN PHYSICAL SPACE
!     =========================== 

!*    SPHERICAL OR CARTESIAN GRID?
!     ----------------------------
      IF (ICASE.EQ.1) THEN

!*      SPHERICAL GRID.
!       ---------------

        IF (IREFRA.NE.2 .AND. IREFRA.NE.3 ) THEN
!*      WITHOUT DEPTH OR/AND CURRENT REFRACTION.
!       ----------------------------------------

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

!*        LOOP OVER DIRECTIONS.
!         ---------------------
          DO K=1,NANG
            IF(COSTH(K).GE.0.) THEN
              JY=1
              JYO=2
            ELSE
              JY=2
              JYO=1
            ENDIF
            IF(SINTH(K).GE.0.) THEN
              JX=1
              JXO=2
              ISIDE=1
              ISDXP=-1
            ELSE
              JX=2
              JXO=1
              ISIDE=0
              ISDXP=1
            ENDIF

!           FIND LOCATION OF GRID BOX ADVECTED AND RESHAPED UPSTREAM
!           DETERMINED THE NEIGHBORING GRID POINTS IN
!           THE CONER TRANSPORT UPSTREAM SCHEME.
            DO IJ=MIJS,MIJL
              IX=IXLG(IJ,IG)
              KY=KXLT(IJ,IG)
              XLON(IJ) = REAL(IX-1)*ZDELLO(KY)

!             FIND GRID POINT ON LATITUDE KY+2*JY-3 (I.E. KY+-1)
!             AND CLOSEST TO CORNER POINT

              KK=KY+2*JY-3
              KKM=MAX(1,MIN(KK,NGY))
              XL(IJ)=XLON(IJ)+ISDXP*ZDELLO(KY)
              XM = XL(IJ)/ZDELLO(KKM)
              IF(ZDELLO(KY).EQ.ZDELLO(KKM) .OR.
     &           IRGG.EQ.0 ) THEN
                IM(IJ) = NINT(XM) + 1
              ELSE
                IF(XM.GE.0.) THEN
                  IM(IJ) = INT(XM) + 1
                ELSE
                  IM(IJ) = INT(XM)
                ENDIF
                IM(IJ) = IM(IJ) + ISIDE
              ENDIF

              KLATN(IJ,K,0)=IJFROMIK(IM(IJ),KK,IG)

!             FIND SECOND CLOSEST
              IM2=IM(IJ)+ISDXP
              KLATN(IJ,K,1)=IJFROMIK(IM2,KK,IG)

!             FIND THIRD CLOSEST
              IM3=IM(IJ)-ISDXP
              KLATN(IJ,K,-1)=IJFROMIK(IM3,KK,IG)

!             GET NEIBOURING GRID POINT ON SAME LATITUDE.
              KLONN(IJ,K)=KLON(IJ,JX)

            ENDDO  ! END LOOP OVER GRID POINTS

!           FIND THE RELATIVE WEIGHT IN
!           THE CONER TRANSPORT UPSTREAM SCHEME.

!*          LOOP OVER FREQUENCIES.
!           ----------------------
            DO M=1,NFRE

!             FIND MEAN GROUP VELOCITY COMPONENTS FOR DIRECTION TH(K)+180
              IF (ISHALLO.NE.1) THEN
!             SHALLOW WATER
                DO IC=1,2
                  DO IJ=MIJS,MIJL
                    IF(LSAMEDEPTH(IJ)) THEN
                      CGX(IJ,IC)=CGR(IJ,M)*COSPHM1(IJ,IG)
                      CGY(IJ,IC)=0.5*CGR(IJ,M)*(1.+DP(IJ,IC))
                    ELSE
                      CGX(IJ,IC)=
     &                   0.5*(CGR(IJ,M)+TCGOND(INDEP(KLON(IJ,IC)),M))
     &                      *COSPHM1(IJ,IG)
!                     IRREGULAR GRID
                      IF(IRGG.EQ.1) THEN
                        CGYP=WLAT(IJ,IC)*TCGOND(INDEP(KLAT(IJ,IC)),M) +
     &                    (1.-WLAT(IJ,IC))*TCGOND(INDEP(KLAT2(IJ,IC)),M)
                      ELSE
!                     REGULAR GRID
                        CGYP=TCGOND(INDEP(KLAT(IJ,IC)),M)
                      ENDIF
                      CGY(IJ,IC)=0.5*(CGR(IJ,M) + DP(IJ,IC)*CGYP)
                    ENDIF
                  ENDDO
                ENDDO
              ELSE
!             DEEP WATER
                DO IC=1,2
                  DO IJ=MIJS,MIJL
                    CGX(IJ,IC)=GOM(M)*COSPHM1(IJ,IG)
                    CGY(IJ,IC)=0.5*GOM(M)*(1.+DP(IJ,IC))
                  ENDDO
                ENDDO
              ENDIF

              DO IJ=MIJS,MIJL
                IX=IXLG(IJ,IG)
                KY=KXLT(IJ,IG)

!               FIND GRID POINT ON LATITUDE KY+2*JY-3 (I.E. KY+-1)
!               AND CLOSEST TO CORNER POINT

                DO IC=1,2
                  DXP=-DELPRO*CGX(IJ,IC)*SINTH(K)*CMTODEG
                  DYP=-DELPRO*CGY(IJ,IC)*COSTH(K)*CMTODEG 

                  ADXP(IC)=ABS(DXP)
                  ADYP(IC)=ABS(DYP)

!                 BASIC CFL CHECKS (IN EACH DIRECTION)
!                 ----------------
                  IF(ADXP(IC).GT.ZDELLO(KY))THEN
                    WRITE (IU06,*) '***********************************'
                    WRITE (IU06,*) '* CTUW:                           *'
                    WRITE (IU06,*) '* CFL VIOLATED IN X DIRECTION.    *'
                    WRITE (IU06,*) '* ADXP SHOULD BE < ZDELLO, BUT    *'
                    WRITE (IU06,*) '* ADXP = ',ADXP(IC),IC
                    WRITE (IU06,*) '* ZDELLO = ',ZDELLO(KY)
                    DTNEW=ZDELLO(KY)*DELPRO/ADXP(IC)
                    WRITE (IU06,*) '* TIME STEP SHOULD REDUCED TO',DTNEW
                    WRITE (IU06,*) '*                                 *'
                    WRITE (IU06,*) '***********************************'
                    LCFLFAIL=.TRUE.
                  ENDIF
                  IF(ADYP(IC).GT.XDELLA)THEN
                    WRITE (IU06,*) '***********************************'
                    WRITE (IU06,*) '* CTUW:                           *'
                    WRITE (IU06,*) '* CFL VIOLATED IN Y DIRECTION.    *'
                    WRITE (IU06,*) '* ADYP SHOULD BE < XDELLA, BUT    *'
                    WRITE (IU06,*) '* ADYP = ',ADYP(IC),IC
                    WRITE (IU06,*) '* XDELLA = ',XDELLA
                    DTNEW=XDELLA*DELPRO/ADYP(IC)
                    WRITE (IU06,*) '* TIME STEP SHOULD REDUCED TO',DTNEW
                    WRITE (IU06,*) '*                                 *'
                    WRITE (IU06,*) '***********************************'
                    LCFLFAIL=.TRUE.
                  ENDIF
                ENDDO

                KK=KY+2*JY-3
                KKM=MAX(1,MIN(KK,NGY))

!               GET ADVECTION WEIGHT FOR ALL NEIGHBOURING GRID POINTS

                DXX=ZDELLO(KY)-ADXP(JXO)
                DYY=XDELLA-ADYP(JYO)

                GRIDAREA =  ZDELLO(KY)*XDELLA

                AREAM=DXX*ADYP(JY)
                AREAC=ADXP(JX)*ADYP(JY)

                XLP=-0.5*ZDELLO(KY)
                XRP= 0.5*ZDELLO(KY)

!               WEIGHTED CONTRIBUTION FROM NORTH-SOUTH
!               (BOX AVERAGE)
                XMID = REAL(IM(IJ)-1)*ZDELLO(KKM)-XLON(IJ)
                XLEFT=XMID-0.5*ZDELLO(KKM)
                XRIGHT=XMID+0.5*ZDELLO(KKM)
                IF(XLP.GE.XLEFT .AND. XRP.LE.XRIGHT) THEN
                  WM1=1.
                ELSE
                  XLEFT=MAX(XLEFT,XLP)
                  XRIGHT=MIN(XRIGHT,XRP)
                  WM1=MIN(1.,ABS(XRIGHT-XLEFT)/ZDELLO(KY))
                ENDIF

!               WEIGHTED CONTRIBUTION FROM CORNER
!               (BOX AVERAGE)
                XMID = REAL(IM(IJ)-1)*ZDELLO(KKM)-XL(IJ)
                XLEFT=XMID-0.5*ZDELLO(KKM)
                XRIGHT=XMID+0.5*ZDELLO(KKM)
                IF(XLP.GE.XLEFT .AND. XRP.LE.XRIGHT) THEN
                  WP1=1.
                ELSE
                  XLEFT=MAX(XLEFT,XLP)
                  XRIGHT=MIN(XRIGHT,XRP)
                  IF(XRIGHT.GT.XLEFT) THEN
                    WP1=MIN(1.,(XRIGHT-XLEFT)/ZDELLO(KY))
                  ELSE
                    WP1=0.
                  ENDIF
                ENDIF

!               CONTRIBUTION FOR KLATN(IJ,K,0)

                WLATN(IJ,K,M,0)=(WM1*AREAM+WP1*AREAC)/GRIDAREA

!               CONTRIBUTION FOR KLATN(IJ,K,1)
                WLATN(IJ,K,M,1)=(1.-WP1)*AREAC/GRIDAREA

!               CONTRIBUTION FOR KLATN(IJ,K,-1)
                WLATN(IJ,K,M,-1)=(1.-WM1)*AREAM/GRIDAREA

!               CONTRIBUTION FOR KLONN(IJ,K)
                WLONN(IJ,K,M)=ADXP(JX)*DYY/GRIDAREA

!               CONTRIBUTION FOR IJ
                SUMWN(IJ,K,M)=(DXX*ADYP(JYO)+DYY*ADXP(JXO)+
     &                             ADXP(JXO)*ADYP(JYO))/GRIDAREA

              ENDDO  ! END LOOP OVER GRID POINTS

            ENDDO  ! END LOOP OVER FREQUENCIES

          ENDDO  ! END LOOP OVER DIRECTIONS

        ELSE
!*      DEPTH AND CURRENT REFRACTION.
!       -----------------------------
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* CTUW:                                  *'
          WRITE (IU06,*) '* CORNER TRANSPORT SCHEME NOT YET READY  *' 
          WRITE (IU06,*) '* FOR DEPTH OR/AND CURRENT REFRACTION !  *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1

        ENDIF

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
          WRITE (IU06,*) '* FOR DEPTH OR/AND CURRENT REFRACTION !  *'
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
        ENDIF


!*      REFRACTION WEIGHTS IN INTEGRATION SCHEME.
!       -----------------------------------------

!*      DEEP WATER OR
!*      SHALLOW WATER (NO DEPTH REFRACTION).
!       ------------------------------------
        IF (ISHALLO.EQ.1 .OR. IREFRA.EQ.0) THEN
          DO M=1,NFRE
            DO IJ=MIJS,MIJL
              DTHP = DRGP(IJ)*CGR(IJ,M)
              DTHM = DRGM(IJ)*CGR(IJ,M)
              WKPMN(IJ,K,M,0)=DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
              WKPMN(IJ,K,M,1)=-DTHP+ABS(DTHP)
              WKPMN(IJ,K,M,-1)=DTHM+ABS(DTHM)
            ENDDO
          ENDDO
        ELSE
!*      SHALLOW WATER AND DEPTH REFRACTION.
!       -----------------------------------
          DO M=1,NFRE
            DO IJ=MIJS,MIJL
              DTHP = DRGP(IJ)*CGR(IJ,M)+SIGSI(IJ,M)*DRDP(IJ)
              DTHM = DRGM(IJ)*CGR(IJ,M)+SIGSI(IJ,M)*DRDM(IJ)
              WKPMN(IJ,K,M,0)=DTHP+ABS(DTHP)-DTHM+ABS(DTHM)
              WKPMN(IJ,K,M,1)=-DTHP+ABS(DTHP)
              WKPMN(IJ,K,M,-1)=DTHM+ABS(DTHM)
            ENDDO
          ENDDO
        ENDIF

      ENDDO  ! END LOOP ON DIRECTIONS

      IF (IREFRA.NE.0 .AND. ISHALLO.NE.1 ) DEALLOCATE(SIGSI)
 

!     CHECK THAT WEIGHTS ARE LESS THAN 1
!     AND COMPUTE THEIR SUM AND CHECK IT IS LESS THAN 1 AS WELL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=MIJS,MIJL
            DO IC=-1,1
              IF(WLATN(IJ,K,M,IC).GT.1.0 .OR. 
     &           WLATN(IJ,K,M,IC).LT.0.0) THEN
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '* CTUW:                           *'
                WRITE (IU06,*) '* CFL VIOLATED IN Y DIRECTION     *'
                WRITE (IU06,*) '* WLATN SHOULD BE < 1 AND > 0, BUT*'
                WRITE (IU06,*) '* WLATN(IJ,K,M,IC)= ',IJ,K,M,IC,
     &                            WLATN(IJ,K,M,IC)
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                LCFLFAIL=.TRUE.
              ENDIF
            ENDDO

            IF(WLONN(IJ,K,M).GT.1.0 .OR. WLONN(IJ,K,M).LT.0.0) THEN
              WRITE (IU06,*) '***********************************'
              WRITE (IU06,*) '* CTUW:                           *'
              WRITE (IU06,*) '* CFL VIOLATED IN X DIRECTION     *'
              WRITE (IU06,*) '* WLONN SHOULD BE < 1 AND > 0, BUT*'
              WRITE (IU06,*) '* WLONN(IJ,K,M)= ',IJ,K,M,
     &                          WLONN(IJ,K,M)
              WRITE (IU06,*) '*                                 *'
              WRITE (IU06,*) '***********************************'
              LCFLFAIL=.TRUE.
            ENDIF

            DO IC=-1,1
              IF(WKPMN(IJ,K,M,IC).GT.1.0 .OR. 
     &           WKPMN(IJ,K,M,IC).LT.0.0) THEN
                WRITE (IU06,*) '***********************************'
                WRITE (IU06,*) '* CTUW:                           *'
                WRITE (IU06,*) '* CFL VIOLATED IN DIRECTION SPACE *'
                WRITE (IU06,*) '* WKPMN SHOULD BE < 1 AND > 0, BUT*'
                WRITE (IU06,*) '* WKPMN(IJ,K,M,IC)= ',IJ,K,M,IC,
     &                            WKPMN(IJ,K,M,IC)
                WRITE (IU06,*) '*                                 *'
                WRITE (IU06,*) '***********************************'
                LCFLFAIL=.TRUE.
              ENDIF
            ENDDO

            SUMWN(IJ,K,M)=SUMWN(IJ,K,M)+WKPMN(IJ,K,M,0)

!           SUM < 1  ?
            IF(SUMWN(IJ,K,M).GT.1.0 .OR. SUMWN(IJ,K,M).LT.0.0) THEN
              WRITE (IU06,*) '***********************************'
              WRITE (IU06,*) '* CTUW:                           *'
              WRITE (IU06,*) '* CFL VIOLATED                    *'
              WRITE (IU06,*) '* SUMW SHOULD BE < 1 AND > 0, BUT*'
              WRITE (IU06,*) '* SUMWN= ',IJ,SUMWN(IJ,K,M)
              WRITE (IU06,*) '*                                 *'
              WRITE (IU06,*) '***********************************'
              LCFLFAIL=.TRUE.
            ENDIF

          ENDDO  ! END LOOP OVER GRID POINTS
        ENDDO  ! END LOOP OVER FREQUENCIES
      ENDDO  ! END LOOP OVER DIRECTIONS

      IF(LCFLFAIL) RETURN


!     CHECK THAT THE KLATN BELONG TO THE PRE-DEFINED MAXIMUM HALO
!     AS DEFINED BY KLAT, KLAT2, AND KCOR. 
      DO IJ=MIJS,MIJL
        LLNOTOK=.TRUE.
        LLNOTNINF=.FALSE.
        OUTER : DO ICC=-1,1
            DO K=1,NANG

              IF(KLATN(IJ,K,ICC).NE.NINF-1) THEN
                LLNOTNINF=.TRUE.
                DO IC=1,2
                  IF(KLATN(IJ,K,ICC).EQ.KLAT(IJ,IC)) THEN
                    LLNOTOK=.FALSE.
                    EXIT OUTER
                  ENDIF
                ENDDO
                DO IC=1,2
                  IF(KLATN(IJ,K,ICC).EQ.KLAT2(IJ,IC)) THEN
                    LLNOTOK=.FALSE.
                    EXIT OUTER
                  ENDIF
                ENDDO
                DO ICL=-1,1
                  DO IC=1,4
                    IF(KLATN(IJ,K,ICC).EQ.KCOR(IJ,IC,ICL)) THEN
                      LLNOTOK=.FALSE.
                      EXIT OUTER
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF

            ENDDO  ! END LOOP OVER DIRECTIONS
        ENDDO  OUTER !  END LOOP ON THREE NEARESTS

        IF(LLNOTOK .AND. LLNOTNINF) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '* CTUW:                                  *'
          WRITE (IU06,*) '* A VALUE OF KLATN IS NOT IN THE HALO    *' 
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* IJ = ', IJ
          DO ICC=-1,1
            WRITE(IU06,*) '* ICC = ',ICC
            DO K=1,NANG
              WRITE(IU06,*) '* K = ',K
              DO IC=1,2
                WRITE(IU06,*) '* KLAT ',
     &                        IC,KLATN(IJ,K,ICC),KLAT(IJ,IC)
              ENDDO
              DO IC=1,2
                WRITE(IU06,*) '* KLAT2 ',
     &                        IC,KLATN(IJ,K,ICC),KLAT2(IJ,IC)
              ENDDO
              DO ICL=-1,1
                DO IC=1,4
                  WRITE(IU06,*)'* KCOR ',
     &                       ICL,IC,KLATN(IJ,K,ICC),KCOR(IJ,IC,ICL)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ENDIF

      ENDDO  ! END LOOP OVER GRID POINTS

!!!!!!INCLUDE THE BLOCKING COEFFCIENTS INTO THE WEIGHTS OF THE
!     SURROUNDING POINTS.

      DO K=1,NANG
        IF(COSTH(K).GE.0.) THEN
          JY=1
        ELSE
          JY=2
        ENDIF
        IF(SINTH(K).GE.0.) THEN
          JX=1
        ELSE
          JX=2
        ENDIF

        IF(JY.EQ.1) THEN
          IF(JX.EQ.1) THEN
            ICR=3
          ELSE
            ICR=2
          ENDIF
        ELSE
          IF(JX.EQ.1) THEN
            ICR=4
          ELSE
            ICR=1
          ENDIF
        ENDIF

        DO IJ=MIJS,MIJL
          IX=IXLG(IJ,IG)
          KY=KXLT(IJ,IG)

!        POINTS ON SURROUNDING LATITUDES 
          DO ICC=-1,1
            IJP=KLATN(IJ,K,ICC)
            IXP=IFROMIJ(IJP)
            KYP=KFROMIJ(IJP)
            IF(IXP.GT.0 .AND. KYP.GT.0) THEN
              XLONP = REAL(IXP-1)*ZDELLO(KYP)

              XDIFF=XLON(IJ)-XLONP
!             NOTE THE BY CONSTRUCT THE PERIODIC POINTS WILL
!             ALWAYS BE A CORNER POINT SO NO NEED TO TEST XDIFF
!             FOR PERIODICITY.
              IF(ABS(XDIFF).LT.ZDELLO(KY)) THEN
                DO M=1,NFRE
                  WLATN(IJ,K,M,ICC)=WLATN(IJ,K,M,ICC)*OBSLAT(IJ,M,JY)
                ENDDO
              ELSE
                DO M=1,NFRE
                  WLATN(IJ,K,M,ICC)=WLATN(IJ,K,M,ICC)*OBSCOR(IJ,M,ICR)
                ENDDO
              ENDIF
            ENDIF

          ENDDO  !  END LOOP ON THREE NEARESTS



!        POINTS ON SURROUNDING LONGITUDE
          DO M=1,NFRE
             WLONN(IJ,K,M)=WLONN(IJ,K,M)*OBSLON(IJ,M,JX)
          ENDDO

        ENDDO  ! END LOOP OVER GRID POINTS
      ENDDO  ! END LOOP OVER DIRECTIONS

! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE CTUW 
