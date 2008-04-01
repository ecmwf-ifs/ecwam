      SUBROUTINE UNSETICE(FL,NSTART,NEND,U10OLD,THWOLD)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      AUGUST 2006 

!*    PURPOSE.
!     --------
!     TRY TO CREATE A JONSWAP SPECTRUM AT POINTS WHICH WERE ICE IN THE
!     PREVIOUS RUN AND WHICH ARE NOW OPEN SEA.

!**   INTERFACE.
!     ----------
!     *CALL* *UNSETICE(FL,NSTART,NEND,U10OLD,THWOLD)
!     *FL*        ARRAY CONTAINING THE SPECTRA CONTRIBUTION ON EACH PE
!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *U10OLD*    WIND SPEED. (used with fetch law to fill empty 
!                 sea points)
!     *THWOLD*    WIND DIRECTION (RADIANS).


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRID  , ONLY : DELPHI
      USE YOWJONS  , ONLY : AJONS    ,BJONS    ,DJONS    ,EJONS
      USE YOWICE   , ONLY : FLMINFR  ,ICEMASK
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,
     &            NBLO     ,CLDOMAIN
      USE YOWPCONS , ONLY : G        ,R        ,EPSMIN
      USE YOWSHAl  , ONLY : DEPTH
      USE YOWTABL  , ONLY : JUMAX    ,DELU
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : LRESTARTED
      USE MPL_MODULE

! ----------------------------------------------------------------------

      INTEGER,DIMENSION(NPROC) :: NSTART, NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD

! ----------------------------------------------------------------------

      REAL,ALLOCATABLE :: ET(:,:), STK(:)
      LOGICAL :: LICE2SEA(NINF:NSUP)

! ----------------------------------------------------------------------

      IF(.NOT.LRESTARTED .AND. CLDOMAIN.NE.'s') THEN
!       TRY TO CREATE A JONSWAP SPECTRUM AT POINTS WHICH WERE ICE IN THE
!       PREVIOUS RUN AND WHICH ARE NOW OPEN SEA.
!       THEY ARE CHARCTERISED BY HAVING NO WAVE ENERGY BESIDE NOISE
!       AND HAVE AN ICEMASK = 1.
        DO IJ=NSTART(IRANK),NEND(IRANK)
          LICE2SEA(IJ)=.TRUE.
        ENDDO
        DO K=1,NANG
          DO M=1,NFRE
            DO IJ=NSTART(IRANK),NEND(IRANK)
               IF(FL(IJ,K,M) .GT. EPSMIN) LICE2SEA(IJ) = .FALSE. 
            ENDDO
          ENDDO
        ENDDO

        ALLOCATE(ET(1,NFRE))
        ALLOCATE(STK(NANG))

        GAMMA=3.000000
        SA=7.000000E-02
        SB=9.000000E-02
        FPMAX=FR(NFRE-1)

        DO IJ=NSTART(IRANK),NEND(IRANK)
          IF(LICE2SEA(IJ) .AND. ICEMASK(IJ,1).EQ.1 ) THEN
            IF(DEPTH(IJ,1).LE.50.) THEN
              FETCH=MIN(DELPHI,50000.)
            ELSEIF(DEPTH(IJ,1).LE.150.) THEN
              FETCH=MIN(2*DELPHI,100000.)
            ELSEIF(DEPTH(IJ,1).LE.250.) THEN
              FETCH=MIN(3*DELPHI,200000.)
            ELSE
              FETCH=MIN(5*DELPHI,250000.)
            ENDIF
            GX = G * FETCH
!           FIND PEAK PERIOD AND ENERGY LEVEL FROM FETCH LAW
!           THE SAME FORMULATION AS IN SUBROUTINE PEAK IS USED.
            IF (U10OLD(IJ,1) .GT. 0.1E-08 .AND.
     &          DEPTH(IJ,1) .GT. 10.) THEN
              U10 = U10OLD(IJ,1)
              GXU = GX/(U10*U10)
              UG = G/U10
              FPK = AJONS * GXU ** DJONS
              FPK = MAX(0.13, FPK)
              FPK = MIN(FPK, FPMAX/UG)
              ALPHJO = BJONS * FPK** EJONS
              ALPHJO = MAX(ALPHJO, 0.0081)
              FPK = FPK*UG
            ELSE
              FPK = 0.
              ALPHJO = 0.
            ENDIF

            CALL JONSWAP(ALPHJO, GAMMA, SA, SB, FPK, 1, 1, ET)

            THESK=THWOLD(IJ,1)

            CALL SPR (NANG, NANG, THESK, TH, STK)

            DO K=1,NANG
              DO M=1,NFRE
                FL(IJ,K,M) = ET(1,M)*STK(K)
              ENDDO
            ENDDO
          ENDIF
        ENDDO

        DO K=1,NANG
          DO M=1,NFRE
            DO IJ=NSTART(IRANK),NEND(IRANK)
              XJ=U10OLD(IJ,1)/DELU
              JU=MIN(JUMAX, MAX(NINT(XJ),1))
              FLLOWEST=FLMINFR(JU,M)*MAX(0.,COS(TH(K)-THWOLD(IJ,1)))**2
              IF(FL(IJ,K,M) .LT. FLLOWEST) FL(IJ,K,M) = FLLOWEST
            ENDDO
          ENDDO
        ENDDO

        DEALLOCATE(ET)
        DEALLOCATE(STK)

        IF (ITEST.GE.2) THEN
         WRITE(IU06,*) ' UNSETICE: SPECTRA INITIALISED OVER OLD SEA ICE'
         WRITE(IU06,*) ' '
        ENDIF

      ENDIF

      RETURN
      END SUBROUTINE UNSETICE
