      SUBROUTINE MSWELL (FL1)
! ----------------------------------------------------------------------

!**** *MSWELL* - MAKES START SWELL FIELDS FOR WAMODEL.

!     J. BIDLOT     ECMWF    

!*    PURPOSE.
!     --------

!       TO GENERATE WAMODEL START FIELDS.

!**   INTERFACE.
!     ----------

!   *CALL* *MSWELL (FL1)
!      *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT 

!     METHOD.
!     -------

!       NLOC SWELL SYSTEMS ARE SPECIFIED (SEE BELOW), CENTERED AT
!       DIFFERENT (XLON0, XLAT0) COORDINATES.
!
!       EACH SWELL SPECTRUM AT (XLON0, XLAT0) IS SPECIFIED AS
!       F(OMEGA,THETA)=F1(OMEGA)*Q(THETA)
!       WHERE
!
!       F1(OMEGA)=(N+1)*m0*(OMEGA_p**N/OMEGA**(N+1))*EXP[-((N+1)/N)*(OMEGA_p/OMEGA)**N]
!       WITH
!          OMEGA=2PI*f ; f: wave frequency
!          m0 the zero spectral moment = Hs**2/(2PI) ; Hs the wave height (to be specified) 
!          OMEGA_p=2PI*f_p  ; f_p the peak frequency 
!          N=5
!
!       Q(THETA)=(8/3PI)*COS(THETA-THETA0)**4 for ABS(THETA-THATA0) >  PI/2 
!       Q(THETA)=0                            for ABS(THETA-THATA0) <= PI/2 
!       WITH
!           THETA the wave direction
!           THETA0 the mean wave direction of the swell system (to be specified)
!
!       THE SPECTRUM AT EACH (XLON0,XLAT0) is SPREAD ON THE SPHERE BY MULTIPLYING
!       F(OMEGA,THETA) BY SPEAD(XLON,XLAT)
!       WHERE
!       SPREAD(XLON,XLAT)=EXP[-AL*SQRT{[(XLON-XLON0)*COS(XLAT0)]**2 + (XLAT-XLAT0]**2}]
!       WITH
!          AL=2*R/L ; R the earth radius and
!          L the correlation distance of the initial disturbance to be specified (in m)
!
!       NOTE THAT XLON,XLAT,XLON0,XLAT0 WILL SPECIFIED IN DEGREE BUT ARE CONVERTED 
!       IN RADIAN FOR CALCULATIONS.
!

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       I. V. LAVRENOV, 2003: WIND-WAVES IN OCEANS, DYNAMICS AND NUMERICAl SIMILATIONS
!       SPRINGER-VERLAG, 376pp.
!       see pages 52-53

!       NONE.

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRID  , ONLY : IGL      ,COSPH
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : ZPI      ,RAD      ,CIRC
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL1 

! ----------------------------------------------------------------------

      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJLOC, KLOC, MLOC
      REAL, ALLOCATABLE, DIMENSION(:) :: H0,THETA0,OMEGAP,XLAT0,XLON0
      REAL, ALLOCATABLE, DIMENSION(:) :: Q0,XL,DISMIN
      REAL, ALLOCATABLE :: FL0(:,:,:)

!----------------------------------------------------------------------

        CQ0=16./(3*ZPI)

        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=NINF-1,NSUP
              FL1(IJ,K,M)=0.
            ENDDO
          ENDDO
        ENDDO

!       TOTAL NUMBER OF SWELL SYSTEMS
        NLOC=4

        ALLOCATE(FL0(NLOC,NANG,NFRE))
        ALLOCATE(Q0(NANG))
        ALLOCATE(H0(NLOC))
        ALLOCATE(THETA0(NLOC))
        ALLOCATE(OMEGAP(NLOC))
        ALLOCATE(XL(NLOC))
        ALLOCATE(XLAT0(NLOC))
        ALLOCATE(XLON0(NLOC))
        ALLOCATE(IJLOC(NLOC))
        ALLOCATE(KLOC(NLOC))
        ALLOCATE(MLOC(NLOC))
        ALLOCATE(DISMIN(NLOC))

!       DEFINE THE SWELL SYSTEMS
        H0(1)=2.0
        THETA0(1)=135.
        OMEGAP(1)=0.3117
        XL(1)=250000.
        XLAT0(1)=47.
        XLON0(1)=165.

        H0(2)=2.0
        THETA0(2)=90.
        OMEGAP(2)=0.3117
        XL(2)=200000.
        XLAT0(2)=-50.
        XLON0(2)=20.

        H0(3)=2.0
        THETA0(3)=180.
        OMEGAP(3)=0.3117
        XL(3)=200000.
        XLAT0(3)=35.
        XLON0(3)=331

        H0(4)=2.0
        THETA0(4)=45.
        OMEGAP(4)=0.3117
        XL(4)=150000.
        XLAT0(4)=52.
        XLON0(4)=329

        NSP=5


        DO ILOC=1,NLOC
          THETA0(ILOC)=RAD*THETA0(ILOC)
        ENDDO

        DO ILOC=1,NLOC
          KLOC(ILOC)=1
          COSDIRMAX=COS(TH(K)-THETA0(1))
          DO K=2,NANG
            COSDIR=COS(TH(K)-THETA0(ILOC))
            IF(COSDIRMAX.LT.COSDIR) THEN
              KLOC(ILOC)=K
              COSDIRMAX=COSDIR
            ENDIF
          ENDDO
          MLOC(ILOC)=1
          OMEGA=ZPI*FR(1)
          OMEGADIFFMIN=ABS(OMEGA-OMEGAP(ILOC))
          DO M=2,NFRE
            OMEGA=ZPI*FR(M)
            OMEGADIFF=ABS(OMEGA-OMEGAP(ILOC))
            IF(OMEGADIFF.LT.OMEGADIFFMIN) THEN
              MLOC(ILOC)=M
              OMEGADIFFMIN=OMEGADIFF
            ENDIF
          ENDDO
        ENDDO

        DO ILOC=1,NLOC
          DO K=1,NANG
            COSDIR=COS(TH(K)-THETA0(ILOC))
            IF(COSDIR.GT.0.) THEN
              Q0(K) = CQ0*COSDIR**4
            ELSE
              Q0(K) = 0. 
            ENDIF
          ENDDO
          E0=H0(ILOC)**2/16.
          CEX=FLOAT(NSP+1)/NSP
          CS0=(NSP+1)*E0*OMEGAP(ILOC)**NSP
          DO M=1,NFRE
            OMEGA=ZPI*FR(M)
            S0=(CS0/OMEGA**(NSP+1))*EXP(-CEX*(OMEGAP(ILOC)/OMEGA)**NSP)
            IF(S0.LT.0.001) S0=0.
            DO K=1,NANG
              FL0(ILOC,K,M)= Q0(K)*S0
            ENDDO
          ENDDO
        ENDDO

        ZR=2*CIRC/ZPI

        DO ILOC=1,NLOC
          DISMIN(ILOC)=360.
        ENDDO

        IG=IGL
        DO IJ=NINF,NSUP
          IX = IXLG(IJ,IG)
          JSN= KXLT(IJ,IG)
          XLON=AMOWEP+(IX-1)*ZDELLO(JSN)
          XLAT=AMOSOP + (JSN-1)*XDELLA

          DO ILOC=1,NLOC
            DPHI=XLAT-XLAT0(ILOC)
            IF(ABS(DPHI).GT.180.) THEN
              DPHI2=(ABS(DPHI)-180)**2
            ELSE
              DPHI2=DPHI**2
            ENDIF
            DIS=SQRT(((XLON-XLON0(ILOC))*COSPH(JSN))**2+DPHI2)
            IF(DIS.LT.DISMIN(ILOC)) THEN
              IJLOC(ILOC)=IJ
              DISMIN(ILOC)=DIS
            ENDIF
            DIS=(ZR/XL(ILOC))*DIS*RAD
            IF(DIS.LT.10.) THEN
             SPRD=EXP(-DIS)
              DO M=1,NFRE
                DO K=1,NANG
                  FL1(IJ,K,M)=FL1(IJ,K,M)+FL0(ILOC,K,M)*SPRD 
                ENDDO
              ENDDO
            ENDIF
          ENDDO

        ENDDO

!?????????????????????????????
! overwrite definition to only input energy at the grid points corresponding to the
! maximum for the prescribed mean direction and peak frequency
!!!        DO M=1,NFRE
!!!          DO K=1,NANG
!!!            DO IJ=NINF-1,NSUP
!!!              FL1(IJ,K,M)=0.
!!!            ENDDO
!!!          ENDDO
!!!        ENDDO
!!!        DO ILOC=1,NLOC
!!!          FL1(IJLOC(ILOC),KLOC(ILOC),MLOC(ILOC))=1.
!!!        ENDDO
!?????????????????????????????
        
        DEALLOCATE(FL0)
        DEALLOCATE(Q0)
        DEALLOCATE(H0)
        DEALLOCATE(THETA0)
        DEALLOCATE(OMEGAP)
        DEALLOCATE(XL)
        DEALLOCATE(XLAT0)
        DEALLOCATE(XLON0)
        DEALLOCATE(IJLOC)
        DEALLOCATE(KLOC)
        DEALLOCATE(MLOC)
        DEALLOCATE(DISMIN)

        IF (ITEST.GE.2) THEN
          WRITE(IU06,*) '    SUB. MSWELL: GENERATED',
     &     ' BLOCK OF SWELL SPECTRA'
          CALL FLUSH (IU06)
        ENDIF

      RETURN
      END SUBROUTINE MSWELL
