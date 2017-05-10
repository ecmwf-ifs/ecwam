      SUBROUTINE TOPOAR (IU01, BATHY)

! ----------------------------------------------------------------------

!**** *TOPOAR* - ARRANGE SUBGRID TOPOGRAPHY.

!     S. HASSELMANN     MPIFM           1/6/86.

!     MODIFIED BY       H. GUNTHER      1/4/90  -  REARANGEMENT OF CODE.
!                       J. BIDLOT    FEB 2003 - 
!                                    INTRODUCE OBSTRUCTION COEFFICIENTS
!                                    MOVED TO MUBUF NOW !
!                       J. BIDLOT    NOV 2005 - 
!                                    INTRODUCE NEW INPUT FORMAT FOR
!                                    BATHYMETRY (the old format can still
!                                    be used).
!*    PURPOSE.
!     --------

!       TO READ IN TOPOGRAPHY ON INPUT GRID AND CONVERT TO OUTPUT GRID.

!**   INTERFACE.
!     ----------

!       *CALL* *TOPOAR (IU01,BATHY)*
!          *IU01*  -  LOGICAL INPUT UNIT OF TOPOGRAPHIC DATA.
!          *BATHY*   -  BATHYMETRY DATA IN OUTPUT GRID.

!     METHOD.
!     -------

!       IF OLD TOPAGRAPHY INPUT:
!       THE TOPOGRAPHY MUST BE ON A REGULAR LATITUDE-LONGITUDE
!       GRID ARRANGED FROM SOUTH TO NORTH, AND FROM WEST TO EAST.
!       IT IS ASSUMED THAT NEGATIVE VALUES ARE SEA DEPTHS (WHICH
!       ARE CONVERTED TO POSITIVE) AND THAT POSITIVE VALUES ARE
!       LAND ELEVATIONS (WHICH ARE CONVERTED TO -999 IDENTIFIERS).

!       THE TOPOGRAPHIC DATA IS READ IN ON THE INPUT GRID AND IT
!       IS STORED ONLY FOR THOSE LATITUDES WITHIN THE REQUESTED
!       GRID. THEN THE TOPOGRAPHIC DATA IS FURTHER RESTRICTED
!       TO LIE WITHIN THE SUBGRID LONGITUDES. IT IS THEN PUT ON
!       THE REQUESTED SUBGRID LAT-LONG RESOLUTION, ALWAYS USING
!       THE NEAREST POINT.
!       FINALLY THE SUBGRID TOPOGRAPHY MAY BE MANUALLY ADJUSTED BY
!       MEANS OF THE CARD INPUT AND A PRINTER OUTPUT IS DONE.

!       FOR NEW TOPOGRAPHY INPUT:
!       THE INPUT IS ALREADY ON THE DESIRED GRID AND THE REDUCTION
!       FACTORS DUE TO OBSTRUCTIONS ARE GIVEN AS WELL
!       !!! THE REDUCTION FACTORS ARE READ AND PROCESSED IN MUBUF!

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *ADJUST*    - CORRECTS LONGITUDE INPUT.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWGRID  , ONLY : NLONRGG
      USE YOWCINP  , ONLY : NOUT     ,XOUTW    ,XOUTS    ,XOUTE    ,
     &            XOUTN    ,NOUTD
      USE YOWMAP   , ONLY : NX       ,NY       ,AMOWEP   ,AMOSOP   ,
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,
     &            LLOBSTRCT
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTHA   ,DEPTHD 
      USE YOWTEST  , ONLY : IU06, ITEST

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IU01
      REAL, DIMENSION(NGX,NGY), INTENT(OUT) :: BATHY

      INTEGER :: I, J, K, JH, L, IX, IAA
      INTEGER :: KLONRGG  
      INTEGER :: MLON, NLATMAX, KMAX, NLAT, N1, N2, LAST
      INTEGER :: NMINADJT
      INTEGER :: ILW, NLON, IH, NJ, JJ, NL, JRGG, KXLO, KAMOEAP
      INTEGER :: ILEN, IPAGE, IA, IE
      INTEGER, ALLOCATABLE :: IDUM(:)
      INTEGER, ALLOCATABLE :: IA2H(:), IA1(:,:)

      REAL :: BATHYMAX, TABLEMAX
      REAL :: XDELA, XDELO, XLAS, XLAN, XLOW, XLOE
      REAL :: XLAT, XLON, XLAG, XLW, XLA, XLO, XLOH, XLOG
      REAL, ALLOCATABLE :: XA2H(:), XA1(:,:)

      CHARACTER(LEN=1), ALLOCATABLE :: AX(:), AXX(:)
      CHARACTER(LEN=4) :: CX
      CHARACTER(LEN=10) :: FORMT
      CHARACTER(LEN=14) :: CHEADER 

      LOGICAL :: LLREALIN

! ----------------------------------------------------------------------

!*    1. READING THE TOPOGRAPHY OF THE INPUT GRID AND STORING THOSE
!*       LATITUDES WITHIN THE OUTPUT SUBGRID AREA.
!        -----------------------------------------

!       MLON - NUMBER OF TOPOGRAPHIC POINTS PER GRID LATITUDE.
!       KMAX - NUMBER OF RECORDS OF INPUT GRID PER LATITUDE.
!       XLAT - LATITUDE OF CURRENT GRID DATA.
!       NLAT - NUMBER OF GRID LATITUDES STORED.
!       XLAG - LATITUDE OF FIRST GRID LATITUDE STORED.

      REWIND (UNIT=IU01)
!     DETERMINE WHICH TYPE OF FILE IS USED
      READ (IU01,'(a14)') CHEADER 
      IF(CHEADER.EQ.'WAM BATHYMETRY') THEN
!       REAL BATHYMETRY INPUT
        LLREALIN=.TRUE.
        WRITE (IU06,*) ' BATHYMETRY FROM REAL INPUT DATA' 
      ELSE
!       INTEGER BATHYMETRY INPUT
        LLREALIN=.FALSE.
        WRITE (IU06,*) ' BATHYMETRY FROM INTEGER INPUT DATA' 
      ENDIF
      REWIND (UNIT=IU01)

      IF(LLREALIN) THEN
        READ (IU01,'(a14)') CHEADER 
        READ (IU01,'(8F13.8)') XDELA, XDELO, XLAS, XLAN, XLOW, XLOE
      ELSE
        READ (IU01,'(8F10.5)') XDELA, XDELO, XLAS, XLAN, XLOW, XLOE
      ENDIF
      CALL ADJUST (XLOW, XLOE)

      WRITE (IU06,'(1H1,'' INPUT GRID''/)')
      WRITE (IU06,'(3X,''RESOLUTION LAT-LON '',2F8.3)') XDELA, XDELO
      WRITE (IU06,'(3X,'' SOUTHERN LAT '','' NORTHERN LAT '',
     &                 '' WESTERN LONG '','' EASTERN LONG'',
     &                 /,2X,4F14.3)') XLAS, XLAN, XLOW, XLOE

      BATHY(:,:) = 999.

      IF(LLOBSTRCT) THEN
!     NEW TREATMENT OF BATHYMETRY INPUT
!     ---------------------------------
        WRITE (IU06,*) '  NEW TREATMENT OF BATHYMETRY INPUT '
        WRITE (IU06,*) ' '
!       TEST DOMAIN
        IF(NINT(10000*XDELA).NE.NINT(10000*XDELLA) .OR.
     &     NINT(10000*XDELO).NE.NINT(10000*XDELLO) .OR.
     &     NINT(10000*XLAS).NE.NINT(10000*AMOSOP) .OR.
     &     NINT(10000*XLAN).NE.NINT(10000*AMONOP) .OR.
     &     NINT(10000*XLOW).NE.NINT(10000*AMOWEP) .OR.
     &     NINT(10000*XLOE).NE.NINT(10000*AMOEAP) ) THEN
          WRITE (IU06,*) ' *******************************************'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR        *'
          WRITE (IU06,*) ' *      ===========================        *'
          WRITE (IU06,*) ' * THE INPUT AND COMPUTED DOMAIN DEFINITION*'
          WRITE (IU06,*) ' * DO NOT AGREE :'
          WRITE (IU06,*) ' * XDELLA ',XDELA,XDELLA
          WRITE (IU06,*) ' * XDELLO ',XDELO,XDELLO
          WRITE (IU06,*) ' * AMOSOP ',XLAS,AMOSOP
          WRITE (IU06,*) ' * AMONOP ',XLAN,AMONOP
          WRITE (IU06,*) ' * AMOWEP ',XLOW,AMOWEP
          WRITE (IU06,*) ' * AMOEAP ',XLOE,AMOEAP
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                 *'
          WRITE (IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF

        DO K=1,NY
!       READ THE NUMBER OF POINTS PER LATITUDE
          READ (IU01,'(I4.4)') KLONRGG  
          IF(KLONRGG.NE.NLONRGG(K)) THEN
            WRITE (IU06,*) ' ******************************************'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR       *'
            WRITE (IU06,*) ' *      ===========================       *'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' * THE INPUT NUMBER OF POINTS PER LATITUDE '
            WRITE (IU06,*) ' * IS NOT THE SAME AS THE COMPUTED ONE !!  '
            WRITE (IU06,*) ' * KLONRGG = ',KLONRGG,K 
            WRITE (IU06,*) ' * NLONRGG = ',NLONRGG(K),K 
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                *'
            WRITE (IU06,*) ' ******************************************'
            CALL ABORT1
          ENDIF
        ENDDO

!       READ THE BATHYMETRY DATA

        ALLOCATE(IDUM(NGX))
        DO K=1,NY
          WRITE(CX,'(I4.4)') NLONRGG(K)
          IF(LLREALIN) THEN
            FORMT='('//CX//'F9.2)'
            READ (IU01,FORMT) (BATHY(IX,K),IX=1,NLONRGG(K))
          ELSE
            FORMT='('//CX//'I4)'
            READ (IU01,FORMT) (IDUM(IX),IX=1,NLONRGG(K))
            DO IX=1,NLONRGG(K)
              BATHY(IX,K)=REAL(IDUM(IX))
            ENDDO
          ENDIF
        ENDDO
        DEALLOCATE(IDUM)
  
!       THE REDUCTION FACTORS ARE NOW READ IN MUBUF (the rest of IU01)

      ELSE
!     OLD TREATMENT OF BATHYMETRY INPUT
!     ---------------------------------
      MLON = NINT((XLOE-XLOW)/XDELO+1.0)

      NLATMAX=NINT((XLAN-XLAS)/XDELA+1.0)

      ALLOCATE(AX(MLON))
      ALLOCATE(IA2H(MLON), IA1(MLON,NLATMAX))
      ALLOCATE(XA2H(MLON), XA1(MLON,NLATMAX))

      KMAX=(MLON+11)/12

      XLAT=XLAS-XDELA
      NLAT=0
1005  CONTINUE
      NLAT=NLAT+1
      IF (NLAT.GT.NLATMAX) THEN
        WRITE (IU06,*) ' ******************************************'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR       *'
        WRITE (IU06,*) ' *      ===========================       *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' * NUMBER OF LATITUDES IN INPUT GRID      *'
        WRITE (IU06,*) ' * EXCEEDS DIMENSION                      *'
        WRITE (IU06,*) ' * DIMENSION IS      NLATMAX = ', NLATMAX
        WRITE (IU06,*) ' * LAST LATITUDE READ IS  XLAT  = ', XLAT
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                *'
        WRITE (IU06,*) ' ******************************************'
        CALL ABORT1
      ENDIF
1010  CONTINUE
      XLAT=XLAT+XDELA
      DO K=1,KMAX
        N1 = 12*(K-1)+1
        N2 = MIN(12*K,MLON)
        IF(LLREALIN) THEN
          READ (IU01, '(12(F9.2))', END=1017)
     &     (XA1(IAA,NLAT),IAA=N1,N2)
        ELSE
          READ (IU01, '(12(I5,A1))', END=1017)
     &     (IA1(IAA,NLAT),AX(IAA),IAA=N1,N2)
        ENDIF
      ENDDO

      IF (XLAT+0.5*XDELA.LE.AMOSOP) GO TO 1010
      IF(NLAT.EQ.1) XLAG=XLAT
      IF(.NOT.LLREALIN) THEN
      DO I=1,MLON
          IF (AX(I).EQ.'E'.AND.IA1(I,NLAT).LT.0)
     &     IA1(I,NLAT)=-IA1(I,NLAT)
          IF (AX(I).EQ.'D'.AND.IA1(I,NLAT).GT.0)
     &     IA1(I,NLAT)=-IA1(I,NLAT)
          IF (AX(I).EQ.'E'.AND.IA1(I,NLAT).EQ.0) IA1(I,NLAT)=1
        ENDDO
      ENDIF
      IF(XLAT+0.5*XDELA.GT.AMONOP) GO TO 1020
      GO TO 1005
1017  CONTINUE
      NLAT=NLAT-1
      WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IU06,*) ' +                                           +'
      WRITE (IU06,*) ' +     WARNING ERROR IN SUB. TOPOAR          +'
      WRITE (IU06,*) ' +     ============================          +'
      WRITE (IU06,*) ' +                                           +'
      WRITE (IU06,*) ' + END OF FILE ON INPUT GRID  UNIT :', IU01
      WRITE (IU06,*) ' + NORTH GRID BOUNDARY CHANGED TO XLAT = ',XLAT
      WRITE (IU06,*) ' +                                           +'
      WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++++'
1020  CONTINUE

! ----------------------------------------------------------------------

!*    2. STORING TOPOGRAPHIC DATA AT LONGITUDES WITHIN SUBGRID AREA.
!        ----------------------------------------------------------

!       ILW      -  INDEX OF NEAREST GRID LONGITUDE EQUIVALENT TO
!                     WESTERN SUBGRID BOUNDARY.
!       NLON     -  NUMBER OF GRID LONGITUDES WITHIN SUBGRID AREA.
!       IH,IH1   -  GRID LONGITUDE NUMBER.


      XLW= MOD(AMOWEP-XLOW+720.,360.)
!!!!      ILW= NINT(XLW/XDELO-0.0001)
      ILW= NINT(XLW/XDELO)

      NLON=NINT((AMOEAP-AMOWEP)/XDELO+1.0)
      DO J=1,NLAT
        IH=ILW
        DO I=1,NLON
          IH=IH+1
          IF(IH.LE.0) IH=IH+MLON
          IF(IH.GT.MLON) IH=IH-MLON
          IF(LLREALIN) THEN
            XA2H(I) =XA1(IH,J)
          ELSE
            IA2H(I) =IA1(IH,J)
          ENDIF
        ENDDO
        DO I=1,NLON
          IF(LLREALIN) THEN
            XA1(I,J)=XA2H(I)
          ELSE
            XA1(I,J)=REAL(IA2H(I))
          ENDIF
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. PUT ON XDELLO BY XDELLA SUBGRID.
!        --------------------------------

!     XLA  - LATITUDE OF OUTPUT GRID.
!     XLAG - LATITUDE OF INPUT GRID.
!     XLO  - LONGITUDE OF OUTPUT GRID.
!     XLOG - LONGITUDE OF INPUT GRID.

      XLA=AMOSOP
      XLAG=XLAG-XDELA
      NJ=0
      JJ=0

!  LOOP THROUGH LATITUDES

      DO J=1,NLAT
        XLAG=XLAG+XDELA
 3010   CONTINUE
        IF(XLA.LT.XLAG-0.5*XDELA.OR.XLA.GE.XLAG+0.5*XDELA) GO TO 3070
        NJ=NJ+1
        IF(NJ.GT.NGY) THEN
          WRITE (IU06,*) ' ******************************************'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR       *'
          WRITE (IU06,*) ' *      ===========================       *'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' * NUMBER OF LATITUDES IN OUTPUT GRID     *'
          WRITE (IU06,*) ' * EXCEEDS DIMENSION.                     *'
          WRITE (IU06,*) ' * DIMENSION IS            NGY = ', NGY
          WRITE (IU06,*) ' *                         NJ  = ', NJ
          WRITE (IU06,*) ' * LAST LATITUDE USED IS   XLA = ', XLA
          WRITE (IU06,*) ' *                         XLAG= ', XLAG
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                *'
          WRITE (IU06,*) ' ******************************************'
          CALL ABORT1
        ENDIF

!  LOOP THROUGH LONGITUDES

        NL=0
        XLOH = XLOW + (REAL(ILW)-1.5000)*XDELO+720.
        XLO=AMOWEP
        IF(XLO .LT.0.) XLO =XLO +360.
        DO I=1,NLON
          XLO = MOD(XLO+720.,360.)
          XLOG = MOD(XLOH + REAL(I)*XDELO,360.)
 3030     CONTINUE
          IF (XLO.LT.XLOG) XLO  = XLO  + 360.

          IF (XLO.LE.XLOG+XDELO) THEN
            NL=NL+1
            IF (NL.GT.NGX) THEN
              WRITE (IU06,*) ' **************************************'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR   *'
              WRITE (IU06,*) ' *      ===========================   *'
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' * NUMBER OF LONG. IN OUTPUT GRID     *'
              WRITE (IU06,*) ' * EXCEEDS DIMENSION.                 *'
              WRITE (IU06,*) ' * DIMENSION IS            NGX = ', NGX
              WRITE (IU06,*) ' * LAST LONGITUDE READ IS  XLO = ', XLO
              WRITE (IU06,*) ' * LONGITUDE INCREMENT IS  ', ZDELLO(J)
              WRITE (IU06,*) ' * INDEX J  IS  ', J
              WRITE (IU06,*) ' * INDEX JRGG  IS  ', JRGG
              WRITE (IU06,*) ' * INDEX I  IS  ', I
              WRITE (IU06,*) ' * NLAT IS  ', NLAT
              WRITE (IU06,*) ' *                                    *'
              WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED            *'
              WRITE (IU06,*) ' **************************************'
              CALL ABORT1
            ENDIF
            BATHY(NL,NJ)=XA1(I,J)
            JRGG = 1 +NINT((XLA-AMOSOP)/XDELLA)
            XLO = AMOWEP + REAL(NL)*ZDELLO(JRGG)
            KXLO=NINT(XLO*100.)
            KAMOEAP=NINT(AMOEAP*100.)
            IF (KXLO.LE.KAMOEAP) THEN
              XLO = MOD(XLO+720.,360.)
              GOTO 3030
            ENDIF
          ENDIF
        ENDDO

        JJ=JJ+1
        XLA=AMOSOP+REAL(JJ)*XDELLA

        IF((XLA-0.5*XDELLA).GT.AMONOP) THEN
          GOTO 3080
        ELSE
          GOTO 3010
        ENDIF
 3070   continue
      ENDDO
 3080 CONTINUE
      IF (NJ.NE.NY .OR. NL.GT.NX) THEN
        WRITE (IU06,*) ' *****************************************'
        WRITE (IU06,*) ' *                                       *'
        WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR      *'
        WRITE (IU06,*) ' *      ===========================      *'
        WRITE (IU06,*) ' *                                       *'
        WRITE (IU06,*) ' * NUMBER OF LONGITUDES OR LATITUDES IN  *'
        WRITE (IU06,*) ' * IS NOT EQUAL TO EXPECTED NUMBER       *'
        WRITE (IU06,*) ' * LATITUDES  FOUND      NJ = ', NJ
        WRITE (IU06,*) ' * LATITUDES  EXPECTED   NY = ', NY
        WRITE (IU06,*) ' * LONGITUDES FOUND      NL = ', NL
        WRITE (IU06,*) ' * LONGITUDES EXPECTED   NX = ', NX
        WRITE (IU06,*) ' *                                       *'
        WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED               *'
        WRITE (IU06,*) ' *****************************************'
        CALL ABORT1
      ENDIF

      ENDIF  ! end on new vs old type bathymetry


! ----------------------------------------------------------------------

!*    4. CONVERT INPUT DEPTH TO MODEL DEPTH
!*       POSITIVE SEA DEPTH IN METRES (-999  FOR LAND).
!        ----------------------------------------------

      DO J=1,NY
        DO I=1,NX
          IF (BATHY(I,J).LT.0.) THEN
            BATHY(I,J) = -BATHY(I,J)
          ELSE
            BATHY(I,J) = -999.
          ENDIF
        ENDDO
      ENDDO

!     CHECK THAT MINIMUM DEPTH USED IN TABLES IS MET. 
      NMINADJT=0
      DO J=1,NY
        DO I=1,NX
          IF (BATHY(I,J).GT.0. .AND. BATHY(I,J).LT. DEPTHA ) THEN
            BATHY(I,J) = DEPTHA 
            NMINADJT=NMINADJT+1
          ENDIF
        ENDDO
      ENDDO
      IF(NMINADJT.GT.0) THEN
        WRITE (IU06,*) ' ' 
        WRITE (IU06,*) ' *******************************************'
        WRITE (IU06,*) ' *                                         *'
        WRITE (IU06,*) ' *      WARNING IN SUB. TOPOAR             *'
        WRITE (IU06,*) ' *      ============================       *'
        WRITE (IU06,*) ' *                                         *'
        WRITE (IU06,*) ' * THE DEPTH AT SOME GRID POINTS WAS RESET *' 
        WRITE (IU06,*) ' * TO THE MINIMUM DEPTH IN TABLES ',DEPTHA
        WRITE (IU06,*) ' * NUMBER OF AFFECTED GRID POINTS: ',NMINADJT
        WRITE (IU06,*) ' *                                         *'
        WRITE (IU06,*) ' *******************************************'
        WRITE (IU06,*) ' ' 
      ENDIF

!     CHECK THAT THE MAXIMUM DEPTH IN TABLES IS SUFFICIENTLY LARGE
      BATHYMAX=999.
      DO J=1,NY
        DO I=1,NX
          BATHYMAX=MIN(BATHY(I,J),BATHYMAX)
        ENDDO
      ENDDO
      TABLEMAX=DEPTHA*DEPTHD**(NDEPTH-1)
      IF(BATHYMAX.GT.TABLEMAX) THEN
        WRITE (IU06,*) ' ******************************************'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *      WARNING ERROR IN SUB. TOPOAR      *'
        WRITE (IU06,*) ' *      ============================      *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *  THE MAXIMUM DEPTH ',BATHYMAX 
        WRITE (IU06,*) ' *  IS LARGER THAN '
        WRITE (IU06,*) ' *  THE MAXIMUM DEPTH IN TABLES ',TABLEMAX 
        WRITE (IU06,*) ' *  ADJUST DEPTHA, DEPTHD, NDEPTH !       *'
        WRITE (IU06,*) ' *  (SEE INPUT NAMELIST)                  *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' ******************************************'
        CALL ABORT1
      ENDIF


! ----------------------------------------------------------------------

!*    5. MANUAL ADJUSTMENT OF TOPOGRAPHY.
!        --------------------------------

      IF (NOUT.NE.0) THEN
        XLAT=AMOSOP-XDELLA
        DO J=1,NY
          XLAT=XLAT+XDELLA
          XLON=AMOWEP-ZDELLO(J)
          IF (XLON.LT.0.) XLON=360.+XLON
          DO I=1,NX
            XLON=XLON+ZDELLO(J)
            IF (XLON.GE.360.) XLON=XLON-360.
            DO JH = 1,NOUT
              IF (XLON.LT.XOUTW(JH)) XLON=XLON+360.
              IF (XLON.GT.XOUTE(JH)) XLON=XLON-360.
              IF (XLON.GE.XOUTW(JH) .AND. XLAT.GE.XOUTS(JH) .AND.
     &         XLON.LE.XOUTE(JH) .AND. XLAT.LE.XOUTN(JH))
     &         BATHY(I,J) = REAL(NOUTD(JH))
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    7. AID TO USERS - SIMPLE PLOT OF GRID.
!        ------------------------------------

      WRITE (IU06,'(''0NUMBER OF LATITUDES IS        NY = '',I5)') NY
      WRITE (IU06,'('' MOST SOUTHERN LATITUDE IS AMOSOP = '',F7.3)')
     &           AMOSOP
      WRITE (IU06,'('' MOST NORTHERN LATITUDE IS AMONOP = '',F7.3)')
     &           AMONOP
      WRITE (IU06,'('' LATITUDE INCREMENT IS     XDELLA = '',F7.3)')
     &           XDELLA
      WRITE (IU06,'(''0MAX NUMBER OF LONGITUDES IS   NX = '',I5)') NX
      WRITE (IU06,'('' MOST WESTERN LONGITUDE IS AMOWEP = '',F7.3)')
     &           AMOWEP
      WRITE (IU06,'('' MOST EASTERN LONGITUDE IS AMOEAP = '',F7.3)')
     &           AMOEAP
      WRITE (IU06,
     &     '('' LONGITUDE INCREMENT AS FUNCTION OF LATITUDE IS'')')
      WRITE (IU06,'(10F8.3)') ZDELLO

      IF (ITEST.GE.3) THEN
        ILEN = 120
        IPAGE = (NX+ILEN-1)/ILEN
        IF (IPAGE.GT.1) THEN
          LAST = (NX-ILEN*(IPAGE-1)+IPAGE-2)/(IPAGE-1)
          IF (LAST.LE.10) THEN
            ILEN = ILEN + 10
            IPAGE = (NX+ILEN-1)/ILEN
          ENDIF
        ENDIF

        ALLOCATE(AXX(NGX))

        DO L=1,IPAGE
          IA = (L-1)*ILEN
          IE = MIN(IA+ILEN,NX)
          IA = IA+1
          WRITE (IU06,'(''0UNBLOCKED GRID               N'',
     &     40X,''PAGE'',I2)') L
          WRITE (IU06,'(''   L = LAND               W -   - E'')')
          WRITE (IU06,'(''   S = SEA                    S'',/)')
          WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
          DO JH =NY,1,-1
            DO I=IA,IE
              IF (BATHY(I,JH).EQ.-999) AXX(I)='L'
              IF (BATHY(I,JH).NE.-999) AXX(I)='S'
            ENDDO
            WRITE (IU06,'(1X,I1,130A1)') MOD(JH,10),(AXX(I),I=IA,IE)
            WRITE (54,'(1X,I1,130A1)') MOD(JH,10),(AXX(I),I=IA,IE)
          ENDDO
          WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
        ENDDO
      ENDIF

      IF(ALLOCATED(AX)) DEALLOCATE(AX)
      IF(ALLOCATED(AXX)) DEALLOCATE(AXX)
      IF(ALLOCATED(IA2H)) DEALLOCATE(IA2H)
      IF(ALLOCATED(XA2H)) DEALLOCATE(XA2H)
      IF(ALLOCATED(IA1)) DEALLOCATE(IA1)
      IF(ALLOCATED(XA1)) DEALLOCATE(XA1)

      RETURN
      END SUBROUTINE TOPOAR
