      SUBROUTINE TOPOAR (IU01, IA2)

! ----------------------------------------------------------------------

!**** *TOPOAR* - ARRANGE SUBGRID TOPOGRAPHY.

!     S. HASSELMANN     MPIFM           1/6/86.

!     MODIFIED BY       H. GUNTHER      1/4/90  -  REARANGEMENT OF CODE.

!*    PURPOSE.
!     --------

!       TO READ IN TOPOGRAPHY ON INPUT GRID AND CONVERT TO OUTPUT GRID.

!**   INTERFACE.
!     ----------

!       *CALL* *TOPOAR (IU01, IA2)*
!          *IU01*  -  LOGICAL INPUT UNIT OF TOPOGRAPHIC DATA.
!          *IA2*   -  TOPOGRAPHIC DATA IN OUTPUT GRID.

!     METHOD.
!     -------

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

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *ADJUST*    - CORRECTS LONGITUDE INPUT.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWCINP  , ONLY : NOUT     ,XOUTW    ,XOUTS    ,XOUTE    ,
     &            XOUTN    ,NOUTD
      USE YOWMAP   , ONLY : NX       ,NY       ,AMOWEP   ,AMOSOP   ,
     &            AMOEAP   ,AMONOP   ,XDELLA   ,ZDELLO
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      DIMENSION IA2(NGX, NGY)

      CHARACTER*1, ALLOCATABLE :: AX(:), AXX(:)
      INTEGER, ALLOCATABLE :: IA2H(:), IA1(:,:)

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
      READ (IU01,'(8F10.5)') XDELA, XDELO, XLAS, XLAN, XLOW, XLOE
      CALL ADJUST (XLOW, XLOE)

      WRITE (IU06,'(1H1,'' INPUT GRID''/)')
      WRITE (IU06,'(3X,''RESOLUTION LAT-LON '',2F8.3)') XDELA, XDELO
      WRITE (IU06,'(3X,'' SOUTHERN LAT '','' NORTHERN LAT '',
     &                 '' WESTERN LONG '','' EASTERN LONG'',
     &                 /,2X,4F14.3)') XLAS, XLAN, XLOW, XLOE

      MLON = NINT((XLOE-XLOW)/XDELO+1.0)

      NLATMAX=NINT((XLAN-XLAS)/XDELA+1.0)

      ALLOCATE(AX(MLON), AXX(NGX))
      ALLOCATE(IA2H(MLON), IA1(MLON,NLATMAX))

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
        READ (IU01, '(12(I5,A1))', END=1017)
     &   (IA1(IAA,NLAT),AX(IAA),IAA=N1,N2)
      ENDDO

      IF (XLAT+0.5*XDELA.LE.AMOSOP) GO TO 1010
      IF(NLAT.EQ.1) XLAG=XLAT
      DO I=1,MLON
        IF (AX(I).EQ.'E'.AND.IA1(I,NLAT).LT.0)
     &   IA1(I,NLAT)=-IA1(I,NLAT)
        IF (AX(I).EQ.'D'.AND.IA1(I,NLAT).GT.0)
     &   IA1(I,NLAT)=-IA1(I,NLAT)
        IF (AX(I).EQ.'E'.AND.IA1(I,NLAT).EQ.0) IA1(I,NLAT)=1
      ENDDO
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
      ILW= NINT(XLW/XDELO-0.0001)
      NLON=NINT((AMOEAP-AMOWEP)/XDELO+1.0)
      DO J=1,NLAT
        IH=ILW
        DO I=1,NLON
          IH=IH+1
          IF(IH.LE.0) IH=IH+MLON
          IF(IH.GT.MLON) IH=IH-MLON
          IA2H(I) =IA1(IH,J)
        ENDDO
        DO I=1,NLON
          IA1(I,J)=IA2H(I)
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
      IA2 = 999

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
          WRITE (IU06,*) ' * LAST LATITUDE USED IS   XLA = ', XLA
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
          IF (XLO.LE.XLOG) XLO  = XLO  + 360.
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
            IA2(NL,NJ)=IA1(I,J)
            JRGG = 1 +NINT((XLA-AMOSOP)/XDELLA)
            XLO = AMOWEP + REAL(NL)*ZDELLO(JRGG)
            IF (XLO.LE.AMOEAP) THEN
              XLO = MOD(XLO+720.,360.)
              GOTO 3030
            ENDIF
          ENDIF
        ENDDO
        XLA=XLA+XDELLA

        IF((XLA-0.5*XDELA).GT.AMONOP) THEN
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

! ----------------------------------------------------------------------

!*    4. CONVERT INPUT DEPTH TO MODEL DEPTH
!*       POSITIVE SEA DEPTH IN METRES (-999  FOR LAND).
!        ----------------------------------------------

      DO J=1,NY
        DO I=1,NX

!PETER CHANGED LE TO LT

          IF (IA2(I,J).LT.0) THEN
            IA2(I,J) = -IA2(I,J)
          ELSE
            IA2(I,J) = -999
          ENDIF
        ENDDO
      ENDDO

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
     &         IA2(I,J) = NOUTD(JH)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    7. AID TO USERS - SIMPLE PLOT OF GRID.
!        ------------------------------------

      WRITE (IU06,'(''0NUMBER OF LATITUDES IS        NY = '',I5)') NY
      WRITE (IU06,'('' MOST SOUTHERN LATITUDE IS AMOSOP = '',F5.1)')
     &           AMOSOP
      WRITE (IU06,'('' MOST NORTHERN LATITUDE IS AMONOP = '',F5.1)')
     &           AMONOP
      WRITE (IU06,'('' LATITUDE INCREMENT IS     XDELLA = '',F5.1)')
     &           XDELLA
      WRITE (IU06,'(''0MAX NUMBER OF LONGITUDES IS   NX = '',I5)') NX
      WRITE (IU06,'('' MOST WESTERN LONGITUDE IS AMOWEP = '',F5.1)')
     &           AMOWEP
      WRITE (IU06,'('' MOST EASTERN LONGITUDE IS AMOEAP = '',F5.1)')
     &           AMOEAP
      WRITE (IU06,
     &     '('' LONGITUDE INCREMENT AS FUNCTION OF LATITUDE IS'')')
      WRITE (IU06,'(10F8.3)') ZDELLO

      ILEN = 120
      IPAGE = (NX+ILEN-1)/ILEN
      IF (IPAGE.GT.1) THEN
        LAST = (NX-ILEN*(IPAGE-1)+IPAGE-2)/(IPAGE-1)
        IF (LAST.LE.10) THEN
          ILEN = ILEN + 10
          IPAGE = (NX+ILEN-1)/ILEN
        ENDIF
      ENDIF
      DO L=1,IPAGE
        IA = (L-1)*ILEN
        IE = MIN(IA+ILEN,NX)
        IA = IA+1
        WRITE (IU06,'(''0UNBLOCKED GRID               N'',
     &   40X,''PAGE'',I2)') L
        WRITE (IU06,'(''   L = LAND               W -   - E'')')
        WRITE (IU06,'(''   S = SEA                    S'',/)')
        WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
        DO JH =NY,1,-1
          DO I=IA,IE
            IF (IA2(I,JH).EQ.-999) AXX(I)='L'
            IF (IA2(I,JH).NE.-999) AXX(I)='S'
          ENDDO
          WRITE (IU06,'(1X,I1,130A1)') MOD(JH,10),(AXX(I),I=IA,IE)
          WRITE (54,'(1X,I1,130A1)') MOD(JH,10),(AXX(I),I=IA,IE)
        ENDDO
        WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
      ENDDO

      DEALLOCATE(AX, AXX)
      DEALLOCATE(IA2H, IA1)

      RETURN
      END SUBROUTINE TOPOAR
