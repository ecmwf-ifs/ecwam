      SUBROUTINE TOPOAR (IU01, IA2)

C ----------------------------------------------------------------------
C
C**** *TOPOAR* - ARRANGE SUBGRID TOPOGRAPHY.
C
C     S. HASSELMANN     MPIFM           1/6/86.
C
C     MODIFIED BY       H. GUNTHER      1/4/90  -  REARANGEMENT OF CODE.
C
C*    PURPOSE.
C     --------
C
C       TO READ IN TOPOGRAPHY ON INPUT GRID AND CONVERT TO OUTPUT GRID.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *TOPOAR (IU01, IA2)*
C          *IU01*  -  LOGICAL INPUT UNIT OF TOPOGRAPHIC DATA.
C          *IA2*   -  TOPOGRAPHIC DATA IN OUTPUT GRID.
C
C     METHOD.
C     -------
C
C       THE TOPOGRAPHY MUST BE ON A REGULAR LATITUDE-LONGITUDE
C       GRID ARRANGED FROM SOUTH TO NORTH, AND FROM WEST TO EAST.
C       IT IS ASSUMED THAT NEGATIVE VALUES ARE SEA DEPTHS (WHICH
C       ARE CONVERTED TO POSITIVE) AND THAT POSITIVE VALUES ARE
C       LAND ELEVATIONS (WHICH ARE CONVERTED TO -999 IDENTIFIERS).
C
C       THE TOPOGRAPHIC DATA IS READ IN ON THE INPUT GRID AND IT
C       IS STORED ONLY FOR THOSE LATITUDES WITHIN THE REQUESTED
C       GRID. THEN THE TOPOGRAPHIC DATA IS FURTHER RESTRICTED
C       TO LIE WITHIN THE SUBGRID LONGITUDES. IT IS THEN PUT ON
C       THE REQUESTED SUBGRID LAT-LONG RESOLUTION, ALWAYS USING
C       THE NEAREST POINT.
C       FINALLY THE SUBGRID TOPOGRAPHY MAY BE MANUALLY ADJUSTED BY
C       MEANS OF THE CARD INPUT AND A PRINTER OUTPUT IS DONE.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT1*     - TERMINATES PROCESSING.
C       *ADJUST*    - CORRECTS LONGITUDE INPUT.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "parameter_topo.h"
C
#include "comcinp.h"
C
#include "commap.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION IA2(NGX, NGY)
      CHARACTER*1 AX(JLONI), AXX(NGX)
      DIMENSION IA2H(JLONI), IA1(JLONI,JLATI)
CAG DEB
      CHARACTER LANDSEA(130)*1
      LOGICAL   LOFORT55
CAG FIN
C
C ----------------------------------------------------------------------
C
C*    1. READING THE TOPOGRAPHY OF THE INPUT GRID AND STORING THOSE
C*       LATITUDES WITHIN THE OUTPUT SUBGRID AREA.
C        -----------------------------------------
C
C       MLON - NUMBER OF TOPOGRAPHIC POINTS PER GRID LATITUDE.
C       KMAX - NUMBER OF RECORDS OF INPUT GRID PER LATITUDE.
C       XLAT - LATITUDE OF CURRENT GRID DATA.
C       NLAT - NUMBER OF GRID LATITUDES STORED.
C       XLAG - LATITUDE OF FIRST GRID LATITUDE STORED.
C
 1000 CONTINUE
      REWIND (UNIT=IU01)
      READ (IU01,'(8F10.5)') XDELA, XDELO, XLAS, XLAN, XLOW, XLOE
      CALL ADJUST (XLOW, XLOE)

      WRITE (IU06,'(1H1,'' INPUT GRID''/)')
      WRITE (IU06,'(3X,''RESOLUTION LAT-LON '',2F8.3)') XDELA, XDELO
      WRITE (IU06,'(3X,'' SOUTHERN LAT '','' NORTHERN LAT '',
     1                 '' WESTERN LONG '','' EASTERN LONG'',
     2                 /,2X,4F14.3)') XLAS, XLAN, XLOW, XLOE

      MLON = NINT((XLOE-XLOW)/XDELO+1.0)
      IF(MLON.GT.JLONI) THEN
         WRITE (IU06,*) ' *********************************************'
         WRITE (IU06,*) ' *                                           *'
         WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR          *'
         WRITE (IU06,*) ' *      ===========================          *'
         WRITE (IU06,*) ' *                                           *'
         WRITE (IU06,*) ' * NUMBER OF LONGITUDES IN INPUT GRID EXCEEDS*'
         WRITE (IU06,*) ' * DIMENSION                                 *'
         WRITE (IU06,*) ' * DIMENSION IS        JLONI = ', JLONI
         WRITE (IU06,*) ' * NUMBER FROM INPUT IS MLON = ', MLON
         WRITE (IU06,*) ' *                                           *'
         WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                   *'
         WRITE (IU06,*) ' *********************************************'
         CALL ABORT1
      ENDIF
      KMAX=(MLON+11)/12
C
      XLAT=XLAS-XDELA
      NLAT=0
1005  CONTINUE
         NLAT=NLAT+1
         IF (NLAT.GT.JLATI) THEN
            WRITE (IU06,*) ' ******************************************'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. TOPOAR       *'
            WRITE (IU06,*) ' *      ===========================       *'
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' * NUMBER OF LATITUDES IN INPUT GRID      *'
            WRITE (IU06,*) ' * EXCEEDS DIMENSION                      *'
            WRITE (IU06,*) ' * DIMENSION IS           JLATI = ', JLATI
            WRITE (IU06,*) ' * LAST LATITUDE READ IS  XLAT  = ', XLAT
            WRITE (IU06,*) ' *                                        *'
            WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                *'
            WRITE (IU06,*) ' ******************************************'
            CALL ABORT1
         ENDIF
1010     CONTINUE
            XLAT=XLAT+XDELA
            DO 1015 K=1,KMAX
               N1 = 12*(K-1)+1
               N2 = MIN(12*K,MLON)
               READ (IU01, '(12(I5,A1))', END=1017)
     1              (IA1(IAA,NLAT),AX(IAA),IAA=N1,N2)
 1015       CONTINUE
C
         IF (XLAT+0.5*XDELA.LE.AMOSOP) GO TO 1010
         IF(NLAT.EQ.1) XLAG=XLAT
         DO 1016 I=1,MLON
             IF (AX(I).EQ.'E'.AND.IA1(I,NLAT).LT.0)
     1                        IA1(I,NLAT)=-IA1(I,NLAT)
             IF (AX(I).EQ.'D'.AND.IA1(I,NLAT).GT.0)
     1                        IA1(I,NLAT)=-IA1(I,NLAT)
             IF (AX(I).EQ.'E'.AND.IA1(I,NLAT).EQ.0) IA1(I,NLAT)=1
 1016    CONTINUE
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
C
C ----------------------------------------------------------------------
C
C*    2. STORING TOPOGRAPHIC DATA AT LONGITUDES WITHIN SUBGRID AREA.
C        ----------------------------------------------------------
C
C       ILW      -  INDEX OF NEAREST GRID LONGITUDE EQUIVALENT TO
C                     WESTERN SUBGRID BOUNDARY.
C       NLON     -  NUMBER OF GRID LONGITUDES WITHIN SUBGRID AREA.
C       IH,IH1   -  GRID LONGITUDE NUMBER.
C
 2000 CONTINUE

      XLW= MOD(AMOWEP-XLOW+720.,360.)
      ILW= NINT(XLW/XDELO-0.0001)
      NLON=NINT((AMOEAP-AMOWEP)/XDELO+1.0)
      DO 2030 J=1,NLAT
         IH=ILW
         DO 2010 I=1,NLON
            IH=IH+1
            IF(IH.LE.0) IH=IH+MLON
            IF(IH.GT.MLON) IH=IH-MLON
            IA2H(I) =IA1(IH,J)
 2010    CONTINUE
         DO 2020 I=1,NLON
            IA1(I,J)=IA2H(I)
 2020    CONTINUE
 2030 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. PUT ON XDELLO BY XDELLA SUBGRID.
C        --------------------------------
C
C     XLA  - LATITUDE OF OUTPUT GRID.
C     XLAG - LATITUDE OF INPUT GRID.
C     XLO  - LONGITUDE OF OUTPUT GRID.
C     XLOG - LONGITUDE OF INPUT GRID.
C
 3000 CONTINUE
      XLA=AMOSOP
      XLAG=XLAG-XDELA
      NJ=0
      IA2 = 999

C  LOOP THROUGH LATITUDES

      DO 3070 J=1,NLAT
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

C  LOOP THROUGH LONGITUDES

        NL=0
        XLOH = XLOW + (REAL(ILW)-1.5000)*XDELO+720.
        XLO=AMOWEP
        IF(XLO .LT.0.) XLO =XLO +360.
        DO 3060 I=1,NLON
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
 3060   CONTINUE
        XLA=XLA+XDELLA
        IF(XLA.GT.AMONOP) THEN
          GOTO 3080
        ELSE
          GOTO 3010
        ENDIF
 3070 CONTINUE
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
C
C ----------------------------------------------------------------------
C
C*    4. CONVERT INPUT DEPTH TO MODEL DEPTH
C*       POSITIVE SEA DEPTH IN METRES (-999  FOR LAND).
C        ----------------------------------------------
C
 4000 CONTINUE
      DO J=1,NY
        DO I=1,NX
C
CPETER CHANGED LE TO LT
C      
          IF (IA2(I,J).LT.0) THEN
            IA2(I,J) = -IA2(I,J)
          ELSE
            IA2(I,J) = -999
          ENDIF
        ENDDO
      ENDDO
C
C ----------------------------------------------------------------------
C
C*    5. MANUAL ADJUSTMENT OF TOPOGRAPHY.
C        --------------------------------
C
 5000 CONTINUE
      IF (NOUT.NE.0) THEN
         XLAT=AMOSOP-XDELLA
         DO 5001 J=1,NY
            XLAT=XLAT+XDELLA
            XLON=AMOWEP-ZDELLO(J)
            IF (XLON.LT.0.) XLON=360.+XLON
            DO 5002 I=1,NX
               XLON=XLON+ZDELLO(J)
               IF (XLON.GE.360.) XLON=XLON-360.
               DO 5003 JH = 1,NOUT
                  IF (XLON.LT.XOUTW(JH)) XLON=XLON+360.
                  IF (XLON.GT.XOUTE(JH)) XLON=XLON-360.
                  IF (XLON.GE.XOUTW(JH) .AND. XLAT.GE.XOUTS(JH) .AND.
     1                XLON.LE.XOUTE(JH) .AND. XLAT.LE.XOUTN(JH))
     2               IA2(I,J) = NOUTD(JH)
 5003          CONTINUE
 5002       CONTINUE
 5001    CONTINUE
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    7. AID TO USERS - SIMPLE PLOT OF GRID.
C        ------------------------------------
C
 7000 CONTINUE
      WRITE (IU06,'(''0NUMBER OF LATITUDES IS        NY = '',I5)') NY
      WRITE (IU06,'('' MOST SOUTHERN LATITUDE IS AMOSOP = '',F5.1)')
     1           AMOSOP
      WRITE (IU06,'('' MOST NORTHERN LATITUDE IS AMONOP = '',F5.1)')
     1           AMONOP
      WRITE (IU06,'('' LATITUDE INCREMENT IS     XDELLA = '',F5.1)')
     1           XDELLA
      WRITE (IU06,'(''0MAX NUMBER OF LONGITUDES IS   NX = '',I5)') NX
      WRITE (IU06,'('' MOST WESTERN LONGITUDE IS AMOWEP = '',F5.1)')
     1           AMOWEP
      WRITE (IU06,'('' MOST EASTERN LONGITUDE IS AMOEAP = '',F5.1)')
     1           AMOEAP
      WRITE (IU06,
     1     '('' LONGITUDE INCREMENT AS FUNCTION OF LATITUDE IS'')')
      WRITE (IU06,'(10F8.3)') ZDELLO

CAG DEB
      INQUIRE ( FILE='fort.55',EXIST=LOFORT55)
      OPEN (FILE='fort.55.hlp', UNIT=54, FORM='formatted')
CAG END
      ILEN = 120
      IPAGE = (NX+ILEN-1)/ILEN
      IF (IPAGE.GT.1) THEN
         LAST = (NX-ILEN*(IPAGE-1)+IPAGE-2)/(IPAGE-1)
         IF (LAST.LE.10) THEN
            ILEN = ILEN + 10
            IPAGE = (NX+ILEN-1)/ILEN
         ENDIF
      ENDIF
      DO 7003 L=1,IPAGE
         IA = (L-1)*ILEN
         IE = MIN(IA+ILEN,NX)
         IA = IA+1
         WRITE (IU06,'(''0UNBLOCKED GRID               N'',
     1                40X,''PAGE'',I2)') L
         WRITE (IU06,'(''   L = LAND               W -   - E'')')
         WRITE (IU06,'(''   S = SEA                    S'',/)')
         WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
         DO 7001 JH =NY,1,-1
CAG DEB
            IF ( LOFORT55) THEN
              IALS=IA-(L-1)*ILEN
              IELS=IE-(L-1)*ILEN
              READ(55,6999) (LANDSEA(I),I=IALS,IELS)
CAG           WRITE(IU06,6999) (LANDSEA(I),I=IALS,IELS)
6999          FORMAT(1X,130A1)
            ENDIF
CAG END
            DO 7002 I=IA,IE
CAG DEB
               IF( LOFORT55 .AND. LANDSEA(I-(L-1)*ILEN).EQ.'L')
     .           IA2(I,JH)=-999
CAG END
               IF (IA2(I,JH).EQ.-999) AXX(I)='L'
               IF (IA2(I,JH).NE.-999) AXX(I)='S'
 7002       CONTINUE
            WRITE (IU06,'(1X,I1,130A1)') MOD(JH,10),(AXX(I),I=IA,IE)
            WRITE (54,'(1X,I1,130A1)') MOD(JH,10),(AXX(I),I=IA,IE)
 7001    CONTINUE
         WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
 7003 CONTINUE

      RETURN
      END
