C ======================================================================
C
      SUBROUTINE READPRE (kuso, kugrid, IREFRA,
     .            NANG, NFRE, NGY, NBLO, NIBLO, MOUTP, NIBLC, NBLC,
     .            NDEPTH, ITAUMAX, JUMAX, IUSTAR, IALPHA,
     .            FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
     .            DELPHI, DELLAM, NLONRGG, SINPH, COSPH, IGL, IJS, IJL2,
     .            IJLS, IJL, IJLT, IXLG, KXLT, NX, NY, IPER,
     .            AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,ZDELLO,
     .            IRGG,IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, AF11,
     .            FKLAP, FKLAP1, FKLAM, FKLAM1, ACL1, ACL2,  CL11,
     .            CL21, DAL1, DAL2, FRH,
     .            BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV, TAUT, DELTAUW,
     .            DELU, TAUHFT, DELUST, DELALP, NGOUT, IGAR, IJAR,
     .            DEPTH, DEPTHA, DEPTHD, TCGOND, TFAK, TSIHKD, U, V)
 
C ----------------------------------------------------------------------
C
C**** *READPRE*  READ GRID OUTPUT FROM PREPROC.
C
C     H. GUNTHER      GKSS/ECMWF     MAY 1990
C     J. BIDLOT       ECMWF          SEPTEMBER 96: read gaussian grid
C
C*    PURPOSE.
C     --------
C
C       INPUT OF PREPROC GRID OUTPUT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *READPRE (kugrid, IREFRA)*
C      I\    *kugrid*      INTEGER   INPUT UNIT OF PREPROC GRID FILE.
C      I\    *IREFRA*    INTEGER   REFRACTION OPTION.
C      I\    *NANG*      INTEGER   NUMBER OF ANGLES.
C      I\    *NFRE*      INTEGER   NUMBER OF FREQUENCIES.
C      I\    *NGY*       INTEGER   NUMBER OF LATITUDES IN GRID.
C      I\    *NBLO*      INTEGER   NUMBER OF BLOCKS.
C      I\    *NIBLO*     INTEGER   NUMBER OF POINTS IN BLOCK.
C      I\    *MOUTP*     INTEGER   MAXIMUM NUMBER OF OUTPUT POINTS.
C      I\    *NIBLC*     INTEGER   = NIBLO IF CURRENT REFRACTION.
C      I\                          = 1     ELSE.
C      I\    *NBLC*      INTEGER   = NIBO  IF CURRENT REFRACTION.
C      I\                          = 1     ELSE.
C      I\    *NDEPTH*    INTEGER   LENGTH OF SHALLOW WATER TABLES.
C      I\    *ITAUMAX*   INTEGER   TABLE DIMENSION.
C      I\    *JUMAX*     INTEGER   TABLE DIMENSION.
C      I\    *IUSTAR*    INTEGER   TABLE DIMENSION.
C      I\    *IALPHA*    INTEGER   TABLE DIMENSION.
C       \O   *FR*        REAL      FREQUENCIES IN HERTZ.
C       \O   *DFIM*      REAL      FREQUENCY INTERVAL*DIRECTION
C                                  INTERVAL.
C       \O   *GOM*       REAL      DEEP WATER GROUP VELOCITIES (M/S).
C       \O   *C*         REAL      DEEP WATER PHASE VELOCITIES (M/S).
C       \O   *DELTH*     REAL      ANGULAR INCREMENT OF SPECTRUM
C                                  (RADIANS).
C       \O   *DELTR*     REAL      DELTH TIMES RADIUS OF EARTH (METRES).
C       \O   *TH*        REAL      DIRECTIONS IN RADIANS.
C       \O   *COSTH*     REAL      COS OF DIRECTION.
C       \O   *SINTH*     REAL      SIN OF DIRECTION.
C       \O   *DELPHI*    REAL      GRID INCREMENT FOR LATITUDE (METRES).
C       \O   *DELLAM*    REAL      GRID INCREMENT FOR LONGITUDE AT
C                                  EQUATOR IN METRES.
C        IT'S DEPENDENT ON THE LATITUDE NOW
C             NLONRGG    INTEGER   NUMBER OF POINT PER LATITUDE
C       \O   *SINPH*     REAL      SIN OF LATITUDE.
C       \O   *COSPH*     REAL      COS OF LATITUDE.
C       \O   *IGL*       INTEGER   NO OF BLOCKS.
C       \O   *IJS*       INTEGER   INDEX OF FIRST POINT OF SECOND LAT.
C       \O   *IJL2*      INTEGER   INDEX OF LAST POINT OF SECOND LAT.
C       \O   *IJLS*      INTEGER   INDEX OF FIRST POINT OF LAT BEFORE
C                                  LAST.
C       \O   *IJL*       INTEGER   INDEX OF LAST POINT OF LAT BEFORE
C                                  LAST.
C       \O   *IJLT*      INTEGER   TOTAL NUMBER OF GRIDPOINTS IN A
C                                  BLOCK.
C       \O   *IXLG*      INTEGER   LONG. GRID INDEX.
C       \O   *KXLT*      INTEGER   LAT. GRID INDEX.
C       \O   *NX*        INTEGER   NUMBER OF LONGITUDES IN GRID.
C       \O   *NY*        INTEGER   NUMBER OF LATITUDES  IN GRID.
C       \O   *IPER*      INTEGER   = 1 IF GRID IS PERIODIC.
C       \O   *AMOWEP*    REAL      MOST WESTERN LONGITUDE IN GRID
C                                  (DEGREE).
C       \O   *AMOSOP*    REAL      MOST SOUTHERN LATITUDE IN GRID
C                                  (DEGREE).
C       \O   *AMOEAP*    REAL      MOST EASTERN LONGITUDE IN GRID
C                                  (DEGREE).
C       \O   *AMONOP*    REAL      MOST NORTHERN LATITUDE IN GRID
C                                  (DEGREE).
C       \O   *XDELLA*    REAL      GRID INCREMENT FOR LATITUDE (DEGREE).
C       \O   *ZDELLO*    REAL      GRID INCREMENT FOR LONGITUDE 
C                                  (DEGREE).
C             IT'S A FUNCTION OF LATITUDE (TO ALLOW FOR GAUSSIAN GRID)
C            *IRGG       INTEGER   GAUSSIAN GRID INDICATOR FLAG
C       \O   *IKP*       INTEGER   FREQUENCY INDEX ARRAY FOR STORING
C                                  ENERGY TRANSFER INCREMENTS INTO
C                                  BINS, WAVE No. 3.
C       \O   *IKP1*      INTEGER   IKP+1.
C       \O   *IKM*       INTEGER   FREQUENCY INDEX ARRAY FOR STORING
C                                  ENERGY TRANSFER INCREMENTS INTO BINS,
C                                  WAVE NO. 4.
C       \O   *IKM1*      INTEGER   IKM+1
C       \O   *K1W*       INTEGER   ANGULAR INDEX ARRAY FOR STORING
C                                  ENERGY TRANSFER INCREMENTS INTO BINS,
C                                  WAVE NO. 3.
C       \O   *K11W*      INTEGER   K1W(.,1)-1, K1W(.,2)+1.
C       \O   *K2W*       INTEGER   ANGULAR INDEX ARRAY FOR STORING
C                                  ENERGY TRANSFER INCREMENTS INTO BINS,
C                                  WAVE NO. 4.
C       \O   *K21W*      INTEGER   K2W(.,1)+1, K2W(.,2)-1.
C       \O   *AF11*      REAL      WEIGHTS FOR DISCRETE APPROXIMATION
C                                  OF NONL TRANSFER (AT PRESENT ONE
C                                  TERM ONLY SET TO 3000). MULTIPLIED
C                                  BY FREQUENCIES **11.
C       \O   *FKLAP*     REAL      WEIGHT IN FREQUENCY GRID FOR
C                                  INTERPOLATION, WAVE NO. 3
C                                  ("1+LAMBDA" TERM).
C       \O   *FKLAP1*    REAL      1-FKLAP.
C       \O   *FKLAM*     REAL      WEIGHT IN FREQUENCY GRID FOR
C                                  INTERPOLATION, WAVE NO. 4
C                                  ("1-LAMBDA" TERM).
C       \O   *ACL1*      REAL      WEIGHT IN ANGULAR GRID FOR
C                                  INTERPOLATION, WAVE NO. 3
C                                  ("1+LAMBDA" TERM).
C       \O   *ACL2*      REAL      WEIGHT IN ANGULAR GRID FOR
C                                  INTERPOLATION, WAVE NO. 4
C                                  ("1-LAMBDA" TERM).
C       \O   *CL11*      REAL      1.-ACL1.
C       \O   *CL21*      REAL      1.-ACL2.
C       \O   *DAL1*      REAL      1./ACL1.
C       \O   *DAL2*      REAL      1./ACL2.
C       \O   *FRH*       REAL      TAIL FREQUENCY RATION **5
C       \O   *BETAMAX*   REAL      PARAMETER FOR WIND INPUT.
C       \O   *ZALP*      REAL      SHIFTS GROWTH CURVE.
C       \O   *ALPHA*     REAL      CHARNOCK CONSTANT.
C       \O   *XKAPPA*    REAL      VON KARMAN CONSTANT.
C       \O   *XNLEV*
C       \O   *TAUT*      REAL      STRESS TABLE.
C       \O   *DELTAUW*   REAL      WAVE STRESS INCREMENT.
C       \O   *DELU*      REAL      WIND INCREMENT.
C       \O   *TAUHFT*    REAL      HIGH FREQUENCY STRESS TABLE.
C       \O   *DELUST*    REAL      USTAR INCREMENT.
C       \O   *DELALP*    REAL      ALPHA INCREMENT.
C       \O   *NGOUT*     INTEGER   NUMBER OF OUTPUT POINTS.
C       \O   *IGAR*      INTEGER   BLOCK NUMBER OF OUTPUT POINTS.
C       \O   *IJAR*      INTEGER   GRIDPOINT NUMBER OF OUTPUT POINT.
C       \O   *DEPTH*     REAL      WATER DEPTH IN METRES.
C       \O   *DEPTHA*    REAL      MINIMUM DEPTH FOR TABLES (METRES).
C       \O   *DEPTHD*    REAL      DEPTH INCREMENT (METRES).
C       \O   *TCGOND*    REAL      SHALLOW WATER GROUP VELOCITY TABLE.
C       \O   *TFAK*      REAL      WAVE NUMBER TABLE.
C       \O   *TSIHKD*    REAL      TABLE FOR OMEGA/SINH(2KD).
C       \O   *U*         REAL      U - COMPONENT OF CURRENT (M/S).
C       \O   *V*         REAL      V - COMPONENT OF CURRENT (M/S).
C
C     METHOD.
C     -------
C
C       UNFORMATED READ FROM UNIT.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      DIMENSION
     .     FR(NFRE), DFIM(NFRE), GOM(NFRE), C(NFRE),
     .     TH(NANG), COSTH(NANG), SINTH(NANG),
     .     SINPH(NGY), COSPH(NGY),
     .     IJS(NBLO), IJL2(NBLO), IJLS(NBLO),
     .     IJL(NBLO), IJLT(NBLO),
     .     IXLG(NIBLO,NBLO), KXLT(NIBLO,NBLO),
     .     IKP(NFRE+4), IKP1(NFRE+4),
     .     IKM(NFRE+4), IKM1(NFRE+4),
     .     K1W(NANG,2), K2W(NANG,2), K11W(NANG,2),
     .     K21W(NANG,2), AF11(NFRE+4), FKLAP(NFRE+4),
     .     FKLAP1(NFRE+4), FKLAM(NFRE+4), FKLAM1(NFRE+4),
     .     FRH(30), IGAR(MOUTP), IJAR(MOUTP),
     .     DEPTH(NIBLO, NBLO),
     .     TCGOND(NDEPTH,NFRE), TFAK(NDEPTH,NFRE),
     .     TSIHKD(NDEPTH,NFRE),
     .     TAUT(0:ITAUMAX,0:JUMAX),
     .     TAUHFT(0:IUSTAR,0:IALPHA),
     .     U(0:NIBLC,NBLC), V(0:NIBLC,NBLC),
     .     ZDELLO(NGY),DELLAM(NGY), NLONRGG(NGY)
C
C ----------------------------------------------------------------------
C
#include <outctrl.h>
      DATA csubna /"readpre"/
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZATIONS.
C        ----------------
  100 CONTINUE
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' ',csubna,': '
      IF (lltest) THEN
        PRINT*, ' kugrid ',  kugrid
        PRINT*, ' IREFRA ',  irefra
        PRINT*, ' NANG   ',  nang
        PRINT*, ' NFRE   ',  nfre
        PRINT*, ' NGY    ',  ngy
        PRINT*, ' NBLO   ',  nblo
        PRINT*, ' NIBLO  ',  niblo
        PRINT*, ' MOUTP  ',  moutp
        PRINT*, ' NIBLC  ',  niblc
        PRINT*, ' NBLC   ',  nblc
        PRINT*, ' NDEPTH ',  ndepth
        PRINT*, ' ITAUMAX',  itaumax
        PRINT*, ' JUMAX  ',  jumax
        PRINT*, ' IUSTAR ',  iustar
        PRINT*, ' IALPHA ',  ialpha
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    1. READ COMMON FREDIR (FREQUENCY DIRECTION GRID).
C        ----------------------------------------------
C
      REWIND (kugrid)
 1000 CONTINUE
      READ(kugrid) FR, DFIM, GOM, C,
     1           DELTH, DELTR, TH, COSTH, SINTH
      IF (lltest) THEN
        PRINT*,' FR    ', fr
        PRINT*,' DFIM  ', dfim
        PRINT*,' GOM   ', gom
        PRINT*,' C     ', c
        PRINT*,' DELTH ', delth
        PRINT*,' DELTR ', deltr
        PRINT*,' TH    ', th
        PRINT*,' COSTH ', costh
        PRINT*,' SINTH ', sinth
      ENDIF
C
C*    2. READ COMMON GRIDPAR (GENERAL GRID ORGANISATION).
C        ------------------------------------------------
C
 2000 CONTINUE
      READ(kugrid) DELPHI, DELLAM, NLONRGG, SINPH, COSPH,
     1           IGL, IJS, IJL2, IJLS, IJL, IJLT
      IF (lltest) THEN
        PRINT*,' DELPHI ', delphi
        PRINT*,' DELLAM ', dellam
        PRINT*,' NLONRGG ', NLONRGG 
        PRINT*,' SINPH  ', sinph
        PRINT*,' COSPH  ', cosph
        PRINT*,' IGL    ', igl
        PRINT*,' IJS    ', ijs
        PRINT*,' IJL2   ', ijl2
        PRINT*,' IJLS   ', ijls
        PRINT*,' IJL    ', ijl
        PRINT*,' IJLT   ', ijlt
      ENDIF
C
C*    3. READ COMMON MAP (LONG. AND LAT. INDICES OF GRID POINTS).
C        --------------------------------------------------------
C
 3000 CONTINUE
      READ(kugrid) IXLG, KXLT, NX, NY, IPER,
     1           AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,ZDELLO,
     2           IRGG
      IF (lltest) THEN
        PRINT*,' IXLG   ', ixlg
        PRINT*,' KXLT   ', kxlt
        PRINT*,' NX     ', nx
        PRINT*,' NY     ', ny
        PRINT*,' IPER   ', iper
        PRINT*,' AMOWEP ', amowep
        PRINT*,' AMOSOP ', amosop
        PRINT*,' AMONOP ', amonop
        PRINT*,' AMOEAP ', amoeap
        PRINT*,' XDELLA ', xdella
        PRINT*,' XDELLO ', xdello
        PRINT*,' ZDELLO ', zdello
      ENDIF
C
C*    4. READ COMMON INDNL (NON-LINEAR INTERACTION).
C       --------------------------------------------
C
 4000 CONTINUE
      READ(kugrid) IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W,
     1           AF11, FKLAP, FKLAP1, FKLAM, FKLAM1,
     2           ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
      IF (lltest) THEN
        PRINT*,' IKP     ', ikp
        PRINT*,' IKP1    ', ikp1
        PRINT*,' IKM     ', ikm
        PRINT*,' IKM1    ', ikm1
        PRINT*,' K1W     ', k1w
        PRINT*,' K2W     ', k2w
        PRINT*,' K11W    ', k11w
        PRINT*,' K21W    ', k21w
        PRINT*,' AF11    ', af11
        PRINT*,' FKLAP   ', fklap
        PRINT*,' FKLAP1  ', fklap1
        PRINT*,' FKLAM   ', fklam
        PRINT*,' FKLAM1  ', fklam1
        PRINT*,' ACL1    ', acl1
        PRINT*,' ACL2    ', acl2
        PRINT*,' CL11    ', cl11
        PRINT*,' CL21    ', cl21
        PRINT*,' DAL1    ', dal1
        PRINT*,' DAL2    ', dal2
        PRINT*,' FRH     ', frh
      ENDIF
C
C*    5. READ COMMON COUPLE.
C        -------------------
C
 5000 CONTINUE
      READ(kugrid) BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV
      IF (lltest) THEN
        PRINT*,' BETAMAX ', betamax
        PRINT*,' ZALP    ', zalp
        PRINT*,' ALPHA   ', alpha
        PRINT*,' XKAPPA  ', xkappa
        PRINT*,' XNLEV   ', xnlev
      ENDIF
C
C*    6. READ COMMON TABLE (INDICES OF OUTPUT POINTS).
C        ---------------------------------------------
C
 6000 CONTINUE
      READ(kugrid) TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP
      IF (lltest) THEN
        PRINT*,' TAUT    ', taut
        PRINT*,' DELTAUW ', deltauw
        PRINT*,' DELU    ', delu
        PRINT*,' TAUHFT  ', tauhft
        PRINT*,' DELUST  ', delust
        PRINT*,' DELALP  ', delalp
      ENDIF
C
C*    7. READ COMMON COUT (INDICES OF OUTPUT POINTS).
C        --------------------------------------------
C
 7000 CONTINUE
      READ(kugrid)  NGOUT, IGAR, IJAR
      IF (lltest) THEN
        PRINT*,' NGOUT ', ngout
        PRINT*,' IGAR  ', igar
        PRINT*,' IJAR  ', ijar
      ENDIF
C
C*    8. READ COMMON SHALLO (DEPTH AND SHALLOW WATER TABLES).
C        ----------------------------------------------------
C
 8000 CONTINUE
      READ(kugrid)  DEPTH, DEPTHA, DEPTHD, TCGOND, TFAK, TSIHKD
      IF (lltest) THEN
        PRINT*,' DEPTH  ', depth
        PRINT*,' DEPTHA ', deptha
        PRINT*,' DEPTHD ', depthd
        PRINT*,' TCGOND ', tcgond
        PRINT*,' TFAK   ', tfak
        PRINT*,' TSIHKD ', tsihkd
      ENDIF
C
C*    9. READ COMMON CURRENT.
C        --------------------
C
 9000 CONTINUE
      IF (IREFRA.EQ.2) THEN
        READ(kugrid)  U, V
        IF (lltest) THEN
          PRINT*, ' U ', u
          PRINT*, ' V ', v
        ENDIF
      ENDIF
      msub = msub - 1
      RETURN
      END
