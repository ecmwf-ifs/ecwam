c--------------------------------------------------------------------
c
      SUBROUTINE expta0 (urtail)
c
c----------------------------------------------------------------------
C
C      PURPOSE.
C      --------
C
C         COMPUTATION OF EXPONENTIAL MEAN VELOCITY SPREAD
C
C
C*     INTERFACE.
C      ----------
C        
C         *CALL* *EXPTA0*(urtail)
c
c         *urtail* - mean orbital vel. of tail (k>0.2) of wave spe.C
C      EXTERNALS.
C      ----------
C         NONE.
C
C      METHOD.
C      -------
C         NONE.
C
C----------------------------------------------------------------------
C
      PARAMETER (KSX = 64, KSY = 64, KSXT = 2*KSX, KSYT = 2*KSY,
     1           KSXP1 = KSX+1, KSXM1 = KSX-1, KSYP1 = KSY+1,
     2           KSYM1 = KSY-1, LOT = KSYT, KSXLOT = 4*KSXT*LOT,
     3           MGLPR = 10, NOR = 13, NFRE = 26, NANG = 12,
     4           NANG1 = NANG+1, NBLO = 20, NIBLO =3000)
C
C-----------------------------------------------------------------------
C
C     VARIABLE      TYPE     PURPOSE
C     --------      ----     -------
C
C      *NANG*      INTEGER  NUMBER OF DIRECT. IN MODEL SPECTRUM.
C      *NFRE*      INTEGER  NUMBER OF FREQUENCIES IN MODEL SPECTRUM.
C      *KSX*       INTEGER  NUMBER OF WAVE NUMBERS IN SAR SPECTRUM(X-D).
C      *KSY*       INTEGER  NUMBER OF WAVE NUMBERS IN SAR SPECTRUM(Y-D).
C      *MGLPR*     INTEGER
C      *NBLO*      INTEGER  NUMBER OF BLOCKS IN WAVE MODEL GRID.
C      *NIBLO*     INTEGER  NUMBER OF INDICES IN A BLOCK.
C      *NOR*       INTEGER  NUMBER OF SERIES EXPANSIONS
C
C-----------------------------------------------------------------------
C
C
C     *COMMON* *AFFT*  EXPRESSIONS COMPUTED AFTER FFTS,DEPENDANT ON
C                      KX,R
C
C-----------------------------------------------------------------------
C
      COMMON/AFFT/ EXETA0(KSXT),SUM3
C
C-----------------------------------------------------------------------
C
C*    VARIABLE    TYPE      PURPOSE.
C     --------    ----      --------
C
C     *EXETAC*    REAL      EXPONENTIAL CUT OFF TERM
C     *SUM3*      REAL
C-----------------------------------------------------------------------
C
C
C      *COMMON* *TKPRE* PARTS OF MAPPING MATRIX INDEPENDANT OF
C                       OF SPECTRUM AND SATELLITE FLIGHT DIRECTION.
C
C
      COMMON/TKPRE/ TKVR(KSXT,KSYT),TKVI(KSXT,KSYT),TKSVABS(KSXT,KSYT),
     1              TK11(KSXT,KSYT),TK12(KSXT,KSYT),TK2(KSYT),
     2              IFLARAR, AMPCOSR, AMPSINR
C----------------------------------------------------------------------
C
C      VARIABLE      TYPE      PURPOSE.
C      --------      ----      --------
C
C      *TKVR*        REAL      - REAL PART OF TKV ,NOT TKVB.
C      *TKVI*        REAL      - IMAG.PART OF TKV.
C      *TKSVABS*     REAL      SQUARE OF ABSOLUTE VALUE OF TKV.
C       *TK1*        REAL      FIRST PART OF HYDRODYNAMIC MTF.
C       *TK2*        REAL      SECOND PART OF HYDRODYNAMIC MTF.
C      *IFLARAR*     INTEGER   = 1 THEORET. RAR TRANSFER FUNCTION
C                              = 0 EMPIRICAL RAR TRANSFER FUNCTION
C      *AMPCOSR*     REAL      AMPLITUDE * COS PHASE FOR RAR TRANSF.FUNC
C      *AMPSINR*     REAL      AMPLITUDE * SIN PHASE FOR RAR TRANSF.FUNC
C
C-----------------------------------------------------------------------
C
C
C      *COMMON* *WAM* - SPECTRUM AND LAT-LONG INFORMATION FOR
C                       MODEL OUTPUT.
C
C
      COMMON /WAM/ FWAM(NFRE,NANG),FWAMSAR(KSXT,KSYT)
     1             ,XLAT(NBLO,NIBLO),XLASP(MGLPR)
     2             ,XLON(NBLO,NIBLO),XLOSP(MGLPR)
C----------------------------------------------------------------------
C
C      VARIABLE      TYPE      PURPOSE.
C      --------      ----      --------
C
C      *FWAM*        REAL      WAM-MODEL SPECTRUM.
C      *FWAMSAR*     REAL      WAM-MODEL SPEC.AFTER TRANSF.INTO FLIGHT D
C      *XLATH*       REAL      MAPS BLOCK INDEX (IG) AND INDEX IN A
C                              BLOCK (IJ) TO LATITUDINAL COORDINATE.
C      *XLASP*       REAL      LATITUDINAL COORDINATE OF SPECTRUM
C                              TO BE COMPARED WITH SAR SPECTRUM.
C      *XLON*        REAL      MAPS BLOCK INDEX (IG) AND INDEX IN A
C                              BLOCK (IJ) TO LONGITUDINAL COORDINATE.
C      *XLOSP*       REAL      LONGITUDINAL COORDINATE OF SPECTRUM
C                              TO BE COMPARED WITH SAR SPECTRUM.
C----------------------------------------------------------------------
C
C
C      *COMMON* *WANU* - WAVENUMBER ARRAYS AND EQUIVALENT FREQUENCY
C                        ARRAYS.
C
      COMMON /WANU/ ANGINC,DELK,FRKX(NFRE,NANG),FRKY(NFRE,NANG),
     1              FR(NFRE),CTHET(NANG),STHET(NANG),THETAN(NANG1),
     2              WANUAB(KSXT,KSYT),SARKX(KSXT),SARKY(KSYT),BETA,
     3              BEKXH(KSXT),BEKXHS(KSXT),BSARKX(KSXT)
C
C-----------------------------------------------------------------------
C
C      VARIABLE      TYPE      PURPOSE.
C      --------      ----      --------
C
C      *ANGINC*     REAL       RADAR INCIDENCE ANGLE
C      *DELK*       REAL       WAVE NUMBER INCREMENT.
C      *FRKX*       REAL       KX COMPUTED FROM (F,THETA).
C      *FRKY*       REAL       KX COMPUTED FROM (F,THETA).
C      *FR*         REAL       MODEL FREQUENCY ARRAY.
C      *THETAN*     REAL       DIRECTION OF WAVENUMBER VECTOR.
C      *CTHET*      REAL       COSINE OF THETA.
C      *STHET*      REAL       SINE OF THETA.
C      *WANUAB*     REAL       ABSOLUTE VALUE OF WAVE NUMBER.
C      *SARKX*      REAL       X-COMPONENT OF WAVE NUMBER.
C      *SARKY*      REAL       Y-COMPONENT OF WAVE NUMBER.
C      *BETA*       REAL       RO/USAT
C      *BEKXH*      REAL       (BETA*KX)**2/2.
C      *BEKXHS*     REAL       BEKXH**2
C      *BSARKX*     REAL       BETA*KX
C
C-----------------------------------------------------------------------
C
C     =====================
C     1. COMPUTATION OF ETA
C     ---------------------
C
 1000 CONTINUE
c
      urtq = urtail*urtail
      SUM3 = 0.
      DO 901    J=1,KSYT+1
      L = J
      IF(J.EQ.KSYT+1) L = 1
      DO 902    I=1,KSXT+1
      K = I
      IF(K.EQ.KSXT+1) K = 1
      SUM3 = SUM3 + FWAMSAR(K,L)*TKSVABS(K,L)*DELK*DELK
  902 CONTINUE
  901 CONTINUE
C
      WRITE(*,*)' ur**2, tail ur**2 ',sum3,urtq
      AMPET0=2*(SUM3 + urtq)

C
      DO 1011 IRX=1,KSXT
      EXETA0(IRX)=EXP(-AMPET0*BEKXH(IRX))
 1011 CONTINUE
C
      RETURN
      END
