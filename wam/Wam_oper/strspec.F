      SUBROUTINE STRSPEC (NFRE, NANG, ML, KL, FR, FL, GAMMA)

C ----------------------------------------------------------------------
C
C**** *STRSPEC* - ROUTINE TO STRETCH A SPECTRUM.
C
C      EVA BAUER      MPI  HAMBURG    MAY 1990.
C      H. GUNTHER     GKSS/ECMWF      JAN 1991  MODIFIED FOR CYCLE_4.
C
C*     PURPOSE.
C      -------
C
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *STRSPEC (NFRE, NANG, ML, KL, FR, FL, GAMMA)*
C         *NFRE*  - FREQUENCY DIMENSION OF ARRAYS.
C         *NANG*  - DIRECTION DIMENSION OF ARRAYS.
C         *ML*    - NUMBER OF FREQUENCIES.
C         *KL*    - NUMBER OF DIRECTIONS.
C         *FR*    - FREQUENCY ARRAY.
C         *FL*    - INPUT AND OUTPUT SPECTRUM.
C         *GAMMA* - STRETCHING PARAMETER.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ---------
C
C       NONE.
C
C     REFERENCES.
C     -----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      DIMENSION FL(NANG,NFRE), FR(NFRE)
      DIMENSION AR1(KL,ML), AR2(ML)
C
C---------------------------------------------------------------------
C
      IF (GAMMA.EQ.1.0) RETURN
C
C*    1. INITIALIZATION.
C        ---------------
C
 1000 CONTINUE
      DO M=1,ML
        DO K=1,KL
         AR1(K,M) = 0.0
        ENDDO
      ENDDO
      ALO = LOG10(1.1)
      GAMS = GAMMA
C
C*    2. DETERMINE ACROSS HOW MANY FREQUENCY BINS THE
C        STRETCHING IS ACTING AND THE STRECHED FREQUENCIES.
C        ---------------------------------------------------
C
 2000 CONTINUE
      INC = IFIX(LOG10(GAMS)/ALO)
      Z = ABS(1.1**INC - GAMS)
      DO 2001 M=1,ML
         AR2(M) = FR(M) * GAMS
 2001 CONTINUE
C
C*    3. STRECH SPECTRUM.
C        ----------------
C
 3000 CONTINUE
      IF (Z.LE.0.001) THEN
C
C*    3.1 SHIFT SPECTRUM IF GAMMA IS A POWER OF 1.1.
C         ------------------------------------------
C
         IF (GAMS.GT.1.0) THEN
C
C*    3.1.1 SHIFT TO LOWER FREQUENIES.
C           --------------------------
C
            DO 3111 M = 1,ML-INC
               MC = M + INC
               DO 3112 K = 1,KL
                  AR1(K,M) = FL(K,MC)
 3112       CONTINUE
 3111       CONTINUE
         ELSE
C
C*    3.1.2 SHIFT TO HIGHER FREQUENCIES.
C           ----------------------------
C
            DO 3121 M = 1-INC,ML
               MC = M + INC
               DO 3122 K = 1,NANG
                  AR1(K,M) = FL(K,MC)
 3122          CONTINUE
 3121       CONTINUE
         ENDIF
      ELSE
C
C*    3.2 SHIFT AND LINEAR INTERPOLATION OF SPECTRAL ENERGY
C*        IF GAMMA IS NOT A POWER OF 1.1.
C*
C*        (SPECTRUM HAS ZERO ENERGY AT FREQUENCY FR(NFRE) )
C         -------------------------------------------------
C
         IF (GAMS.GT.1.0) THEN
C
C*    3.2.1 SHIFT TO LOWER FREQUENCIES.
C           ---------------------------
C
            DO 3211 M=1,ML-INC-1
               IFR = IFIX(ALOG10(AR2(M)/FR(1))/ALO+1.)
               IFRP1 = IFR+1
               MC = M + INC
               ADIF = (FR(IFRP1)-AR2(M)) / (FR(IFRP1)-FR(IFR))
               BDIF = 1. - ADIF
               DO 3212 K=1,KL
                  AR1(K,M) = ADIF * FL(K,MC) + BDIF * FL(K,MC+1)
 3212          CONTINUE
 3211       CONTINUE
         ELSE
C
C*     3.2.2 SHIFT TO HIGHER FREQUENCIES.
C            ----------------------------
C
            DO 3221 M = 2-INC,ML
               IFR = IFIX(ALOG10(AR2(M)/FR(1))/ALO+1.)
               IFRP1 = IFR+1
               MC = M + INC - 1
               ADIF = (FR(IFRP1)-AR2(M)) / (FR(IFRP1)-FR(IFR))
               BDIF = 1. - ADIF
               DO 3222 K=1,KL
                  AR1(K,M) = ADIF * FL(K,MC) + BDIF * FL(K,MC+1)
 3222          CONTINUE
 3221       CONTINUE
         ENDIF
      ENDIF
C
C*    4. COPY NEW TO OLD.
C        ----------------
C
 4000 CONTINUE

      DO M=1,ML
        DO K=1,KL
         FL(K,M) = AR1(K,M)
        ENDDO
      ENDDO

      RETURN
      END
