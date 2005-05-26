      SUBROUTINE STRSPEC (NFRE, NANG, ML, KL, FR, FL, GAMMA)

! ----------------------------------------------------------------------

!**** *STRSPEC* - ROUTINE TO STRETCH A SPECTRUM.

!      EVA BAUER      MPI  HAMBURG    MAY 1990.
!      H. GUNTHER     GKSS/ECMWF      JAN 1991  MODIFIED FOR CYCLE_4.

!*     PURPOSE.
!      -------


!**   INTERFACE.
!     ----------

!       *CALL* *STRSPEC (NFRE, NANG, ML, KL, FR, FL, GAMMA)*
!         *NFRE*  - FREQUENCY DIMENSION OF ARRAYS.
!         *NANG*  - DIRECTION DIMENSION OF ARRAYS.
!         *ML*    - NUMBER OF FREQUENCIES.
!         *KL*    - NUMBER OF DIRECTIONS.
!         *FR*    - FREQUENCY ARRAY.
!         *FL*    - INPUT AND OUTPUT SPECTRUM.
!         *GAMMA* - STRETCHING PARAMETER.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ---------

!       NONE.

!     REFERENCES.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      DIMENSION FL(NANG,NFRE), FR(NFRE)
      DIMENSION AR1(KL,ML), AR2(ML)

!---------------------------------------------------------------------

      IF (GAMMA.EQ.1.0) RETURN

!*    1. INITIALIZATION.
!        ---------------

      DO M=1,ML
        DO K=1,KL
          AR1(K,M) = 0.0
        ENDDO
      ENDDO
      ALO = LOG10(1.1)
      GAMS = GAMMA

!*    2. DETERMINE ACROSS HOW MANY FREQUENCY BINS THE
!        STRETCHING IS ACTING AND THE STRECHED FREQUENCIES.
!        ---------------------------------------------------

      INC = INT(LOG10(GAMS)/ALO)
      Z = ABS(1.1**INC - GAMS)
      DO M=1,ML
        AR2(M) = FR(M) * GAMS
      ENDDO

!*    3. STRECH SPECTRUM.
!        ----------------

      IF (Z.LE.0.001) THEN

!*    3.1 SHIFT SPECTRUM IF GAMMA IS A POWER OF 1.1.
!         ------------------------------------------

        IF (GAMS.GT.1.0) THEN

!*    3.1.1 SHIFT TO LOWER FREQUENIES.
!           --------------------------

          DO M = 1,ML-INC
            MC = M + INC
            DO K = 1,KL
              AR1(K,M) = FL(K,MC)
            ENDDO
          ENDDO
        ELSE

!*    3.1.2 SHIFT TO HIGHER FREQUENCIES.
!           ----------------------------

          DO M = 1-INC,ML
            MC = M + INC
            DO K = 1,NANG
              AR1(K,M) = FL(K,MC)
            ENDDO
          ENDDO
        ENDIF
      ELSE

!*    3.2 SHIFT AND LINEAR INTERPOLATION OF SPECTRAL ENERGY
!*        IF GAMMA IS NOT A POWER OF 1.1.
!*
!*        (SPECTRUM HAS ZERO ENERGY AT FREQUENCY FR(NFRE) )
!         -------------------------------------------------

        IF (GAMS.GT.1.0) THEN

!*    3.2.1 SHIFT TO LOWER FREQUENCIES.
!           ---------------------------

          DO M=1,ML-INC-1
            IFR = INT(LOG10(AR2(M)/FR(1))/ALO+1.)
            IFRP1 = IFR+1
            MC = M + INC
            ADIF = (FR(IFRP1)-AR2(M)) / (FR(IFRP1)-FR(IFR))
            BDIF = 1. - ADIF
            DO K=1,KL
              AR1(K,M) = ADIF * FL(K,MC) + BDIF * FL(K,MC+1)
            ENDDO
          ENDDO
        ELSE

!*     3.2.2 SHIFT TO HIGHER FREQUENCIES.
!            ----------------------------

          DO M = 2-INC,ML
            IFR = INT(LOG10(AR2(M)/FR(1))/ALO+1.)
            IFRP1 = IFR+1
            MC = M + INC - 1
            ADIF = (FR(IFRP1)-AR2(M)) / (FR(IFRP1)-FR(IFR))
            BDIF = 1. - ADIF
            DO K=1,KL
              AR1(K,M) = ADIF * FL(K,MC) + BDIF * FL(K,MC+1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!*    4. COPY NEW TO OLD.
!        ----------------


      DO M=1,ML
        DO K=1,KL
          FL(K,M) = AR1(K,M)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE STRSPEC




