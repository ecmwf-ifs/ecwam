      SUBROUTINE MUBUF (IA1, IOBSLAT, IOBSLON, IG, ML, IU08, IU18,IFORM)

! ----------------------------------------------------------------------

!**** *MUBUF* - ROUTINE TO ARRANGE YOWMON UBUF FOR ONE BLOCK.

!     H.GUNTHER            ECMWF       04/04/1990
!     J. BIDLOT            ECMWF       APRIL 2000: add second closest
!                                                  grid points.

!*    PURPOSE.
!     -------

!       TO ARRANGE NEIGHBOUR GRID POINT INDICES FOR A BLOCK

!**   INTERFACE.
!     ----------

!       *CALL* *MUBUF (IA1, IOBSLAT, IOBSLON,IG,ML,IU08, IU18, IFORM)*
!          *IA1*     - TOPOGRAPHIC DATA.
!          *IOBSLAT* - REDUCTION FACTOR DUE TO OBSTRUCTIONS IN
!                      THE MERIDIONAL DIRECTION (0.1%).
!          *IOBSLON* - REDUCTION FACTOR DUE TO OBSTRUCTIONS IN
!                      THE ZONAL DIRECTION (0.1%).
!          *IG*      - BLOCK NUMBER.
!          *ML*      - NUMBER OF FREQUENCIES.
!          *IU08*    - LOGICAL UNIT FOR OUTPUT OF GRID BLOCKING
!                      COMMON UBUF (UNFORMATED)
!          *IU18*    - LOGICAL UNIT FOR OUTPUT OF GRID BLOCKING
!                      COMMON UBUF (FORMATED)
!          *IFORM*   - OUTPUT FORMAT OPTION = 1 UNFORMATED
!                                           = 2 FORMATED
!                                           OTHERWISE BOTH

!     METHOD.
!     -------

!       THE INDICES OF THE NEXT POINTS ON LAT. AND LONG. ARE
!       COMPUTED. ZERO INDICATES A LAND POINT IS NEIGHBOUR.
!       THE FINAL COMMON UBUF IS WRITTEN OUT.

!     EXTERNALS.
!     ----------

!       *OUTUBUF*   - WRITE OUT COMMON UBUF.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWGRID  , ONLY : NLONRGG  ,IJLT
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,NY       ,IPER     ,
     &            ZDELLO   ,IRGG
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2    ,WADV

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IG,IU08,IU18,IFORM,IJ,IJP,I,K,IP,IH,IMIN,IPLUS,IS,ML,M
      INTEGER ::  IA1(NGX, NGY)
      INTEGER, DIMENSION(NGX, NGY, ML, 2) :: IOBSLAT, IOBSLON
      INTEGER, DIMENSION(NIBLO,ML,2) :: KOBSLON, KOBSLAT
      REAL :: XMIN, XPLUS, D0, D1, D2, D3, D4, D5, D6

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ---------------

      ALLOCATE(KLAT(NIBLO,2))
      ALLOCATE(KLON(NIBLO,2))
      ALLOCATE(KLAT2(NIBLO,2))
      ALLOCATE(WADV(NIBLO,2))


!*    1. INITIALISE ARRAYS.
!        ------------------

      DO IJP=1,2
        DO IJ=1,NIBLO
           KLAT(IJ,IJP) = 0
           KLAT2(IJ,IJP) = 0
           KLON(IJ,IJP) = 0
           WADV(IJ,IJP)=1.
        ENDDO
        DO M=1,ML
          DO IJ=1,NIBLO
             KOBSLAT(IJ,M,IJP)=1000
             KOBSLON(IJ,M,IJP)=1000
          ENDDO
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    2. COMPUTE INDICES OF NEIGHBOUR SEA POINTS.
!        ----------------------------------------


!*    2.1 LONGITUDE NEIGHBOURS.
!         ---------------------

      DO IP = 1,IJLT(IG)
        I = IXLG(IP,IG)
        K = KXLT(IP,IG)
        IF (I.GT.1) THEN
          IF (IA1(I-1,K).GT.-990) KLON(IP,1) = IP-1
        ELSE
          IF (IPER.EQ.1 .AND. IA1(NLONRGG(K),K).GT.-990) THEN
            KLON(IP,1) = IP
            DO IH=2,NLONRGG(K)
              IF (IA1(IH,K).GT.-990) KLON(IP,1) = KLON(IP,1)+1
            ENDDO
          ENDIF
        ENDIF
        IF (I.LT.NLONRGG(K)) THEN
          IF (IA1(I+1,K).GT.-990) KLON(IP,2) = IP+1
        ELSE
          IF (IPER.EQ.1 .AND. IA1(1,K).GT.-990) THEN
            KLON(IP,2) = IP
            DO IH=NLONRGG(K)-1,1,-1
              IF (IA1(IH,K).GT.-990) KLON(IP,2) = KLON(IP,2)-1
            ENDDO
          ENDIF
        ENDIF
      ENDDO

!*    2.2 LATITUDE NEIGHBOURS.
!         --------------------

      DO IP = 1,IJLT(IG)
        I = IXLG(IP,IG)
        K = KXLT(IP,IG)
        IF (K.GT.1) THEN
          XMIN = REAL(I-1)*ZDELLO(K)/ZDELLO(K-1)
          IMIN = NINT(XMIN) + 1

!         CLOSEST GRID POINT
          IF (IA1(IMIN,K-1).GT.-990) THEN
            DO IH = IP,1,-1
              IF (IXLG(IH,IG).EQ.IMIN .AND.
     &            KXLT(IH,IG).EQ.K-1) EXIT
            ENDDO
            KLAT(IP,1) = IH
          ENDIF

!         SECOND CLOSEST GRID POINT
          IF(IRGG.EQ.1) THEN
            IF(XMIN.LE.FLOAT(IMIN-1)) THEN
              IF(IMIN.LE.1) THEN
                IMIN=1
              ELSE
                IMIN=IMIN-1
              ENDIF
            ELSE
              IF(IMIN.GE.NLONRGG(K-1)) THEN
                IMIN=NLONRGG(K-1)
              ELSE
                IMIN=IMIN+1
              ENDIF
            ENDIF

            IF (IA1(IMIN,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN .AND.
     1              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
              KLAT2(IP,1) = IH
            ENDIF
          ELSE
            KLAT2(IP,1) = KLAT(IP,1)
          ENDIF

        ENDIF

        IF (K.LT.NY) THEN
          XPLUS = REAL(I-1)*ZDELLO(K)/ZDELLO(K+1)
          IPLUS = NINT(XPLUS) + 1
          IF (IA1(IPLUS,K+1).GT.-990) THEN
            DO IH = IP,IJLT(IG)
              IF (IXLG(IH,IG).EQ.IPLUS .AND.
     &            KXLT(IH,IG).EQ.K+1) EXIT
            ENDDO
            KLAT(IP,2) = IH
          ENDIF

          IF(IRGG.EQ.1) THEN
            IF(XPLUS.LE.FLOAT(IPLUS-1)) THEN
              IF(IPLUS.LE.1) THEN
                IPLUS=1
              ELSE
                IPLUS=IPLUS-1
              ENDIF
            ELSE
              IF(IPLUS.GE.NLONRGG(K+1)) THEN
                IPLUS=NLONRGG(K+1)
              ELSE
                IPLUS=IPLUS+1
              ENDIF
            ENDIF

            IF (IA1(IPLUS,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS .AND.
     1              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KLAT2(IP,2) = MIN(IH,IJLT(IG))
            ENDIF
          ELSE
            KLAT2(IP,2) = KLAT(IP,2)
          ENDIF
        ENDIF

      ENDDO

!     COMPUTE THE WEIGHT FOR THE ADVECTION SCHEME ON AN IRREGULAR GRID

      IF(IRGG.EQ.1) THEN
        DO IP = 1,IJLT(IG)
          I = IXLG(IP,IG)
          K = KXLT(IP,IG)
          D0 = FLOAT(I-1)*ZDELLO(K)
          D3=D0-0.5*ZDELLO(K)
          D5=D0+0.5*ZDELLO(K)

          IF (K.GT.1) THEN
            XMIN = D0/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1

            D4=FLOAT(IMIN-1)*ZDELLO(K-1)-0.5*ZDELLO(K-1)
            D6=FLOAT(IMIN-1)*ZDELLO(K-1)+0.5*ZDELLO(K-1)

            IF(D0.LE.FLOAT(IMIN-1)*ZDELLO(K-1)) THEN

              IF(D4.LE.D3.OR.D6.LE.D5) THEN
                 WADV(IP,1) = 1.
              ELSE
                 D2=D4-D3
                 D1=ZDELLO(K)-D2
                 WADV(IP,1)=MIN(1.,D1/ZDELLO(K))
              ENDIF
            ELSE

              IF(D4.GE.D3.OR.D6.GE.D5) THEN
                 WADV(IP,1) = 1.
              ELSE
                 D2=D5-D6
                 D1=ZDELLO(K)-D2
                 WADV(IP,1)=MIN(1.,D1/ZDELLO(K))
              ENDIF
            ENDIF

          ENDIF

          IF (K.LT.NY) THEN
            XPLUS = D0/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1

            D4=FLOAT(IPLUS-1)*ZDELLO(K+1)-0.5*ZDELLO(K+1)
            D6=FLOAT(IPLUS-1)*ZDELLO(K+1)+0.5*ZDELLO(K+1)

            IF(D0.LE.FLOAT(IPLUS-1)*ZDELLO(K+1)) THEN

              IF(D4.LE.D3.OR.D6.LE.D5) THEN
                 WADV(IP,2) = 1.
              ELSE
                 D2=D4-D3
                 D1=ZDELLO(K)-D2
                 WADV(IP,2)=MIN(1.,D1/ZDELLO(K))
              ENDIF
            ELSE

              IF(D4.GE.D3.OR.D6.GE.D5) THEN
                 WADV(IP,2) = 1.
              ELSE
                 D2=D5-D6
                 D1=ZDELLO(K)-D2
                 WADV(IP,2)=MIN(1.,D1/ZDELLO(K))
              ENDIF
            ENDIF

          ENDIF

        ENDDO
      ENDIF

!     DETERMINE THE REAL REDUCTION FACTORS TO BE USED IN THE PROPAGATION

      DO IS=1,2
        DO M=1,ML
          DO IP = 1,IJLT(IG)
            I = IXLG(IP,IG)
            K = KXLT(IP,IG)
            KOBSLAT(IP,M,IS)=IOBSLAT(I,K,M,IS)
            KOBSLON(IP,M,IS)=IOBSLON(I,K,M,IS)
          ENDDO
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. OUTPUT OF COMMON UBUF FOR ONE BLOCK.
!        ------------------------------------
      CALL OUTUBUF (IU08, IU18, IFORM,NIBLO,ML,KOBSLAT,KOBSLON)

      RETURN
      END SUBROUTINE MUBUF
