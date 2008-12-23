      SUBROUTINE MUBUF (IU01, BATHY, IG, ML, IU08, IU18,IFORM)
! ----------------------------------------------------------------------

!**** *MUBUF* - ROUTINE TO ARRANGE YOWMON UBUF FOR ONE BLOCK.

!     H.GUNTHER            ECMWF       04/04/1990
!     J. BIDLOT            ECMWF       APRIL 2000: add second closest
!                                                  grid points.
!     J. BIDLOT            ECMWF       CLOSEST AND SECOND CLOSEST GRID
!                                      POINT FOR THE ROTATED CELL.
!                               IN ORDER TO SAVE MEMORY THE OBSTRUCTION
!                               COEFFICIENTS ARE READ IN AND PROCESSED
!                               SEQUENTIALLY.

!*    PURPOSE.
!     -------

!       TO ARRANGE NEIGHBOUR GRID POINT INDICES FOR A BLOCK

!**   INTERFACE.
!     ----------

!       *CALL* *MUBUF (BATHY,IG,ML,IU08, IU18, IFORM)*
!          *IU01*  -  LOGICAL INPUT UNIT OF TOPOGRAPHIC DATA.
!          *BATHY*     -  BATHYMETRY DATA.
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
     &            XDELLA   ,ZDELLO   ,IRGG     ,LLOBSTRCT
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2    ,KCOR     ,
     &                      KRLAT    ,KRLAT2   ,KRLON    ,KRLON2   ,
     &                      WLAT     ,WRLAT    ,WRLON

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IU01
      INTEGER :: NFREMAX, IX
      INTEGER :: IG,IU08,IU18,IFORM,IJ,IJP,I,K,IP,IH,IS,ML,M
      INTEGER :: IMIN,IPLUS, IMIN2, IPLUS2, IMIN3, IPLUS3
      INTEGER :: IC, ICP
      REAL ::  BATHY(NGX, NGY)
      INTEGER, DIMENSION(NGX, NGY) :: IDUM
      INTEGER, DIMENSION(NIBLO,ML,2) :: KOBSLON, KOBSLAT
      INTEGER, DIMENSION(NIBLO,ML,4) :: KOBSCOR
      INTEGER, DIMENSION(NIBLO,ML,2) :: KOBSRLON, KOBSRLAT
      REAL :: XMIN, XPLUS, D0, D1, D2, D3, D4, D5, D6
      REAL :: XLON, XL, XP, TWOXDELLA
      CHARACTER*4 :: CX
      CHARACTER*10 :: FORMAT
      LOGICAL :: LLABORT

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ---------------

      ALLOCATE(KLAT(NIBLO,2))
      ALLOCATE(KLON(NIBLO,2))
      ALLOCATE(KLAT2(NIBLO,2))
      ALLOCATE(KCOR(NIBLO,4,-1:1))
      ALLOCATE(KRLAT(NIBLO,2))
      ALLOCATE(KRLON(NIBLO,2))
      ALLOCATE(KRLAT2(NIBLO,2))
      ALLOCATE(KRLON2(NIBLO,2))
      ALLOCATE(WLAT(NIBLO,2))
      ALLOCATE(WRLAT(NIBLO,2))
      ALLOCATE(WRLON(NIBLO,2))


!*    1. INITIALISE ARRAYS.
!        ------------------

      TWOXDELLA=2*XDELLA

      DO IJP=1,2
        DO IJ=1,NIBLO
           KLAT(IJ,IJP) = 0
           KLAT2(IJ,IJP) = 0
           KLON(IJ,IJP) = 0
           KRLAT(IJ,IJP) = 0
           KRLAT2(IJ,IJP) = 0
           KRLON(IJ,IJP) = 0
           KRLON2(IJ,IJP) = 0
           WLAT(IJ,IJP)=1.
           WRLAT(IJ,IJP)=1.
           WRLON(IJ,IJP)=1.
        ENDDO
        DO M=1,ML
          DO IJ=1,NIBLO
             KOBSLAT(IJ,M,IJP)=1000
             KOBSLON(IJ,M,IJP)=1000
             KOBSRLAT(IJ,M,IJP)=1000
             KOBSRLON(IJ,M,IJP)=1000
          ENDDO
        ENDDO
      ENDDO
      DO ICP=-1,1
        DO IC=1,4
          DO IJ=1,NIBLO
             KCOR(IJ,IC,ICP) = 0
          ENDDO
        ENDDO
      ENDDO
      DO IC=1,4
        DO M=1,ML
          DO IJ=1,NIBLO
             KOBSCOR(IJ,M,IC)=1000
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
          IF (BATHY(I-1,K).GT.-990.) KLON(IP,1) = IP-1
        ELSE
          IF (IPER.EQ.1 .AND. BATHY(NLONRGG(K),K).GT.-990.) THEN
            KLON(IP,1) = IP
            DO IH=2,NLONRGG(K)
              IF (BATHY(IH,K).GT.-990.) KLON(IP,1) = KLON(IP,1)+1
            ENDDO
          ENDIF
        ENDIF
        IF (I.LT.NLONRGG(K)) THEN
          IF (BATHY(I+1,K).GT.-990.) KLON(IP,2) = IP+1
        ELSE
          IF (IPER.EQ.1 .AND. BATHY(1,K).GT.-990.) THEN
            KLON(IP,2) = IP
            DO IH=NLONRGG(K)-1,1,-1
              IF (BATHY(IH,K).GT.-990.) KLON(IP,2) = KLON(IP,2)-1
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
          IF (BATHY(IMIN,K-1).GT.-990.) THEN
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
                IMIN2=1
              ELSE
                IMIN2=IMIN-1
              ENDIF
            ELSE
              IF(IMIN.GE.NLONRGG(K-1)) THEN
                IMIN2=NLONRGG(K-1)
              ELSE
                IMIN2=IMIN+1
              ENDIF
            ENDIF

            IF (BATHY(IMIN2,K-1).GT.-990.) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN2 .AND.
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
          IF (BATHY(IPLUS,K+1).GT.-990.) THEN
            DO IH = IP,IJLT(IG)
              IF (IXLG(IH,IG).EQ.IPLUS .AND.
     &            KXLT(IH,IG).EQ.K+1) EXIT
            ENDDO
            KLAT(IP,2) = IH
          ENDIF

          IF(IRGG.EQ.1) THEN
            IF(XPLUS.LE.FLOAT(IPLUS-1)) THEN
              IF(IPLUS.LE.1) THEN
                IPLUS2=1
              ELSE
                IPLUS2=IPLUS-1
              ENDIF
            ELSE
              IF(IPLUS.GE.NLONRGG(K+1)) THEN
                IPLUS2=NLONRGG(K+1)
              ELSE
                IPLUS2=IPLUS+1
              ENDIF
            ENDIF

            IF (BATHY(IPLUS2,K+1).GT.-990.) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS2 .AND.
     1              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KLAT2(IP,2) = MIN(IH,IJLT(IG))
            ENDIF
          ELSE
            KLAT2(IP,2) = KLAT(IP,2)
          ENDIF
        ENDIF

      ENDDO

!     COMPUTE THE CORNER GRID POINT AND THE CLOSEST GRID POINT ON EITHER
!     SIDE (ON A GIVEN LATITUDE)
!     ------------------------------------------------------------------
      DO IP = 1,IJLT(IG)
        I = IXLG(IP,IG)
        K = KXLT(IP,IG)
        XLON = REAL(I-1)*ZDELLO(K)

        IF (K.GT.1) THEN
!         CLOSEST GRID POINT IN SW GRID CORNER POINT
          XL=XLON-ZDELLO(K)
          XMIN = XL/ZDELLO(K-1)
          IMIN = NINT(XMIN) + 1
          IF (IPER.EQ.1 .AND. IMIN .LT. 1)  THEN
            IMIN = IMIN + NLONRGG(K-1)
            XMIN = XMIN + FLOAT(NLONRGG(K-1))
          ENDIF
          IF (IMIN .GE. 1)  THEN
            IF (BATHY(IMIN,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN .AND.
     &              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
             KCOR(IP,3,0) = IH
            END IF
          ENDIF

!         SECOND CLOSEST GRID POINT IN SW GRID CORNER POINT
          IF (IMIN .GE. 1)  THEN
            IF(XMIN.LE.FLOAT(IMIN-1)) THEN
              IF(IMIN.LE.1) THEN
!!test                IMIN2=1
                IMIN2=NLONRGG(K-1)
              ELSE
                IMIN2=IMIN-1
              ENDIF
            ELSE
              IF(IMIN.GE.NLONRGG(K-1)) THEN
!!test                IMIN2=NLONRGG(K-1)
                IMIN2=1
              ELSE
                IMIN2=IMIN+1
              ENDIF
            ENDIF

            IF (BATHY(IMIN2,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN2 .AND.
     1              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
              KCOR(IP,3,1) = IH
            ENDIF
          ENDIF

!         THIRD CLOSEST GRID POINT IN SW GRID CORNER POINT
          IF (IMIN .GE. 1)  THEN
            IF(XMIN.LE.FLOAT(IMIN-1)) THEN
              IF(IMIN.GE.NLONRGG(K-1)) THEN
                IMIN3=1
              ELSE
                IMIN3=IMIN+1
              ENDIF
            ELSE
              IF(IMIN.LE.1) THEN
                IMIN3=NLONRGG(K-1)
              ELSE
                IMIN3=IMIN-1
              ENDIF
            ENDIF

            IF (BATHY(IMIN3,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN3 .AND.
     1              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
              KCOR(IP,3,-1) = IH
            ENDIF
          ELSE IF (IMIN.GE.0) THEN
            IMIN3=IMIN+1
            IF (BATHY(IMIN3,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN3 .AND.
     1              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
              KCOR(IP,3,-1) = IH
            ENDIF
          ENDIF

!         CLOSEST GRID POINT IN SE GRID CORNER
          XL=XLON+ZDELLO(K)
          XMIN = XL/ZDELLO(K-1)
          IMIN = NINT(XMIN) + 1
          IF (IPER.EQ.1 .AND. IMIN .GT. NLONRGG(K-1))  THEN
            IMIN = IMIN - NLONRGG(K-1)
            XMIN = XMIN - FLOAT(NLONRGG(K-1))
          ENDIF
          IF (IMIN .LE. NLONRGG(K-1))  THEN
            IF (BATHY(IMIN,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN .AND.
     &              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
             KCOR(IP,2,0) = IH
            END IF
          ENDIF

!         SECOND CLOSEST GRID POINT IN SE GRID CORNER
          IF (IMIN .LE. NLONRGG(K-1))  THEN
            IF(XMIN.LE.FLOAT(IMIN-1)) THEN
              IF(IMIN.LE.1) THEN
!!test                IMIN2=1
                IMIN2=NLONRGG(K-1)
              ELSE
                IMIN2=IMIN-1
              ENDIF
            ELSE
              IF(IMIN.GE.NLONRGG(K-1)) THEN
!!test                IMIN2=NLONRGG(K-1)
                IMIN2=1
              ELSE
                IMIN2=IMIN+1
              ENDIF
            ENDIF

            IF (BATHY(IMIN2,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN2 .AND.
     1              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
              KCOR(IP,2,1) = IH
            ENDIF
          ENDIF

!         THIRD CLOSEST GRID POINT IN SE GRID CORNER
          IF (IMIN .LE. NLONRGG(K-1))  THEN
            IF(XMIN.LE.FLOAT(IMIN-1)) THEN
              IF(IMIN.GE.NLONRGG(K-1)) THEN
                IMIN3=1
              ELSE
                IMIN3=IMIN+1
              ENDIF
            ELSE
              IF(IMIN.LE.1) THEN
                IMIN3=NLONRGG(K-1)
              ELSE
                IMIN3=IMIN-1
              ENDIF
            ENDIF

            IF (BATHY(IMIN3,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN3 .AND.
     1              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
              KCOR(IP,2,-1) = IH
            ENDIF
          ELSE IF (IMIN-1 .LE. NLONRGG(K-1))  THEN
            IMIN3=IMIN-1
            IF (BATHY(IMIN3,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN3 .AND.
     1              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
              KCOR(IP,2,-1) = IH
            ENDIF
          ENDIF


        ENDIF ! END K > 1


        IF (K.LT.NY) THEN
!         CLOSEST GRID POINT IN NW GRID CORNER 
          XL=XLON-ZDELLO(K)
          XPLUS = XL/ZDELLO(K+1)
          IPLUS = NINT(XPLUS) + 1
          IF (IPER.EQ.1 .AND. IPLUS .LT. 1)  THEN
            IPLUS = NLONRGG(K+1) + IPLUS 
            XPLUS = XPLUS + FLOAT(NLONRGG(K+1))
          ENDIF

          IF (IPLUS .GE. 1)  THEN
            IF (BATHY(IPLUS,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS .AND.
     &              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,4,0) = IH
            ENDIF
          ENDIF

!         SECOND CLOSEST GRID POINT IN NW GRID CORNER 
          IF (IPLUS .GE. 1)  THEN
            IF(XPLUS.LE.FLOAT(IPLUS-1)) THEN
              IF(IPLUS.LE.1) THEN
!!test                IPLUS2=1
                IPLUS2=NLONRGG(K+1)
              ELSE
                IPLUS2=IPLUS-1
              ENDIF
            ELSE
              IF(IPLUS.GE.NLONRGG(K+1)) THEN
!!test                IPLUS2=NLONRGG(K+1)
                IPLUS2=1
              ELSE
                IPLUS2=IPLUS+1
              ENDIF
            ENDIF
 
            IF (BATHY(IPLUS2,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS2 .AND.
     1              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,4,1) = MIN(IH,IJLT(IG))
            ENDIF
          ENDIF

!         THIRD CLOSEST GRID POINT IN NW GRID CORNER 
          IF (IPLUS .GE. 1)  THEN
            IF(XPLUS.LE.FLOAT(IPLUS-1)) THEN
              IF(IPLUS.GE.NLONRGG(K+1)) THEN
                IPLUS3=1
              ELSE
                IPLUS3=IPLUS+1
              ENDIF
            ELSE
              IF(IPLUS.LE.1) THEN
                IPLUS3=NLONRGG(K+1)
              ELSE
                IPLUS3=IPLUS-1
              ENDIF
            ENDIF
 
            IF (BATHY(IPLUS3,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS3 .AND.
     1              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,4,-1) = MIN(IH,IJLT(IG))
            ENDIF
          ELSE IF (IPLUS .GE. 0)  THEN
            IPLUS3=IPLUS+1
            IF (BATHY(IPLUS3,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS3 .AND.
     1              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,4,-1) = MIN(IH,IJLT(IG))
            ENDIF
          ENDIF

!         CLOSEST GRID POINT IN NE GRID CORNER
          XL=XLON+ZDELLO(K)
          XPLUS = XL/ZDELLO(K+1)
          IPLUS = NINT(XPLUS) + 1
          IF (IPER.EQ.1 .AND. IPLUS .GT. NLONRGG(K+1))  THEN
            IPLUS = IPLUS - NLONRGG(K+1)
            XPLUS = XPLUS - FLOAT(NLONRGG(K+1))
          ENDIF

          IF (IPLUS .LE. NLONRGG(K+1))  THEN
            IF (BATHY(IPLUS,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS .AND.
     &              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,1,0) = IH
            ENDIF
          ENDIF

!         SECOND CLOSEST GRID POINT IN NE GRID CORNER 
          IF (IPLUS .LE. NLONRGG(K+1))  THEN
            IF(XPLUS.LE.FLOAT(IPLUS-1)) THEN
              IF(IPLUS.LE.1) THEN
!!test                IPLUS2=1
                IPLUS2=NLONRGG(K+1)
              ELSE
                IPLUS2=IPLUS-1
              ENDIF
            ELSE
              IF(IPLUS.GE.NLONRGG(K+1)) THEN
!!test                IPLUS2=NLONRGG(K+1)
                IPLUS2=1
              ELSE
                IPLUS2=IPLUS+1
              ENDIF
            ENDIF

            IF (BATHY(IPLUS2,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS2 .AND.
     1              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,1,1) = MIN(IH,IJLT(IG))
            ENDIF
          ENDIF

!         THIRD CLOSEST GRID POINT IN NE GRID CORNER 
          IF (IPLUS .LE. NLONRGG(K+1))  THEN
            IF(XPLUS.LE.FLOAT(IPLUS-1)) THEN
              IF(IPLUS.GE.NLONRGG(K+1)) THEN
                IPLUS3=1
              ELSE
                IPLUS3=IPLUS+1
              ENDIF
            ELSE
              IF(IPLUS.LE.1) THEN
                IPLUS3=NLONRGG(K+1)
              ELSE
                IPLUS3=IPLUS-1
              ENDIF
            ENDIF

            IF (BATHY(IPLUS3,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS3 .AND.
     1              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,1,-1) = MIN(IH,IJLT(IG))
            ENDIF
          ELSE IF (IPLUS-1 .LE. NLONRGG(K+1))  THEN
            IPLUS3=IPLUS-1
            IF (BATHY(IPLUS3,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS3 .AND.
     1              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,1,-1) = MIN(IH,IJLT(IG))
            ENDIF
          ENDIF

        ENDIF  ! END K < NY

      ENDDO


!     COMPUTE THE CLOSEST GRID POINTS IN THE SW-NE and SE-NW DIRECTIONS
!     (I.E. GOING AT 45 DEGREE)
!     -----------------------------------------------------------------
      DO IP = 1,IJLT(IG)
        I = IXLG(IP,IG)
        K = KXLT(IP,IG)
        XLON = REAL(I-1)*ZDELLO(K)

        IF (K.GT.1) THEN
!         CLOSEST GRID POINT IN SW CORNER
          XL=XLON-XDELLA
          XMIN = XL/ZDELLO(K-1)
          IMIN = NINT(XMIN) + 1
          IF (IPER.EQ.1 .AND. IMIN .LT. 1)  THEN
            IMIN = IMIN + NLONRGG(K-1)
            XMIN = XMIN + FLOAT(NLONRGG(K-1))
          ENDIF
          IF (IMIN .GE. 1)  THEN
            IF (BATHY(IMIN,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN .AND.
     &              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
             KRLON(IP,1) = IH
            END IF
          ENDIF

!         SECOND CLOSEST GRID POINT IN SW CORNER
          IF(IRGG.EQ.1) THEN
            IF (IMIN .GE. 1)  THEN
              IF(XMIN.LE.FLOAT(IMIN-1)) THEN
                IF(IMIN.LE.1) THEN
                  IMIN2=1
                ELSE
                  IMIN2=IMIN-1
                ENDIF
              ELSE
                IF(IMIN.GE.NLONRGG(K-1)) THEN
                  IMIN2=NLONRGG(K-1)
                ELSE
                  IMIN2=IMIN+1
                ENDIF
              ENDIF

              IF (BATHY(IMIN2,K-1).GT.-990) THEN
                DO IH = IP,1,-1
                  IF (IXLG(IH,IG).EQ.IMIN2 .AND.
     1                KXLT(IH,IG).EQ.K-1) EXIT
                ENDDO
                KRLON2(IP,1) = IH
              ENDIF
            ENDIF
          ELSE
            KRLON2(IP,1) = KRLON(IP,1)
          ENDIF

!         CLOSEST GRID POINT IN SE CORNER
          XL=XLON+XDELLA
          XMIN = XL/ZDELLO(K-1)
          IMIN = NINT(XMIN) + 1
          IF (IPER.EQ.1 .AND. IMIN .GT. NLONRGG(K-1))  THEN
            IMIN = IMIN - NLONRGG(K-1)
            XMIN = XMIN - FLOAT(NLONRGG(K-1))
          ENDIF
          IF (IMIN .LE. NLONRGG(K-1))  THEN
            IF (BATHY(IMIN,K-1).GT.-990) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN .AND.
     &              KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
             KRLAT(IP,1) = IH
            END IF
          ENDIF

!         SECOND CLOSEST GRID POINT IN SE CORNER
          IF(IRGG.EQ.1) THEN
            IF (IMIN .LE. NLONRGG(K-1))  THEN
              IF(XMIN.LE.FLOAT(IMIN-1)) THEN
                IF(IMIN.LE.1) THEN
                  IMIN2=1
                ELSE
                  IMIN2=IMIN-1
                ENDIF
              ELSE
                IF(IMIN.GE.NLONRGG(K-1)) THEN
                  IMIN2=NLONRGG(K-1)
                ELSE
                  IMIN2=IMIN+1
                ENDIF
              ENDIF

              IF (BATHY(IMIN2,K-1).GT.-990) THEN
                DO IH = IP,1,-1
                  IF (IXLG(IH,IG).EQ.IMIN2 .AND.
     1                KXLT(IH,IG).EQ.K-1) EXIT
                ENDDO
                KRLAT2(IP,1) = IH
              ENDIF
            ENDIF
          ELSE
            KRLAT2(IP,1) = KRLAT(IP,1)
          ENDIF
        ENDIF


        IF (K.LT.NY) THEN
!         CLOSEST GRID POINT IN NW CORNER 
          XL=XLON-XDELLA
          XPLUS = XL/ZDELLO(K+1)
          IPLUS = NINT(XPLUS) + 1
          IF (IPER.EQ.1 .AND. IPLUS .LT. 1)  THEN
            IPLUS = NLONRGG(K+1) + IPLUS 
            XPLUS = XPLUS + FLOAT(NLONRGG(K+1))
          ENDIF

          IF (IPLUS .GE. 1)  THEN
            IF (BATHY(IPLUS,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS .AND.
     &              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KRLAT(IP,2) = IH
            ENDIF
          ENDIF

!         SECOND CLOSEST GRID POINT IN NW CORNER 
          IF(IRGG.EQ.1) THEN
            IF (IPLUS .GE. 1)  THEN
              IF(XPLUS.LE.FLOAT(IPLUS-1)) THEN
                IF(IPLUS.LE.1) THEN
                  IPLUS2=1
                ELSE
                  IPLUS2=IPLUS-1
                ENDIF
              ELSE
                IF(IPLUS.GE.NLONRGG(K+1)) THEN
                  IPLUS2=NLONRGG(K+1)
                ELSE
                  IPLUS2=IPLUS+1
                ENDIF
              ENDIF
 
              IF (BATHY(IPLUS2,K+1).GT.-990) THEN
                DO IH = IP,IJLT(IG)
                  IF (IXLG(IH,IG).EQ.IPLUS2 .AND.
     1                KXLT(IH,IG).EQ.K+1) EXIT
                ENDDO
                KRLAT2(IP,2) = MIN(IH,IJLT(IG))
              ENDIF
            ENDIF
          ELSE
            KRLAT2(IP,2) = KRLAT(IP,2)
          ENDIF

!         CLOSEST GRID POINT IN NE CORNER
          XL=XLON+XDELLA
          XPLUS = XL/ZDELLO(K+1)
          IPLUS = NINT(XPLUS) + 1
          IF (IPER.EQ.1 .AND. IPLUS .GT. NLONRGG(K+1))  THEN
            IPLUS = IPLUS - NLONRGG(K+1)
            XPLUS = XPLUS - FLOAT(NLONRGG(K+1))
          ENDIF

          IF (IPLUS .LE. NLONRGG(K+1))  THEN
            IF (BATHY(IPLUS,K+1).GT.-990) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS .AND.
     &              KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KRLON(IP,2) = IH
            ENDIF
          ENDIF

!         SECOND CLOSEST GRID POINT IN NE CORNER 
          IF(IRGG.EQ.1) THEN
            IF (IPLUS .LE. NLONRGG(K+1))  THEN
              IF(XPLUS.LE.FLOAT(IPLUS-1)) THEN
                IF(IPLUS.LE.1) THEN
                  IPLUS2=1
                ELSE
                  IPLUS2=IPLUS-1
                ENDIF
              ELSE
                IF(IPLUS.GE.NLONRGG(K+1)) THEN
                  IPLUS2=NLONRGG(K+1)
                ELSE
                  IPLUS2=IPLUS+1
                ENDIF
              ENDIF

              IF (BATHY(IPLUS2,K+1).GT.-990) THEN
                DO IH = IP,IJLT(IG)
                  IF (IXLG(IH,IG).EQ.IPLUS2 .AND.
     1                KXLT(IH,IG).EQ.K+1) EXIT
                ENDDO
                KRLON2(IP,2) = MIN(IH,IJLT(IG))
              ENDIF
            ENDIF
          ELSE
            KRLON2(IP,2) = KRLON(IP,2)
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
!           SOUTHERN POINT 
            XMIN = D0/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            XP=FLOAT(IMIN-1)*ZDELLO(K-1)
            D4=XP-0.5*ZDELLO(K-1)
            D6=XP+0.5*ZDELLO(K-1)
            IF(D0.LE.XP) THEN
              IF(D4.LE.D3.OR.D6.LE.D5) THEN
                 WLAT(IP,1) = 1.
              ELSE
                 D2=D4-D3
                 D1=ZDELLO(K)-D2
                 WLAT(IP,1)=MIN(1.,D1/ZDELLO(K))
              ENDIF
            ELSE
              IF(D4.GE.D3.OR.D6.GE.D5) THEN
                 WLAT(IP,1) = 1.
              ELSE
                 D2=D5-D6
                 D1=ZDELLO(K)-D2
                 WLAT(IP,1)=MIN(1.,D1/ZDELLO(K))
              ENDIF
            ENDIF

!           SOUTH-WEST CORNER 
            XL=D0-XDELLA
            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            XP=FLOAT(IMIN-1)*ZDELLO(K-1)
            D4=XP-0.5*ZDELLO(K-1)
            D6=XP+0.5*ZDELLO(K-1)
            IF(XL.LE.XP) THEN
              IF(D4.LE.D0-TWOXDELLA) THEN
                 WRLON(IP,1)=1.
              ELSE IF(D6.LE.D0) THEN
                 D2=XP-XL
                 D1=ZDELLO(K-1)-D2
                 WRLON(IP,1)=MIN(1.,D1/ZDELLO(K-1))
              ELSE
                 D2=D4-(D0-TWOXDELLA)
                 D1=TWOXDELLA-D2
                 WRLON(IP,1)=MIN(1.,D1/TWOXDELLA)
              ENDIF
            ELSE
              IF(D6.GE.D0) THEN
                 WRLON(IP,1)=1.
              ELSE IF(D4.GE.D0-TWOXDELLA) THEN
                 D2=XL-XP
                 D1=ZDELLO(K-1)-D2
                 WRLON(IP,1)=MIN(1.,D1/ZDELLO(K-1))
              ELSE
                 D2=D0-D6
                 D1=TWOXDELLA-D2
                 WRLON(IP,1)=MIN(1.,D1/TWOXDELLA)
              ENDIF
            ENDIF

!           SOUTH-EAST CORNER 
            XL=D0+XDELLA
            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            XP=FLOAT(IMIN-1)*ZDELLO(K-1)
            D4=XP-0.5*ZDELLO(K-1)
            D6=XP+0.5*ZDELLO(K-1)
            IF(XL.LE.XP) THEN
              IF(D4.LE.D0) THEN
                 WRLAT(IP,1)=1.
              ELSE IF(D6.LE.D0+TWOXDELLA) THEN
                 D2=XP-XL
                 D1=ZDELLO(K-1)-D2
                 WRLAT(IP,1)=MIN(1.,D1/ZDELLO(K-1))
              ELSE
                 D2=D4-D0
                 D1=TWOXDELLA-D2
                 WRLAT(IP,1)=MIN(1.,D1/TWOXDELLA)
              ENDIF
            ELSE
              IF(D6.GE.D0+TWOXDELLA) THEN
                 WRLAT(IP,1)=1.
              ELSE IF(D4.GE.D0) THEN
                 D2=XL-XP
                 D1=ZDELLO(K-1)-D2
                 WRLAT(IP,1)=MIN(1.,D1/ZDELLO(K-1))
              ELSE
                 D2=D0+TWOXDELLA-D6
                 D1=TWOXDELLA-D2
                 WRLAT(IP,1)=MIN(1.,D1/TWOXDELLA)
              ENDIF
            ENDIF

          ENDIF


          IF (K.LT.NY) THEN
!           NORTHERN POINT 
            XPLUS = D0/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            XP=FLOAT(IPLUS-1)*ZDELLO(K+1)
            D4=XP-0.5*ZDELLO(K+1)
            D6=XP+0.5*ZDELLO(K+1)
            IF(D0.LE.XP) THEN
              IF(D4.LE.D3.OR.D6.LE.D5) THEN
                 WLAT(IP,2) = 1.
              ELSE
                 D2=D4-D3
                 D1=ZDELLO(K)-D2
                 WLAT(IP,2)=MIN(1.,D1/ZDELLO(K))
              ENDIF
            ELSE
              IF(D4.GE.D3.OR.D6.GE.D5) THEN
                 WLAT(IP,2) = 1.
              ELSE
                 D2=D5-D6
                 D1=ZDELLO(K)-D2
                 WLAT(IP,2)=MIN(1.,D1/ZDELLO(K))
              ENDIF
            ENDIF

!           NORTH-WEST CORNER 
            XL=D0-XDELLA
            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            XP=FLOAT(IPLUS-1)*ZDELLO(K+1)
            D4=XP-0.5*ZDELLO(K+1)
            D6=XP+0.5*ZDELLO(K+1)
            IF(XL.LE.XP) THEN
              IF(D4.LE.D0-TWOXDELLA) THEN
                 WRLAT(IP,2)=1.
              ELSE IF(D6.LE.D0) THEN
                 D2=XP-XL
                 D1=ZDELLO(K+1)-D2
                 WRLAT(IP,2)=MIN(1.,D1/ZDELLO(K+1))
              ELSE
                 D2=D4-(D0-TWOXDELLA)
                 D1=TWOXDELLA-D2
                 WRLAT(IP,2)=MIN(1.,D1/TWOXDELLA)
              ENDIF
            ELSE
              IF(D6.GE.D0) THEN
                 WRLAT(IP,2)=1.
              ELSE IF(D4.GE.D0-TWOXDELLA) THEN
                 D2=XL-XP
                 D1=ZDELLO(K+1)-D2
                 WRLAT(IP,2)=MIN(1.,D1/ZDELLO(K+1))
              ELSE
                 D2=D0-D6
                 D1=TWOXDELLA-D2
                 WRLAT(IP,2)=MIN(1.,D1/TWOXDELLA)
              ENDIF
            ENDIF

!           NORTH-EAST CORNER 
            XL=D0+XDELLA
            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            XP=FLOAT(IPLUS-1)*ZDELLO(K+1)
            D4=XP-0.5*ZDELLO(K+1)
            D6=XP+0.5*ZDELLO(K+1)
            IF(XL.LE.XP) THEN
              IF(D4.LE.D0) THEN
                 WRLON(IP,2)=1.
              ELSE IF(D6.LE.D0+TWOXDELLA) THEN
                 D2=XP-XL
                 D1=ZDELLO(K+1)-D2
                 WRLON(IP,2)=MIN(1.,D1/ZDELLO(K+1))
              ELSE
                 D2=D4-D0
                 D1=TWOXDELLA-D2
                 WRLON(IP,2)=MIN(1.,D1/TWOXDELLA)
              ENDIF
            ELSE
              IF(D6.GE.D0+TWOXDELLA) THEN
                 WRLON(IP,2)=1.
              ELSE IF(D4.GE.D0) THEN
                 D2=XL-XP
                 D1=ZDELLO(K+1)-D2
                 WRLON(IP,2)=MIN(1.,D1/ZDELLO(K+1))
              ELSE
                 D2=D0+TWOXDELLA-D6
                 D1=TWOXDELLA-D2
                 WRLON(IP,2)=MIN(1.,D1/TWOXDELLA)
              ENDIF
            ENDIF

          ENDIF

        ENDDO
      ENDIF

      IF(IRGG.EQ.1) THEN
        LLABORT=.FALSE.
        DO IS=1,2
          DO IP = 1,IJLT(IG)
            IF(WLAT(IP,IS).LT.0. .OR. WLAT(IP,IS).GT.1.) THEN
              WRITE(IU06,*) ' WLAT < 0 or > 1 !!!! ',IP,IS,WLAT(IP,IS)
              LLABORT=.TRUE.
            ENDIF
            IF(WRLON(IP,IS).LT.0. .OR. WRLON(IP,IS).GT.1.) THEN
              WRITE(IU06,*) ' WRLON < 0 or > 1 !!!! ',IP,IS,WRLON(IP,IS)
              LLABORT=.TRUE.
            ENDIF
            IF(WRLAT(IP,IS).LT.0. .OR. WRLAT(IP,IS).GT.1.) THEN
              WRITE(IU06,*) ' WRLAT < 0 or > 1 !!!! ',IP,IS,WRLAT(IP,IS)
              LLABORT=.TRUE.
            ENDIF
          ENDDO
        ENDDO
        IF(LLABORT) THEN
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) ' PROBLEM IN MUBUF  (see above) !!!! '
          WRITE(IU06,*) ' ' 
          CALL ABORT1
        ENDIF
      ENDIF

      IF(LLOBSTRCT) THEN
!     DETERMINE THE REAL REDUCTION FACTORS TO BE USED IN THE PROPAGATION

        READ (IU01,'(I4)') NFREMAX

        IF(NFREMAX.NE.ML ) THEN
          WRITE (IU06,*) ' *******************************************'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *      FATAL  ERROR IN SUB. MUBUF         *'
          WRITE (IU06,*) ' * NFREMAX MUST BE = ML                    *'
          WRITE (IU06,*) ' * NFREMAX = ',NFREMAX
          WRITE (IU06,*) ' * ML      = ',ML 
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' * PROGRAM WILL BE ABORTED                 *'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF

        DO M=1,NFREMAX
          DO IS=1,2
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT) (IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KOBSLAT(IP,M,IS)=IDUM(I,K)
            ENDDO
          ENDDO

          DO IS=1,2
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT) (IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KOBSLON(IP,M,IS)=IDUM(I,K)
            ENDDO
          ENDDO

          DO IS=1,2
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT)(IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KOBSRLAT(IP,M,IS)=IDUM(I,K)
            ENDDO
          ENDDO

          DO IS=1,2
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT)(IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KOBSRLON(IP,M,IS)=IDUM(I,K)
            ENDDO
          ENDDO

          DO IS=1,4
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT) (IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KOBSCOR(IP,M,IS)=IDUM(I,K)
            ENDDO
          ENDDO

        ENDDO ! END LOOP ON NFREMAX

      ELSE

        DO IS=1,2
          DO M=1,ML
            DO IP = 1,IJLT(IG)
              KOBSLAT(IP,M,IS)=1000
              KOBSLON(IP,M,IS)=1000
              KOBSRLAT(IP,M,IS)=1000
              KOBSRLON(IP,M,IS)=1000
            ENDDO
          ENDDO
        ENDDO

        DO IS=1,4
          DO M=1,ML
            DO IP = 1,IJLT(IG)
              KOBSCOR(IP,M,IS)=1000
            ENDDO
          ENDDO
        ENDDO

      ENDIF

! ----------------------------------------------------------------------

!*    3. OUTPUT OF COMMON UBUF FOR ONE BLOCK.
!        ------------------------------------
      CALL OUTUBUF (IU08, IU18, IFORM,NIBLO,ML,
     &              KOBSLAT,KOBSLON,KOBSCOR,KOBSRLAT,KOBSRLON)

      RETURN
      END SUBROUTINE MUBUF
