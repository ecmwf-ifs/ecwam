      SUBROUTINE MUBUF (IU01, BATHY, IG, ML, IU08, NPROPAGS)
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
!                          ECMWF       MODIFIED TO COMPUTE AND WRITE ALL
!                                      ARRAYS SEQUENTIALLY.

!*    PURPOSE.
!     -------

!       TO ARRANGE NEIGHBOUR GRID POINT INDICES FOR A BLOCK

!**   INTERFACE.
!     ----------

!       *CALL* *MUBUF (BATHY,IG,ML,IU08)*
!          *IU01*  -  LOGICAL INPUT UNIT OF TOPOGRAPHIC DATA.
!          *BATHY*     -  BATHYMETRY DATA.
!          *IG*      - BLOCK NUMBER.
!          *ML*      - NUMBER OF FREQUENCIES.
!          *IU08*    - LOGICAL UNITS FOR OUTPUT OF GRID BLOCKING
!                      COMMON UBUF (UNFORMATED)

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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWGRID  , ONLY : NLONRGG  ,IJLT
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,NY       ,IPER     ,    &
     &            XDELLA   ,ZDELLO   ,IRGG     ,LLOBSTRCT
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KCOR     ,              &
     &                      KRLAT    ,KRLON    ,                        &
     &                      WLAT     ,WCOR     ,WRLAT    ,WRLON

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM) :: IU01
      INTEGER(KIND=JWIM) :: NPROPAGS
      INTEGER(KIND=JWIM) :: NFREMAX, IX
      INTEGER(KIND=JWIM) :: IG, IJ, IJP, I, K, IP, IH, IS, ML, M
      INTEGER(KIND=JWIM) :: IMIN, IPLUS, IMIN2, IPLUS2, IMIN3, IPLUS3
      INTEGER(KIND=JWIM) :: IC, ICP, ICL, ICR
      INTEGER(KIND=JWIM) :: IU08(0:NPROPAGS)
      INTEGER(KIND=JWIM), DIMENSION(NGX, NGY) :: IDUM
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: KDUM
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: KOBSLON, KOBSLAT
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: KOBSCOR
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:,:) :: KOBSRLON, KOBSRLAT

      REAL(KIND=JWRB) :: XMIN, XPLUS, D0, D1, D2, D3, D4, D5, D6
      REAL(KIND=JWRB) :: XLON, XL, XP, TWOXDELLA
      REAL(KIND=JWRB) :: XLL, XPL, XLR, XPR
      REAL(KIND=JWRB) ::  BATHY(NGX, NGY)

      CHARACTER(LEN=4) :: CX
      CHARACTER(LEN=10) :: FORMAT

      LOGICAL :: LLABORT

! ----------------------------------------------------------------------

      TWOXDELLA=2.0_JWRB*XDELLA

!*    2. COMPUTE INDICES OF NEIGHBOUR SEA POINTS.
!        ----------------------------------------

!*    2.1 LATITUDE NEIGHBOURS (KLAT)
!         --------------------------

      ALLOCATE(KLAT(NIBLO,2,2))
      DO ICL=1,2
        DO IJP=1,2
          DO IJ=1,NIBLO
             KLAT(IJ,IJP,ICL) = 0
          ENDDO
        ENDDO
      ENDDO

      DO IP = 1,IJLT(IG)
        I = IXLG(IP,IG)
        K = KXLT(IP,IG)
        IF (K.GT.1) THEN
          XMIN = REAL(I-1)*ZDELLO(K)/ZDELLO(K-1)
          IMIN = NINT(XMIN) + 1

!         CLOSEST GRID POINT
          IF (BATHY(IMIN,K-1).GT.-990.0_JWRB) THEN
            DO IH = IP,1,-1
              IF (IXLG(IH,IG).EQ.IMIN .AND. KXLT(IH,IG).EQ.K-1) EXIT 
            ENDDO
            KLAT(IP,1,1) = IH
          ENDIF

!         SECOND CLOSEST GRID POINT
          IF (IRGG.EQ.1) THEN
            IF (XMIN.LE.FLOAT(IMIN-1)) THEN
              IF (IMIN.LE.1) THEN
                IMIN2=1
              ELSE
                IMIN2=IMIN-1
              ENDIF
            ELSE
              IF (IMIN.GE.NLONRGG(K-1)) THEN
                IMIN2=NLONRGG(K-1)
              ELSE
                IMIN2=IMIN+1
              ENDIF
            ENDIF

            IF (BATHY(IMIN2,K-1).GT.-990.0_JWRB) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN2 .AND. KXLT(IH,IG).EQ.K-1) EXIT 
              ENDDO
              KLAT(IP,1,2) = IH
            ENDIF
          ELSE
            KLAT(IP,1,2) = KLAT(IP,1,1)
          ENDIF

        ENDIF

        IF (K.LT.NY) THEN
          XPLUS = REAL(I-1)*ZDELLO(K)/ZDELLO(K+1)
          IPLUS = NINT(XPLUS) + 1
          IF (BATHY(IPLUS,K+1).GT.-990.0_JWRB) THEN
            DO IH = IP,IJLT(IG)
              IF (IXLG(IH,IG).EQ.IPLUS .AND. KXLT(IH,IG).EQ.K+1) EXIT 
            ENDDO
            KLAT(IP,2,1) = IH
          ENDIF

          IF (IRGG.EQ.1) THEN
            IF (XPLUS.LE.FLOAT(IPLUS-1)) THEN
              IF (IPLUS.LE.1) THEN
                IPLUS2=1
              ELSE
                IPLUS2=IPLUS-1
              ENDIF
            ELSE
              IF (IPLUS.GE.NLONRGG(K+1)) THEN
                IPLUS2=NLONRGG(K+1)
              ELSE
                IPLUS2=IPLUS+1
              ENDIF
            ENDIF

            IF (BATHY(IPLUS2,K+1).GT.-990.0_JWRB) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS2 .AND.  KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KLAT(IP,2,2) = MIN(IH,IJLT(IG))
            ENDIF
          ELSE
            KLAT(IP,2,2) = KLAT(IP,2,1)
          ENDIF
        ENDIF

      ENDDO

      DO IJP=0,NPROPAGS
        WRITE (IU08(IJP)) KLAT
      ENDDO

      DEALLOCATE(KLAT)

!*    2.2 LONGITUDE NEIGHBOURS (KLON)
!         ---------------------------

      ALLOCATE(KLON(NIBLO,2))
      DO IJP=1,2
        DO IJ=1,NIBLO
           KLON(IJ,IJP) = 0
        ENDDO
      ENDDO

      DO IP = 1,IJLT(IG)
        I = IXLG(IP,IG)
        K = KXLT(IP,IG)
        IF (I.GT.1) THEN
          IF (BATHY(I-1,K).GT.-990.0_JWRB) KLON(IP,1) = IP-1
        ELSE
          IF (IPER.EQ.1 .AND. BATHY(NLONRGG(K),K).GT.-990.0_JWRB) THEN
            KLON(IP,1) = IP
            DO IH=2,NLONRGG(K)
              IF (BATHY(IH,K).GT.-990.0_JWRB) KLON(IP,1) = KLON(IP,1)+1
            ENDDO
          ENDIF
        ENDIF
        IF (I.LT.NLONRGG(K)) THEN
          IF (BATHY(I+1,K).GT.-990.0_JWRB) KLON(IP,2) = IP+1
        ELSE
          IF (IPER.EQ.1 .AND. BATHY(1,K).GT.-990.0_JWRB) THEN
            KLON(IP,2) = IP
            DO IH=NLONRGG(K)-1,1,-1
              IF (BATHY(IH,K).GT.-990.0_JWRB) KLON(IP,2) = KLON(IP,2)-1
            ENDDO
          ENDIF
        ENDIF
      ENDDO

      DO IJP=0,NPROPAGS
        WRITE (IU08(IJP)) KLON
      ENDDO

      DEALLOCATE(KLON)


!     COMPUTE THE CORNER GRID POINT AND THE CLOSEST GRID POINT ON EITHER
!     SIDE (ON A GIVEN LATITUDE) (KCOR)
!     ------------------------------------------------------------------

      ALLOCATE(KCOR(NIBLO,4,2))
      DO IC=1,2
        DO ICR=1,4
          DO IJ=1,NIBLO
             KCOR(IJ,ICR,IC) = 0
          ENDDO
        ENDDO
      ENDDO

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
            IF (BATHY(IMIN,K-1).GT.-990.0_JWRB) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN .AND. KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
             KCOR(IP,3,1) = IH
            END IF
          ENDIF

!         SECOND CLOSEST GRID POINT IN SW GRID CORNER POINT
          IF (IMIN .GE. 1)  THEN
            IF (XMIN.LE.FLOAT(IMIN-1)) THEN
              IMIN2=IMIN-1
              IF (IMIN.LE.1) THEN
!!test                IMIN2=1
                IMIN2=NLONRGG(K-1)
              ELSE
                IMIN2=IMIN-1
              ENDIF
            ELSE
               IMIN2=IMIN+1
              IF (IMIN.GE.NLONRGG(K-1)) THEN
!!test                IMIN2=NLONRGG(K-1)
                IMIN2=1
              ELSE
                IMIN2=IMIN+1
              ENDIF
            ENDIF

            IF (BATHY(IMIN2,K-1).GT.-990.0_JWRB) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN2 .AND. KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
              KCOR(IP,3,2) = IH
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
            IF (BATHY(IMIN,K-1).GT.-990.0_JWRB) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN .AND. KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
             KCOR(IP,2,1) = IH
            END IF
          ENDIF

!         SECOND CLOSEST GRID POINT IN SE GRID CORNER
          IF (IMIN .LE. NLONRGG(K-1))  THEN
            IF (XMIN.LE.FLOAT(IMIN-1)) THEN
               IMIN2=IMIN-1
              IF (IMIN.LE.1) THEN
!!test                IMIN2=1
                IMIN2=NLONRGG(K-1)
              ELSE
                IMIN2=IMIN-1
              ENDIF
            ELSE
               IMIN2=IMIN+1
              IF (IMIN.GE.NLONRGG(K-1)) THEN
!!test                IMIN2=NLONRGG(K-1)
                IMIN2=1
              ELSE
                IMIN2=IMIN+1
              ENDIF
            ENDIF

            IF (BATHY(IMIN2,K-1).GT.-990.0_JWRB) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN2 .AND. KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
              KCOR(IP,2,2) = IH
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
            IF (BATHY(IPLUS,K+1).GT.-990.0_JWRB) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS .AND. KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,4,1) = IH
            ENDIF
          ENDIF

!         SECOND CLOSEST GRID POINT IN NW GRID CORNER 
          IF (IPLUS .GE. 1)  THEN
            IF (XPLUS.LE.FLOAT(IPLUS-1)) THEN
               IPLUS2=IPLUS-1
              IF (IPLUS.LE.1) THEN
!!test                IPLUS2=1
                IPLUS2=NLONRGG(K+1)
              ELSE
                IPLUS2=IPLUS-1
              ENDIF
            ELSE
               IPLUS2=IPLUS+1
              IF (IPLUS.GE.NLONRGG(K+1)) THEN
!!test                IPLUS2=NLONRGG(K+1)
                IPLUS2=1
              ELSE
                IPLUS2=IPLUS+1
              ENDIF
            ENDIF
 
            IF (BATHY(IPLUS2,K+1).GT.-990.0_JWRB) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS2 .AND. KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,4,2) = MIN(IH,IJLT(IG))
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
            IF (BATHY(IPLUS,K+1).GT.-990.0_JWRB) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS .AND. KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,1,1) = IH
            ENDIF
          ENDIF

!         SECOND CLOSEST GRID POINT IN NE GRID CORNER 
          IF (IPLUS .LE. NLONRGG(K+1))  THEN
            IF (XPLUS.LE.FLOAT(IPLUS-1)) THEN
               IPLUS2=IPLUS-1
              IF (IPLUS.LE.1) THEN
!!test                IPLUS2=1
                IPLUS2=NLONRGG(K+1)
              ELSE
                IPLUS2=IPLUS-1
              ENDIF
            ELSE
               IPLUS2=IPLUS+1
              IF (IPLUS.GE.NLONRGG(K+1)) THEN
!!test                IPLUS2=NLONRGG(K+1)
                IPLUS2=1
              ELSE
                IPLUS2=IPLUS+1
              ENDIF
            ENDIF

            IF (BATHY(IPLUS2,K+1).GT.-990.0_JWRB) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS2 .AND. KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KCOR(IP,1,2) = MIN(IH,IJLT(IG))
            ENDIF
          ENDIF

        ENDIF  ! END K < NY

      ENDDO


      WRITE (IU08(2)) KCOR

      DEALLOCATE(KCOR)


!     COMPUTE THE CLOSEST GRID POINTS IN THE SW-NE and SE-NW DIRECTIONS
!     (I.E. GOING AT 45 DEGREE) (KRLAT, KRLON)
!     -----------------------------------------------------------------

      ALLOCATE(KRLAT(NIBLO,2,2))
      ALLOCATE(KRLON(NIBLO,2,2))
      DO ICL=1,2
        DO IJP=1,2
          DO IJ=1,NIBLO
             KRLAT(IJ,IJP,ICL) = 0
             KRLON(IJ,IJP,ICL) = 0
          ENDDO
        ENDDO
      ENDDO

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
            IF (BATHY(IMIN,K-1).GT.-990.0_JWRB) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN .AND. KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
             KRLON(IP,1,1) = IH
            END IF
          ENDIF

!         SECOND CLOSEST GRID POINT IN SW CORNER
          IF (IRGG.EQ.1) THEN
            IF (IMIN .GE. 1)  THEN
              IF (XMIN.LE.FLOAT(IMIN-1)) THEN
                IF (IMIN.LE.1) THEN
                  IMIN2=1
                ELSE
                  IMIN2=IMIN-1
                ENDIF
              ELSE
                IF (IMIN.GE.NLONRGG(K-1)) THEN
                  IMIN2=NLONRGG(K-1)
                ELSE
                  IMIN2=IMIN+1
                ENDIF
              ENDIF

              IF (BATHY(IMIN2,K-1).GT.-990.0_JWRB) THEN
                DO IH = IP,1,-1
                  IF (IXLG(IH,IG).EQ.IMIN2 .AND. KXLT(IH,IG).EQ.K-1) EXIT
                ENDDO
                KRLON(IP,1,2) = IH
              ENDIF
            ENDIF
          ELSE
            KRLON(IP,1,2) = KRLON(IP,1,1)
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
            IF (BATHY(IMIN,K-1).GT.-990.0_JWRB) THEN
              DO IH = IP,1,-1
                IF (IXLG(IH,IG).EQ.IMIN .AND. KXLT(IH,IG).EQ.K-1) EXIT
              ENDDO
             KRLAT(IP,1,1) = IH
            END IF
          ENDIF

!         SECOND CLOSEST GRID POINT IN SE CORNER
          IF (IRGG.EQ.1) THEN
            IF (IMIN .LE. NLONRGG(K-1))  THEN
              IF (XMIN.LE.FLOAT(IMIN-1)) THEN
                IF (IMIN.LE.1) THEN
                  IMIN2=1
                ELSE
                  IMIN2=IMIN-1
                ENDIF
              ELSE
                IF (IMIN.GE.NLONRGG(K-1)) THEN
                  IMIN2=NLONRGG(K-1)
                ELSE
                  IMIN2=IMIN+1
                ENDIF
              ENDIF

              IF (BATHY(IMIN2,K-1).GT.-990.0_JWRB) THEN
                DO IH = IP,1,-1
                  IF (IXLG(IH,IG).EQ.IMIN2 .AND. KXLT(IH,IG).EQ.K-1) EXIT
                ENDDO
                KRLAT(IP,1,2) = IH
              ENDIF
            ENDIF
          ELSE
            KRLAT(IP,1,2) = KRLAT(IP,1,1)
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
            IF (BATHY(IPLUS,K+1).GT.-990.0_JWRB) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS .AND. KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KRLAT(IP,2,1) = IH
            ENDIF
          ENDIF

!         SECOND CLOSEST GRID POINT IN NW CORNER 
          IF (IRGG.EQ.1) THEN
            IF (IPLUS .GE. 1)  THEN
              IF (XPLUS.LE.FLOAT(IPLUS-1)) THEN
                IF (IPLUS.LE.1) THEN
                  IPLUS2=1
                ELSE
                  IPLUS2=IPLUS-1
                ENDIF
              ELSE
                IF (IPLUS.GE.NLONRGG(K+1)) THEN
                  IPLUS2=NLONRGG(K+1)
                ELSE
                  IPLUS2=IPLUS+1
                ENDIF
              ENDIF
 
              IF (BATHY(IPLUS2,K+1).GT.-990.0_JWRB) THEN
                DO IH = IP,IJLT(IG)
                  IF (IXLG(IH,IG).EQ.IPLUS2 .AND. KXLT(IH,IG).EQ.K+1) EXIT
                ENDDO
                KRLAT(IP,2,2) = MIN(IH,IJLT(IG))
              ENDIF
            ENDIF
          ELSE
            KRLAT(IP,2,2) = KRLAT(IP,2,1)
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
            IF (BATHY(IPLUS,K+1).GT.-990.0_JWRB) THEN
              DO IH = IP,IJLT(IG)
                IF (IXLG(IH,IG).EQ.IPLUS .AND. KXLT(IH,IG).EQ.K+1) EXIT
              ENDDO
              KRLON(IP,2,1) = IH
            ENDIF
          ENDIF

!         SECOND CLOSEST GRID POINT IN NE CORNER 
          IF (IRGG.EQ.1) THEN
            IF (IPLUS .LE. NLONRGG(K+1))  THEN
              IF (XPLUS.LE.FLOAT(IPLUS-1)) THEN
                IF (IPLUS.LE.1) THEN
                  IPLUS2=1
                ELSE
                  IPLUS2=IPLUS-1
                ENDIF
              ELSE
                IF (IPLUS.GE.NLONRGG(K+1)) THEN
                  IPLUS2=NLONRGG(K+1)
                ELSE
                  IPLUS2=IPLUS+1
                ENDIF
              ENDIF

              IF (BATHY(IPLUS2,K+1).GT.-990.0_JWRB) THEN
                DO IH = IP,IJLT(IG)
                  IF (IXLG(IH,IG).EQ.IPLUS2 .AND. KXLT(IH,IG).EQ.K+1) EXIT
                ENDDO
                KRLON(IP,2,2) = MIN(IH,IJLT(IG))
              ENDIF
            ENDIF
          ELSE
            KRLON(IP,2,2) = KRLON(IP,2,1)
          ENDIF

        ENDIF

      ENDDO

      WRITE (IU08(1)) KRLAT
      WRITE (IU08(1)) KRLON

      DEALLOCATE(KRLAT)
      DEALLOCATE(KRLON)
 

!     COMPUTE THE WEIGHT FOR THE ADVECTION SCHEME ON AN IRREGULAR GRID
!     ()

      ALLOCATE(WLAT(NIBLO,2))
      ALLOCATE(WRLAT(NIBLO,2))
      ALLOCATE(WRLON(NIBLO,2))
      ALLOCATE(WCOR(NIBLO,4))

      DO IJP=1,2
        DO IJ=1,NIBLO
           WLAT(IJ,IJP)=1.0_JWRB
           WRLAT(IJ,IJP)=1.0_JWRB
           WRLON(IJ,IJP)=1.0_JWRB
        ENDDO
      ENDDO
      DO ICR=1,4
        DO IJ=1,NIBLO
           WCOR(IJ,ICR)=1.0_JWRB
        ENDDO
      ENDDO

      IF (IRGG.EQ.1) THEN
        DO IP = 1,IJLT(IG)
          I = IXLG(IP,IG)
          K = KXLT(IP,IG)
          D0 = FLOAT(I-1)*ZDELLO(K)
          D3=D0-0.5_JWRB*ZDELLO(K)
          D5=D0+0.5_JWRB*ZDELLO(K)

          IF (K.GT.1) THEN
!           SOUTHERN POINT 
            XMIN = D0/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            XP=FLOAT(IMIN-1)*ZDELLO(K-1)
            D4=XP-0.5_JWRB*ZDELLO(K-1)
            D6=XP+0.5_JWRB*ZDELLO(K-1)
            IF (D0.LE.XP) THEN
              IF (D4.LE.D3 .OR. D6.LE.D5) THEN
                 WLAT(IP,1) = 1.0_JWRB
              ELSE
                 D2=D4-D3
                 D1=ZDELLO(K)-D2
                 WLAT(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K))
              ENDIF
            ELSE
              IF (D4.GE.D3 .OR. D6.GE.D5) THEN
                 WLAT(IP,1) = 1.0_JWRB
              ELSE
                 D2=D5-D6
                 D1=ZDELLO(K)-D2
                 WLAT(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K))
              ENDIF
            ENDIF

!           SOUTH-WEST CORNER 
!           FOR DUAL ROTATED SCHEME
            XL=D0-XDELLA
            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            XP=FLOAT(IMIN-1)*ZDELLO(K-1)
            D4=XP-0.5_JWRB*ZDELLO(K-1)
            D6=XP+0.5_JWRB*ZDELLO(K-1)
            IF (XL.LE.XP) THEN
              IF (D4.LE.D0-TWOXDELLA) THEN
                 WRLON(IP,1)=1.0_JWRB
              ELSEIF (D6.LE.D0) THEN
                 D2=XP-XL
                 D1=ZDELLO(K-1)-D2
                 WRLON(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K-1))
              ELSE
                 D2=D4-(D0-TWOXDELLA)
                 D1=TWOXDELLA-D2
                 WRLON(IP,1)=MIN(1.0_JWRB,D1/TWOXDELLA)
              ENDIF
            ELSE
              IF (D6.GE.D0) THEN
                 WRLON(IP,1)=1.0_JWRB
              ELSE IF (D4.GE.D0-TWOXDELLA) THEN
                 D2=XL-XP
                 D1=ZDELLO(K-1)-D2
                 WRLON(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K-1))
              ELSE
                 D2=D0-D6
                 D1=TWOXDELLA-D2
                 WRLON(IP,1)=MIN(1.0_JWRB,D1/TWOXDELLA)
              ENDIF
            ENDIF

!           FOR CORNER TRANSPORT SCHEME
!           SOUTH-WEST CORNER 
            XL=D0-ZDELLO(K)
            XLL=XL-0.5_JWRB*ZDELLO(K)
            XLR=XL+0.5_JWRB*ZDELLO(K)

            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            XP=FLOAT(IMIN-1)*ZDELLO(K-1)
            XPL=XP-0.5_JWRB*ZDELLO(K-1)
            XPR=XP+0.5_JWRB*ZDELLO(K-1)

            IF (XPL.GT.XLL .AND. XPR.LT.XLR) THEN
              D1=ZDELLO(K)
            ELSE
              D1=MIN(XLR,XPR)-MAX(XLL,XPL)
            ENDIF
            WCOR(IP,3)=MIN(1.0_JWRB,D1/ZDELLO(K))

!           SOUTH-EAST CORNER 
!           FOR DUAL ROTATED SCHEME
            XL=D0+XDELLA
            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            XP=FLOAT(IMIN-1)*ZDELLO(K-1)
            D4=XP-0.5_JWRB*ZDELLO(K-1)
            D6=XP+0.5_JWRB*ZDELLO(K-1)
            IF (XL.LE.XP) THEN
              IF (D4.LE.D0) THEN
                 WRLAT(IP,1)=1.0_JWRB
              ELSE IF (D6.LE.D0+TWOXDELLA) THEN
                 D2=XP-XL
                 D1=ZDELLO(K-1)-D2
                 WRLAT(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K-1))
              ELSE
                 D2=D4-D0
                 D1=TWOXDELLA-D2
                 WRLAT(IP,1)=MIN(1.0_JWRB,D1/TWOXDELLA)
              ENDIF
            ELSE
              IF (D6.GE.D0+TWOXDELLA) THEN
                 WRLAT(IP,1)=1.0_JWRB
              ELSE IF (D4.GE.D0) THEN
                 D2=XL-XP
                 D1=ZDELLO(K-1)-D2
                 WRLAT(IP,1)=MIN(1.0_JWRB,D1/ZDELLO(K-1))
              ELSE
                 D2=D0+TWOXDELLA-D6
                 D1=TWOXDELLA-D2
                 WRLAT(IP,1)=MIN(1.0_JWRB,D1/TWOXDELLA)
              ENDIF
            ENDIF

!           FOR CORNER TRANSPORT SCHEME
!           SOUTH-EAST CORNER 
            XL=D0+ZDELLO(K)
            XLL=XL-0.5_JWRB*ZDELLO(K)
            XLR=XL+0.5_JWRB*ZDELLO(K)

            XMIN = XL/ZDELLO(K-1)
            IMIN = NINT(XMIN) + 1
            XP=FLOAT(IMIN-1)*ZDELLO(K-1)
            XPL=XP-0.5_JWRB*ZDELLO(K-1)
            XPR=XP+0.5_JWRB*ZDELLO(K-1)

            IF (XPL.GT.XLL .AND. XPR.LT.XLR) THEN
              D1=ZDELLO(K)
            ELSE
              D1=MIN(XLR,XPR)-MAX(XLL,XPL)
            ENDIF
            WCOR(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K))

          ENDIF

          IF (K.LT.NY) THEN
!           NORTHERN POINT 
            XPLUS = D0/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            XP=FLOAT(IPLUS-1)*ZDELLO(K+1)
            D4=XP-0.5_JWRB*ZDELLO(K+1)
            D6=XP+0.5_JWRB*ZDELLO(K+1)
            IF (D0.LE.XP) THEN
              IF (D4.LE.D3.OR.D6.LE.D5) THEN
                 WLAT(IP,2) = 1.0_JWRB
              ELSE
                 D2=D4-D3
                 D1=ZDELLO(K)-D2
                 WLAT(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K))
              ENDIF
            ELSE
              IF (D4.GE.D3 .OR. D6.GE.D5) THEN
                 WLAT(IP,2) = 1.0_JWRB
              ELSE
                 D2=D5-D6
                 D1=ZDELLO(K)-D2
                 WLAT(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K))
              ENDIF
            ENDIF

!           NORTH-WEST CORNER 
!           FOR DUAL ROTATED SCHEME
            XL=D0-XDELLA
            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            XP=FLOAT(IPLUS-1)*ZDELLO(K+1)
            D4=XP-0.5_JWRB*ZDELLO(K+1)
            D6=XP+0.5_JWRB*ZDELLO(K+1)
            IF (XL.LE.XP) THEN
              IF (D4.LE.D0-TWOXDELLA) THEN
                 WRLAT(IP,2)=1.0_JWRB
              ELSE IF (D6.LE.D0) THEN
                 D2=XP-XL
                 D1=ZDELLO(K+1)-D2
                 WRLAT(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K+1))
              ELSE
                 D2=D4-(D0-TWOXDELLA)
                 D1=TWOXDELLA-D2
                 WRLAT(IP,2)=MIN(1.0_JWRB,D1/TWOXDELLA)
              ENDIF
            ELSE
              IF (D6.GE.D0) THEN
                 WRLAT(IP,2)=1.0_JWRB
              ELSE IF (D4.GE.D0-TWOXDELLA) THEN
                 D2=XL-XP
                 D1=ZDELLO(K+1)-D2
                 WRLAT(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K+1))
              ELSE
                 D2=D0-D6
                 D1=TWOXDELLA-D2
                 WRLAT(IP,2)=MIN(1.0_JWRB,D1/TWOXDELLA)
              ENDIF
            ENDIF

!           FOR CORNER TRANSPORT SCHEME
!           NORTH-WEST CORNER 
            XL=D0-ZDELLO(K)
            XLL=XL-0.5_JWRB*ZDELLO(K)
            XLR=XL+0.5_JWRB*ZDELLO(K)

            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            XP=FLOAT(IPLUS-1)*ZDELLO(K+1)
            XPL=XP-0.5_JWRB*ZDELLO(K+1)
            XPR=XP+0.5_JWRB*ZDELLO(K+1)

            IF (XPL.GT.XLL .AND. XPR.LT.XLR) THEN
              D1=ZDELLO(K)
            ELSE
              D1=MIN(XLR,XPR)-MAX(XLL,XPL)
            ENDIF
            WCOR(IP,4)=MIN(1.0_JWRB,D1/ZDELLO(K))


!           NORTH-EAST CORNER 
!           FOR DUAL ROTATED SCHEME
            XL=D0+XDELLA
            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            XP=FLOAT(IPLUS-1)*ZDELLO(K+1)
            D4=XP-0.5_JWRB*ZDELLO(K+1)
            D6=XP+0.5_JWRB*ZDELLO(K+1)
            IF (XL.LE.XP) THEN
              IF (D4.LE.D0) THEN
                 WRLON(IP,2)=1.0_JWRB
              ELSE IF (D6.LE.D0+TWOXDELLA) THEN
                 D2=XP-XL
                 D1=ZDELLO(K+1)-D2
                 WRLON(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K+1))
              ELSE
                 D2=D4-D0
                 D1=TWOXDELLA-D2
                 WRLON(IP,2)=MIN(1.0_JWRB,D1/TWOXDELLA)
              ENDIF
            ELSE
              IF (D6.GE.D0+TWOXDELLA) THEN
                 WRLON(IP,2)=1.0_JWRB
              ELSE IF (D4.GE.D0) THEN
                 D2=XL-XP
                 D1=ZDELLO(K+1)-D2
                 WRLON(IP,2)=MIN(1.0_JWRB,D1/ZDELLO(K+1))
              ELSE
                 D2=D0+TWOXDELLA-D6
                 D1=TWOXDELLA-D2
                 WRLON(IP,2)=MIN(1.0_JWRB,D1/TWOXDELLA)
              ENDIF
            ENDIF

!           FOR CORNER TRANSPORT SCHEME
!           NORTH-EAST CORNER 
            XL=D0+ZDELLO(K)
            XLL=XL-0.5_JWRB*ZDELLO(K)
            XLR=XL+0.5_JWRB*ZDELLO(K)

            XPLUS = XL/ZDELLO(K+1)
            IPLUS = NINT(XPLUS) + 1
            XP=FLOAT(IPLUS-1)*ZDELLO(K+1)
            XPL=XP-0.5_JWRB*ZDELLO(K+1)
            XPR=XP+0.5_JWRB*ZDELLO(K+1)

            IF (XPL.GT.XLL .AND. XPR.LT.XLR) THEN
              D1=ZDELLO(K)
            ELSE
              D1=MIN(XLR,XPR)-MAX(XLL,XPL)
            ENDIF
            WCOR(IP,1)=MIN(1.,D1/ZDELLO(K))

          ENDIF

        ENDDO
      ENDIF

      IF (IRGG.EQ.1) THEN
        LLABORT=.FALSE.
        DO IS=1,2
          DO IP = 1,IJLT(IG)
            IF (WLAT(IP,IS).LT.0.0_JWRB .OR.                            &
     &          WLAT(IP,IS).GT.1.0_JWRB     ) THEN
              WRITE(IU06,*) ' WLAT < 0 or > 1 !!!! ',IP,IS,WLAT(IP,IS)
              LLABORT=.TRUE.
            ENDIF
            IF (WRLON(IP,IS).LT.0.0_JWRB .OR.                           &
     &          WRLON(IP,IS).GT.1.0_JWRB     ) THEN
              WRITE(IU06,*) ' WRLON < 0 or > 1 !!!! ',IP,IS,WRLON(IP,IS)
              LLABORT=.TRUE.
            ENDIF
            IF (WRLAT(IP,IS).LT.0.0_JWRB .OR.                           &
     &          WRLAT(IP,IS).GT.1.0_JWRB     ) THEN
              WRITE(IU06,*) ' WRLAT < 0 or > 1 !!!! ',IP,IS,WRLAT(IP,IS)
              LLABORT=.TRUE.
            ENDIF
          ENDDO
        ENDDO
        DO IS=1,4
          DO IP = 1,IJLT(IG)
            IF (WCOR(IP,IS).LT.0.0_JWRB .OR.                            &
     &          WCOR(IP,IS).GT.1.0_JWRB     ) THEN
              WRITE(IU06,*) ' WCOR < 0 or > 1 !!!! ',IP,IS,WCOR(IP,IS)
              LLABORT=.TRUE.
            ENDIF
          ENDDO
        ENDDO
        IF (LLABORT) THEN
          WRITE(IU06,*) ' ' 
          WRITE(IU06,*) ' PROBLEM IN MUBUF  (see above) !!!! '
          WRITE(IU06,*) ' ' 
          CALL ABORT1
        ENDIF
      ENDIF

      DO IJP=0,NPROPAGS
        WRITE (IU08(IJP)) WLAT 
      ENDDO

      WRITE (IU08(1)) WRLAT
      WRITE (IU08(1)) WRLON

      WRITE (IU08(2)) WCOR 

      DEALLOCATE(WLAT)
      DEALLOCATE(WRLAT)
      DEALLOCATE(WRLON)
      DEALLOCATE(WCOR)


!     DETERMINE THE REAL REDUCTION FACTORS TO BE USED IN THE PROPAGATION

!     CHECK INOUT IS CONSISTENT WITH CURRENT SETUP
      IF (LLOBSTRCT) THEN
        READ (IU01,'(I4)') NFREMAX
        IF (NFREMAX.NE.ML ) THEN
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
      ENDIF

      ALLOCATE(KDUM(NIBLO))

      DO M=1,ML  ! loop over frequencies

!     KOBSLAT
        DO IS=1,2
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT) (IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          DO IJP=0,NPROPAGS
            WRITE(IU08(IJP)) KDUM 
          ENDDO
        ENDDO

!     KOBSLON
        DO IS=1,2
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT) (IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          DO IJP=0,NPROPAGS
            WRITE(IU08(IJP)) KDUM 
          ENDDO
        ENDDO

!     KOBSRLAT
        DO IS=1,2
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT)(IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          WRITE(IU08(1)) KDUM 
        ENDDO

!     KOBSRLON
        DO IS=1,2
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT)(IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          WRITE(IU08(1)) KDUM 
        ENDDO

!     KOBSCOR
        DO IS=1,4
          DO IJ=1,NIBLO
            KDUM(IJ)=1000
          ENDDO
          IF (LLOBSTRCT) THEN
            DO K=1,NY
              WRITE(CX,'(I4.4)') NLONRGG(K)
              FORMAT='('//CX//'I4)'
              READ (IU01,FORMAT) (IDUM(IX,K),IX=1,NLONRGG(K))
            ENDDO
            DO IP = 1,IJLT(IG)
              I = IXLG(IP,IG)
              K = KXLT(IP,IG)
              KDUM(IP)=IDUM(I,K)
            ENDDO
          ENDIF
          WRITE(IU08(2)) KDUM 
        ENDDO

      ENDDO  ! end loop over frequencies

      DEALLOCATE(KDUM)

      END SUBROUTINE MUBUF
