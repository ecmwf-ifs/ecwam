      SUBROUTINE MUBUF (IA1, IG, IU08, IU18, IFORM)

! ----------------------------------------------------------------------

!**** *MUBUF* - ROUTINE TO ARRANGE YOWMON UBUF FOR ONE BLOCK.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       TO ARRANGE NEIGHBOUR GRID POINT INDICES FOR A BLOCK

!**   INTERFACE.
!     ----------

!       *CALL* *MUBUF (IA1, IG, IU08, IU18, IFORM)*
!          *IA1*     - TOPOGRAPHIC DATA.
!          *IG*      - BLOCK NUMBER.
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
     &            ZDELLO
      USE YOWUBUF  , ONLY : KLAT     ,KLON

! ----------------------------------------------------------------------

      DIMENSION IA1(NGX, NGY)

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ---------------

      ALLOCATE(KLAT(NIBLO,2))
      ALLOCATE(KLON(NIBLO,2))

!*    1. INITIALISE ARRAYS.
!        ------------------

      DO IJ=1,NIBLO
        KLAT(IJ,1) = 0
        KLAT(IJ,2) = 0
        KLON(IJ,1) = 0
        KLON(IJ,2) = 0
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
          IF (IA1(IMIN,K-1).GT.-990) THEN
            DO IH = IP,1,-1
              IF (IXLG(IH,IG).EQ.IMIN .AND.
     &         KXLT(IH,IG).EQ.K-1) GOTO 2202
            ENDDO
 2202       CONTINUE
            KLAT(IP,1) = IH
          ENDIF
        ENDIF
        IF (K.LT.NY) THEN
          XPLUS = REAL(I-1)*ZDELLO(K)/ZDELLO(K+1)
          IPLUS = NINT(XPLUS) + 1
          IF (IA1(IPLUS,K+1).GT.-990) THEN
            DO IH = IP,IJLT(IG)
              IF (IXLG(IH,IG).EQ.IPLUS .AND.
     &         KXLT(IH,IG).EQ.K+1) GOTO 2204
            ENDDO
 2204       CONTINUE
            KLAT(IP,2) = IH
          ENDIF
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    3. OUTPUT OF COMMON UBUF FOR ONE BLOCK.
!        ------------------------------------

      CALL OUTUBUF (IU08, IU18, IFORM)

      RETURN
      END SUBROUTINE MUBUF
