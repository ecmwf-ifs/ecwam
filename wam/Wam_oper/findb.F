      SUBROUTINE FINDB (NDIM, NBOUN, BLATB, BLNGB, IGARB, IJARB)

! ----------------------------------------------------------------------

!**** *FINDB* - FIND BLOCK AND GRID POINT NUMBER.

!     R. PORTZ     MPI         15/01/1991

!*    PURPOSE.
!     -------

!       FIND BLOCK AND GRID POINT NUMBER FOR A GIVEN ARRAY
!       OF LONGITUDES AND LATITUDES.

!**   INTERFACE.
!     ----------

!       *CALL* *FINDB (NDIM, NBOUN, BLATB, BLNGB, IGARB, IJARB)*
!          *NDIM*    - DIMENSION OF ARRAYS.
!          *NBOUN*   - NUMBER OF POINTS IN ARRAYS.
!          *BLATB*   - INPUT LATITUDES.
!          *BLNGB*   - INPUT LONGITUDES.
!          *IGARB*   - OUTPUT BLOCK NUMBERS.
!          *IJARB*   - OUTPUT SEA POINT NUMBERS.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IPER     ,AMOWEP   ,
     &            AMOSOP   ,XDELLA   ,ZDELLO
      USE YOWPARAM , ONLY : NGY

! ----------------------------------------------------------------------

      DIMENSION BLATB(NDIM), BLNGB(NDIM), IGARB(NDIM), IJARB(NDIM)

! ----------------------------------------------------------------------

!*    1. LOOP OVER INPUT LATITUDES, LONGITUDES.
!        --------------------------------------

      DO IO = 1,NBOUN

!*    1.1 COMPUTE GRID MATRIX INDICES.
!         ----------------------------

        IOLT = NINT((BLATB(IO)-AMOSOP)/XDELLA+1.0)
        IF(IOLT.LT.1.OR.IOLT.GT.NGY) GOTO 1002
        ALONG = MOD(BLNGB(IO)-AMOWEP+720.,360.)
        IOLG = NINT(ALONG/ZDELLO(IOLT)+1.0)
        IF (IOLG.EQ.(NLONRGG(IOLT)+1) .AND. IPER.EQ.1) IOLG = 1
        IF(IOLG.LT.1.OR.IOLG.GT.NLONRGG(IOLT)) GOTO 1002

!*    1.2 SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
!         -----------------------------------------

        DO IG=1,IGL
          DO IJ = IJS(IG),IJL(IG)
            IF (IXLG(IJ,IG).EQ.IOLG .AND. KXLT(IJ,IG).EQ.IOLT) THEN
              IGARB(IO) = IG
              IJARB(IO) = IJ
              GOTO 1001
            ENDIF
          ENDDO
        ENDDO
 1002   IGARB(IO) = 0
        IJARB(IO) = 0
 1001   continue
      ENDDO

      RETURN
      END SUBROUTINE FINDB
