      SUBROUTINE FINDB (NDIM, NBOUN, BLATB, BLNGB, IGARB, IJARB,
     &                  NPROC, NSTART,NEND, IRANK)

! ----------------------------------------------------------------------

!**** *FINDB* - FIND BLOCK AND GRID POINT NUMBER.

!     R. PORTZ     MPI         15/01/1991
!     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING

!*    PURPOSE.
!     -------

!       FIND BLOCK AND GRID POINT NUMBER FOR A GIVEN ARRAY
!       OF LONGITUDES AND LATITUDES.
!       IN CASE OF MESSAGE PASSING ARCHITECHTURE, THE SEARCH IS ONLY
!       DONE ON THE SUBREGION WHICH IS DEFINED ON THAT PE.

!**   INTERFACE.
!     ----------

!       *CALL* *FINDB (NDIM, NBOUN, BLATB, BLNGB, IGARB, IJARB)*
!          *NDIM*    - DIMENSION OF ARRAYS.
!          *NBOUN*   - NUMBER OF POINTS IN ARRAYS.
!          *BLATB*   - INPUT LATITUDES.
!          *BLNGB*   - INPUT LONGITUDES.
!          *IGARB*   - OUTPUT BLOCK NUMBERS.
!          *IJARB*   - OUTPUT SEA POINT NUMBERS.
!          *NPROC*     DIMENSION OF NSTART AND NEND
!          *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!          *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!          *IRANK*     PE INDEX


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
     &            AMOSOP   ,AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS

      USE YOWPARAM , ONLY : NGX

! ----------------------------------------------------------------------

      DIMENSION BLATB(NDIM), BLNGB(NDIM), IGARB(NDIM), IJARB(NDIM)
      INTEGER,DIMENSION(NPROC) :: NSTART,NEND

! ----------------------------------------------------------------------

!*    1. LOOP OVER INPUT LATITUDES, LONGITUDES.
!        --------------------------------------

      IF(LMESSPASS) THEN
        ILATS=KXLT(NSTART(IRANK),1)-1
        ILATN=KXLT(NEND(IRANK),1)+1
      ELSE
        ILATS=1
        ILATN=NGX
      ENDIF

      DO IO = 1,NBOUN

!*    1.1 COMPUTE GRID MATRIX INDICES.
!         ----------------------------

        IOLT = NINT((BLATB(IO)-AMOSOP)/XDELLA+1.0)
        ALONG = MOD(BLNGB(IO)-AMOWEP+720.,360.)
        IF(IOLT.LT.ILATS.OR.IOLT.GT.ILATN ) THEN
          ZDEL=XDELLO
        ELSE
          ZDEL=ZDELLO(IOLT)
        ENDIF
        IOLG = NINT(ALONG/ZDEL+1.0)
        IF (IOLG.EQ.(NLONRGG(IOLT)+1) .AND. IPER.EQ.1) IOLG = 1

        IF(IOLG.LT.1.OR.IOLG.GT.NLONRGG(IOLT).OR.
     &     IOLT.LT.ILATS.OR.IOLT.GT.ILATN ) THEN
          IGARB(IO) = 0
          IJARB(IO) = 0
        ELSE
!*    1.2 SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
!         -----------------------------------------
          IF(LMESSPASS) THEN
            IGARB(IO) = 0
            IJARB(IO) = 0
            DO IJ = NSTART(IRANK),NEND(IRANK)
              IF (IXLG(IJ,1).EQ.IOLG .AND. KXLT(IJ,1).EQ.IOLT) THEN
                IGARB(IO) = 1
                IJARB(IO) = IJ
                EXIT
              ENDIF
            ENDDO
          ELSE
            IGARB(IO) = 0
            IJARB(IO) = 0
            DO IG=1,IGL
              DO IJ = IJS(IG),IJL(IG)
                IF (IXLG(IJ,IG).EQ.IOLG .AND. KXLT(IJ,IG).EQ.IOLT) THEN
                  IGARB(IO) = IG
                  IJARB(IO) = IJ
                  EXIT 
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE FINDB
