      SUBROUTINE FILLBL (FL, IG, IU13)

! ----------------------------------------------------------------------

!**** *FILLBL* - ADDS FIRST AND LAST LATITUDE TO BLOCK
!****            AND SAVES SECOND TO LAST LAT FOR NEXT BLOCK.

!     HEINZ GUNTHER         ECMWF         SEPTEMBER 1989

!*    PURPOSE.
!     --------

!       THIS ROUTINE ADDS TO A BLOCK OF SPECTRA THE SEA POINTS IN
!       THE FIRST AND LAST LATITUDE. THE INFORMATION IS READ FROM
!       THE INTERMEDIATE STORAGE (SEE SUB SPLITBL) FOR THE LAST
!       LATITUDE AND TAKEN FROM COMMON OVER FOR THE FIRST LATITUDE.
!       THE YOWMON IS FILLED WITH THE SECOND TO LAST LATITUDE TO
!       PASS THE INFORMATION TO THE NEXT BLOCK.

!**   INTERFACE.
!     ----------

!       *CALL* *FILLBL (FL, IG, IU13)*
!          *FL*      REAL     BLOCK OF SPECTRA.
!          *IG*      INTEGER  BLOCK NUMBER.
!          *IU13*    INTEGER  INPUT UNIT OF LAST LATITUDE.

!     EXTERNALS.
!     ----------

!       NONE.

!     METHOD.
!     -------

!          UNFORMATED READ FROM SPECIFIED UNITS OR COPY.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------


      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NOVER
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJLS     ,IJL      ,IJLT
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWOVER  , ONLY : OV



! ----------------------------------------------------------------------

      DIMENSION FL(NINF-1:NSUP,NANG,NFRE)

! ----------------------------------------------------------------------
!     0. ALLOCATE

      IF(.NOT.ALLOCATED(OV)) ALLOCATE(OV(NOVER,NANG,NFRE))

!*    1. ADD LATITUDE BEFORE LAST OF PREVIOUS BLOCK INTO FIRST
!        LATITUDE OF BLOCK.
!        ------------------------------------------------------

      IF(IG.NE.1) THEN
        IJA = 1
        IJE = IJS(IG)-1
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJA,IJE
              FL(IJ,K,M) = OV(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2. ADD SECOND LATITUDE OF NEXT BLOCK INTO LAST LATITUDE OF BLOCK.
!        --------------------------------------------------------------

      IF(IG.NE.IGL) THEN
        IJA = IJL(IG)+1
        IJE = IJLT(IG)
        READ(IU13) (((FL(IJ,K,M),IJ=IJA,IJE),K=1,NANG),M=1,NFRE)
      ENDIF

! ----------------------------------------------------------------------

!*    3. LATITUDE BEFORE LAST OF BLOCK IS SAVED FOR FIRST
!        LATITUDE OF NEXT BLOCK.
!        -------------------------------------------------

      IF (IG.NE.IGL) THEN
        IJA = IJLS(IG)
        IJH = 1-IJA
        IJE = IJL(IG)
        DO M=1,NFRE
          DO K=1,NANG
            DO IJ=IJA,IJE
              OV(IJ+IJH,K,M) = FL(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      RETURN
      END SUBROUTINE FILLBL
