      SUBROUTINE MGRID (BATHY)

! ----------------------------------------------------------------------

!**** *MGRID* - ROUTINE TO ARRANGE WAMODEL GRID.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       TO ARRANGE WAMODEL GRID FOR A GIVEN AREA AND COMPUTE VARIOUS
!       MODEL CONSTANTS.

!**   INTERFACE.
!     ----------

!       *CALL* *MGRID (BATHY)*
!          *BATHY*     - BATHYMETRY /DATA OF PART

!     METHOD.
!     -------

!       THE NUMBER OF SEA POINTS PER LATITUDE IS COUNTED AND MODEL
!       BLOCKS OF MAXIMUM LENGTH OF NIBLO ARE CONSTRUCTED.

!     EXTERNALS.
!     ----------

!       *MBLOCK*    - SUB. TO GENERATE A BLOCK.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO    ,NBLO     ,
     &            NOVER    ,NIBL1
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : NY

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL, INTENT(INOUT) :: BATHY(NGX,NGY)

      INTEGER :: K, I, IL, KA, KE
      INTEGER, DIMENSION(NGY) :: IPP

! ----------------------------------------------------------------------


!*    1. COUNT NUMBER OF SEA POINTS PER LATITUDE.
!        ----------------------------------------

      DO K=1,NY
        IPP(K) = 0
        DO I=1,NLONRGG(K)
          IF (BATHY(I,K).GT.-990.) THEN
            IPP(K) = IPP(K) + 1
          ENDIF
        ENDDO
      ENDDO

      IF(NIBLO.LE.0) THEN
        NIBLO=0
        DO K=1,NY
          NIBLO=NIBLO+IPP(K)
        ENDDO
        IF(NBLO.GT.1)  NIBLO=NIBLO/NBLO+1
      ENDIF

! ----------------------------------------------------------------------

!*    2. MAKE BLOCKS.
!        ------------

      IL = 0
      KA = 1
      DO K = 1,NY
        IL = IL + IPP(K)
        IF (IL.GT.NIBLO) THEN
          KE = K-1
          CALL MBLOCK (BATHY, KA, KE, IPP)
          KA = KE-1
          IL = IPP(KA)+IPP(KE)+IPP(KE+1)
        ENDIF
      ENDDO
      CALL MBLOCK (BATHY, KA, NY, IPP)

      IF(NBLO.GT.1) THEN
        NOVER=1
        DO K=1,NY
          NOVER=MAX(NOVER,IPP(K)) 
        ENDDO
        NIBL1=NIBLO
      ELSE
        NOVER=1
        NIBL1=1
      ENDIF

      RETURN
      END SUBROUTINE MGRID
