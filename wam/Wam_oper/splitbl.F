      SUBROUTINE SPLITBL (FL, IJA, IJE, IU14)

! ----------------------------------------------------------------------

!**** *SPLITBL*   SAVES SECOND LATITUDE OF BLOCK IG
!****             FOR BLOCK IG-1 FOR COMPUTATION OF NEXT TIME STEP.

!     H. GUNTHER       GKSS/ECMWF         SEPTEMBER 1989

!*    PURPOSE.
!     --------

!       INTERMEDIATE STORAGE OF WAVE INFORMATION FOR USE
!       AT NEXT TIME STEP (SEE SUB FILLBL).

!**   INTERFACE.
!     ----------

!       *CALL* *SPLITBL (FL, IJA, IJE, IU14)*
!          *FL*      REAL ARRAY   BLOCK OF SPECTRA
!          *IJA*     INTEGER      FIRST GRID POINT TO BE SAVED.
!          *IJE*     INTEGER      LAST  GRID POINT TO BE SAVED.
!          *IU14*    INTEGER      UNIT FOR OUTPUT OF SECOND LATITUDE.

!     EXTERNALS.
!     ----------

!       NONE.

!     METHOD.
!     -------

!       UNFORMATED WRITE TO SPECIFIED UNIT.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------


      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWMPP   , ONLY : NINF     ,NSUP



! ----------------------------------------------------------------------

      REAL,DIMENSION (NINF-1:NSUP,NANG,NFRE) :: FL

! ----------------------------------------------------------------------

!*    1. OUTPUT OF SECOND LATITUDE OF BLOCK FOR LAST
!        LATITUDE OF PREVIOUS BLOCK.
!        ---------------------------------------------

      WRITE(IU14) (((FL(IJ,K,M),IJ=IJA,IJE),K=1,NANG),M=1,NFRE)

      RETURN
      END SUBROUTINE SPLITBL
