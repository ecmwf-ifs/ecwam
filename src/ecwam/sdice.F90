! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE SDICE (KIJS, KIJL, FL1, FLD, SL,             &
     &                  INDEP, WAVNUM, CGROUP,                &
     &                  CICV,CITH)
! ----------------------------------------------------------------------

!**** *SDICE* - COMPUTATION OF SEA ICE ATTENUATION


!     LOTFI AOUF       METEO FRANCE 2023


!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!       *CALL* *SDICE (KIJS, KIJL, FL1, FLD,SL,*
!                      CICVR,CITH)*
!          *KIJS*   - INDEX OF FIRST GRIDPOINT
!          *KIJL*   - INDEX OF LAST GRIDPOINT
!          *FL1*    - SPECTRUM.
!          *FLD*    - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *SL*     - TOTAL SOURCE FUNCTION ARRAY
!          *INDEP*  - DEPTH INDEX
!          *WAVNUM* - WAVE NUMBER
!          *CGROUP* - GROUP SPEED
!          *CICV*   - SEA ICE COVER
!          *CITH*   - SEA ICE THICKNESS

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!     Jie Yu * , W. Erick Rogers, David W. Wang, 2022 
!     DOI:10.1016/j.coldregions.2022.103582

! ----------------------------------------------------------------------

! TODO: sort out modules/params needed
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR      ,GOM     ,TH    ,ZPIFR
!      USE YOWICE   , ONLY : CICOVER  ,CITHICK ! TODO: make ice info available like so
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWPCONS , ONLY : G       ,EPSMIN  ,ZPI

      USE YOWTEST  , ONLY : IU06

      USE YOMHOOK  , ONLY : LHOOK   ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL

      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(IN) :: FL1
      REAL(KIND=JWRB), DIMENSION(KIJL,NANG,NFRE), INTENT(INOUT) :: FLD, SL
      INTEGER(KIND=JWIM), DIMENSION(KIJL), INTENT(IN) :: INDEP
      REAL(KIND=JWRB), DIMENSION(KIJL,NFRE), INTENT(IN) :: WAVNUM, CGROUP
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: CICV
      REAL(KIND=JWRB), DIMENSION(KIJL), INTENT(IN) :: CITH

      REAL(KIND=JWRB), DIMENSION(NFRE)    :: XK2

      INTEGER(KIND=JWIM) :: IMODEL                       !! DAMPING MODEL: 1=FIT TO TEMPELFJORD DATA, 2=Jie Yu 2022

      INTEGER(KIND=JWIM) :: IJ, K, M
      REAL(KIND=JWRB)    :: EWH, X, XK, ZPI4G2
      REAL(KIND=JWRB)    :: ALP, ALP1, ALP2              !! ALP=SPATIAL ATTENUATION RATE OF ENERGY?
      REAL(KIND=JWRB)    :: BETA, TEMP
      REAL(KIND=JWRB)    :: CDICE, DICE, D1, D2, RIND
      REAL(KIND=JWRB)    :: HICEMAX, HICEMIN, HICE
      REAL(KIND=JWRB)    :: MUSHFRAC, HBDY               !! FRACTION OF ICE THICKNESS WHICH IS MUSHY
      
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('SDICE',0,ZHOOK_HANDLE)

      IMODEL = 2
!      IF (ITEST.GE.2) WRITE (IU06,*)'IMODEL =',IMODEL
      DICE=1.0_JWRB
      HICEMAX=4.0_JWRB
      HICEMIN=0.1_JWRB
      IF (IMODEL.EQ.1) THEN
         ! Ice damping with best fit from Tempelfjorden observations
         WRITE (IU06,*)'Ice damping with best fit from Tempelfjorde obs'
         CDICE=0.0656_JWRB
      ELSE IF (IMODEL.EQ.2) THEN
         WRITE (IU06,*)'Ice damping based on: '
         WRITE (IU06,*)'  Jie Yu, W. Erik Rogers, David W. Wang 2022'
         CDICE=0.1274_JWRB*( ZPI/SQRT(G) )**(4.5_JWRB)
      END IF
      
      DO M = 1,NFRE
         DO K = 1,NANG
            DO IJ = KIJS,KIJL
               IF (IMODEL.EQ.1) THEN
                  ! XK = TFAK(INDEP(IJ),M)                ! OLD
                  ! ALP = CDICE*CITH(IJ)*XK**2            ! OLD 
                  ALP = CDICE*CITH(IJ)*WAVNUM(IJ,M)**2    ! NEW
               ELSE IF (IMODEL.EQ.2) THEN
                  ALP = 2._JWRB*CDICE*(CITH(IJ)**(1.25_JWRB))*(FR(M)**(4.5_JWRB))
               END IF
               BETA=1._JWRB-CICV(IJ)
               ! TEMP = -CICV(IJ)*ALP*TCGOND(INDEP(IJ),M) ! OLD
               TEMP = -CICV(IJ)*ALP*CGROUP(IJ,M)         ! NEW
               SL(IJ,K,M)  = BETA*SL(IJ,K,M)  + FL1(IJ,K,M)*TEMP
               FLD(IJ,K,M) = BETA*FLD(IJ,K,M) + TEMP

            END DO
         END DO
      END DO
      
      ! Deep water case only, this worked
         
     !    ZPI4G2=ZPI**4/G**2
     !    DO M = 1,NFRE
     !       XK2(M)=ZPI4G2*FR(M)**4
     !       DO K = 1,NANG
     !          DO IJ = KIJS,KIJL
               
     !             IF (IMODEL.EQ.1) THEN
     !                ALP = CDICE*CITH(IJ)*XK2(M)
     !             ELSE IF (IMODEL.EQ.2) THEN
     !   !               write(*,*) 'test CIT', CITHICK(IJ)
     !                ALP = 2._JWRB*CDICE*(CITH(IJ)**(1.25_JWRB))*(FR(M)**(4.5_JWRB))
     !             END IF
               
     !             BETA=1._JWRB-CICV(IJ)
     !             TEMP = -CICV(IJ)*ALP*GOM(M) 
     !             SL(IJ,K,M)  = BETA*SL(IJ,K,M)  + FL1(IJ,K,M)*TEMP
     !             FLD(IJ,K,M) = BETA*FLD(IJ,K,M) + TEMP
               
     !          END DO
     !       END DO
     !    END DO 
      
      ! END IF

      
      IF (LHOOK) CALL DR_HOOK('SDICE',1,ZHOOK_HANDLE)

      END SUBROUTINE SDICE
