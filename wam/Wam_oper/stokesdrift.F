      SUBROUTINE STOKESDRIFT(F3, IJS, IJL, USTOKES, VSTOKES)
 
!
!***  *STOKESDRIFT*   DETERMINES THE STOKES DRIFT
!
!     PETER JANSSEN MARCH 2009
!
!     PURPOSE.
!     --------
!
!              DETERMINATION OF STOKES DRIFT VECTOR
!
!     INTERFACE.
!     ----------
!              *CALL*  *STOKESDRIFT(F3,IJS,IJL,USTOKES,VSTOKES)*
!
!                       INPUT:
!                            *F3*    - 2-D SPECTRUM
!                            *IJS*   - FIRST GRIDPOINT              
!                            *IJL*   - LAST GRIDPOINT              
!
!                       OUTPUT: 
!                            *USTOKES*   - U-COMPONENT STOKES DRIFT
!                            *VSTOKES*   - V-COMPONENT STOKES DRIFT
!
!     METHOD.
!     -------
!              DETERMINE U- AND V-COMPONENT OF STOKES DRIFT FOLLOWING
!              K.E. KENYON, J.G.R., 74, 6991-6994
!
!     EXTERNALS.
!     ----------
!              NONE
!
!
!-----------------------------------------------------------------------
!     MODULES:
!     --------

      USE YOWPCONS , ONLY : G        ,ZPI
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DELTH    ,TH       ,
     &                      FRATIO   ,COSTH    ,SINTH  
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : TFAC_ST  ,INDEP
       
! ----------------------------------------------------------------------
 
      IMPLICIT NONE

      INTEGER :: IJ, M, K, IJS, IJL
      INTEGER :: NFRE_ODD

      REAL :: CONST, FAC, FAC1, FAC2, FAC3
      REAL :: XLOGFRATIO
      REAL :: DFIM_SIM(NFRE)
      REAL :: F3(IJS:IJL,NANG,NFRE)
      REAL, DIMENSION(IJS:IJL) ::  USTOKES,VSTOKES 
 
!***  1. DETERMINE STOKE DRIFT VECTOR.
!     --------------------------------

      NFRE_ODD=NFRE-1+MOD(NFRE,2)
       
      CONST = 2.*ZPI**3/G*FR(NFRE_ODD)**4


!!!! debile test:
      XLOGFRATIO=LOG(FRATIO)
      DFIM_SIM(1)=DELTH*XLOGFRATIO*FR(1)/3.
      DO M=2,NFRE_ODD-1,2
        DFIM_SIM(M)=4.*DELTH*XLOGFRATIO*FR(M)/3.
        DFIM_SIM(M+1)=2.*DELTH*XLOGFRATIO*FR(M+1)/3.
      ENDDO
      DFIM_SIM(NFRE_ODD)=DELTH*XLOGFRATIO*FR(NFRE_ODD)/3.
      
 
!***  1.1 PERFORM INTEGRATION.
!     ------------------------
 
      DO IJ = IJS,IJL
         USTOKES(IJ) = 0.
         VSTOKES(IJ) = 0.
      ENDDO

      DO M=1,NFRE_ODD
         FAC = DFIM_SIM(M)
         DO K=1,NANG
            FAC1 = FAC*SINTH(K)
            FAC2 = FAC*COSTH(K)
            DO IJ = IJS,IJL
               FAC3 = TFAC_ST(INDEP(IJ),M)*F3(IJ,K,M)
               USTOKES(IJ) = USTOKES(IJ)+FAC1*FAC3
               VSTOKES(IJ) = VSTOKES(IJ)+FAC2*FAC3
            ENDDO
         ENDDO
      ENDDO
 
!***  1.2 ADD CONTRIBUTION OF UNRESOLVED WAVES.
!     -----------------------------------------
 
      DO K=1,NANG
         FAC1 = DELTH*SINTH(K)
         FAC2 = DELTH*COSTH(K)
         DO IJ = IJS,IJL
            FAC3 = CONST*F3(IJ,K,NFRE_ODD)
            USTOKES(IJ) = USTOKES(IJ)+FAC1*FAC3
            VSTOKES(IJ) = VSTOKES(IJ)+FAC2*FAC3
         ENDDO
      ENDDO

      END SUBROUTINE STOKESDRIFT
