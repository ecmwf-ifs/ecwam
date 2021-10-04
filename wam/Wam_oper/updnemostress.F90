SUBROUTINE UPDNEMOSTRESS


!****  *UPDNEMOSTRESS* - UPDATE FIELDS PASSED TO NEMO WITH WAVE INFORMATION

!      KRISTIAN MOGENSEN ECMWF    JULY 2013                    

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PASSES WAVE INFORMATION THROUGH TO
!          NEMO VIA THE NEMO SINGLE EXECUTABLE COUPLING INTERFACE

!          !!!!! APPLIES TO ACCUMULATED FIELDS !!!!!
!          !!!!! SEE UPDNEMOFILDS FOR NON-ACCUMULATED FIELDS !!!!!

!*    INTERFACE.                                                        
!     ----------                                                        


!          THE NEMO COUPLING IN nemo/coupled/src/nemointerface

!     METHOD.                                                           
!     -------                                                           

!          PARALLEL INTERPOLATION BASED ON PREDETERMINED WEIGHTS

!     EXTERNALS.                                                        
!     ----------                                                        

!          NEMOGCMCOUP_WAM_TAU_UPDATE  -  UPDATE WAM STRESS IN NEMO

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWNEMOP    , ONLY : NEMODP

! MODULES NEED FOR GRID DEFINTION      
      USE YOWGRID  , ONLY : IJS, IJL
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! DATA FOR NEMO 
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOCOUDEBUG, NEMONTAU
      USE YOWNEMOFLDS, ONLY : WAM2NEMO
      USE YOWSTAT  , ONLY : CDTPRO

! -------------------------------------------------------------------   

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IJ, NPOINTS, IC

      REAL(KIND=JWRB) :: ZNEMONTAUM1
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=NEMODP), DIMENSION(IJL-IJS+1) :: NEMOTAUX, NEMOTAUY, NEMOWSWAVE, NEMOPHIF

! -------------------------------------------------------------------   
IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',0,ZHOOK_HANDLE)


      IF (LWNEMOCOU .AND. ALLOCATED(WAM2NEMO) ) THEN

        NPOINTS = IJL-IJS+1

        IF (NEMONTAU > 0) THEN
          ZNEMONTAUM1= 1.0_JWRB/NEMONTAU
          DO IC = 1, NPOINTS
            IJ = IJS+IC-1
            NEMOTAUX(IC)  = WAM2NEMO(IJ)%NEMOTAUX*ZNEMONTAUM1
            NEMOTAUY(IC) = WAM2NEMO(IJ)%NEMOTAUY*ZNEMONTAUM1
            NEMOWSWAVE(IC) = WAM2NEMO(IJ)%NEMOWSWAVE*ZNEMONTAUM1
            NEMOPHIF(IC) = WAM2NEMO(IJ)%NEMOPHIF*ZNEMONTAUM1
          ENDDO
        ELSE
          DO IC = 1, NPOINTS
            NEMOTAUX(IC) = 0.0_NEMODP
            NEMOTAUY(IC) = 0.0_NEMODP
            NEMOWSWAVE(IC) = 0.0_NEMODP
            NEMOPHIF(IC) = 0.0_NEMODP
          ENDDO
        ENDIF

#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_UPDATE_STRESS( IRANK-1, NPROC, MPL_COMM,   &
     &      NPOINTS, NEMOTAUX, NEMOTAUY, NEMOWSWAVE, NEMOPHIF, &
     &      CDTPRO, LWNEMOCOUDEBUG )
#endif
        ! INITIALIZE STRESS FOR ACCUMULATION

        NEMONTAU = 0
        WAM2NEMO(IJS:IJL)%NEMOTAUX = 0.0_NEMODP
        WAM2NEMO(IJS:IJL)%NEMOTAUY = 0.0_NEMODP
        WAM2NEMO(IJS:IJL)%NEMOWSWAVE = 0.0_NEMODP
        WAM2NEMO(IJS:IJL)%NEMOPHIF = 0.0_NEMODP

      ENDIF

IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',1,ZHOOK_HANDLE)

END SUBROUTINE UPDNEMOSTRESS
