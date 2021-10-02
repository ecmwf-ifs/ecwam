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

      INTEGER(KIND=JWIM) :: IJ

      REAL(KIND=JWRB) :: TAU, ZNEMONTAUM1
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! -------------------------------------------------------------------   
IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',0,ZHOOK_HANDLE)

ASSOCIATE(NEMOTAUX  => WAM2NEMO%NEMOTAUX,  &
 &        NEMOTAUY => WAM2NEMO%NEMOTAUY,  &
 &        NEMOWSWAVE => WAM2NEMO%NEMOWSWAVE,  &
 &        NEMOPHIF => WAM2NEMO%NEMOPHIF)

      IF (LWNEMOCOU) THEN

        IF (NEMONTAU > 0) THEN
          ZNEMONTAUM1= 1.0_JWRB/NEMONTAU
          NEMOTAUX(IJS:IJL)=NEMOTAUX(IJS:IJL)*ZNEMONTAUM1
          NEMOTAUY(IJS:IJL)=NEMOTAUY(IJS:IJL)*ZNEMONTAUM1
          NEMOWSWAVE(IJS:IJL)=NEMOWSWAVE(IJS:IJL)*ZNEMONTAUM1
          NEMOPHIF(IJS:IJL)=NEMOPHIF(IJS:IJL)*ZNEMONTAUM1
        ELSE
          NEMOTAUX(IJS:IJL) = 0.0_NEMODP
          NEMOTAUY(IJS:IJL) = 0.0_NEMODP
          NEMOWSWAVE(IJS:IJL) = 0.0_NEMODP
          NEMOPHIF(IJS:IJL) = 0.0_NEMODP
        ENDIF

#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_UPDATE_STRESS( IRANK-1, NPROC, MPL_COMM,   &
     &      IJL-IJS+1, NEMOTAUX, NEMOTAUY, NEMOWSWAVE, NEMOPHIF, &
     &      CDTPRO, LWNEMOCOUDEBUG )
#endif
        ! INITIALIZE STRESS FOR ACCUMULATION

        NEMONTAU = 0
        NEMOTAUX(:) = 0.0_NEMODP
        NEMOTAUY(:) = 0.0_NEMODP
        NEMOWSWAVE(:) = 0.0_NEMODP
        NEMOPHIF(:) = 0.0_NEMODP

      ENDIF

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',1,ZHOOK_HANDLE)

END SUBROUTINE UPDNEMOSTRESS
