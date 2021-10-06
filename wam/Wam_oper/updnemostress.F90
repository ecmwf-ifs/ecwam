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
      USE PARKIND1    , ONLY : JPRO

! MODULES NEED FOR GRID DEFINTION      
      USE YOWGRID  , ONLY : IJS, IJL
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! DATA FOR NEMO 
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOCOUDEBUG, KCOUSTEP,        &
     &                      NEMOTAUX, NEMOTAUY, NEMONEW10,              &
     &                      NEMOPHIF, NEMONTAU
      USE YOWSTAT  , ONLY : CDTPRO

! -------------------------------------------------------------------   

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IG, IJ

      REAL(KIND=JWRB) :: TAU
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


! -------------------------------------------------------------------   
      IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',0,ZHOOK_HANDLE)

      IF (LWNEMOCOU .AND. ALLOCATED(NEMOTAUX)) THEN

        IG=1

        IF (NEMONTAU.GT.0) THEN
          NEMOTAUX(IJS(IG):IJL(IG))=NEMOTAUX(IJS(IG):IJL(IG))/NEMONTAU
          NEMOTAUY(IJS(IG):IJL(IG))=NEMOTAUY(IJS(IG):IJL(IG))/NEMONTAU
          NEMONEW10(IJS(IG):IJL(IG))=NEMONEW10(IJS(IG):IJL(IG))/NEMONTAU
          NEMOPHIF(IJS(IG):IJL(IG))=NEMOPHIF(IJS(IG):IJL(IG))/NEMONTAU
        ELSE

          NEMOTAUX(IJS(IG):IJL(IG)) = 0.0_JPRO
          NEMOTAUY(IJS(IG):IJL(IG)) = 0.0_JPRO
          NEMONEW10(IJS(IG):IJL(IG)) = 0.0_JPRO
          NEMOPHIF(IJS(IG):IJL(IG)) = 0.0_JPRO
        ENDIF
#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_UPDATE_STRESS( IRANK-1, NPROC, MPL_COMM,   &
     &      IJL(IG)-IJS(IG)+1, NEMOTAUX, NEMOTAUY, NEMONEW10, NEMOPHIF, &
     &      CDTPRO, LWNEMOCOUDEBUG )
#endif
        ! INITIALIZE STRESS FOR ACCUMULATION

        NEMONTAU = 0
        NEMOTAUX(:) = 0.0_JPRO
        NEMOTAUY(:) = 0.0_JPRO
        NEMONEW10(:) = 0.0_JPRO
        NEMOPHIF(:) = 0.0_JPRO

      ENDIF

      IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',1,ZHOOK_HANDLE)

      END SUBROUTINE UPDNEMOSTRESS
