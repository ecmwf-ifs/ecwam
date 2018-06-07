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

! MODULES NEED FOR GRID DEFINTION      
      USE YOWGRID  , ONLY : IJS, IJL
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
! CONSTANTS
      USE YOWPCONS , ONLY : EPSUS, EPSU10
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! INFORMATION FO OPTIONAL DEBUGGING
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOCOUDEBUG, KCOUSTEP,
     &                      NEMOTAUX, NEMOTAUY, NEMONEW10, 
     &                      NEMOPHIF, NEMONTAU
      USE YOWSTAT  , ONLY : CDATEF, CDTPRO

! -------------------------------------------------------------------   

      IMPLICIT NONE

      INTEGER :: IG, IJ

      REAL :: TAU
      REAL :: ZHOOK_HANDLE


! -------------------------------------------------------------------   
#ifdef ECMWF 
      IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',0,ZHOOK_HANDLE)
#endif

      IF (LWNEMOCOU) THEN

        IG=1

        IF (NEMONTAU.GT.0) THEN
          NEMOTAUX(IJS(IG):IJL(IG))=NEMOTAUX(IJS(IG):IJL(IG))/NEMONTAU
          NEMOTAUY(IJS(IG):IJL(IG))=NEMOTAUY(IJS(IG):IJL(IG))/NEMONTAU
          NEMONEW10(IJS(IG):IJL(IG))=NEMONEW10(IJS(IG):IJL(IG))/NEMONTAU
          NEMOPHIF(IJS(IG):IJL(IG))=NEMOPHIF(IJS(IG):IJL(IG))/NEMONTAU
        ELSE
!!!debile
        if(.not. allocated(NEMOTAUX)) then
          write(*,*) 'updnemostress, we have a serious problem !!!!'
        endif

          NEMOTAUX(IJS(IG):IJL(IG)) = 0.0
          NEMOTAUY(IJS(IG):IJL(IG)) = 0.0
          NEMONEW10(IJS(IG):IJL(IG)) = 0.0 
          NEMOPHIF(IJS(IG):IJL(IG)) = 0.0 
        ENDIF
#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_UPDATE_STRESS( IRANK-1, NPROC, MPL_COMM,
     &      IJL(IG)-IJS(IG)+1, NEMOTAUX, NEMOTAUY, NEMONEW10, NEMOPHIF,
     &      CDTPRO, LWNEMOCOUDEBUG )
#endif
        ! INITIALIZE STRESS FOR ACCUMULATION

        NEMONTAU = 0
        NEMOTAUX(:) = 0.0
        NEMOTAUY(:) = 0.0
        NEMONEW10(:) = 0.0
        NEMOPHIF(:) = 0.0

      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('UPDNEMOSTRESS',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE UPDNEMOSTRESS
