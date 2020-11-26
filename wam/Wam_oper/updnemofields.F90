      SUBROUTINE UPDNEMOFIELDS


!****  *UPDNEMOFIELDS* - UPDATE FIELDS IN NEMO WITH WAVE INFORMATION

!      KRISTIAN MOGENSEN ECMWF    NOVEMBER 2012                    

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PASSES WAVE INFORMATION THROUGH TO
!          NEMO VIA THE NEMO SINGLE EXECUTABLE COUPLING INTERFACE

!          !!!!! APPLIES TO NON-ACCUMULATED FIELDS !!!!!
!          !!!!! SEE UPDNEMOSTRESS FOR ACCUMULATED FIELDS !!!!!

!*    INTERFACE.                                                        
!     ----------                                                        

!          THE NEMO COUPLING IN nemo/coupled/src/nemointerface

!     METHOD.                                                           
!     -------                                                           

!          PARALLEL INTERPOLATION BASED ON PREDETERMINED WEIGHTS

!     EXTERNALS.                                                        
!     ----------                                                        

!          NEMOGCMCOUP_WAM_UPDATE  -  UPDATE WAM FIELDS IN NEMO

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! MODULES NEED FOR GRID DEFINTION      
      USE YOWGRID  , ONLY : IJS, IJL
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
! COUPLING INFORMATION
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOCOUDEBUG,                  &
     &                      LWNEMOCOUSTK, LWNEMOCOUSTRN, KCOUSTEP,      &
     &                      NSWH, NMWP, NPHIEPS, NTAUOC,                &
     &                      NEMOSTRN, NEMOUSTOKES, NEMOVSTOKES
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
! INFORMATION FOR OPTIONAL DEBUGGING
      USE YOWSTAT  , ONLY : CDTPRO

! -------------------------------------------------------------------   

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IG

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
! -------------------------------------------------------------------   

      IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',0,ZHOOK_HANDLE)

      IF (LWNEMOCOU) THEN

        IG=1

#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_UPDATE( IRANK-1, NPROC, MPL_COMM,          &
     &       IJL(IG)-IJS(IG)+1,                                         &
     &       NSWH(IJS(IG):IJL(IG)), NMWP(IJS(IG):IJL(IG)),              &
     &       NPHIEPS(IJS(IG):IJL(IG)), NTAUOC(IJS(IG):IJL(IG)),         &
     &       NEMOSTRN(IJS(IG):IJL(IG)),                                 &
     &       NEMOUSTOKES(IJS(IG):IJL(IG)), NEMOVSTOKES(IJS(IG):IJL(IG)),&
     &       CDTPRO, LWNEMOCOUDEBUG )
#endif


      ENDIF

      IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',1,ZHOOK_HANDLE)

      END SUBROUTINE UPDNEMOFIELDS
