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
     &                      LWNEMOCOUSTK, LWNEMOCOUSTRN, KCOUSTEP
      USE YOWNEMOP   , ONLY : NEMODP
      USE YOWNEMOFLDS, ONLY : NEMOSTRN, NEMOUSTOKES, NEMOVSTOKES
! DATA TO BE SEND
      USE YOWMEAN  , ONLY : NSWH, NMWP, NPHIEPS, NTAUOC
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! INFORMATION FOR OPTIONAL DEBUGGING
      USE YOWSTAT  , ONLY : CDATEF, CDTPRO

! -------------------------------------------------------------------   

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IG, IJ

      REAL(KIND=JWRB) :: TAU
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=NEMODP), ALLOCATABLE,DIMENSION(:) :: ZDUMU, ZDUMV, ZDUMSTRN

! -------------------------------------------------------------------   

      IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',0,ZHOOK_HANDLE)

      IF (LWNEMOCOU) THEN

        IG=1

! MEAN SQUARE WAVE STRAIN IN SEA ICE
      
        ALLOCATE(ZDUMSTRN(IJS(IG):IJL(IG)))
        IF (LWNEMOCOUSTRN) THEN
           ZDUMSTRN(IJS(IG):IJL(IG)) = NEMOSTRN(IJS(IG):IJL(IG))
        ELSE
           ZDUMSTRN(IJS(IG):IJL(IG)) = 0.0_NEMODP
        ENDIF
  
! STOKES DRIFT

        ALLOCATE(ZDUMU(IJS(IG):IJL(IG)),ZDUMV(IJS(IG):IJL(IG)))
        IF (LWNEMOCOUSTK) THEN
           ZDUMU(IJS(IG):IJL(IG)) = NEMOUSTOKES(IJS(IG):IJL(IG))
           ZDUMV(IJS(IG):IJL(IG)) = NEMOVSTOKES(IJS(IG):IJL(IG))
        ELSE
           ZDUMU(IJS(IG):IJL(IG)) = 0.0_NEMODP
           ZDUMV(IJS(IG):IJL(IG)) = 0.0_NEMODP
        ENDIF

#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_UPDATE( IRANK-1, NPROC, MPL_COMM,          &
     &       IJL(IG)-IJS(IG)+1,                                         &
     &       NSWH(IJS(IG):IJL(IG)), NMWP(IJS(IG):IJL(IG)),              &
     &       NPHIEPS(IJS(IG):IJL(IG)), NTAUOC(IJS(IG):IJL(IG)),         &
     &       ZDUMSTRN(IJS(IG):IJL(IG)),                                 &
     &       ZDUMU(IJS(IG):IJL(IG)), ZDUMV(IJS(IG):IJL(IG)),            &
     &       CDTPRO, LWNEMOCOUDEBUG )
#endif

        DEALLOCATE(ZDUMSTRN,ZDUMU,ZDUMV)

      ENDIF

      IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',1,ZHOOK_HANDLE)

      END SUBROUTINE UPDNEMOFIELDS
