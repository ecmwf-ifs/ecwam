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
      USE YOWCOUP  , ONLY : LWNEMOCOU, LWNEMOCOUDEBUG
      USE YOWNEMOFLDS, ONLY : WAM2NEMO 
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! INFORMATION FOR OPTIONAL DEBUGGING
      USE YOWSTAT  , ONLY : CDTPRO

! -------------------------------------------------------------------   

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IJ, NPOINTS, IC
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJL-IJS+1) :: NSWH, NMWP, NPHIEPS, NTAUOC
      REAL(KIND=JWRB), DIMENSION(IJL-IJS+1) :: NEMOSTRN, NEMOUSTOKES, NEMOVSTOKES

! -------------------------------------------------------------------   

IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',0,ZHOOK_HANDLE)


      NPOINTS = IJL-IJS+1

      IF (LWNEMOCOU) THEN

         DO IC = 1, NPOINTS
            IJ = IJS+IC-1
            NSWH(IC)  = WAM2NEMO(IJ)%NSWH
            NMWP(IC)  = WAM2NEMO(IJ)%NMWP
            NPHIEPS(IC) = WAM2NEMO(IJ)%NPHIEPS
            NTAUOC(IC) = WAM2NEMO(IJ)%NTAUOC
            NEMOSTRN(IC) = WAM2NEMO(IJ)%NEMOSTRN
            NEMOUSTOKES(IC) = WAM2NEMO(IJ)%NEMOUSTOKES
            NEMOVSTOKES(IC) = WAM2NEMO(IJ)%NEMOVSTOKES
         ENDDO

#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_UPDATE( IRANK-1, NPROC, MPL_COMM,             &
     &                               NPOINTS                               &
     &                               NSWH, NMWP, NPHIEPS, NTAUOC,          &
     &                               NEMOSTRN, NEMOUSTOKES, NEMOVSTOKES,   &
     &                               CDTPRO, LWNEMOCOUDEBUG )
#endif

      ENDIF

IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',1,ZHOOK_HANDLE)

END SUBROUTINE UPDNEMOFIELDS
