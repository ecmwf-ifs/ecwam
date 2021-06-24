      SUBROUTINE RECVNEMOFIELDS(LREST,LINIT)


!****  *RECVNEMOFIELDS* - UPDATE FIELDS WAVE FIELDS WITH NEMO INFORMATION

!      KRISTIAN MOGENSEN ECMWF    MARCH 2013

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PASSES NEMO INFORMATION THROUGH TO
!          WAM VIA THE NEMO SINGLE EXECUTABLE COUPLING INTERFACE

!*    INTERFACE.                                                        
!     ----------                                                        


!     METHOD.                                                           
!     -------                                                           

!          PARALLEL INTERPOLATION BASED ON PREDETERMINED WEIGHTS

!     EXTERNALS.                                                        
!     ----------                                                        

!          NEMOGCMCOUP_WAM_GET  -  UPDATE NEMO FIELDS IN WAM

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! MODULES NEED FOR GRID DEFINTION      
      USE YOWGRID  , ONLY : IJS, IJL
! MODULES NEEDED FOR LAKE MASK HANDLING
      USE YOWMAP   , ONLY : IFROMIJ  ,JFROMIJ
      USE YOWWIND  , ONLY : FIELDG   ,LLNEWCURR 
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
      USE MPL_MODULE
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! COUPLING INFORMATION
      USE YOWCOUP  , ONLY : LWCOU, LWNEMOCOUCIC, LWNEMOCOUCIT, LWNEMOCOUCUR
! ICE AND CURRENT INFORMATION 
      USE YOWICE   , ONLY : CITHICK, CICOVER
      USE YOWCURR  , ONLY : U        ,V        , CURRENT_MAX
! OUTPUT FORTRAN UNIT
      USE YOWTEST  , ONLY : IU06
! NEMO FIELDS ON WAVE GRID
      USE YOWNEMOFLDS,ONLY: NEMOSST, NEMOCICOVER, NEMOCITHICK,          &
     &                      NEMOUCUR, NEMOVCUR, LNEMOCITHICK,           &
     &                      LNEMOICEREST

! -------------------------------------------------------------------   

      IMPLICIT NONE
#include "abort1.intfb.h"
      LOGICAL, INTENT(IN) :: LREST ! RESTART SO UPDATE FROM RESTART VALUES
      LOGICAL, INTENT(IN) :: LINIT ! UPDATE CICOVER, CITHICK, And U and V CURRENTS AT INITIAL TIME
                                   ! IF NEEDED.
      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      INTEGER(KIND=JWIM) :: IX, IY, IJ 
! -------------------------------------------------------------------   

      IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',0,ZHOOK_HANDLE)

      ! IF WE ARE IN RESTART JUST COPY THE RESTART FILE INFO.
      IF (LREST) THEN
        NEMOCICOVER(IJS:IJL)=CICOVER(IJS:IJL)
        NEMOCITHICK(IJS:IJL)=CITHICK(IJS:IJL)
        IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
          NEMOUCUR(IJS:IJL)=U(IJS:IJL)
          NEMOVCUR(IJS:IJL)=V(IJS:IJL)
        ELSE
          NEMOUCUR(IJS:IJL)=0.0_JWRB
          NEMOVCUR(IJS:IJL)=0.0_JWRB
        ENDIF
        LNEMOICEREST=.TRUE.
        LNEMOCITHICK=.TRUE.
      ELSE
#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_GET( IRANK-1, NPROC, MPL_COMM,             &
     &                            IJL-IJS+1,                    &
     &                            NEMOSST(IJS:IJL),             &
     &                            NEMOCICOVER(IJS:IJL),         &
     &                            NEMOCITHICK(IJS:IJL),         &
     &                            NEMOUCUR(IJS:IJL),            &
     &                            NEMOVCUR(IJS:IJL),            &
     &                            LNEMOCITHICK )

        LLNEWCURR=.TRUE. 

#endif
        LNEMOICEREST=.FALSE.

      ENDIF

!     UPDATE CICOVER, CITHICK AND U and V CURRENTS AT INITIAL TIME ONLY !!!!
      IF (LINIT) THEN

        WRITE(IU06,*)' RECVNEMOFIELDS: INITIALISE OCEAN FIELDS'

        IF(LWCOU) THEN
          IF (LWNEMOCOUCIC) THEN
            DO IJ = IJS,IJL
              IX = IFROMIJ(IJ)
              IY = JFROMIJ(IJ)
!            if lake cover = 0, we assume open ocean point, then get sea ice directly from NEMO
              IF (FIELDG(IX,IY)%LKFR .LE. 0.0_JWRB ) THEN
                CICOVER(IJ)=NEMOCICOVER(IJ)
              ELSE
!            get ice information from atmopsheric model
                CICOVER(IJ)=FIELDG(IX,IY)%CIFR 
              ENDIF
            ENDDO
          ENDIF

          IF (LWNEMOCOUCIT) THEN
            DO IJ = IJS,IJL
              IX = IFROMIJ(IJ)
              IY = JFROMIJ(IJ)
!            if lake cover = 0, we assume open ocean point, then get sea ice thickness directly from NEMO
              IF (FIELDG(IX,IY)%LKFR .LE. 0.0_JWRB ) THEN
                CITHICK(IJ)=NEMOCICOVER(IJ)*NEMOCITHICK(IJ)
              ELSE
                CICOVER(IJ)=0.5_JWRB*NEMOCICOVER(IJ)
              ENDIF
            ENDDO
          ENDIF

          IF (LWNEMOCOUCUR) THEN
            IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
              DO IJ = IJS,IJL
                IX = IFROMIJ(IJ)
                IY = JFROMIJ(IJ)
!              if lake cover = 0, we assume open ocean point, then get currents directly from NEMO
                IF (FIELDG(IX,IY)%LKFR .LE. 0.0_JWRB ) THEN
                  U(IJ) = SIGN(MIN(ABS(NEMOUCUR(IJ)),CURRENT_MAX),NEMOUCUR(IJ))
                  V(IJ) = SIGN(MIN(ABS(NEMOVCUR(IJ)),CURRENT_MAX),NEMOVCUR(IJ))
                ELSE
                  U(IJ)=0.0_JWRB
                  V(IJ)=0.0_JWRB
                ENDIF
              ENDDO

            ENDIF
          ENDIF

        ELSE

          IF (LWNEMOCOUCIC) THEN
            CICOVER(IJS:IJL)=NEMOCICOVER(IJS:IJL)
          ENDIF
          IF (LWNEMOCOUCIT) THEN
            CITHICK(IJS:IJL)=NEMOCICOVER(IJS:IJL)*     &
     &                                  NEMOCITHICK(IJS:IJL)
          ENDIF
          IF (LWNEMOCOUCUR) THEN
             IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
               U(IJS:IJL)=NEMOUCUR(IJS:IJL)
               V(IJS:IJL)=NEMOVCUR(IJS:IJL)
             ENDIF
          ENDIF
        ENDIF

      ENDIF

      IF (LWNEMOCOUCIT.AND.(.NOT.LNEMOCITHICK)) THEN
        WRITE(IU06,*) ' --------------------------------'
        WRITE(IU06,*) ' LWNEMOCOUCIT ONLY MAKES SENSES  '
        WRITE(IU06,*) ' IF LIM IS ACTIVATED IN NEMO     '
        WRITE(IU06,*) ' --------------------------------'
        CALL ABORT1
      ENDIF

      IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',1,ZHOOK_HANDLE)

      END SUBROUTINE RECVNEMOFIELDS
