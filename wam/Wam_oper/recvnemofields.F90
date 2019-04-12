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
      USE YOWCURR  , ONLY : U        ,V
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
      INTEGER(KIND=JWIM) :: IG
      INTEGER(KIND=JWIM) :: IX, IY, IJ 
! -------------------------------------------------------------------   

      IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',0,ZHOOK_HANDLE)

      IG=1
      ! IF WE ARE IN RESTART JUST COPY THE RESTART FILE INFO.
      IF (LREST) THEN
        NEMOCICOVER(IJS(IG):IJL(IG))=CICOVER(IJS(IG):IJL(IG),IG)
        NEMOCITHICK(IJS(IG):IJL(IG))=CITHICK(IJS(IG):IJL(IG),IG)
        IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
          NEMOUCUR(IJS(IG):IJL(IG))=U(IJS(IG):IJL(IG),IG)
          NEMOVCUR(IJS(IG):IJL(IG))=V(IJS(IG):IJL(IG),IG)
        ELSE
          NEMOUCUR(IJS(IG):IJL(IG))=0.0_JWRB
          NEMOVCUR(IJS(IG):IJL(IG))=0.0_JWRB
        ENDIF
        LNEMOICEREST=.TRUE.
        LNEMOCITHICK=.TRUE.
      ELSE
#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_GET( IRANK-1, NPROC, MPL_COMM,             &
     &                            IJL(IG)-IJS(IG)+1,                    &
     &                            NEMOSST(IJS(IG):IJL(IG)),             &
     &                            NEMOCICOVER(IJS(IG):IJL(IG)),         &
     &                            NEMOCITHICK(IJS(IG):IJL(IG)),         &
     &                            NEMOUCUR(IJS(IG):IJL(IG)),            &
     &                            NEMOVCUR(IJS(IG):IJL(IG)),            &
     &                            LNEMOCITHICK )

        LLNEWCURR=.TRUE. 
#endif
        LNEMOICEREST=.FALSE.

!       MASK CURRENTS IF SEA ICE FROM NEMO
        IF (LWNEMOCOUCUR) THEN
          IF (LWNEMOCOUCIC) THEN
            DO IJ = IJS(IG),IJL(IG)
              IF(NEMOCICOVER(IJ,IG).GT.0.0_JWRB) THEN
                NEMOUCUR(IJ,IG)=0.0_JWRB
                NEMOVCUR(IJ,IG)=0.0_JWRB
              ENDIF
            ENDDO
          ENDIF 
        ENDIF 

      ENDIF

!     UPDATE CICOVER, CITHICK AND U and V CURRENTS AT INITIAL TIME ONLY !!!!
      IF (LINIT) THEN

        WRITE(IU06,*)' RECVNEMOFIELDS: INITIALISE OCEAN FIELDS'

        IF(LWCOU) THEN
          IF (LWNEMOCOUCIC) THEN
            DO IJ = IJS(IG),IJL(IG)
              IX = IFROMIJ(IJ,IG)
              IY = JFROMIJ(IJ,IG)
!            if lake cover = 0, we assume open ocean point, then get sea ice directly from NEMO
              IF (FIELDG(IX,IY)%LKFR .LE. 0.0_JWRB ) THEN
                CICOVER(IJ,IG)=NEMOCICOVER(IJ)
              ELSE
!            get ice information from atmopsheric model
                CICOVER(IJ,IG)=FIELDG(IX,IY)%CIFR 
              ENDIF
            ENDDO
          ENDIF

          IF (LWNEMOCOUCIT) THEN
            DO IJ = IJS(IG),IJL(IG)
              IX = IFROMIJ(IJ,IG)
              IY = JFROMIJ(IJ,IG)
!            if lake cover = 0, we assume open ocean point, then get sea ice thickness directly from NEMO
              IF (FIELDG(IX,IY)%LKFR .LE. 0.0_JWRB ) THEN
                CITHICK(IJ,IG)=NEMOCICOVER(IJ)*NEMOCITHICK(IJ)
              ELSE
                CICOVER(IJ,IG)=0.5_JWRB*NEMOCICOVER(IJ)
              ENDIF
            ENDDO
          ENDIF

          IF (LWNEMOCOUCUR) THEN
            IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
              DO IJ = IJS(IG),IJL(IG)
                IX = IFROMIJ(IJ,IG)
                IY = JFROMIJ(IJ,IG)
!              if lake cover = 0, we assume open ocean point, then get currents directly from NEMO
                IF (FIELDG(IX,IY)%LKFR .LE. 0.0_JWRB ) THEN
                  U(IJ,IG)=NEMOUCUR(IJ)
                  V(IJ,IG)=NEMOVCUR(IJ)
                ELSE
                  U(IJ,IG)=0.0_JWRB
                  V(IJ,IG)=0.0_JWRB
                ENDIF
              ENDDO

            ENDIF
          ENDIF

        ELSE

          IF (LWNEMOCOUCIC) THEN
            CICOVER(IJS(IG):IJL(IG),IG)=NEMOCICOVER(IJS(IG):IJL(IG))
          ENDIF
          IF (LWNEMOCOUCIT) THEN
            CITHICK(IJS(IG):IJL(IG),IG)=NEMOCICOVER(IJS(IG):IJL(IG))*     &
     &                                  NEMOCITHICK(IJS(IG):IJL(IG))
          ENDIF
          IF (LWNEMOCOUCUR) THEN
             IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
               U(IJS(IG):IJL(IG),IG)=NEMOUCUR(IJS(IG):IJL(IG))
                V(IJS(IG):IJL(IG),IG)=NEMOVCUR(IJS(IG):IJL(IG))
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
