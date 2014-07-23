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

! MODULES NEED FOR GRID DEFINTION      
      USE YOWGRID  , ONLY : IJS, IJL
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
      USE MPL_MODULE
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! COUPLING INFORMATION
      USE YOWCOUP  , ONLY : LWNEMOCOUCIC, LWNEMOCOUCIT, LWNEMOCOUCUR
! ICE AND CURRENT INFORMATION 
      USE YOWICE   , ONLY : CITHICK, CICOVER
      USE YOWCURR  , ONLY : U        ,V
! OUTPUT FORTRAN UNIT
      USE YOWTEST  , ONLY : IU06
! NEMO FIELDS ON WAVE GRID
      USE YOWNEMOFLDS,ONLY: NEMOSST, NEMOCICOVER, NEMOCITHICK, 
     &                      NEMOUCUR, NEMOVCUR, LNEMOCITHICK,
     &                      LNEMOICEREST
      IMPLICIT NONE
      LOGICAL :: LREST ! RESTART SO UPDATE FROM RESTART VALUES
      LOGICAL :: LINIT ! UPDATE CICOVER and CITHICK AT INITIAL TIME
                       ! FOR OUTPUT
      REAL :: ZHOOK_HANDLE
      INTEGER :: IG,IP,JSTORE

      IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',0,ZHOOK_HANDLE)

      ! Temporary debugging.
      WRITE(IU06,*)'WAM: Receiving NEMO files. LREST IS ',LREST

      IG=1
      ! IF WE ARE IN RESTART JUST COPY THE RESTART FILE INFO.
      IF (LREST) THEN
        NEMOCICOVER(IJS(IG):IJL(IG))=CICOVER(IJS(IG):IJL(IG),IG)
        NEMOCITHICK(IJS(IG):IJL(IG))=CITHICK(IJS(IG):IJL(IG),IG)
        IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
          NEMOUCUR(IJS(IG):IJL(IG))=U(IJS(IG):IJL(IG),IG)
          NEMOVCUR(IJS(IG):IJL(IG))=V(IJS(IG):IJL(IG),IG)
        ELSE
          NEMOUCUR(IJS(IG):IJL(IG))=0.0
          NEMOVCUR(IJS(IG):IJL(IG))=0.0
        ENDIF
        LNEMOICEREST=.TRUE.
        LNEMOCITHICK=.TRUE.
      ELSE
        CALL NEMOGCMCOUP_WAM_GET( IRANK-1, NPROC, MPL_COMM, 
     &                            IJL(IG)-IJS(IG)+1,
     &                            NEMOSST(IJS(IG):IJL(IG)),
     &                            NEMOCICOVER(IJS(IG):IJL(IG)),
     &                            NEMOCITHICK(IJS(IG):IJL(IG)),
     &                            NEMOUCUR(IJS(IG):IJL(IG)),
     &                            NEMOVCUR(IJS(IG):IJL(IG)), 
     &                            LNEMOCITHICK )
        LNEMOICEREST=.FALSE.
      ENDIF

      WRITE(IU06,*)'ICECOVER:',SUM(NEMOCICOVER(IJS(IG):IJL(IG)))
      WRITE(IU06,*)'ICETHICK:',SUM(NEMOCITHICK(IJS(IG):IJL(IG)))

      IF (LINIT) THEN
        WRITE(IU06,*)'RECVNEMOFIELDS INIT'
        IF (LWNEMOCOUCIC) THEN
          CICOVER(IJS(IG):IJL(IG),IG)=NEMOCICOVER(IJS(IG):IJL(IG))
        ENDIF
        IF (LWNEMOCOUCIT) THEN
          CITHICK(IJS(IG):IJL(IG),IG)=NEMOCICOVER(IJS(IG):IJL(IG))*
     &                                NEMOCITHICK(IJS(IG):IJL(IG))
        ENDIF
        IF (LWNEMOCOUCUR) THEN
           IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
              U(IJS(IG):IJL(IG),IG)=NEMOUCUR(IJS(IG):IJL(IG))
              U(IJS(IG):IJL(IG),IG)=NEMOVCUR(IJS(IG):IJL(IG))
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

      END
