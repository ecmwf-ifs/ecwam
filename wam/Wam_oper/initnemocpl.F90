      SUBROUTINE INITNEMOCPL(LRECV)

!****  *INITNEMOCPL* - INITIALIZE THE COUPLING OF WAM TO NEMO

!      KRISTIAN MOGENSEN ECMWF    NOVEMBER 2012                    

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PREPARES FOR THE COUPLING OF WAM
!          TO NEMO BY DEFINING A GLOBAL GRID INDEX NEEDED TO
!          CALL NEMOGCMCOUP_WAM_COUPINIT         

!*    INTERFACE.                                                        
!     ----------                                                        

!          THE NEMO COUPLING IN nemo/coupled/src/nemointerface

!     METHOD.                                                           
!     -------                                                           

!          CREATE A GLOBAL INDEX AND NEMOGCMCOUP_WAM_COUPINIT

!     EXTERNALS.                                                        
!     ----------                                                        

!          NEMOGCMCOUP_WAM_COUPINIT  -  INITIALIZES INTERNAL
!                      DATA STRUCTURES FOR THE COUPLING

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

! MODULES NEED FOR KIND DEFINITION
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWNEMOP    , ONLY : NEMODP

! MODULES NEED FOR GRID DEFINTION      
      USE YOWPARAM , ONLY : NGX, NGY, NBLO
      USE YOWGRID  , ONLY : IJS, IJL, NLONRGG
      USE YOWMAP   , ONLY : IXLG, KXLT
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
! WAVE FIELDS TO NEMO 
      USE YOWCOUP  , ONLY : NPHIEPS, NTAUOC, NSWH, NMWP,                &
     &                      NEMOTAUX, NEMOTAUY, NEMONEW10,              &
     &                      NEMOPHIF, NEMONTAU, NEMOSTRN, NEMOUSTOKES,  &
     &                      NEMOVSTOKES
! FOR PRINTING
      USE YOWTEST  , ONLY : IU06
! NEMO FIELDS ON WAVE GRID
      USE YOWNEMOFLDS,ONLY: NEMOSST, NEMOCICOVER, NEMOCITHICK,          &
     &                      NEMOUCUR, NEMOVCUR
! DR. HOOK
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! -------------------------------------------------------------------   

      IMPLICIT NONE
#include "abort1.intfb.h"

! RECEIVE DATA FROM NEMO
      LOGICAL, INTENT(IN) :: LRECV
! GLOBAL INDEX DEFINITON.      
      INTEGER(KIND=JWIM) :: IGLOBAL(NGX,NGY)
! INDICES
      INTEGER(KIND=JWIM) :: IG, IJ, IX, JSN, I, J, K
! TEMPORARY ARRAYS FOR LOCAL MASK AND GLOBAL INDICES
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: NLOCMSK, NGLOIND
! MISC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! -------------------------------------------------------------------   

      IF (LHOOK) CALL DR_HOOK('INITNEMOCPL',0,ZHOOK_HANDLE)

      WRITE(IU06,*)
      WRITE(IU06,*)' **************************************************'
      WRITE(IU06,*)' * INITIALIZING THE WAM <-> NEMO COUPLING         *'
      WRITE(IU06,*)' **************************************************'
      WRITE(IU06,*)
      CALL FLUSH(IU06)

      IF (NBLO/=1) THEN  
         WRITE(IU06,*) ' ------------------------------------------'
         WRITE(IU06,*) ' NEMO COUPLING ASSUMES THAT NBLO==1        '
         WRITE(IU06,*) ' ------------------------------------------'
         CALL ABORT1
      ENDIF

! CONSTRUCT GLOBAL INDICES GRID.
! THE CONVENTION IS NORTH TO SOUTH AND WEST TO EAST 
! STARTING AT 0 DEG LONGITUDE

      IGLOBAL(:,:)=-1
      K=0
      DO J=1,NGY
         DO I=1,NLONRGG(NGY-J+1)
            K=K+1
            IGLOBAL(I,NGY-J+1)=K
         ENDDO
      ENDDO

! SETUP GLOBAL INDICES 1D ARRAY + LOCAL MASK
! SINCE WE ONLY HAVE OCEAN POINTS THE LOCAL MASK IS ONE FOR ALL POINTS

      IG=1
      ALLOCATE( NLOCMSK(IJS(IG):IJL(IG)), NGLOIND(IJS(IG):IJL(IG)) )
      DO IJ=IJS(IG),IJL(IG)
         IX          = IXLG(IJ,IG)
         JSN         = KXLT(IJ,IG)
         NLOCMSK(IJ) = 1
         NGLOIND(IJ) = IGLOBAL(IX,JSN)
      ENDDO

! CHECK FOR GLOBAL INDICES LESS THAN ZERO
      IF (MINVAL(NGLOIND)<1) THEN
         WRITE(0,*)'GLOBAL INDEX LESS THAN 1 FOUND!!!'
         DO IJ=IJS(IG),IJL(IG)
            IX          = IXLG(IJ,IG)
            JSN         = KXLT(IJ,IG)
            WRITE(0,*) IJ, IX, JSN, NLONRGG(JSN), NGLOIND(IJ)
         ENDDO
         CALL ABORT1
      ENDIF         

! ALLOCATE ACCUMULATED U,V STRESS
      NEMONTAU     = 0
      ALLOCATE( NEMOTAUX(IJS(IG):IJL(IG)) )
      NEMOTAUX(:)  = 0.0_NEMODP
      ALLOCATE( NEMOTAUY(IJS(IG):IJL(IG)) )
      NEMOTAUY(:)  = 0.0_NEMODP
      ALLOCATE( NEMONEW10(IJS(IG):IJL(IG)) )
      NEMONEW10(:) = 0.0_NEMODP
      ALLOCATE( NEMOPHIF(IJS(IG):IJL(IG)) )
      NEMOPHIF(:)  = 0.0_NEMODP

! INTEGRATED PARAMETERS

      ALLOCATE(NSWH(IJS(IG):IJL(IG)))
      NSWH(:) = 0.0_NEMODP
      ALLOCATE(NMWP(IJS(IG):IJL(IG)))
      NMWP(:) = 0.0_NEMODP
      ALLOCATE(NPHIEPS(IJS(IG):IJL(IG)))
      NPHIEPS(:) = 0.0_NEMODP
      ALLOCATE(NTAUOC(IJS(IG):IJL(IG)))
      NTAUOC(:) = 0.0_NEMODP

! CALL COUPLING INTERFACE.
! THE TASKS RANK CONVENTION IS 0 TO NPROC-1 SO WE SUBTRACT ONE FOR IRANK
! SINCE WE ONLY HAVE LOCAL OCEAN POINTS WE SPECIFY GLOBAL NUMBER OF
! POINTS AS THE SUM OF NLOGRGG

#ifdef WITH_NEMO
      CALL NEMOGCMCOUP_WAM_COUPINIT( IRANK-1, NPROC, MPL_COMM,          &
     &     IJL(IG)-IJS(IG)+1, SUM(NLONRGG(1:NGY)),                      &
     &     NLOCMSK(IJS(IG):IJL(IG)), NGLOIND(IJS(IG):IJL(IG)), 0)
#endif


! DEALLOCATE DATA NOT NEEDED ANYMORE

      DEALLOCATE( NLOCMSK, NGLOIND )
#ifdef ECMWF

! ALLOCATE RECEIVE DATA IF NEEDED

      IF (LRECV) THEN
         IG=1
         ALLOCATE(NEMOSST(IJS(IG):IJL(IG)),                             &
     &            NEMOCICOVER(IJS(IG):IJL(IG)),                         &
     &            NEMOCITHICK(IJS(IG):IJL(IG)),                         &
     &            NEMOUCUR(IJS(IG):IJL(IG)),                            &
     &            NEMOVCUR(IJS(IG):IJL(IG)))
      ENDIF

! ALLOCATE SEND DATA

      IG=1
      ALLOCATE(NEMOSTRN(IJS(IG):IJL(IG)))
      NEMOSTRN(:) = 0.0_NEMODP
      ALLOCATE(NEMOUSTOKES(IJS(IG):IJL(IG)))
      NEMOUSTOKES(:) = 0.0_NEMODP
      ALLOCATE(NEMOVSTOKES(IJS(IG):IJL(IG)))
      NEMOVSTOKES(:) = 0.0_NEMODP
      
#endif

      IF (LHOOK) CALL DR_HOOK('INITNEMOCPL',1,ZHOOK_HANDLE)

      END SUBROUTINE INITNEMOCPL
