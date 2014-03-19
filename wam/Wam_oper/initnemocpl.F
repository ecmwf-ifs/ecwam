      SUBROUTINE INITNEMOCPL

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

! MODULES NEED FOR GRID DEFINTION      
      USE YOWPARAM , ONLY : NGX, NGY
      USE YOWGRID  , ONLY : IJS, IJL, NLONRGG
      USE YOWMAP   , ONLY : IXLG, KXLT
     &                    , AMOWEP,AMOSOP,XDELLA,ZDELLO
      ! AMOWEP ETC IS FOR INITIAL TESTING AND SHOULD BE REMOVED!!
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_MODULE, ONLY : MPL_COMM
! FOR PRINTING
      USE YOWTEST  , ONLY : IU06
! DR. HOOK
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

      IMPLICIT NONE
! GLOBAL INDEX DEFINITON.      
      INTEGER :: IGLOBAL(NGX,NGY)
! INDICES
      INTEGER :: IG,IJ,IX,JSN,I,J,K
! TEMPORARY ARRAYS FOR LOCAL MASK AND GLOBAL INDICES
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NLOCMSK, NGLOIND
! MISC
      REAL :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('INITNEMOCPL',0,ZHOOK_HANDLE)

      WRITE(IU06,*)
      WRITE(IU06,*)' **************************************************'
      WRITE(IU06,*)' * INITIALIZING THE WAM <-> NEMO COUPLING         *'
      WRITE(IU06,*)' **************************************************'
      WRITE(IU06,*)
      CALL FLUSH(IU06)

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

! CALL COUPLING INTERFACE.
! THE TASKS RANK CONVENTION IS 0 TO NPROC-1 SO WE SUBTRACT ONE FOR IRANK
! SINCE WE ONLY HAVE LOCAL OCEAN POINTS WE SPECIFY GLOBAL NUMBER OF
! POINTS AS THE SUM OF NLOGRGG
      
      CALL NEMOGCMCOUP_WAM_COUPINIT( IRANK-1, NPROC, MPL_COMM, 
     &     IJL(IG)-IJS(IG)+1, SUM(NLONRGG(1:NGY)),
     &     NLOCMSK(IJS(IG):IJL(IG)), NGLOIND(IJS(IG):IJL(IG)), 0)


! DEALLOCATE DATA NOT NEEDED ANYMORE

      DEALLOCATE( NLOCMSK, NGLOIND )
      
      IF (LHOOK) CALL DR_HOOK('INITNEMOCPL',1,ZHOOK_HANDLE)

      END

