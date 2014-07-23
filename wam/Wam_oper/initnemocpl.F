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

! MODULES NEED FOR GRID DEFINTION      
      USE YOWPARAM , ONLY : NGX, NGY, NBLO
      USE YOWGRID  , ONLY : IJS, IJL, NLONRGG
      USE YOWMAP   , ONLY : IXLG, KXLT
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
! ACCUMULATED STRESSES
      USE YOWCOUP  , ONLY : NEMOTAUX, NEMOTAUY, NEMONEW10, 
     &                      NEMOPHIF, NEMONTAU
! INTEGRATED PARAMETERS
      USE YOWMEAN  , ONLY : NPHIEPS, NTAUOC, NSWH, NMWP
! FOR PRINTING
      USE YOWTEST  , ONLY : IU06
! NEMO FIELDS ON WAVE GRID
      USE YOWNEMOFLDS,ONLY: NEMOSST, NEMOCICOVER, NEMOCITHICK, 
     &                      NEMOUCUR, NEMOVCUR, NEMOSTRN, NEMOUSTOKES,
     &                      NEMOVSTOKES      
! DR. HOOK
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

      IMPLICIT NONE
! RECEIVE DATA FROM NEMO
      LOGICAL :: LRECV
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
      ALLOCATE( NEMOTAUX(IJS(IG):IJL(IG)) )
      ALLOCATE( NEMOTAUY(IJS(IG):IJL(IG)) )
      ALLOCATE( NEMONEW10(IJS(IG):IJL(IG)) )
      ALLOCATE( NEMOPHIF(IJS(IG):IJL(IG)) )
      NEMOTAUX(:)  = 0.0
      NEMOTAUY(:)  = 0.0
      NEMONEW10(:) = 0.0
      NEMOPHIF(:) = 0.0
      NEMONTAU     = 0

! INTEGRATED PARAMETERS

      ALLOCATE(NPHIEPS(IJS(IG):IJL(IG)))
      ALLOCATE(NTAUOC(IJS(IG):IJL(IG)))
      ALLOCATE(NSWH(IJS(IG):IJL(IG)))
      ALLOCATE(NMWP(IJS(IG):IJL(IG)))
      

! CALL COUPLING INTERFACE.
! THE TASKS RANK CONVENTION IS 0 TO NPROC-1 SO WE SUBTRACT ONE FOR IRANK
! SINCE WE ONLY HAVE LOCAL OCEAN POINTS WE SPECIFY GLOBAL NUMBER OF
! POINTS AS THE SUM OF NLOGRGG
      
      CALL NEMOGCMCOUP_WAM_COUPINIT( IRANK-1, NPROC, MPL_COMM, 
     &     IJL(IG)-IJS(IG)+1, SUM(NLONRGG(1:NGY)),
     &     NLOCMSK(IJS(IG):IJL(IG)), NGLOIND(IJS(IG):IJL(IG)), 0)


! DEALLOCATE DATA NOT NEEDED ANYMORE

      DEALLOCATE( NLOCMSK, NGLOIND )

! ALLOCATE RECEIVE DATA IF NEEDED

      IF (LRECV) THEN
         IG=1
         ALLOCATE(NEMOSST(IJS(IG):IJL(IG)),
     &            NEMOCICOVER(IJS(IG):IJL(IG)),
     &            NEMOCITHICK(IJS(IG):IJL(IG)), 
     &            NEMOUCUR(IJS(IG):IJL(IG)),
     &            NEMOVCUR(IJS(IG):IJL(IG)))
      ENDIF

! ALLOCATE SEND DATA

      IG=1
      ALLOCATE(NEMOSTRN(IJS(IG):IJL(IG)))
      ALLOCATE(NEMOUSTOKES(IJS(IG):IJL(IG)))
      ALLOCATE(NEMOVSTOKES(IJS(IG):IJL(IG)))
      
      IF (LHOOK) CALL DR_HOOK('INITNEMOCPL',1,ZHOOK_HANDLE)


      END

