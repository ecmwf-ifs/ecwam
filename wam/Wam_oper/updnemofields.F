      SUBROUTINE UPDNEMOFIELDS(LWSTOKES, USTOKES, VSTOKES)


!****  *UPDNEMOFIELDS* - UPDATE FIELDS IN NEMO WITH WAVE INFORMATION

!      KRISTIAN MOGENSEN ECMWF    NOVEMBER 2012                    

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PASSES WAVE INFORMATION THROUGH TO
!          NEMO VIA THE NEMO SINGLE EXECUTABLE COUPLING INTERFACE

!*    INTERFACE.                                                        
!     ----------                                                        

!          LWSTOKES, USTOKES, VSTOKES : STOKES VELOCITIES

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

! MODULES NEED FOR GRID DEFINTION      
      USE YOWGRID  , ONLY : IJS, IJL
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_MODULE, ONLY : MPL_COMM
! DATA TO BE SEND
      USE YOWMEAN  , ONLY : NEMEAN, NFMEAN, PHIEPS, TAUOC
      USE YOWSPEC  , ONLY : USNEW, U10NEW
! CONSTANTS
      USE YOWPCONS , ONLY : EPSUS, EPSU10
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
! INFORMATION FO OPTIONAL DEBUGGING
      USE YOWCOUP  , ONLY : LWNEMOCOUDEBUG, KCOUSTEP
      USE YOWSTAT  , ONLY : CDATEF, CDTPRO
      IMPLICIT NONE
      LOGICAL :: LWSTOKES
      REAL, DIMENSION(IJS(1):IJL(1)) :: USTOKES, VSTOKES
      REAL, ALLOCATABLE, DIMENSION(:) :: MWP, CD, ZDUMU, ZDUMV
      INTEGER :: IG, IFCST, IJ
      REAL :: TAU
      REAL :: ZHOOK_HANDLE
 
      IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',0,ZHOOK_HANDLE)

      IG=1
      if (LWNEMOCOUDEBUG) THEN
         CALL DIFDATE(CDATEF,CDTPRO,IFCST)
      ELSE
         IFCST=0
      ENDIF

! MEAN WAVE PERIOD

      ALLOCATE(MWP(IJS(IG):IJL(IG)))
      DO IJ = IJS(IG),IJL(IG)
         IF (NFMEAN(IJ)/=0) THEN
            MWP(IJ) = 1.0/NFMEAN(IJ)
         ELSE
            MWP(IJ) = 0.0
         ENDIF
      ENDDO

! COMPUTE CD FOR TESTING PURPOSES

      ALLOCATE(CD(IJS(IG):IJL(IG)))
      DO IJ = IJS(IG),IJL(IG)
         TAU = MAX(USNEW(IJ)**2,EPSUS)
         CD(IJ) = TAU/MAX(U10NEW(IJ)**2,EPSU10)
      ENDDO

! STOKES DRIFT
      ALLOCATE(ZDUMU(IJS(IG):IJL(IG)),ZDUMV(IJS(IG):IJL(IG)))
      IF (LWSTOKES) THEN
         ZDUMU(IJS(IG):IJL(IG)) = USTOKES(IJS(IG):IJL(IG))
         ZDUMV(IJS(IG):IJL(IG)) = VSTOKES(IJS(IG):IJL(IG))
      ELSE
         ZDUMU(IJS(IG):IJL(IG)) = 0.0
         ZDUMV(IJS(IG):IJL(IG)) = 0.0
      ENDIF

      CALL NEMOGCMCOUP_WAM_UPDATE( IRANK-1, NPROC, MPL_COMM,
     &     IJL(IG)-IJS(IG)+1, 4.0*SQRT(NEMEAN(IJS(IG):IJL(IG))),
     &     MWP, PHIEPS(IJS(IG):IJL(IG)), TAUOC(IJS(IG):IJL(IG)), 
     &     CD(IJS(IG):IJL(IG)), U10NEW(IJS(IG):IJL(IG)),
     &     ZDUMU(IJS(IG):IJL(IG)), ZDUMV(IJS(IG):IJL(IG)),
     &     IFCST, LWNEMOCOUDEBUG )

      DEALLOCATE(MWP,CD,ZDUMU,ZDUMV)

      IF (LHOOK) CALL DR_HOOK('UPDNEMOFIELDS',1,ZHOOK_HANDLE)

      END
