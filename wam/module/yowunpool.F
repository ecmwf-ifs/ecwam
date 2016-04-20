!**********************************************************************
!*                                                                    *
!**********************************************************************
      MODULE YOWUNPOOL 
      USE YOWPARAM , ONLY : MDC => NANG,
     &                      MSC => NFRE 

         LOGICAL :: LLUNSTR
         LOGICAL :: LPREPROC
         LOGICAL :: LLUNBINOUT ! controls output into erghs.bin, ergtm.bin and ergtest.bin

         REAL*8,  PARAMETER     :: ZERO = 0.d0
         REAL*8,  PARAMETER     :: ONE = 1.d0
         REAL*8,  PARAMETER     :: TWO = 2.d0
         REAL*8,  PARAMETER     :: ONEHALF = 0.5d0
         REAL*8, PARAMETER      :: PI        = 3.141592653589793d0
         REAL*8, PARAMETER      :: PI2       = 2.d0*PI
         REAL*8, PARAMETER      :: INVPI     = 1.d0/PI
         REAL*8, PARAMETER      :: INVPI2    = 1.d0/PI2

         REAL*8, PARAMETER      :: REARTH    = 6378137.d0 ! WGS84
         REAL*8, PARAMETER      :: DEGRAD    = PI/180.d0
         REAL*8, PARAMETER      :: RADDEG    = 180.d0/PI

         REAL*8, PARAMETER      :: THR      = TINY(1.)
         REAL*8, PARAMETER      :: INVTHR   = 1.d0/TINY(1.)

         REAL*8, PARAMETER      :: THR8      = TINY(1.0d0)
         REAL*8, PARAMETER      :: INVTHR8   = 1.d0/TINY(1.0d0)

         REAL*8, PARAMETER      :: SMALL     = 1000.d0 * THR8
         REAL*8, PARAMETER      :: LARGE     = 1.d0/SMALL

         REAL*8,  PARAMETER     :: ONESIXTH = 1.d0/6.d0
         REAL*8,  PARAMETER     :: ONETHIRD = 1.0d0/3.0d0
         REAL*8,  PARAMETER     :: TWOTHIRD = 2.d0/3.d0

         integer :: OUT_METHOD = 2
         integer MaxLen
         integer NIBLO_FD, NIBLO_OUT
         
         
!         INTEGER                :: MNP ! Number of nodes ...
!         INTEGER                :: MNE ! Number of elements ...
         INTEGER                :: MAXMNECON ! max. number of connected elements on all nodes ...

         !INTEGER, ALLOCATABLE   :: INE(:,:) ! element connection table ..
         INTEGER, ALLOCATABLE   :: CCON(:)  ! number of connected elements per node ....
         INTEGER, ALLOCATABLE   :: IE_CELL(:) ! element pointer in an ordered way. points do the needed element when looping over all nodes and all connectes elements 
         INTEGER, ALLOCATABLE   :: POS_CELL(:) ! position of the above element in the connection table ...
         INTEGER, ALLOCATABLE   :: POS_CELL2(:,:) ! position of the above element in the connection table ...
         INTEGER, ALLOCATABLE   :: IE_CELL2(:,:) ! element pointer in an ordered way. points do the needed element when looping over all nodes and all connectes elements
         INTEGER, ALLOCATABLE   :: CELLVERTEX(:,:,:) ! Temp array for pointer creation
         INTEGER, ALLOCATABLE   :: ITER_EXP(:,:) ! needed amount of iterations to fullfill stability for each direction and each freq. 
         INTEGER, ALLOCATABLE   :: ITER_EXPD(:) ! needed amount of iterations to fullfill stability for each freq.
         INTEGER                :: ITER_MAXA ! max. value of the above array
         INTEGER, ALLOCATABLE   :: IE_OUTPTS(:,:) ! ELEMENT NUMBER CORRESPONDING TO ALL STRUCTURED GRID OUTPUT POINTS

!         REAL*8, ALLOCATABLE    :: XP(:), YP(:) ! x.-y. coordinates 

         REAL*8,  ALLOCATABLE   :: SI(:) ! median dual patch area for each node 
         REAL*8,  ALLOCATABLE   :: TRIA(:) ! triangle area of each element 
         REAL*8,  ALLOCATABLE   :: IEN(:,:) ! normal vectors for each element in x and y direction, inward pointing normalized by the edge length 

         REAL*8, ALLOCATABLE    :: CURTXY(:,:) ! current vel. x and y direction for each node ...
         REAL*8, ALLOCATABLE    :: CFLCXY(:,:) ! cfl nunmber of each direction and freq. 

         INTEGER                :: IVECTOR 

         INTEGER, SAVE          :: ITER_MAX

         LOGICAL                :: LCALC = .TRUE.  ! if true do cfl number 
         LOGICAL                :: LCFL = .FALSE.  ! later use for writing cfl numbers for 
         LOGICAL                :: LCUR = .FALSE.  ! prepared for taking currnets into account 
         LOGICAL                :: LADVTEST = .FALSE. ! switch for testing advection ... 
         LOGICAL                :: LSPHE = .TRUE. ! spherical coordinates 
         LOGICAL                :: LVECTOR = .FALSE. 

!         REAL*8, ALLOCATABLE    :: WBACOLD(:,:,:)
!         REAL*8, ALLOCATABLE    :: DSPEC(:,:,:)

         REAL*8, ALLOCATABLE    :: CG(:,:), WK(:,:) ! group vel. and wave numbers ...

         TYPE FILEDEF
            CHARACTER(LEN=40)   :: FNAME ! unWAM files names ... 
            INTEGER             :: FHNDL ! unWAM file handles ...
         END TYPE

         TYPE (FILEDEF)         :: BND,WAV,STAT,DBG,WAMOUT,GRID,XFN_HS,
     &                             XFN_TM, XFN_TEST
      END MODULE YOWUNPOOL
!**********************************************************************
!*                                                                    *
!**********************************************************************
