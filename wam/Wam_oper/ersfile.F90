      SUBROUTINE ERSFILE (IU91, IU92, LLOUTPUT)
                                                                        
! -------------------------------------------------------------------   

!**** *ERSFILE* - FILE MANAGEMENT OF ERS FILES.                      

!     H. GUNTHER               ECMWF         JULY 1991                  

!     PURPOSE.                                                          
!     --------                                                          

!       CONTROL INPUT AND OUTPUT FILES OF ERS SPECTRAL DATA.          


!*** INTERFACE.                                                         
!    ----------                                                         
!      *CALL* *ERSFILE (IU91, IU92 LLOUTPUT)*                            
!         *IU91*   - INPUT UNIT OF LCATIOB FILE.                        
!         *IU92*   - OUTPUT UNIT FOR SPECTRA.                           
!         *LLOUTPUT* - IF TRUE OUTPUT IS DONE.

!     EXTERNALS.                                                        
!     ----------                                                        

!        NONE.                                                          


!     METHOD.                                                           
!     -------                                                           

!       NONE.                                                           

!     REFERENCES.                                                       
!     -----------                                                       

!        NONE.                                                          

!     MODIFICATIONS.                                                    
!     --------------                                                    

!        B. HANSEN    *ECMWF*      FEB 93'                              
!           ACCESS TO ECFILE DELETED. INPUT FILES COLYYMMDDhhmm HAVE    
!           TO BE ONLINE ALREADY OTHERWISE NO OUTPUT OF WAVE SPECTRA    
!           CLOSE TO POSITIONS OF ERS1-SAR-SPECTRA.                     
!           WAVE SPECTRA WRITTEN TO UNIT 71. FILES ARE OPEND FOLLOWING  
!           THE KNOWN NAMING KONVENTION COSYYMMDDhhmm. THESE FILES      
!           ARE CLOSED (STATUS=KEEP) AND ARE TO BE SAVED AFTER THE      
!           MODEL RUN.                                                  

!     J. BIDLOT       *ECMWF*      JUNE 1996 MIGRATION TO FUJITSU

!---------------------------------------------------------------------- 

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOER  , ONLY : NERS     ,CDTERS   ,IDELERS  ,IERS     ,    &
     &            IJERS    ,IGERS
      USE YOWSTAT  , ONLY : CDATEF   ,CDTPRO   ,CDATEA   ,MARSTYPE
      USE YOWTEST  , ONLY : IU06
      USE YOWMPP   , ONLY : IRANK

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "incdate.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU91, IU92
      LOGICAL, INTENT(IN) :: LLOUTPUT

      INTEGER(KIND=JWIM) :: I, JC
      INTEGER(KIND=JWIM) :: IERR, IOS, NERSOLD
      INTEGER(KIND=JWIM), ALLOCATABLE :: IDUM(:)

      REAL(KIND=JWRB) :: ZLATIS, ZLONGS, ZLATIW, ZLONGW
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      CHARACTER(LEN= 2) :: CLH       ! HELPER STRING TO STORE HOUR.            
      CHARACTER(LEN=15) :: FNAME
      CHARACTER(LEN=12) :: CL_DATE
      CHARACTER(LEN=14) :: CL_TIME

      LOGICAL, SAVE :: FRSTIME
      DATA FRSTIME / .TRUE./

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('ERSFILE',0,ZHOOK_HANDLE)

!*    1.  SAVE SPECTRA FILE.                                            
!         ------------------                                            

      IF (LLOUTPUT) THEN

        IF (ALLOCATED(IJERS)) DEALLOCATE(IJERS)
        IF (ALLOCATED(IGERS)) DEALLOCATE(IGERS)

!       DO NOT DEAL WITH IU91 AND IU92 IF NOT ON PE 1
        IF (IRANK == 1) THEN
          CLOSE (IU92,STATUS='KEEP')
        ENDIF
        IERS = 0
        IF (LHOOK) CALL DR_HOOK('ERSFILE',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

! ----------------------------------------------------------------------

!*    2. INITIAL DATES IF FIRST CALL.                                   
!        ----------------------------                                   

      IF (FRSTIME) THEN
        NERS = 1000
        IERS = 0
        IDELERS = 21600
        CDTERS  = CDATEA
        CLH=CDATEA(9:10)
        IF (CLH == '09' .OR. CLH == '15' .OR.                           &
     &      CLH == '21' .OR. CLH == '03'     ) THEN
          CALL INCDATE (CDTERS, -10800)
        ENDIF
        FRSTIME = .FALSE.
      ENDIF                                                             

      DO WHILE (CDTERS.LT.CDTPRO)
        CALL INCDATE (CDTERS, IDELERS)
      ENDDO

! ----------------------------------------------------------------------

!*    3. FETCH LOCATION FILE.                                           
!        --------------------                                           

!     DO NOT DEAL WITH IU91 AND IU92 IF NOT ON PE 1
      IF (IRANK /= 1) THEN
        IF (LHOOK) CALL DR_HOOK('ERSFILE',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

      IF (.NOT.ALLOCATED(IJERS)) ALLOCATE(IJERS(NERS))
      IF (.NOT.ALLOCATED(IGERS)) ALLOCATE(IGERS(NERS))

!     DO NOT DEAL WITH IU91 AND IU92 IF IT IS NOT COLLOCATION TIME 
      IF (CDTPRO /= CDTERS) THEN
        IERS=0
        IF (LHOOK) CALL DR_HOOK('ERSFILE',1,ZHOOK_HANDLE)
        RETURN
      ENDIF

      IERS = 0
      IF (CDTERS <= CDATEF .OR. MARSTYPE == 'an' .OR.                   &
     &    MARSTYPE == 'fg' .OR. MARSTYPE == '4v'      ) THEN
        WRITE(FNAME,'(''COL'',A12)') CDTERS(1:12)
        OPEN (UNIT=IU91, FILE=FNAME, FORM='UNFORMATTED', STATUS='OLD',  &
     &   IOSTAT=IERR)
        IF (IERR /= 0) THEN
          IF (LHOOK) CALL DR_HOOK('ERSFILE',1,ZHOOK_HANDLE)
          RETURN
        ELSE
          READ(IU91)CL_DATE
          IF (IRANK == 1) THEN
            WRITE (IU06,*) '                                     '
            WRITE (IU06,*) '      SUB. ERSFILE: FILE ',                 &
     &       FNAME, ' OPENED    '
          ENDIF

          WRITE(FNAME,'(''COS'',A12)') CDTERS(1:12)
          OPEN (UNIT=IU92, FILE=FNAME(1:LEN_TRIM(FNAME))//MARSTYPE,     &
     &     FORM='UNFORMATTED', STATUS='UNKNOWN')
          WRITE(IU92) CDTPRO
        ENDIF

 3001   CONTINUE
        IERS = IERS + 1

        IF (IERS > NERS) THEN
          IF (IRANK == 1) THEN
            WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++'
            WRITE (IU06,*) ' +                                     +'
            WRITE (IU06,*) ' +        WARNING SUB ERSFILE.         +'
            WRITE (IU06,*) ' +     ==========================      +'
            WRITE (IU06,*) ' + NUMBER OF OUTPUT POINTS             +'
            WRITE (IU06,*) ' + EXCEEDS DIMENSION NERS = ', NERS
            WRITE (IU06,*) ' + THE DIMENSION WILL BE DOUBLED       +'
            WRITE (IU06,*) ' + IN ORDER TO CONTNUE                 +'
            WRITE (IU06,*) ' +                                     +'
            WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++'
          ENDIF
          NERSOLD=NERS
          NERS=NERS*2
          ALLOCATE(IDUM(NERSOLD))
          IDUM=IJERS
          DEALLOCATE(IJERS)
          ALLOCATE(IJERS(NERS))
          DO JC=1,NERSOLD
            IJERS(JC)=IDUM(JC)
          ENDDO

          IDUM=IGERS
          DEALLOCATE(IGERS)
          ALLOCATE(IGERS(NERS))
          DO JC=1,NERSOLD
            IGERS(JC)=IDUM(JC)
          ENDDO
        
          DEALLOCATE(IDUM)

        ENDIF

        READ (IU91, END=3002, ERR=3002, IOSTAT=IOS) CL_TIME,            &
     &   ZLATIS, ZLONGS, ZLATIW, ZLONGW, IGERS(IERS), IJERS(IERS)
        IF (IGERS(IERS) == 0) THEN
          IERS = IERS-1
          GOTO 3001
        ENDIF
        IF (IERS > 1) THEN
          DO I=1,IERS-1
            IF (IGERS(I) == IGERS(IERS) .AND.                           &
     &       IJERS(I).EQ.IJERS(IERS)) THEN
              IERS = IERS-1
              GOTO 3001
            ENDIF
          ENDDO
        ENDIF
        GOTO 3001
 3002   CONTINUE
        IERS = IERS - 1
        CLOSE (IU91, STATUS='KEEP')
      ENDIF

      IF (LHOOK) CALL DR_HOOK('ERSFILE',1,ZHOOK_HANDLE)

      END SUBROUTINE ERSFILE
