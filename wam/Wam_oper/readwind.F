      SUBROUTINE READWIND(CDTWIR,IUGB,IVGB,ILEN,lWNDFILE,LCLOSEWND)

!***  *READWIND* - PROGRAM TO REFORMAT WIND FIELDS FROM MARS 
!                 EITHER FROM A FILE OR FROM GRIB DATA IN IUGB, IVGB

!    PETER JANSSEN      ECMWF        FEBRUARY 1987                


!    MODIFIED BY:                                                 
!    ------------                                                       

!    LIANA ZAMBRESKY     GKSS/ECMWF   JULY 1988                  

!        1. PASSING OF ARRAYS AS FORMAL PARMETERS TO SAVE MEMORY
!        2. PROCESSING WINDS FOR ONLY ONE WIND TIME STEP       

!    PEDRO VITERBO       ECMWF        OCTOBER 1988            

!        DEALING WITH MARS PACKED DATA.                      

!    H. GUNTHER      ECMWF/GKSS       OCTOBER 1991          

!        GRIB EDITION 1.                                   

!     J. DOYLE    NRL/ECMWF     OCT. 1996
!                 DECODES GRID WIND FIELDS WHICH ARE PASSED AS ARRAYS
!                 WHEN MODEL IS COUPLED TO ATMOPHERIC MODEL
!                 SEE LWNDFILE

!    B. HANSEN    ECMWF 1997
!                 RENAME SUBROUTINE TO MATCH THE NAME OF THE FILE.
!                 USE INCLUDE FILES FOR ALL COMMON BLOCKS.

!    J. BIDLOT    ECMWF 1998 RENAMED READWND_MARS AS READWIND

!     PURPOSE                                                       
!     -------                                                      

!     *READWIND*
!         READS TARGET FILES FROM THE MARS ARCHIVE AND INFERS FROM THIS 
!             A) THE DATE OF THE WIND FIELD                             
!             B) THE DEFINITION OF THE MARS GRID                        
!             C) THE U AND V COMPONENT OF THE WIND FIELD                
!             D) THE TYPE OF WIND FIELD AND                             
!             E) THE UNPACKED WINDS                                     

!     INTERFACE                                                         
!     ---------                                                         

!     *CALL* *READWIND (CDTWIR, IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)*
!        *CDTWIR* - DATE/TIME OF THE DATA READ                          
!        *IUGB*   - ARRAY CONTAINING THE GRIB U COMPONENT
!        *IVGB*   - ARRAY CONTAINING THE GRIB V COMPONENT
!        *ILEN*   - SIZE OF IUGB AND IVGB
!        *LWNDFILE*  IF TRUE DATA ARE ACCESSED FROM THE INPUT FILE
!                    sfcwind
!                    ELSE THE DATA WILL BE SUPPLIED IN IUGB AND IVGB
!        *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND THE UNIT 
!                    RESET

!     EXTERNALS                                                         
!     ---------                                                         

!     *ABORT1*         TERMINATES PROCESSING
!     *ADJUST*         ADJUSTS EAST AND WEST BOUNDARIES
!     *GRIBEX*         UNPACKS MARS WINDS                               
!     *GRPRS1*         PRINTS PRODUCT DEFINITION BLOCK                  
!     *GRPRS2*         PRINTS GRID  DEFINITION BLOCK                    
!     *INCDATE*        INCREMENTS DATE                                  
!     *PBCLOSE*
!     *PBGRIB*
!     *PBOPEN*

! --------------------------------------------------------------------- 

      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,  
     &            XDELLA
      USE YOWPARAM , ONLY : NGX      ,NGY      ,SWAMPWIND,CLDOMAIN
      USE YOWSTAT  , ONLY : CDATEA   ,IDELWI 
      USE YOWWNDG  , ONLY : DLAM     ,DPHI     ,RLATS    ,RLATN    ,
     &            RLONL    ,RLONR    ,KCOL     ,KROW     ,ICODE    ,
     &            IWPER    ,ICOORD
      USE YOWWIND  , ONLY : NC       ,NR       ,UWND     ,VWND  

! --------------------------------------------------------------------- 

      INTEGER, ALLOCATABLE ::  KSEC1(:), KSEC2(:), KSEC3(:), KSEC4(:)
      INTEGER, ALLOCATABLE :: KSEC0(:), KGRIB(:)
      REAL, ALLOCATABLE :: PSEC2(:), PSEC3(:), PSEC4(:)

      INTEGER IUGB(ILEN), IVGB(ILEN)
      CHARACTER HOPER*1, FILNM*24
      CHARACTER*12 CDTWIR
      LOGICAL FRSTIME, LWNDFILE, LCLOSEWND

      DATA HOPER, FRSTIME / 'R', .TRUE. /
      SAVE IUNIT1      

! --------------------------------------------------------------------  

!     for the swamp case create your own wind fields
      IF(CLDOMAIN.EQ.'s') goto 8000
 
      ALLOCATE(KSEC0(2), KSEC1(1024), KSEC2(1024))
      ALLOCATE(KSEC3(2), KSEC4(512))
      ALLOCATE(PSEC2(512), PSEC3(2))

      IF(LWNDFILE) THEN

!       WIND FIELDS ARE SUPPLIED FROM INPUT FILE sfcwindin
!       AS IS THE CASE IN THE STAND ALONE VERSION OR AT INITIAL INPUT
!       TIME IN COUPLED MODEL.

!       1.0 INITIALISE VALUES. 
!           ------------------ 

        IF(CLDOMAIN.EQ.'g' .AND. XDELLA.LE.0.5) THEN
          NBIT=130000
        ELSEIF(CLDOMAIN.EQ.'m' .AND. XDELLA.LE.0.25) THEN
          NBIT=130000
        ELSE
          NBIT=60000
        ENDIF

        ICOORD = 1

        IF (FRSTIME) THEN
          FILNM='sfcwindin'
          CALL PBOPEN(IUNIT1,FILNM,'r',KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '****************************************'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBOPEN      *'
            WRITE (IU06,*) '*  IN READWIND                         *'
            IF(KRET.EQ.-1)
     &       WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILNM
            IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILNM
            IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '****************************************'
            CALL ABORT1
          ENDIF

!         TRYING TO FIGURE OUT THE SIZE OF EACH FIELD

          FRSTIME = .FALSE.
        ENDIF  


! --------------------------------------------------------------------  

!*      2.0 READ MARS U and V WIND COMPONENTS  IN GRIB CODE FORMAT.
!           -------------------------------------------------------

        ISIZE=NBIT

2000  CONTINUE   

2001    KPLENG=ISIZE*KIND(ISIZE)
        IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))

        CALL PBGRIB(IUNIT1,KGRIB,KPLENG,KLEN,KRET)

        IF (KRET.EQ.-1) THEN
          WRITE(IU06,*) ' *******************************************' 
          WRITE(IU06,*) ' *                                         *' 
          WRITE(IU06,*) ' *       FATAL ERROR IN SUB. READWIND      *' 
          WRITE(IU06,*) ' *       ===========================       *' 
          WRITE(IU06,*) ' * PBGRIB EOF ON FILE ',FILNM                  
          WRITE(IU06,*) ' *       NO MORE WINDS                     *' 
          WRITE(IU06,*) ' *                                         *' 
          WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS      *' 
          WRITE(IU06,*) ' *                                         *' 
          WRITE(IU06,*) ' *******************************************' 
          CALL ABORT1
        ELSE IF (KRET.EQ.-2) THEN 
          WRITE(IU06,*) ' *******************************************' 
          WRITE(IU06,*) ' *                                         *' 
          WRITE(IU06,*) ' *       FATAL ERROR IN SUB. READWIND      *' 
          WRITE(IU06,*) ' *       ===========================       *' 
          WRITE(IU06,*) ' * PBGRIB ERROR ON FILE ',FILNM              
          WRITE(IU06,*) ' *                                         *' 
          WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS      *' 
          WRITE(IU06,*) ' *                                         *' 
          WRITE(IU06,*) ' *******************************************' 
          CALL ABORT1
        ELSE IF (KRET.EQ.-3) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' ***** WARNING IN MODEL FOR WAVES *****'
          WRITE (IU06,*) ' SUB. READWIND : '
          WRITE (IU06,*) ' SIZE  OF KGRIB IS NOT BIG ENOUGH.'
          WRITE (IU06,*) ' IT WAS ', ISIZE

          DEALLOCATE(KGRIB)
          ISIZE=(KLEN+KIND(ISIZE)-1)/KIND(ISIZE)

          WRITE (IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
          WRITE (IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
          WRITE (IU06,*) ' AND THE FIELD READ WITH THE NEW SIZE'
          WRITE (IU06,*) ' IF THIS PROBLEM OCCURS TOO OFTEN'
          WRITE (IU06,*) ' MODIFY THE VALUE OF NBIT IN SOURCE'
          WRITE (IU06,*) ' ***** WARNING ****** WARNING *****'
          WRITE (IU06,*) ' '
          CALL FLUSH(IU06)

!         RESET THE FILE POINTER TO READ FIELD AGAIN
          KOFFSET=-KLEN
          CALL PBSEEK(IUNIT1,KOFFSET,1,KRET)
          IF (KRET.EQ.-1) THEN
            WRITE (IU06,*) ' SUB. READWIND : '
            WRITE (IU06,*) ' PBSEEK REACHED EOF ON UNIT ',IUNIT1
            CALL ABORT1
          ELSEIF (KRET .EQ. -2) THEN
            WRITE (IU06,*) ' SUB. READWIND  : '
            WRITE (IU06,*) ' FILE HANDLING ERROR ON UNIT ',IUNIT1
            CALL ABORT1
          ENDIF

          GOTO 2001
        ENDIF

        IF (ITEST.GT.1) THEN
          WRITE(IU06,*) ' SUB. READWIND - READ FROM ',FILNM            
        ENDIF              

! ----------------------------------------------------------------------

!*      3.0 UNPACK MARS FIELDS.
!           -------------------

        IF(.NOT.ALLOCATED(PSEC4)) THEN

!         GET SIZE OF THE UNPACKED DATA
          ALLOCATE(PSEC4(1))
          KLENG = KLEN
          KLENP = 1
          KRET  = 1
!         CALL GRSDBG(1)
          CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &                 PSEC4, KLENP, KGRIB, KLENG , KWORD, 'I', KRET)
          IF(KRET.NE.0) THEN 
            WRITE(IU06,*) ' **************************************' 
            WRITE(IU06,*) ' *                                    *' 
            WRITE(IU06,*) ' * FATAL ERROR IN SUB. READWIND       *' 
            WRITE(IU06,*) ' * ===========================        *' 
            WRITE(IU06,*) ' * WHILE READING HEADERS FROM FILE    *' 
            WRITE(IU06,*) ' * GRIBEX ERROR - KRET = ',KRET
            WRITE(IU06,*) ' *                                    *' 
            WRITE(IU06,*) ' *  PROGRAM ABORTS   PROGRAM ABORTS   *' 
            WRITE(IU06,*) ' *                                    *' 
            WRITE(IU06,*) ' **************************************' 
            CALL ABORT1
          ENDIF

          IF(KSEC2(17).EQ.1) THEN
!           the size of the decoded data array will be bigger than
!           the size given by KSEC2(2) and KSEC2(3) BECAUSE we are
!           decoding the grib data with option HOPER=R which means
!           that the data will be produced on a regular lat lon grid.
!           Nevertheless, the size of the data array can be estimated 
!           as follows
            NC=KSEC2(3)*2
            NR=KSEC2(3)+1
          ELSE
            NC=KSEC2(2)
            NR=KSEC2(3)
          ENDIF

          DEALLOCATE(PSEC4)

          ALLOCATE(PSEC4(NC*NR))

          IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
          ALLOCATE(UWND(NC,NR))                                
          IF(ALLOCATED(VWND)) DEALLOCATE(VWND)
          ALLOCATE(VWND(NC,NR))                                

        ENDIF

!       3.0.2 GET THE DATA


        KLENG = KLEN 
        KLENP = NC*NR
        KRET  = 1
!       CALL GRSDBG(1)
        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &               PSEC4, KLENP, KGRIB, KLENG , KWORD, HOPER, KRET)
        IF(KRET.NE.0) THEN 
          WRITE(IU06,*) ' **************************************' 
          WRITE(IU06,*) ' *                                    *' 
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. READWIND       *' 
          WRITE(IU06,*) ' * ===========================        *' 
          WRITE(IU06,*) ' * WHILE READING FROM FILE            *' 
          WRITE(IU06,*) ' * GRIBEX ERROR - KRET = ',KRET
          WRITE(IU06,*) ' *                                    *' 
          WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS *' 
          WRITE(IU06,*) ' *                                    *' 
          WRITE(IU06,*) ' **************************************' 
          CALL ABORT1
        ENDIF


!*      3.1 PRODUCT DEFINITION BLOCK.
!           -------------------------

        IF (ITEST.GT.2)  CALL GRPRS1 (KSEC0, KSEC1)

!*      3.1.1 DETERMINE DATE OF WIND FIELD FROM KSEC1(10)-KSEC1(14).

!       START DATE

        IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
        WRITE(CDTWIR,'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     &                               KSEC1(13),KSEC1(14)

!       DETERMINE TIME UNIT

        IF (KSEC1(15).EQ.0) THEN
          JCONS = 60
        ELSEIF (KSEC1(15).EQ.1) THEN
          JCONS = 3600 
        ELSEIF (KSEC1(15).EQ.2) THEN 
          JCONS = 86400 
        ELSE 
          WRITE(IU06,*) ' ********************************************' 
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' *        FATAL ERROR IN SUB. READWIND      *' 
          WRITE(IU06,*) ' *        ============================      *' 
          WRITE(IU06,*) ' * SUSPICIOUS TIME UNIT IN MARS BLOCK 1     *' 
          WRITE(IU06,*) ' * TIME UNIT IS KSEC1(15) = ', KSEC1(15)       
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS     *' 
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' ********************************************' 
          CALL ABORT1 
        ENDIF 

!       DETERMINE WIND FIELD DATE.

        CALL INCDATE(CDTWIR,JCONS*KSEC1(16))

!*      3.1.2 DETERMINE CODE FOR WIND FIELD TYPE:
!             165 = U10 FIELD , 180 = STRESS FIELD.
!             ------------------------------------ 

        IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.166) THEN
          ICODE = 3
        ELSEIF (KSEC1(6).EQ.180 .OR. KSEC1(6).EQ.181) THEN
          ICODE = 2
        ELSE 
          WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++' 
          WRITE(IU06,*) ' +                                          +' 
          WRITE(IU06,*) ' +      WARNING ERROR IN SUB. READWIND      +' 
          WRITE(IU06,*) ' +      ==============================      +' 
          WRITE(IU06,*) ' + SUSPICIOUS WIND FIELD TYPE IN MARS BLOCK +' 
          WRITE(IU06,*) ' + TYPE IS KSEC1(6) = ', KSEC1(6)              
          WRITE(IU06,*) ' +                                          +' 
          WRITE(IU06,*) ' +    PROGRAM CONTINUES                     +' 
          WRITE(IU06,*) ' +                                          +' 
          WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++' 
          GOTO 2000
        ENDIF

!*      3.2 GRID DEFINITION BLOCK.
!           ---------------------

        IF (ITEST.GT.2) CALL GRPRS2 (KSEC0, KSEC2, PSEC2)

!*      3.2.1 DETERMINE GRID PARAMETERS.                
!             -------------------------- 

        KCOL = KSEC2(2)
        KROW = KSEC2(3) 
        IF ((KCOL.GT.NC).OR.(KROW.GT.NR)) THEN
          WRITE(IU06,*) ' ********************************************' 
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' *        FATAL ERROR IN SUB. READWIND      *' 
          WRITE(IU06,*) ' *        ============================      *' 
          WRITE(IU06,*) ' * DIMENSION DO NOT FIT                     *' 
          WRITE(IU06,*) ' * REQUIRED IS     KCOL = ', KCOL              
          WRITE(IU06,*) ' * DIMENSION IS      NC = ', NC                
          WRITE(IU06,*) ' * REQUIRED IS     KROW = ', KROW              
          WRITE(IU06,*) ' * DIMENSION IS      NR = ', NR                
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS     *' 
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' ********************************************' 
          CALL ABORT1
        ENDIF                               
        RLATN = REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &           0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)
        RLONL = REAL(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &           0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)
        RLATS = REAL(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &           0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)
        RLONR = REAL(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &           0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)
        CALL ADJUST (RLONL, RLONR)
        DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
        DPHI  = (RLATN-RLATS)/REAL(KROW-1)
        IWPER = 0
        IF (MOD(RLONR-RLONL+DLAM+720., 360.).EQ.0.) IWPER = 1


! --------------------------------------------------------------------  

!*      4. REARRANGE WIND DATA.
!          --------------------

        NR1 = KROW+1
        IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.180) THEN
          L = 0  
          DO J=1,KROW                                               
            DO I=1,KCOL                                               
              L = L+1
              UWND(I,NR1-J)=PSEC4(L)
            ENDDO
          ENDDO
!C!!!     go get the other wind component
          GOTO 2000
        ELSEIF (KSEC1(6).EQ.166 .OR. KSEC1(6).EQ.181) THEN
          L = 0                                           
          DO J=1,KROW                                               
            DO I=1,KCOL                                               
              L = L+1
              VWND(I,NR1-J)=PSEC4(L)
            ENDDO
          ENDDO
        ENDIF                      

! --------------------------------------------------------------------  

!*      5. WRITE TEST OUTPUT.
!          ------------------

        IF (ITEST.GT.2) THEN
          WRITE(IU06,*) ' READWIND -  WIND FIELD FOR THE CDTWIR = ',
     &     CDTWIR
          WRITE(IU06,66) ((UWND(I,J),I=1,24),J=1,5)
          WRITE(IU06,66) ((VWND(I,J),I=1,24),J=1,5)
 66       FORMAT(1X,24F5.2)
        ENDIF

! --------------------------------------------------------------------  
                                                                        
        DEALLOCATE(KGRIB)

        IF (LCLOSEWND) THEN
          CALL PBCLOSE(IUNIT1,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '****************************************'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBCLOSE     *'
            WRITE (IU06,*) '*  IN READWIND                         *'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '****************************************'
            CALL ABORT1
          ENDIF
          FRSTIME = .TRUE.
          IF (ITEST.GT.1)
     &     WRITE(IU06,*) ' SUB. READWIND - PBCLOSE ON ',FILNM
        ENDIF  


      ELSE


!*      6.0 READ ATMOS. MODEL GRIB FIELDS FOR COUPLED MODEL
!           -----------------------------------------------

        IWPER = 1
        ICOORD = 1

        WIND: DO IVAR=1,2

          IF(IVAR.EQ.1) THEN

!           GET SIZE OF THE UNPACKED DATA
            ALLOCATE(PSEC4(1))
            KLENP = 1
            KRET  = 1
!           CALL GRSDBG(1)
            CALL GRIBEX (KSEC0, KSEC1, KSEC2,PSEC2,KSEC3,PSEC3,KSEC4,
     &                   PSEC4, KLENP, IUGB,  ILEN , KWORD,'I', KRET)
            IF(KRET.NE.0) THEN 
              WRITE(IU06,*) ' **************************************' 
              WRITE(IU06,*) ' *                                    *' 
              WRITE(IU06,*) ' * FATAL ERROR IN SUB. READWIND       *' 
              WRITE(IU06,*) ' * ===========================        *' 
              WRITE(IU06,*) ' * WHILE READING HEADERS FROM ARRAY   *' 
              WRITE(IU06,*) ' * GRIBEX ERROR - KRET = ',KRET
              WRITE(IU06,*) ' *                                    *' 
              WRITE(IU06,*) ' *  PROGRAM ABORTS   PROGRAM ABORTS   *' 
              WRITE(IU06,*) ' *                                    *' 
              WRITE(IU06,*) ' **************************************' 
              CALL ABORT1
            ENDIF

            IF(KSEC2(17).EQ.1) THEN
              NC=KSEC2(3)*2
              NR=KSEC2(3)+1
            ELSE
              NC=KSEC2(2)
              NR=KSEC2(3)
            ENDIF

            DEALLOCATE(PSEC4)

            ALLOCATE(PSEC4(NC*NR))

            IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
            ALLOCATE(UWND(NC,NR))                                
            IF(ALLOCATED(VWND)) DEALLOCATE(VWND)
            ALLOCATE(VWND(NC,NR))                                

          ENDIF

          KLENP = NC*NR
          KRET  = 1
!         CALL GRSDBG(1)
          IF(IVAR.EQ.1) THEN
            CALL GRIBEX(KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &                  PSEC4, KLENP, IUGB,  ILEN , KWORD, HOPER, KRET)
          ELSE
            CALL GRIBEX(KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &                   PSEC4, KLENP, IVGB,  ILEN , KWORD, HOPER, KRET)
          ENDIF
          IF(KRET.NE.0) THEN 
            WRITE(IU06,*) ' **************************************' 
            WRITE(IU06,*) ' *                                    *' 
            WRITE(IU06,*) ' * FATAL  ERROR IN SUB. READWIND      *' 
            WRITE(IU06,*) ' * ============================       *' 
            WRITE(IU06,*) ' * WHILE DECODING FROM ARRAY          *' 
            WRITE(IU06,*) ' * GRIBEX ERROR -  KRET = ',KRET
            WRITE(IU06,*) ' *                                    *' 
            WRITE(IU06,*) ' *  PROGRAM ABORTS   PROGRAM ABORTS   *' 
            WRITE(IU06,*) ' *                                    *' 
            WRITE(IU06,*) ' **************************************' 
            CALL ABORT1
          ENDIF

!*      6.1 PRODUCT DEFINITION BLOCK.
!           -------------------------

          IF (ITEST.GT.2)  CALL GRPRS1 (KSEC0, KSEC1)

!*        6.1.1 DETERMINE DATE OF WIND FIELD FROM KSEC1(10)-KSEC1(14).

!         START DATE

          IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
          WRITE(CDTWIR,'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     &                                 KSEC1(13),KSEC1(14)

!         DETERMINE TIME UNIT

          IF (KSEC1(15).EQ.0) THEN
            JCONS = 60
          ELSEIF (KSEC1(15).EQ.1) THEN
            JCONS = 3600
          ELSEIF (KSEC1(15).EQ.2) THEN
            JCONS = 86400
          ELSE
            WRITE(IU06,*) ' *******************************************'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *        FATAL ERROR IN SUB. READWIND     *'
            WRITE(IU06,*) ' *        ============================     *'
            WRITE(IU06,*) ' * SUSPICIOUS TIME UNIT IN MARS BLOCK 1    *'
            WRITE(IU06,*) ' * TIME UNIT IS KSEC1(15) = ', KSEC1(15)
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS    *'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *******************************************'
            CALL ABORT1
          ENDIF

!         DETERMINE WIND FIELD DATE.

          CALL INCDATE(CDTWIR,JCONS*KSEC1(16))


!*        6.1.2 DETERMINE CODE FOR WIND FIELD TYPE:
!               165 = U10 FIELD , 180 = STRESS FIELD.
!               ------------------------------------

          IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.166) THEN
            ICODE = 3
          ELSEIF (KSEC1(6).EQ.180 .OR. KSEC1(6).EQ.181) THEN
            ICODE = 2
          ELSE
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +      WARNING ERROR IN SUB. READWIND     +'
            WRITE(IU06,*) ' +      ==============================     +'
            WRITE(IU06,*) ' + SUSPICIOUS WIND FIELD TYPE IN MARS BLOCK+'
            WRITE(IU06,*) ' + TYPE IS KSEC1(6) = ', KSEC1(6)
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +    PROGRAM CONTINUES                    +'
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF

!*      6.2 GRID DEFINITION BLOCK.
!           ---------------------

          IF (ITEST.GT.2) CALL GRPRS2 (KSEC0, KSEC2, PSEC2)

!*        6.2.1 DETERMINE GRID PARAMETERS.
!               --------------------------

          KCOL = KSEC2(2)
          KROW = KSEC2(3)
          IF ((KCOL.GT.NC).OR.(KROW.GT.NR)) THEN
            WRITE(IU06,*) ' ******************************************'
            WRITE(IU06,*) ' *                                        *'
            WRITE(IU06,*) ' *        FATAL ERROR IN SUB. READWIND    *'
            WRITE(IU06,*) ' *        ============================    *'
            WRITE(IU06,*) ' * DIMENSION DO NOT FIT                   *'
            WRITE(IU06,*) ' * REQUIRED IS     KCOL = ', KCOL
            WRITE(IU06,*) ' * DIMENSION IS      NC = ', NC
            WRITE(IU06,*) ' * REQUIRED IS     KROW = ', KROW
            WRITE(IU06,*) ' * DIMENSION IS      NR = ', NR
            WRITE(IU06,*) ' *                                        *'
            WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS   *'
            WRITE(IU06,*) ' *                                        *'
            WRITE(IU06,*) ' ******************************************'
            CALL ABORT1
          ENDIF
          RLATN = REAL(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &             0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)
          RLONL = REAL(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &             0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)
          RLATS = REAL(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &             0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)
          RLONR = REAL(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &             0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)
          CALL ADJUST (RLONL, RLONR)

          IF(KCOL.NE.1) THEN
            DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
          ELSE
            DLAM  = 0.
          ENDIF
          IF(KROW.NE.1) THEN
            DPHI  = (RLATN-RLATS)/REAL(KROW-1)
          ELSE
            DPHI  = 0.
          ENDIF

! --------------------------------------------------------------------

!*      7. REARRANGE WIND DATA.
!          --------------------

          NR1 = KROW+1
          IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.180) THEN
            L = 0
            DO J=1,KROW
              DO I=1,KCOL
                L = L+1
                UWND(I,NR1-J)=PSEC4(L)
              ENDDO
            ENDDO
          ELSEIF (KSEC1(6).EQ.166 .OR. KSEC1(6).EQ.181) THEN
            L = 0
            DO J=1,KROW
              DO I=1,KCOL
                L = L+1
                VWND(I,NR1-J)=PSEC4(L)
              ENDDO
            ENDDO
          ENDIF

        ENDDO WIND


      ENDIF

      IF(ALLOCATED(KSEC0)) DEALLOCATE(KSEC0)
      IF(ALLOCATED(KSEC1)) DEALLOCATE(KSEC1)
      IF(ALLOCATED(KSEC2)) DEALLOCATE(KSEC2)
      IF(ALLOCATED(KSEC3)) DEALLOCATE(KSEC3)
      IF(ALLOCATED(KSEC4)) DEALLOCATE(KSEC4)

      IF(ALLOCATED(PSEC2)) DEALLOCATE(PSEC2)
      IF(ALLOCATED(PSEC3)) DEALLOCATE(PSEC3)
      IF(ALLOCATED(PSEC4)) DEALLOCATE(PSEC4)

      RETURN


!     SWAMP CASE : CREATE YOUR OWN WIND FIELD

8000  CONTINUE
      WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',
     &                SWAMPWIND 

!     THE DIMENSIONS ARE HARD CODED !!
      NC=NGX
      NR=NGY

      IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
      ALLOCATE(UWND(NC,NR))                                
      IF(ALLOCATED(VWND)) DEALLOCATE(VWND)
      ALLOCATE(VWND(NC,NR))                                

      IF (FRSTIME) THEN
         RLATS = AMOSOP
         RLATN = AMONOP 
         RLONL = AMOWEP 
         RLONR = AMOEAP 
         ICODE = 3
         CDTWIR = CDATEA
         CALL INCDATE(CDTWIR,-IDELWI)
         KCOL = NC
         KROW = NR
         IWPER = 0
         ICOORD = 1
         IF(KCOL.NE.1) THEN
           DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
         ELSE
           DLAM  = 0.
         ENDIF
         IF(KROW.NE.1) THEN
           DPHI  = (RLATN-RLATS)/REAL(KROW-1)
         ELSE
           DPHI  = 0.
         ENDIF
         FRSTIME = .FALSE.
      ENDIF

      CALL INCDATE(CDTWIR,IDELWI)

      DO J=1,NR
        DO I=1,NC
          UWND(I,J)=0.0
        ENDDO
      ENDDO
      DO J=1,NR
        DO I=1,NC
         VWND(I,J)=SWAMPWIND
        ENDDO
      ENDDO

!*      WRITE TEST OUTPUT.
!       ------------------

        IF (ITEST.GT.2) THEN
          WRITE(IU06,*) ' READWIND -  WIND FIELD FOR THE CDTWIR = ',
     &     CDTWIR
          WRITE(IU06,67) ((UWND(I,J),I=1,24),J=1,5)
          WRITE(IU06,67) ((VWND(I,J),I=1,24),J=1,5)
 67       FORMAT(1X,24F5.2)
        ENDIF

      RETURN


      END SUBROUTINE READWIND
