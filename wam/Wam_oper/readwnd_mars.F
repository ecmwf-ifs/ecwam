      SUBROUTINE READWND_MARS (CDTWIR, UWND, VWND, NC, NR,
     1                         IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)
C                                                                       
C***  *READWND* - PROGRAM TO REFORMAT WIND FIELDS FROM MARS             
C                                                                       
C    PETER JANSSEN       ECMWF        FEBRUARY 1987                     
C                                                                       
C    MODIFIED BY:                                                       
C    ------------                                                       
C                                                                       
C    LIANA ZAMBRESKY     GKSS/ECMWF   JULY 1988                         
C                                                                       
C        1. PASSING OF ARRAYS AS FORMAL PARMETERS TO SAVE MEMORY        
C        2. PROCESSING WINDS FOR ONLY ONE WIND TIME STEP                
C                                                                       
C    PEDRO VITERBO       ECMWF        OCTOBER 1988                      
C                                                                       
C        DEALING WITH MARS PACKED DATA.                                 
C                                                                       
C    H. GUNTHER      ECMWF/GKSS       OCTOBER 1991                      
C                                                                       
C        GRIB EDITION 1.                                                
C                                                                       
C    B. HANSEN    ECMWF 1997
C                 RENAME SUBROUTINE TO MATCH THE NAME OF THE FILE.
C                 USE INCLUDE FILES FOR ALL COMMON BLOCKS.
C
C     PURPOSE                                                           
C     -------                                                           
C                                                                       
C     *READWND_MARS*
C         READS TARGET FILES FROM THE MARS ARCHIVE AND INFERS FROM THIS 
C             A) THE DATE OF THE WIND FIELD                             
C             B) THE DEFINITION OF THE MARS GRID                        
C             C) THE U AND V COMPONENT OF THE WIND FIELD                
C             D) THE TYPE OF WIND FIELD AND                             
C             E) THE UNPACKED WINDS                                     
C                                                                       
C     INTERFACE                                                         
C     ---------                                                         
C                                                                       
C     *CALL* *READWND_MARS (CDTWIR, UWND, VWND, NC, NR)*
C    1                      IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)
C        *CDTWIR* - DATE/TIME OF THE DATA READ                          
C        *UWND*   - HORIZONTAL WIND COMPONENT                           
C        *VWND*   - VERTICAL WIND COMPONENT                             
C        *NC*     - NUMBER OF COLUMNS IN INPUT WIND ARRAY (DIMENSION)   
C        *NR*     - NUMBER OF ROWS    IN INPUT WIND ARRAY (DIMENSION)   
C        *IUGB*   - ARRAY CONTAINING THE GRIB U COMPONENT
C        *IVGB*   - ARRAY CONTAINING THE GRIB V COMPONENT
C        *ILEN*   - SIZE OF IUGB AND IVGB
C        *LWNDFILE*  IF TRUE DATA ARE ACCESSED FROM THE INPUT FILE sfcwindin
C                    ELSE THE DATA WILL BE SUPPLIED IN IUGB AND IVGB
C        *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND THE UNIT RESET
C                                                                       
C     UNIT             DESCRIPTION                                      
C     ----             -----------                                      
CT63                                                                    
C     TAPE63           WINDFIELDS FROM MARS ARCHIVE AS WRITTEN BY       
C                      PROG GETINPT.                                    
CT63  TAPE39           WINDFIELD FROM CT63 MODEL. (THESE CHANGES ARE    
CT63                   COMMENTED OUT BY CT63.)                          
C                                                                       
C     EXTERNALS                                                         
C     ---------                                                         
C                                                                       
C     *ABORT1*         TERMINATES PROCESSING
C     *GRIBEX*         UNPACKS MARS WINDS                               
C     *GRPRS1*         PRINTS PRODUCT DEFINITION BLOCK                  
C     *GRPRS2*         PRINTS GRID  DEFINITION BLOCK                    
C     *INCDATE*        INCREMENTS DATE                                  
C     *PBCLOSE*
C     *PBGRIB*
C     *PBOPEN*
C                                                                       
C --------------------------------------------------------------------- 
C                                                                       
#include "comtest.h"
C
#include "comwndg.h"
C                                                                       

C --------------------------------------------------------------------- 
C                                                                       
#if   region=='g' && resolution==50
      PARAMETER (NBIT = 130000)
#elif region=='m' && resolution==25
      PARAMETER (NBIT = 130000)
#else
      PARAMETER (NBIT = 60000)
#endif
CT106 PARAMETER (NBIT = 13000)
CT63  PARAMETER (NBIT = 14000)
      DIMENSION UWND(NC,NR), VWND(NC,NR)                                
      INTEGER   KSEC1(60), KSEC2(22+NR), KSEC3(2), KSEC4(42), KSEC0(2)
      INTEGER, ALLOCATABLE :: KGRIB(:)
C ---------------------------------------------------------------------
      INTEGER IUGB(ILEN), IVGB(ILEN)
C ---------------------------------------------------------------------
      DIMENSION PSEC2(10), PSEC3(2), PSEC4(NC*NR)                       
      CHARACTER HOPER*1, FILNM*24
      CHARACTER*12 CDTWIR
      LOGICAL FRSTIME, LWNDFILE, LCLOSEWND
      DATA IFIL, HOPER, FRSTIME / 63, 'R', .TRUE. /
      SAVE IUNIT1      
C                                                                       
C --------------------------------------------------------------------  
C                                                                       

      IF(LWNDFILE) THEN


C       1.0 INITIALISE VALUES. 
C           ------------------ 
C                                                                       
        IPER = 1 
        ICOORD = 1
C                                                      
        IF (FRSTIME) THEN
           FILNM='sfcwindin'
           CALL PBOPEN(IUNIT1,FILNM,'r',KRET)
           IF(KRET.LT.0) THEN
             WRITE (IU06,*) '****************************************'
             WRITE (IU06,*) '*                                      *'
             WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBOPEN      *'
             WRITE (IU06,*) '*  IN READWND_MARS                     *'
             IF(KRET.EQ.-1)
     1       WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILNM
             IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILNM
             IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
             WRITE (IU06,*) '*                                      *'
             WRITE (IU06,*) '****************************************'
             CALL ABORT1
           ENDIF
           FRSTIME = .FALSE.
        END IF  
C
C                                                                       
C --------------------------------------------------------------------  
C                                                                       
C*      2.0 READ MARS WINDS IN GRIB CODE FORMAT.
C           ------------------------------------
C                                                                       
        ISIZE=NBIT
C             
 2000 CONTINUE   
C                                                        
2001    KPLENG=ISIZE*KIND(ISIZE)
        IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))
C
        CALL PBGRIB(IUNIT1,KGRIB,KPLENG,KLEN,KRET)
C  
        IF (KRET.EQ.-1) GO TO 5000 
        IF (KRET.EQ.-2) GO TO 5100
        IF (KRET.EQ.-3) THEN
          WRITE (IU06,*) ' '
          WRITE (IU06,*) ' ***** WARNING ****** WARNING *****'
          WRITE (IU06,*) ' SUB. READWND_MARS : '
          WRITE (IU06,*) ' SIZE  OF KGRIB IS NOT BIG ENOUGH.'
          WRITE (IU06,*) ' IT WAS ', ISIZE
C
          DEALLOCATE(KGRIB)
          ISIZE=(KLEN+KIND(ISIZE)-1)/KIND(ISIZE)
C
          WRITE (IU06,*) ' IT SHOULD AT LEAST BE ', ISIZE
          WRITE (IU06,*) ' THE SIZE WAS RESET AUTOMATICALLY'
          WRITE (IU06,*) ' AND THE FIELD READ WITH THE NEW SIZE'
          WRITE (IU06,*) ' IF THIS PROBLEM OCCURS TOO OFTEN'
          WRITE (IU06,*) ' MODIFY THE VALUE OF NBIT IN SOURCE'
          WRITE (IU06,*) ' ***** WARNING ****** WARNING *****'
          WRITE (IU06,*) ' '
          CALL FLUSH(IU06)
C
C         RESET THE FILE POINTER TO READ FIELD AGAIN
          KOFFSET=-KLEN
          CALL PBSEEK(IUNIT1,KOFFSET,1,KRET)
          IF (KRET.EQ.-1) THEN
            WRITE (IU06,*) ' SUB. READWND_MARS : '
            WRITE (IU06,*) ' PBSEEK REACHED EOF ON UNIT ',IUNIT1
            CALL ABORT1
          ELSEIF (KRET .EQ. -2) THEN
            WRITE (IU06,*) ' SUB. READWND_MARS  : '
            WRITE (IU06,*) ' FILE HANDLING ERROR ON UNIT ',IUNIT1
            CALL ABORT1
          ENDIF
C
          GOTO 2001
        ENDIF

        IF (ITEST.GT.1) THEN
           WRITE(IU06,*) ' SUB. READWND - READ FROM ',FILNM            
        ENDIF              
 2100   CONTINUE
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
C*      3.0 UNPACK MARS FIELDS.
C           -------------------
C                                                                       
        KLENG = KLEN 
        KLENP = NC*NR
        KRET  = 1
C       CALL GRSDBG(1)
        CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     1               PSEC4, KLENP, KGRIB, KLENG , KWORD, HOPER, KRET)
        IF(KRET.NE.0) GOTO 5200

C                                                                       
C*      3.1 PRODUCT DEFINITION BLOCK.
C           -------------------------
C                                                                       
        IF (ITEST.GT.2)  CALL GRPRS1 (KSEC0, KSEC1)
C                                                                       
C*      3.1.1 DETERMINE DATE OF WIND FIELD FROM KSEC1(10)-KSEC1(14).
C                                                                       
C       START DATE
C                                                                       
        IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
        WRITE(CDTWIR,'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     1                               KSEC1(13),KSEC1(14)
C                                                                       
C       DETERMINE TIME UNIT
C                                                                       
        IF (KSEC1(15).EQ.0) THEN
          JCONS = 60
        ELSE IF (KSEC1(15).EQ.1) THEN
          JCONS = 3600 
        ELSE IF (KSEC1(15).EQ.2) THEN 
          JCONS = 86400 
        ELSE 
          WRITE(IU06,*) ' ********************************************' 
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND      *' 
          WRITE(IU06,*) ' *         ===========================      *' 
          WRITE(IU06,*) ' * SUSPICIOUS TIME UNIT IN MARS BLOCK 1     *' 
          WRITE(IU06,*) ' * TIME UNIT IS KSEC1(15) = ', KSEC1(15)       
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS     *' 
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' ********************************************' 
          CALL ABORT1 
        ENDIF 
C                                                                       
C       DETERMINE WIND FIELD DATE.
C                                                                       
        CALL INCDATE(CDTWIR,JCONS*KSEC1(16))
C                                                                       
C*      3.1.2 DETERMINE CODE FOR WIND FIELD TYPE:
C             165 = U10 FIELD , 180 = STRESS FIELD.
C             ------------------------------------ 
C                                                                       
        IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.166) THEN
           ICODE = 3
        ELSE IF (KSEC1(6).EQ.180 .OR. KSEC1(6).EQ.181) THEN
           ICODE = 2
        ELSE 
          WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++' 
          WRITE(IU06,*) ' +                                          +' 
          WRITE(IU06,*) ' +       WARNING ERROR IN SUB. READWND      +' 
          WRITE(IU06,*) ' +       =============================      +' 
          WRITE(IU06,*) ' + SUSPICIOUS WIND FIELD TYPE IN MARS BLOCK +' 
          WRITE(IU06,*) ' + TYPE IS KSEC1(6) = ', KSEC1(6)              
          WRITE(IU06,*) ' +                                          +' 
          WRITE(IU06,*) ' +    PROGRAM CONTINUES                     +' 
          WRITE(IU06,*) ' +                                          +' 
          WRITE(IU06,*) ' ++++++++++++++++++++++++++++++++++++++++++++' 
          GOTO 2000
        ENDIF
C                                                                       
C*      3.2 GRID DEFINITION BLOCK.
C           ---------------------
C                                                                       
        IF (ITEST.GT.2) CALL GRPRS2 (KSEC0, KSEC2, PSEC2)
C                                                                       
C*      3.2.1 DETERMINE GRID PARAMETERS.                
C             -------------------------- 
C                                                                       
        KCOL = KSEC2(2)
        KROW = KSEC2(3) 
        IF ((KCOL.GT.NC).OR.(KROW.GT.NR)) THEN
          WRITE(IU06,*) ' ********************************************' 
          WRITE(IU06,*) ' *                                          *' 
          WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND      *' 
          WRITE(IU06,*) ' *         ===========================      *' 
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
        END IF                               
        RLATN = FLOAT(KSEC2(4)/1000)+0.1*(MOD(KSEC2(4),1000)/100)+
     &           0.01*(MOD(KSEC2(4),100)/10)+0.001*MOD(KSEC2(4),10)
        RLONL = FLOAT(KSEC2(5)/1000)+0.1*(MOD(KSEC2(5),1000)/100)+
     &           0.01*(MOD(KSEC2(5),100)/10)+0.001*MOD(KSEC2(5),10)
        RLATS = FLOAT(KSEC2(7)/1000)+0.1*(MOD(KSEC2(7),1000)/100)+
     &           0.01*(MOD(KSEC2(7),100)/10)+0.001*MOD(KSEC2(7),10)
        RLONR = FLOAT(KSEC2(8)/1000)+0.1*(MOD(KSEC2(8),1000)/100)+
     &           0.01*(MOD(KSEC2(8),100)/10)+0.001*MOD(KSEC2(8),10)

#if region=='m'
        IF (RLONL.GT.RLONR)RLONL=RLONL-360.
#endif
        DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
        DPHI  = (RLATN-RLATS)/REAL(KROW-1)
C                                                                       
C --------------------------------------------------------------------  
C                                                                       
C*      4. REARRANGE WIND DATA.
C          --------------------
C                                                                       
        NR1 = KROW+1
        IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.180) THEN
          L = 0  
          DO J=1,KROW                                               
            DO I=1,KCOL                                               
              L = L+1
              UWND(I,NR1-J)=PSEC4(L)
            ENDDO
          ENDDO
CC!!!     go get the other wind component
          GOTO 2000
        ELSE IF (KSEC1(6).EQ.166 .OR. KSEC1(6).EQ.181) THEN
          L = 0                                           
          DO J=1,KROW                                               
            DO I=1,KCOL                                               
              L = L+1
              VWND(I,NR1-J)=PSEC4(L)
            ENDDO
          ENDDO
        ENDIF                      
C                                                                       
C --------------------------------------------------------------------  
C                                                                       
C*      5. WRITE TEST OUTPUT.
C          ------------------
C                                                                       
        IF (ITEST.GT.2) THEN
           WRITE(IU06,*) ' READWND -  WIND FIELD FOR THE CDTWIR = ',
     1                    CDTWIR
           WRITE(IU06,66) ((UWND(I,J),I=1,24),J=1,5)
           WRITE(IU06,66) ((VWND(I,J),I=1,24),J=1,5)
 66        FORMAT(1X,24F5.2)
        ENDIF
C                                                                       
C --------------------------------------------------------------------  
                                                                        
        DEALLOCATE(KGRIB)
C
        IF (LCLOSEWND) THEN
           CALL PBCLOSE(IUNIT1,KRET)
           IF(KRET.LT.0) THEN
             WRITE (IU06,*) '****************************************'
             WRITE (IU06,*) '*                                      *'
             WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBCLOSE     *'
             WRITE (IU06,*) '*  IN READWND_MARS                     *'
             WRITE (IU06,*) '*                                      *'
             WRITE (IU06,*) '****************************************'
             CALL ABORT1
           ENDIF
           FRSTIME = .TRUE.
           IF (ITEST.GT.1)
     &        WRITE(IU06,*) ' SUB. READWND - PBCLOSE ON ',FILNM            
        ENDIF  

        RETURN
C                                                                       
C --------------------------------------------------------------------  
C                                                                       
C*      6. ERROR EXITS
C          -----------
C                                                                       
 5000   CONTINUE
        WRITE(IU06,*) ' *********************************************' 
        WRITE(IU06,*) ' *                                           *' 
        WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND       *' 
        WRITE(IU06,*) ' *         ===========================       *' 
        WRITE(IU06,*) ' * PBGRIB EOF ON FILE ',FILNM                  
        WRITE(IU06,*) ' *         NO MORE WINDS                     *' 
        WRITE(IU06,*) ' *                                           *' 
        WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS      *' 
        WRITE(IU06,*) ' *                                           *' 
        WRITE(IU06,*) ' *********************************************' 
        CALL ABORT1                                                     
                                                                        
 5100   CONTINUE
        WRITE(IU06,*) ' *********************************************' 
        WRITE(IU06,*) ' *                                           *' 
        WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND       *' 
        WRITE(IU06,*) ' *         ===========================       *' 
        WRITE(IU06,*) ' * PBGRIB ERROR ON FILE ',FILNM              
        WRITE(IU06,*) ' *                                           *' 
        WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS      *' 
        WRITE(IU06,*) ' *                                           *' 
        WRITE(IU06,*) ' *********************************************' 
        CALL ABORT1                                                     
 5200   CONTINUE
        WRITE(IU06,*) ' *********************************************' 
        WRITE(IU06,*) ' *                                           *' 
        WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND       *' 
        WRITE(IU06,*) ' *         ===========================       *' 
        WRITE(IU06,*) ' * GRIBEX ERROR - KRET = ',KRET                 
        WRITE(IU06,*) ' *                                           *' 
        WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS      *' 
        WRITE(IU06,*) ' *                                           *' 
        WRITE(IU06,*) ' *********************************************' 
        CALL ABORT1                                                     
C
C
      ELSE
C
C*      3.0 READ ATMOS. MODEL GRIB FILES FOR COUPLED MODEL
C           ----------------------------------------------
C
        DO 4000 IVAR=1,2
          KLENP = NC*NR
          KRET  = 1
C         CALL GRSDBG(1)
          IF(IVAR.EQ.1) THEN
            CALL GRIBEX(KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     1                  PSEC4, KLENP, IUGB,  ILEN , KWORD, HOPER, KRET)
          ELSE
            CALL GRIBEX(KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     1                   PSEC4, KLENP, IVGB,  ILEN , KWORD, HOPER, KRET)
          ENDIF
          IF(KRET.NE.0) GOTO 6200
C
C*      3.1 PRODUCT DEFINITION BLOCK.
C           -------------------------
C
          IF (ITEST.GT.2)  CALL GRPRS1 (KSEC0, KSEC1)
C
C*        3.1.1 DETERMINE DATE OF WIND FIELD FROM KSEC1(10)-KSEC1(14).
C
C         START DATE
C
          IYYYY=(KSEC1(21)-1)*100+KSEC1(10)
          WRITE(CDTWIR,'(I4.4,4I2.2)') IYYYY,KSEC1(11),KSEC1(12),
     1                                 KSEC1(13),KSEC1(14)
C
C         DETERMINE TIME UNIT
C
          IF (KSEC1(15).EQ.0) THEN
            JCONS = 60
          ELSE IF (KSEC1(15).EQ.1) THEN
            JCONS = 3600
          ELSE IF (KSEC1(15).EQ.2) THEN
           JCONS = 86400
          ELSE
            WRITE(IU06,*) ' *******************************************'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND     *'
            WRITE(IU06,*) ' *         ===========================     *'
            WRITE(IU06,*) ' * SUSPICIOUS TIME UNIT IN MARS BLOCK 1    *'
            WRITE(IU06,*) ' * TIME UNIT IS KSEC1(15) = ', KSEC1(15)
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS    *'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *******************************************'
            CALL ABORT1
          ENDIF
C
C         DETERMINE WIND FIELD DATE.
C
          CALL INCDATE(CDTWIR,JCONS*KSEC1(16))
C
C
C*        3.1.2 DETERMINE CODE FOR WIND FIELD TYPE:
C               165 = U10 FIELD , 180 = STRESS FIELD.
C               ------------------------------------
C
          IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.166) THEN
            ICODE = 3
          ELSE IF (KSEC1(6).EQ.180 .OR. KSEC1(6).EQ.181) THEN
            ICODE = 2
          ELSE
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +       WARNING ERROR IN SUB. READWND     +'
            WRITE(IU06,*) ' +       =============================     +'
            WRITE(IU06,*) ' + SUSPICIOUS WIND FIELD TYPE IN MARS BLOCK+'
            WRITE(IU06,*) ' + TYPE IS KSEC1(6) = ', KSEC1(6)
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +    PROGRAM CONTINUES                    +'
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF
C
C*      3.2 GRID DEFINITION BLOCK.
C           ---------------------
C
          IF (ITEST.GT.2) CALL GRPRS2 (KSEC0, KSEC2, PSEC2)
C
C*        3.2.1 DETERMINE GRID PARAMETERS.
C               --------------------------
C
          KCOL = KSEC2(2)
          KROW = KSEC2(3)
          IF ((KCOL.GT.NC).OR.(KROW.GT.NR)) THEN
            WRITE(IU06,*) ' ******************************************'
            WRITE(IU06,*) ' *                                        *'
            WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND    *'
            WRITE(IU06,*) ' *         ===========================    *'
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
          END IF
          RLATN=.001*FLOAT(KSEC2(4))
          RLONL=.001*FLOAT(KSEC2(5))
          RLATS=.001*FLOAT(KSEC2(7))
          RLONR=.001*FLOAT(KSEC2(8))
C
          DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
          DPHI  = (RLATN-RLATS)/REAL(KROW-1)
C
C --------------------------------------------------------------------
C
C*      4. REARRANGE WIND DATA.
C          --------------------
C
          NR1 = KROW+1
          IF (KSEC1(6).EQ.165 .OR. KSEC1(6).EQ.180) THEN
            L = 0
            DO 4100 J=1,KROW
            DO 4100 I=1,KCOL
              L = L+1
              UWND(I,NR1-J)=PSEC4(L)
 4100       CONTINUE
          ELSE IF (KSEC1(6).EQ.166 .OR. KSEC1(6).EQ.181) THEN
           L = 0
           DO 4200 J=1,KROW
           DO 4200 I=1,KCOL
             L = L+1
             VWND(I,NR1-J)=PSEC4(L)
 4200      CONTINUE
          ENDIF
C
 4000   CONTINUE
C
        RETURN
C
 6200   CONTINUE
        WRITE(IU06,*) ' *********************************************'
        WRITE(IU06,*) ' *                                           *'
        WRITE(IU06,*) ' *         FATAL ERROR IN SUB. READWND       *'
        WRITE(IU06,*) ' *         ===========================       *'
        WRITE(IU06,*) ' * GRIBEX ERROR - KRET = ',KRET
        WRITE(IU06,*) ' *                                           *'
        WRITE(IU06,*) ' *    PROGRAM ABORTS     PROGRAM ABORTS      *'
        WRITE(IU06,*) ' *                                           *'
        WRITE(IU06,*) ' *********************************************'
        CALL ABORT1
C
      ENDIF
      RETURN
C
      END                                                               
