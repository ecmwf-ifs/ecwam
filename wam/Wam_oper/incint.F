      FUNCTION INCINT(KINT,KINC)
C
C**** *INCINT* - Increment integer under multi-tasking lock.
C
C     Purpose.
C     --------
C           Increment KINT with KINC under multi-tasking lock. If
C       multi-tasking not requested just increment KINT.         
C**   Interface.                                                 
C     ----------                                                 
C        *ILAT =INCINT(KINT,KINC)*                               
C                                                                
C        Explicit arguments :  KINT - integer to increment       
C        --------------------  KINC - increment                  
C                                                                
C                                                                
C        Implicit arguments :   NVEROU in YOMMUL                 
C        --------------------                                    
C                                                                
C     Method.                                                    
C     -------                                                    
C        See documentation                                       
C                                                                
C     Externals.                                                 
C     ----------                                                 
C                   ISELFADD- Cray multi-tasking routine         
C                                                                
C     Reference.                                                 
C     ----------                                                 
C        ECMWF Research Department documentation of the IFS      
C                                                                
C     Author.                                                    
C     -------                                                    
C        Mats Hamrud and Philippe Courtier  *ECMWF*              
C                                                                
C     Modifications.                                             
C     --------------                                             
C        Original : 89-04-01                                     
C     ------------------------------------------------------------------
      IMPLICIT LOGICAL (L)                                              
C                                                                       
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C*       1.    INCREMENT INTEGER.                                       
C              ------------------                                       
C                                                                       
 100  CONTINUE                                                          
C                                                                       
      INCINT=ISELFADD(KINT,KINC)+KINC                                   
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      RETURN                                                            
      END                                                               
