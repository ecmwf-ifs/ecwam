      SUBROUTINE MAKEBLO (BLOCK, GRID, IJS, IJL, IG)                    
                                                                        
C**** *MAKEBLO* -  TRANSFER VARIABLES FROM GRID TO BLOCK.               
C                                                                       
C     PIERO LIONELLO      ECMWF     JUNE 1990                           
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       TRANSFER VARIABLES FROM GRID TO BLOCK.                          
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *MAKEBLO (BLOCK, GRID, IJS, IJL, IG)*                    
C                                                                       
C         *BLOCK*   REAL      BLOCKED DATA.                             
C         *GRID*    REAL      GRIDDED DATA.                             
C         *IJS*     INTEGER   INDEX OF FIRST POINT IN BLOCK.            
C         *IJL*     INTEGER   INDEX OF LAST POINT IN BLOCK.             
C         *IG*      INTEGER   BLOCK NUMBER.                             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
#include "param.h"
C                                                                       
#include "commap.h"
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      DIMENSION BLOCK(IJS:IJL), GRID(NGX,NGY)                             
C                                                                       
C ----------------------------------------------------------------------
C                                                                       
      DO 1000 IJ = IJS,IJL                                              
         IX = IXLG(IJ,IG)                                               
         JY = NGY-KXLT(IJ,IG) + 1                                       
         BLOCK(IJ) = GRID(IX,JY)                                        
 1000 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
