      SUBROUTINE MAKEBLO (BLOCK, GRID, IJS, IJL, IG)                    
                                                                        
!**** *MAKEBLO* -  TRANSFER VARIABLES FROM GRID TO BLOCK.               

!     PIERO LIONELLO      ECMWF     JUNE 1990                           

!     PURPOSE.                                                          
!     --------                                                          

!       TRANSFER VARIABLES FROM GRID TO BLOCK.                          

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *MAKEBLO (BLOCK, GRID, IJS, IJL, IG)*                    

!         *BLOCK*   REAL      BLOCKED DATA.                             
!         *GRID*    REAL      GRIDDED DATA.                             
!         *IJS*     INTEGER   INDEX OF FIRST POINT IN BLOCK.            
!         *IJL*     INTEGER   INDEX OF LAST POINT IN BLOCK.             
!         *IG*      INTEGER   BLOCK NUMBER.                             

! ----------------------------------------------------------------------


      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWMAP   , ONLY : IXLG     ,KXLT



! ----------------------------------------------------------------------

      DIMENSION BLOCK(IJS:IJL), GRID(NGX,NGY)

! ----------------------------------------------------------------------

      DO IJ = IJS,IJL                                              
        IX = IXLG(IJ,IG)                                               
        JY = NGY-KXLT(IJ,IG) + 1                                       
        BLOCK(IJ) = GRID(IX,JY)                                        
      ENDDO
                                                                        
      RETURN                                                            
      END SUBROUTINE MAKEBLO
