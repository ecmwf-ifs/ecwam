      FUNCTION INCINT(KINT,KINC)

!**** *INCINT* - Increment integer under multi-tasking lock.

!     Purpose.
!     --------
!           Increment KINT with KINC under multi-tasking lock. If
!       multi-tasking not requested just increment KINT.         
!**   Interface.                                                 
!     ----------                                                 
!        *ILAT =INCINT(KINT,KINC)*                               

!        Explicit arguments :  KINT - integer to increment       
!        --------------------  KINC - increment                  


!        Implicit arguments :   NVEROU in YOMMUL                 
!        --------------------                                    

!     Method.                                                    
!     -------                                                    
!        See documentation                                       

!     Externals.                                                 
!     ----------                                                 
!                   ISELFADD- Cray multi-tasking routine         

!     Reference.                                                 
!     ----------                                                 
!        ECMWF Research Department documentation of the IFS      

!     Author.                                                    
!     -------                                                    
!        Mats Hamrud and Philippe Courtier  *ECMWF*              

!     Modifications.                                             
!     --------------                                             
!        Original : 89-04-01                                     
!     ------------------------------------------------------------------
      IMPLICIT LOGICAL (L)                                              


!     ------------------------------------------------------------------

!*       1.    INCREMENT INTEGER.                                       
!              ------------------                                       


      INCINT=ISELFADD(KINT,KINC)+KINC                                   

!     ------------------------------------------------------------------

      RETURN                                                            
      END FUNCTION INCINT




