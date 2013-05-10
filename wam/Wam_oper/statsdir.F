      SUBROUTINE STATSDIR(NTOT, X, Y,
     &                    XMEAN, YMEAN, STDEVX, STDEVY,      
     &                    DMAX, BIAS, RMSE, STDEVD, SLOPES)
                                                                        
! ----------------------------------------------------------------------

!**** *STATSDIR* - EVALUATE STATISTICS FOR DIRECTION EXPRESSED IN DEGREE



!*    PURPOSE.                                                          
!     --------                                                          

!       THIS SUB COMPUTES THE FOLLOWING STATISTICS                      

!       1.  MEAN AND STANDARD DEVIATION OF BOTH X AND Y VALUES          
!       2.  MAXIMUM VALUE OF ABSOLUTE DIFFERNCE                         
!       3.  MEAN VALUE OF DIFFERENCE (BIAS)                           
!       4.  ROOT-MEAN-SQUARE ERROR (RMSE)                               
!       5.  STANDARD DEVIATION OF DIFFERENCES
!       6.  SYMMETRIC SLOPE.

!***  INTERFACE.                                                        
!     ----------                                                        

!        *CALL STATSDIR(NTOT, X, Y,
!                       XMEAN, YMEAN, STDEVX, STDEVY,      
!                       DMAX, BIAS, RMSE, STDEVD, SLOPES)

!*    VARIABLE.    TYPE.     PURPOSE.                                   
!     ---------    -------   --------                                   
!       *NTOT*     INTEGER   NUMBER OF POINTS.                          
!       *X*        REAL      INDEPENDENT VARIABLE
!       *Y*        REAL      DEPENDENT VARIABLE
!       *XMEAN*    REAL      MEAN VALUE OF INDEPENDENT VARIABLE         
!                            GIVEN FOR DIRECTION IN [0,360.[
!       *YMEAN*    REAL      MEAN VALUE OF DEPENDENT VARIABLE           
!                            GIVEN FOR DIRECTION IN [0,360.[
!       *STDEVX*   REAL      STANDARD DEVIATION OF INDEPENDENT VARIABLE 
!       *STDEVY*   REAL      STANDARD DEVIATION OF DEPENDENT VARIABLE   
!       *DMAX*     REAL      MAX ABSOLUTE VALUE OF DIFFERENCE           
!       *BIAS*     REAL      MEAN VALUE OF DIFFERENCES                  
!       *RMSE*     REAL      ROOT-MEAN-SQUARE ERROR                     
!       *STDEVD*   REAL      STANDARD DEVIATION OF DIFFERENCE           
!       *SLOPES*   REAL      SYMMETRIC SLOPE 

!     METHOD.                                                           
!     -------                                                           

!       THE COMPUTATION OF THESE STATISTICS IS STRAIGHTFORWARD.         


!     EXTERNALS.                                                        
!     ----------                                                        

!       NONE.                                                           

!     REFERENCES.                                                       
!     -----------                                                       

! ----------------------------------------------------------------------

      REAL :: X(NTOT), Y(NTOT)

!*    1.0 INITIALISE VALUES.                                            
!         ------------------                                            

      XMEAN  = 0.
      YMEAN  = 0.
      STDEVX = 0.
      STDEVY = 0.
      BIAS   = 0.
      RMSE   = 0.
      STDEVD = 0.
      SLOPES = 0.
      DMAX = 0.

      PI=4*ATAN(1.)
      RAD=PI/180.
      DEG=1./RAD

! ----------------------------------------------------------------------

!*    2. ENOUGH DATA?                                                   
!        ------------                                                   

      IF (NTOT.LT.1) RETURN                                             

      XPT  = REAL(NTOT)                                                 
      XPTD = 1./XPT                                                     

! ----------------------------------------------------------------------

!*    3. COMPUTE MEANS
!        -------------

      SUMXC=0.
      SUMXS=0.
      SUMYC=0.
      SUMYS=0.

      DO IC=1,NTOT
        XRAD=RAD*X(IC)
        SUMXC=SUMXC+COS(XRAD)
        SUMXS=SUMXS+SIN(XRAD)
        YRAD=RAD*Y(IC)
        SUMYC=SUMYC+COS(YRAD)
        SUMYS=SUMYS+SIN(YRAD)
      ENDDO
      XMEAN=DEG*ATAN2(SUMXS,SUMXC)
      YMEAN=DEG*ATAN2(SUMYS,SUMYC)

      BIAS = YMEAN-XMEAN                                             
      IF(BIAS.GT.180.) THEN
         BIAS=BIAS-360.
      ELSE IF(BIAS.LT.-180.) THEN
         BIAS=BIAS+360.
      ENDIF


!*    3. COMPUTE RMSE STDEVX STDEVY STDEVD SLOPES
!        ----------------------------------------
                                                                        
      IF (NTOT.GT.1) THEN
        SUMXX2=0.
        SUMYY2=0.
        SUMYX2=0.
        SUMYYXX2=0.
        DO IC=1,NTOT
          XDIFF=X(IC)-XMEAN
          XDIFF=XDIFF-INT(XDIFF/180.)*180.
          SUMXX2=SUMXX2+XDIFF**2
          YDIFF=Y(IC)-YMEAN
          YDIFF=YDIFF-INT(YDIFF/180.)*180.
          SUMYY2=SUMYY2+YDIFF**2
          YXDIFF=Y(IC)-X(IC)
          YXDIFF=YXDIFF-INT(YXDIFF/180.)*180.
          SUMYX2=SUMYX2+YXDIFF**2
          DMAX=MAX(DMAX,ABS(YXDIFF))
          YYXXDIFF=YDIFF-XDIFF
          YYXXDIFF=YYXXDIFF-INT(YYXXDIFF/180.)*180.
          SUMYYXX2=SUMYYXX2+YYXXDIFF**2
        ENDDO
        STDEVX = SQRT(SUMXX2/(XPT-1.))
        STDEVY = SQRT(SUMYY2/(XPT-1.))
        RMSE = SQRT(SUMYX2*XPTD)
        STDEVD = SQRT(SUMYYXX2/(XPT-1.))
        IF(SUMXX2.GT.0.) THEN
          SLOPES = SQRT(SUMYY2/SUMXX2)
        ENDIF

      ENDIF                                                             

!     INSURE MEANS ARE IN [0,360.]
      IF(XMEAN.LT.0.) XMEAN=XMEAN+360.
      IF(YMEAN.LT.0.) YMEAN=YMEAN+360.

      RETURN                                                            
      END SUBROUTINE STATSDIR
