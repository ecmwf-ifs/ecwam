      SUBROUTINE ANALYSE (XM, NOBS, XME, XMO, XOI, IMEAS, JMEAS, I, J,
     &                    DIST, NDIM2) 

! ----------------------------------------------------------------------

!**** *ANALYSE* - COMPUTE FIELD AT (I,J) BY OPTIMUM INTERPOLATION.      

!     PIERO LIONELLO      ECMWF     JUNE 1990                           

!     PURPOSE.                                                          
!     --------                                                          

!       COMPUTE THE FIELD IN THE POINT (I,J) BY OPTIMUM INTERPOLATION.  
!       ANALYSIS IS DONE FOR ONE GRID POINT.                            

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *ANALYSE (XM, NOBS, XME, XMO, XOI, IMEAS, JMEAS, I, J,
!                        DIST, NDIM2)*  

!         *XM*      REAL      INVERSE OF (P+M) : P=MODEL CORRELATION    
!                             MATRIX , M=MEASUREMENTS CORRELATION MATRIX
!         *NOBS*    INTEGER   NUMBER OF OBSERVATION IN THE              
!                             BOX AROUND POINT (I,J)                    
!         *XME*     REAL      FIELD FROM  MEASUREMENT                   
!         *XMO*     REAL      FIELD FROM MODEL                          
!         *XOI*     REAL      RESULT OF O.I.                            
!         *IMEAS*   INTEGER   X-INDICES OF MEASUREMENTS                 
!         *JMEAS*   INTEGER   Y-INDICES OF MEASUREMENTS                 
!         *I*       INTEGER   POINT WHERE O.I.FIELD IS PRODUCED         
!         *J*       INTEGER   POINT WHERE O.I.FIELD IS PRODUCED
!         *DIST*    REAL      CUTOFF LENGTH OF THE SPREADING FUNCTION
!         *NDIM2*   INTEGER   DIMENSION OF IMEAS AND JMEAS

! ----------------------------------------------------------------------

      USE YOWGRID  , ONLY : SINPH    ,COSPH
      USE YOWMAP   , ONLY : ZDELLO
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWPCONS , ONLY : RAD

! ----------------------------------------------------------------------

      DIMENSION XME(NGX,NGY), XMO(NGX,NGY), IMEAS(NDIM2), JMEAS(NDIM2)  
      DIMENSION XM(NOBS,NOBS)                                           

! ----------------------------------------------------------------------

      DO L=1,NOBS-1                                                
        DO K=L+1,NOBS                                                
          XM(L,K) = XM(K,L)
        ENDDO
      ENDDO

      XNEW = 0.                                                         
      DO IOBS=1,NOBS 
        I1 = I
        I2 = IMEAS(IOBS)
        J1 = NGY - J + 1
        J2 = NGY - JMEAS(IOBS) + 1

        IF((I1.EQ.I2).AND.(J1.EQ.J2)) THEN
          DOBS = 0.
          W = 1.
        ELSE
          DELLON  = REAL(I1-1)*ZDELLO(J1)
     &     - REAL(I2-1)*ZDELLO(J2)
          COSLON  = COS(DELLON*RAD)
          DOBS2   = COSLON*COSPH(J1)*COSPH(J2) +
     &     SINPH(J1)*SINPH(J2)
          DOBS    = ACOS(DOBS2)
          W       = EXP(-DOBS/DIST)
        ENDIF

        DO JOBS=1,NOBS                                            
          IX = IMEAS(JOBS)                                            
          JX = JMEAS(JOBS)                                            
          XNEW = XNEW + W*XM(IOBS,JOBS)*(XME(IX,JX)-XMO(IX,JX))       
        ENDDO
      ENDDO
                                                                        
      XOI = XMO(I,J) + XNEW                                             

      RETURN                                                            
      END SUBROUTINE ANALYSE
