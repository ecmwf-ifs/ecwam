      SUBROUTINE OIFIELD (XME, XMO, XOI, SIGOBS, SIGMOD, 
     .                    LMAX, NDIM1, NDIM2)            
C                                                                       
C**** *OIFIELD* - OPTIMUM INTERPOLATION.                                
C                                                                       
C     P.LIONELLO     ECMWF       APRIL 1990                             
C     J. BIDLOT      ECMWF       SEPTEMBER 96  MESSAGE PASSING
C                                                                       
C     PURPOSE.                                                          
C     --------                                                          
C                                                                       
C       TO PRODUCE A MAP OF THE FIELD X, MERGING MEASUREMENT            
C       AND MODEL , BY OPTIMUM INTERPOLATION. THE ARRAY XOI             
C       AT THE END OF THE SUBROUTINE CONTAINS THE VALUES                
C       TO BE USED TO ANALYSE THE SPECTRA, HAVING NEGATIVE RETURN       
C       CODES WHERE O.I. PRODUCED NO RESULTS.                           
C       THE FIRST GUESS FIELD XMO IS NOT MODIFIED                       
C       IN THIS SUBROUTINE.                                             
C                                                                       
C**   INTERFACE.                                                        
C     ----------                                                        
C                                                                       
C       *CALL* *OIFIELD (XME, XMO, XOI, SIGOBS, SIGMOD)*      
C                                                                       
C         *XME*     REAL     FIELD FROM MEASUREMENTS.                   
C         *XMO*     REAL     FIELD FROM MODEL (FIRST GUESS).            
C         *XOI*     REAL     FIELD FROM O.I.                            
C         *SIGOBS*  REAL     MEASUREMENT SCATTER.                       
C         *SIGMOD*  REAL     MODEL SCATTER.                             
C         *NGX*     INTEGER  FIRST FIELD DIMENSION.                     
C         *NGY*     INTEGER  SECOND FIELD DIMENSION.                    
C                                                                       
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "parcons.h"
C
#include "comgrid.h"
C
#include "commap.h"
C                                                                       
      DIMENSION XMO(NGX,NGY), XME(NGX,NGY), XOI(NGX,NGY)                
      DIMENSION IMEAS(NDIM2), JMEAS(NDIM2), XM(NDIM2,NDIM2),            
     1          D(NDIM2,NDIM2), P(NDIM2,NDIM2), V(NDIM2)                
      DIMENSION LLON(NGY)
C                                                                       
C   LMAX = MAXIMUM NUMBER OF GRID STEPS TO SPREAD MEASURED HS AND USTAR 
C                                                                       
C ----------------------------------------------------------------------
C    
C
C*    1. INITIALISE SEARCH DISTANCE AND CORRELATION LENGTH.
C     -----------------------------------------------------
C
      DIST  = 300000./R
      DIST2 = 2.*DIST
                                                                    
C       THE MEASUREMENTS ARE STORED IN THE O.I. FIELD.                  
C       NEGATIVE VALUES WILL REMAIN WHERE O.I. WILL PRODUCE NO RESULTS  
C                                                                       
      DO 1001 J=1,NGY                                                   
         DO 1002 I=1,NGX                                                
            XOI(I,J) = XME(I,J)                                         
 1002    CONTINUE                                                       
 1001 CONTINUE                                                          
C                                                                       
C       DETERMINE OBSERVATION CORRELATION MATRIX  D                     
C                                                                       
      DO 1010 J=1,NDIM2                                                 
         DO 1011 I=1,NDIM2                                              
            D(I,J) = 0.                                                 
 1011    CONTINUE                                                       
 1010 CONTINUE                                                          
C                                                                       
      DO 1020 J=1,NDIM2                                                 
         D(J,J) = (SIGOBS/SIGMOD)**2                                    
         V(J) = 0.                                                      
 1020 CONTINUE
C
      DO 1030 JSN=1,NGY
         SIN2 = SINPH(JSN)**2
         COS2 = COSPH(JSN)**2
         COSLON = (COS(DIST2)-SIN2)/COS2
         DELLON = ACOS(COSLON)
         LLON(JSN) = NINT(DEG*DELLON/ZDELLO(JSN))
         LLON(JSN) = MIN(LMAX,LLON(JSN))
 1030 CONTINUE         
C
         LLAT = NINT(DEG*DIST2/XDELLA)
         LLAT = MIN(LMAX,LLAT)
C
C                                                                       
C       LOOP OVER GRID POINTS.                                          
C       ----------------------                                          
C                                                                       
      DO 2000 IG = 1,IGL
         DO 2001 IJ = IJS(IG),IJL(IG)
                                                                        
            I = IXLG(IJ,IG)
            J = NGY-KXLT(IJ,IG) + 1
C        
         IF(XMO(I,J).LT.0.01) GOTO 2001                                                               
            JSN = NGY-J+1                                                   
            NOBS=0                                                      
            DO 1050 IOBS=1,NDIM2                                        
               IMEAS(IOBS) = 0                                          
               JMEAS(IOBS) = 0                                          
 1050       CONTINUE                                                    
                                                                        
C      FIND NUMBER AND POSITION OF OBSERVATIONS IN A BOX                
C      OF SIZE 2*LMAX+1 AROUND ANALYSIS POINT I,J.                      
                                                                        
            J1 = J-LLAT                                                 
            J1 = MAX(J1,1)                                              
            J2 = J+LLAT
            J2 = MIN(J2,NGY)                                            
            DO 1100 JX=J1,J2
               JXSN = NGY-JX+1

               NEWI=NINT((I-1)*ZDELLO(JSN)/ZDELLO(JXSN))+1
               I1 = MAX(1,NEWI-LLON(JXSN))
               I2 = MIN(NLONRGG(JXSN),NEWI+LLON(JXSN))

               DO 1110 IX = I1,I2
                  IF (XME(IX,JX).GT.0.) THEN
                     IF (XMO(IX,JX).GE.0.01) THEN
                       NOBS = NOBS+1
                       IMEAS(NOBS) = IX
                       JMEAS(NOBS) = JX
                     ENDIF
                  ENDIF
1110           CONTINUE

               IF(IPER.EQ.1) THEN
                  I1 = 1
                  I2 = NEWI+LLON(JXSN)-NLONRGG(JXSN)
                  DO 1120 IX = I1,I2
                    IF (XME(IX,JX).GT.0.) THEN
                       IF (XMO(IX,JX).GE.0.01) THEN
                         NOBS = NOBS+1
                         IMEAS(NOBS) = IX
                         JMEAS(NOBS) = JX
                       ENDIF
                     ENDIF
1120              CONTINUE

                 I1 = NEWI-LLON(JXSN)+NLONRGG(JXSN)
                 I2 = NLONRGG(JXSN)
                 DO 1130 IX = I1,I2
                    IF (XME(IX,JX).GT.0.) THEN
                       IF (XMO(IX,JX).GE.0.01) THEN
                         NOBS = NOBS+1
                         IMEAS(NOBS) = IX
                         JMEAS(NOBS) = JX
                       ENDIF
                    ENDIF
1130             CONTINUE
               ENDIF


1100        CONTINUE

C                                                                       
            IF (NOBS.EQ.0) GOTO 2001                                    
                                                                        
C     DETERMINE MODEL CORRELATION MATRIX AT OBSERVATION                 
C     POINTS.                                                           
                                                                        
            DO 1150 IOBS=1,NOBS                                         
               DO 1160 JOBS=IOBS+1,NOBS                                      
                  I1 = IMEAS(IOBS)
                  I2 = IMEAS(JOBS)
                  J1 = NGY - JMEAS(IOBS) + 1
                  J2 = NGY - JMEAS(JOBS) + 1
C                  
                  DELLON = FLOAT(I1-1)*ZDELLO(J1)
     1                     - FLOAT(I2-1)*ZDELLO(J2)
                  COSLON = COS(DELLON*RAD)
                  DOBS2  = COSLON*COSPH(J1)*COSPH(J2) +
     1                          SINPH(J1)*SINPH(J2)
                  DOBS   = ACOS(DOBS2)
C
                  P(IOBS,JOBS) = EXP(-DOBS/DIST)
                  P(JOBS,IOBS) = P(IOBS,JOBS) 
 1160          CONTINUE                                                 
                  P(IOBS,IOBS) = 1.
 1150       CONTINUE                                                    
                                                                        
C   DETERMINE M=P+D                                                     
                                                                        
            CALL SETMAT(XM,NOBS,P,D,NDIM2)                              
                                                                        
C   INVERSE MATRIX M                                                    
                                                                        
            COND=0.                                                     
            EPS=0.                                                      
                                                                        
            CALL SYMINV(XM,NOBS,NOBS,COND,V,EPS)                        
                                                                        
C   PRODUCE FIELD AT POINT I,J ACCORDING TO OPTIMUM INTERPOLATION       
                                                                        
            CALL ANALYSE (XM, NOBS, XME, XMO, XOI(I,J),
     1                    IMEAS, JMEAS, I, J, DIST, LMAX, NDIM1, NDIM2)
 2001    CONTINUE                                                       
 2000 CONTINUE                                                          
                                                                        
      RETURN                                                            
      END                                                               
