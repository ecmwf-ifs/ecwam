      SUBROUTINE OIFIELD (XME, XMO, XOI, SIGOBS, LWAVE, SIGMOD, 
     &                    LMAX, NDIM2)            

!**** *OIFIELD* - OPTIMUM INTERPOLATION.                                

!     P.LIONELLO     ECMWF       APRIL 1990                             
!     J. BIDLOT      ECMWF       SEPTEMBER 96  MESSAGE PASSING

!     PURPOSE.                                                          
!     --------                                                          

!       TO PRODUCE A MAP OF THE FIELD X, MERGING MEASUREMENT            
!       AND MODEL , BY OPTIMUM INTERPOLATION. THE ARRAY XOI             
!       AT THE END OF THE SUBROUTINE CONTAINS THE VALUES                
!       TO BE USED TO ANALYSE THE SPECTRA, HAVING NEGATIVE RETURN       
!       CODES WHERE O.I. PRODUCED NO RESULTS.                           
!       THE FIRST GUESS FIELD XMO IS NOT MODIFIED                       
!       IN THIS SUBROUTINE.                                             

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *OIFIELD (XME, XMO, XOI, SIGOBS, SIGMOD)*      

!         *XME*     REAL     FIELD FROM MEASUREMENTS.                   
!         *XMO*     REAL     FIELD FROM MODEL (FIRST GUESS).            
!         *XOI*     REAL     FIELD FROM O.I.                            
!         *SIGOBS*  REAL     MEASUREMENT SCATTER.                       
!         *LWAVE*   LOGICAL  TRUE IF ASSIMILATION FOR WAVE HEIGHTS WITH
!                            ADJUSTED OBSERVATIONAL ERRORS.
!         *SIGMOD*  REAL     MODEL SCATTER.                             

! ----------------------------------------------------------------------

      USE YOWGRID  , ONLY : SINPH    ,COSPH    ,NLONRGG  ,IGL      ,
     &            IJS      ,IJL
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IPER     ,XDELLA   ,
     &            ZDELLO
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWPCONS , ONLY : RAD      ,DEG      ,R

! ----------------------------------------------------------------------

      DIMENSION XMO(NGX,NGY), XME(NGX,NGY), XOI(NGX,NGY)                
      INTEGER, DIMENSION(NDIM2) :: IMEAS, JMEAS
      REAL :: SIGMAOBS(NDIM2)
      REAL, ALLOCATABLE, DIMENSION(:) :: V
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XM, D, P 

      DIMENSION LLON(NGY)
 
      LOGICAL LWAVE

!   LMAX = MAXIMUM NUMBER OF GRID STEPS TO SPREAD MEASURED HS AND USTAR 

! ----------------------------------------------------------------------


!*    1. INITIALISE SEARCH DISTANCE AND CORRELATION LENGTH.
!     -----------------------------------------------------

      IF(XDELLA.LT.0.5) THEN
        DIST  = 150000./R
      ELSE
        DIST  = 300000./R
      ENDIF
      DIST2 = 2.*DIST

!     DEFINE COEFFICIENT FOR THE ALT DATA ERROR
      IF(LWAVE) THEN
        A=SIGOBS
        B=-LOG(0.001)
        HSMIN=1.
      ELSE
        A=0.
        B=0.
        HSMIN=1.
      ENDIF

                                                                    
!       THE MEASUREMENTS ARE STORED IN THE O.I. FIELD.                  
!       NEGATIVE VALUES WILL REMAIN WHERE O.I. WILL PRODUCE NO RESULTS  

      DO J=1,NGY                                                   
        DO I=1,NGX                                                
          XOI(I,J) = XME(I,J)                                         
        ENDDO
      ENDDO

      DO JSN=1,NGY
        SIN2 = SINPH(JSN)**2
        COS2 = COSPH(JSN)**2
        COSLON = (COS(DIST2)-SIN2)/COS2
        DELLON = ACOS(COSLON)
        LLON(JSN) = NINT(DEG*DELLON/ZDELLO(JSN))
        LLON(JSN) = MIN(LMAX,LLON(JSN))
      ENDDO

      LLAT = NINT(DEG*DIST2/XDELLA)
      LLAT = MIN(LMAX,LLAT)

!       LOOP OVER GRID POINTS.                                          
!       ----------------------                                          

      DO IG = 1,IGL
        DO IJ = IJS(IG),IJL(IG)
                                                                        
          I = IXLG(IJ,IG)
          J = NGY-KXLT(IJ,IG) + 1

          IF(XMO(I,J).LT.0.01) GOTO 2001
          JSN = NGY-J+1
          NOBS=0
          DO IOBS=1,NDIM2
            IMEAS(IOBS) = 0
            JMEAS(IOBS) = 0
          ENDDO
                                                                        
!      FIND NUMBER AND POSITION OF OBSERVATIONS IN A BOX                
!      OF SIZE 2*LMAX+1 AROUND ANALYSIS POINT I,J.                      
                                                                        
          J1 = J-LLAT                                                 
          J1 = MAX(J1,1)                                              
          J2 = J+LLAT
          J2 = MIN(J2,NGY)                                            
          DO JX=J1,J2
            JXSN = NGY-JX+1

            NEWI=NINT((I-1)*ZDELLO(JSN)/ZDELLO(JXSN))+1
            I1 = MAX(1,NEWI-LLON(JXSN))
            I2 = MIN(NLONRGG(JXSN),NEWI+LLON(JXSN))

            DO IX = I1,I2
              IF (XME(IX,JX).GT.0.) THEN
                IF (XMO(IX,JX).GE.0.01) THEN
                  NOBS = NOBS+1
                  SIGMAOBS(NOBS)=SIGOBS+A*EXP(-B*(XME(IX,JX)-HSMIN))
                  IMEAS(NOBS) = IX
                  JMEAS(NOBS) = JX
                ENDIF
              ENDIF
            ENDDO

            IF(IPER.EQ.1) THEN
              I1 = 1
              I2 = NEWI+LLON(JXSN)-NLONRGG(JXSN)
              DO IX = I1,I2
                IF (XME(IX,JX).GT.0.) THEN
                  IF (XMO(IX,JX).GE.0.01) THEN
                    NOBS = NOBS+1
                    SIGMAOBS(NOBS)=SIGOBS+A*EXP(-B*(XME(IX,JX)-HSMIN))
                    IMEAS(NOBS) = IX
                    JMEAS(NOBS) = JX
                  ENDIF
                ENDIF
              ENDDO

              I1 = NEWI-LLON(JXSN)+NLONRGG(JXSN)
              I2 = NLONRGG(JXSN)
              DO IX = I1,I2
                IF (XME(IX,JX).GT.0.) THEN
                  IF (XMO(IX,JX).GE.0.01) THEN
                    NOBS = NOBS+1
                    SIGMAOBS(NOBS)=SIGOBS+A*EXP(-B*(XME(IX,JX)-HSMIN))
                    IMEAS(NOBS) = IX
                    JMEAS(NOBS) = JX
                  ENDIF
                ENDIF
              ENDDO
            ENDIF


          ENDDO


          IF (NOBS.EQ.0) GOTO 2001                                    

          ALLOCATE(D(NOBS,NOBS))            
          ALLOCATE(P(NOBS,NOBS))            
          ALLOCATE(XM(NOBS,NOBS))            

!         DETERMINE OBSERVATION CORRELATION MATRIX  D                     

          D = 0.                                                 

          DO IOBS=1,NOBS                                         
            D(IOBS,IOBS) = (SIGMAOBS(IOBS)/SIGMOD)**2
          ENDDO

                                                                        
!         DETERMINE MODEL CORRELATION MATRIX AT OBSERVATION                 
!         POINTS.                                                           
                                                                        
          DO IOBS=1,NOBS                                         
            DO JOBS=IOBS+1,NOBS                                      
              I1 = IMEAS(IOBS)
              I2 = IMEAS(JOBS)
              J1 = NGY - JMEAS(IOBS) + 1
              J2 = NGY - JMEAS(JOBS) + 1

              DELLON = REAL(I1-1)*ZDELLO(J1)
     &         - REAL(I2-1)*ZDELLO(J2)
              COSLON = COS(DELLON*RAD)
              DOBS2  = COSLON*COSPH(J1)*COSPH(J2) +
     &         SINPH(J1)*SINPH(J2)
              DOBS   = ACOS(DOBS2)

              P(IOBS,JOBS) = EXP(-DOBS/DIST)
              P(JOBS,IOBS) = P(IOBS,JOBS) 
            ENDDO
            P(IOBS,IOBS) = 1.
          ENDDO
                                                                        
          XM=P+D

          DEALLOCATE(D)            
          DEALLOCATE(P)            

!         INVERSE MATRIX XM

          COND=0.
          EPS=0. 
          ALLOCATE(V(NOBS))            
          V = 0.

          CALL SYMINV(XM,NOBS,NOBS,COND,V,EPS)

          DEALLOCATE(V)            

!         PRODUCE FIELD AT POINT I,J ACCORDING TO OPTIMUM INTERPOLATION
                                                                        
          CALL ANALYSE (XM, NOBS, XME, XMO, XOI(I,J),
     &                  IMEAS, JMEAS, I, J, DIST, NDIM2)

          DEALLOCATE(XM)            

 2001     CONTINUE

        ENDDO
      ENDDO
                                                                        
      RETURN                                                            
      END SUBROUTINE OIFIELD
