      SUBROUTINE FLDINTER (IU06, ITEST, NGPTOTG, NC, NR, FLDIN,
     &                     NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                     RLONRGG, FLDOUT,  
     &                     LLINTERPOL,
     &                     DK1M, DII1M, DIIP1M, KKM, IIM, IIPM)
! ----------------------------------------------------------------------    

!***  *FLDINTER* - INTERPOLATION OF ATMOSPHERIC FIELDS OVER
!                  THE WAVE MODEL GRID USING BILINEAR
!                  INTERPOLATION.

!      S. ABDALLA   ECMWF   DECEMBER 2000. 

!     PURPOSE.                                                          
!     --------                                                          

!     IT INTERPOLATES AN ATMOSPHERIC DATA FIELD TO PRODUCE A FIELD 
!     ON THE WAVE MODEL GRID USING BILINEAR INTERPOLATION.
!     IT MAKES USE OF THE INTERPOLATION COEFFICIENTS PRODUCED BY 
!     SUBROUTINE "INITIALINT" WHICH SHOULD BE CALLED BEFORE ANY
!     CALL TO "FLDINTER".
!     IT ONLY WORKS IF THE INPUT AND OUTPUT GRIDS ARE LATITUDE/LONGITUDE
!     (REGULAR OR IRREGULAR) OR GAUSSIAN GRIDS (FULL OR REDUCED) !!!!

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL* *FLDINTER* (IU06, ITEST, NGPTOTG, NC, NR, FLDIN,
!    &                    NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
!    &                    AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
!    &                    RLONRGG, FLDOUT, LLINTERPOL,
!    &                    DK1M, DII1M, DIIP1M, KKM, IIM, IIPM)
!
!        *IU06*   - OUTPUT UNIT.
!        *ITEST*  - TEST MESSAGE LEVEL.
!        ATMOSPHERIC MODEL GRID AND FIELD (INPUT):
!        *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!        *RLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE (INPUT FIELD)
!        *FLDIN*  - INPUT ATMOSPHERIC FIELD WHICH WILL BE INTERPOLATED
!        WAVE MODEL GRID SPECIFICATION (INPUT):
!        *NGX*    - NUMBER OF COLUMNS IN ARRAY FIELD USED.              
!        *NGY*    - NUMBER OF ROWS    IN ARRAY FIELD USED.              
!        *KRGG*   - GRID DEFINITION PARAMETER (0=REGULAR, 1=IRREGULAR)
!        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE (OUTPUT FIELD)
!        *XDELLA* - GRID POINT SPACING BETWEEN LATITUDES.
!        *ZDELLO* - GRID POINT SPACING PER LATITUDES. 
!        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).           
!        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).           
!        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).           
!        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).           
!        *IPERIODIC* - SPECIFIES IF WAM GRID IS PERIODIC (=1) OR NOT (=0)
!        INTERPOLATION COEFFICIENTS & INDECES (INPUT):
!        *LLINTERPOL* - FLAG (TRUE=DO INTERPOLATION, FALSE=ASSIGN ONLY)
!        *DK1M*   - COEFFICIENT
!        *DII1M*  - COEFFICIENT
!        *DIIP1M* - COEFFICIENT
!        *KKM*    - INDEX
!        *IIM*    - INDEX
!        *IIPM*   - INDEX
!        OUTPUT:
!        *FLDOUT* - INTERPOLATED FIELD

!     EXTERNALS.                                                        
!     ----------                                                        

! ----------------------------------------------------------------------

      INTEGER KLONRGG(NGY)

      INTEGER NGPTOTG, NC, NR
      INTEGER RLONRGG(NR)
      REAL FLDIN(NGPTOTG)
      LOGICAL LLINTERPOL
      
      REAL FLDOUT(NGX,NGY)                            
      REAL ZDELLO(NGY)
      REAL, ALLOCATABLE :: WORK(:,:)

      INTEGER  KKM(NGY), IIM(NGX,NGY), IIPM(NGX,NGY)
      REAL DK1M(NGY), DII1M(NGX,NGY), DIIP1M(NGX,NGY)
! ----------------------------------------------------------------------

      IF(.NOT.LLINTERPOL) THEN

!       REARRANGE DATA FIELD.
!       --------------------

        L = 0                                                          
        DO K=1,NGY                                                
          JSN=NGY-K+1
          DO I=1,KLONRGG(JSN)
            L = L+1                                                     
            FLDOUT(I,K) = FLDIN(L) 
          ENDDO
        ENDDO

      ELSE

!       INTERPOLATE TO WAVE MODEL GRID
!       ------------------------------

        IF (IPERIODIC.EQ.1) THEN
          ALLOCATE(WORK(0:NC+1, NR))
        ELSE
          ALLOCATE(WORK(NC,NR))
        ENDIF

        L = 0                                                          
        DO KK=1,NR                                                
          JSN=NR-KK+1
          DO II=1,RLONRGG(JSN)
            L = L+1                                                     
            WORK(II,KK) = FLDIN(L)                                       
          ENDDO
        ENDDO
        IF(NGPTOTG.NE.L)THEN
          WRITE(IU06,*)' ######################################'
          WRITE(IU06,*)' ####  WARNING   WARNING   WARNING  ###'
          WRITE(IU06,*)' ######################################'
          WRITE(IU06,*)' ##  NGPTOTG  IS  NOT  EQUAL TO  L  ###'
          WRITE(IU06,*)' ##  NGPTOTG = ', NGPTOTG 
          WRITE(IU06,*)' ##        L = ', L 
          WRITE(IU06,*)' ######################################'
          WRITE(IU06,*)' ##  NC = ', NC, '   NR = ', NR
          WRITE(IU06,*)' ##  SUM RLONRGG = ',SUM(RLONRGG(1:NR))
          WRITE(IU06,*)' ##  RLONRGG = ', RLONRGG(1:NR)
          WRITE(IU06,*)' ######################################'
          WRITE(IU06,*)' ######################################'
          WRITE(IU06,*)' ######################################'
          CALL FLUSH(IU06)
        ENDIF
        IF(IPERIODIC.EQ.1) THEN
          DO KK=1,NR                                                
          JSN=NR-KK+1
          WORK(0,KK)= WORK(RLONRGG(JSN),KK)
          WORK(RLONRGG(JSN)+1,KK)= WORK(1,KK)
          ENDDO
        ENDIF

        DO K=1,NGY                                                

          JSN=NGY-K+1

          KK = KKM(K)
          KSN= NR-KK+1
          KK1= MIN(KK+1,NR)
          KSN1=NR-KK1+1
          DK1= DK1M(K) 
          DK2=1.-DK1


          DO I=1,KLONRGG(JSN)

            II = IIM(I,K)
            II1 = MIN(II+1,RLONRGG(KSN)+IPERIODIC)
            DII1=DII1M(I,K) 
            DII2=1.-DII1

            IIP = IIPM(I,K) 
            IIP1 = MIN(IIP+1,RLONRGG(KSN1)+IPERIODIC)
            DIIP1= DIIP1M(I,K) 
            DIIP2=1.-DIIP1

            FLDOUT(I,K)=DK2*(DII2*WORK(II,KK)+DII1*WORK(II1,KK)) +
     &             DK1*(DIIP2*WORK(IIP,KK1)+DIIP1*WORK(IIP1,KK1))
          ENDDO
        ENDDO

        DEALLOCATE(WORK)
      ENDIF

      RETURN
      END SUBROUTINE FLDINTER

