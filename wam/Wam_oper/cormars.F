      SUBROUTINE CORMARS (FELD1, ID1, ID2, NGX, NGY, IPAR)              

! ----------------------------------------------------------------------

!**** *CORMARS* - CORRECTS MARS WAVE DATA.                              

!     HEINZ GUNTHER    ECMWF/GKSS    JANUARY 1990                       

!     PURPOSE                                                           
!     -------                                                           

!       TO CORRECT WAVE DATA WHICH ARE STORED IN MARS.                  

!     METHOD                                                            
!     ------                                                            

!        THE DATA ARE REORDERED AND MISSING VALUES ARE INTEPOLATED.     

!**   INTERFACE                                                         
!     ---------                                                         

!       *CALL* *CORMARS (FELD1, ID1, ID2, NGX, NGY, IPAR)*              

!          FELD1   REAL      INPUT/OUTPUT DATA FIELD.                   
!          ID1     INTEGER   FIRST ARRAY DIMENSION                      
!          ID2     INTEGER   SECOND ARRAY DIMENSION                     
!          NGX     INTEGER   FIRST ARRAY DIMENSION USED.                
!          NGY     INTEGER   SECOND ARRAY DIMENSION USED.               
!          IPAR    INTEGER   MARS PARAMETER OF DATA .                   

!     EXTERNALS                                                         
!     ---------                                                         

!          NONE                                                         

!     REFERENCES                                                        
!     ----------                                                        

!          NONE                                                         

! ----------------------------------------------------------------------

      DIMENSION FELD1(ID1,ID2), FELD2(ID1,ID2)                          

! ----------------------------------------------------------------------

      PARAMETER (PI=3.141592654, DEG=180./PI, RAD = PI/180.)            
      DIMENSION IN(121), IE(7), IL(16), KL(16)                          
                                                                        
      DATA IN /  2,  3,  4,  5,  6,  7,  7,  9, 10, 11, 12, 12, 13, 15, 
     &          16, 17, 18, 19, 20, 21, 22, 23, 23, 24, 25, 26, 28, 29, 
     &          30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 
     &          44, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 
     &          58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 
     &          72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 
     &          86, 87, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 
     &          99,100,101,102,103,104,105,107,108,109,110,111,112,113, 
     &         114,115,116,117,118,119,120, 121, 1/                     
                                                                        
      DATA IE / 1, 4, 8, 14, 27, 53, 106 /                              
                                                                        
      DATA IL / 1,  1,  1,  4,  4, 14, 14, 14, 27, 27, 27, 53,          
     &         53, 53, 106, 106 /                                       
                                                                        
      DATA KL / 7, 10, 23,  4, 25, 26, 29, 30, 20, 21, 22,  6,          
     &          7,  8,   5,  26 /                                       

! ----------------------------------------------------------------------

!*    1. CHECK DIMENSION.                                               
!        ----------------                                               

      IF (NGX.NE.121) THEN                                              
        WRITE (6,*) 'NGX IS NOT 121, NGX = ', NGX                      
        CALL ABORT1                                                     
      ENDIF                                                             

! ----------------------------------------------------------------------

!*    2. REARANGE DATA.                                                 
!        --------------                                                 

      DO  K=1,NGY 
        DO I=1,NGX                 
          FELD2(IN(I),k) = FELD1(I,K)
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    3. INTERPOLATE MISSING VALUES.                                    
!        ---------------------------                                    

      DO K=1,NGY                                                   
        DO I=1,7                                                     
          IE1 = IE(I)-1
          IF (IE1.EQ.0) IE1 = 120
          IE2 = IE(I)+1
          FELD2(IE(I),K) = 0.
          IF (IPAR.EQ.230) THEN
            CO = COS(FELD2(IE1,K)*RAD)+COS(FELD2(IE2,K)*RAD)
            SI = SIN(FELD2(IE1,K)*RAD)+SIN(FELD2(IE2,K)*RAD)
            IF (SI.NE.0..OR.CO.NE.0.) FELD2(IE(I),K) = ATAN2(SI,CO)*DEG
          ELSE
            FELD2(IE(I),K) = 0.5*(FELD2(IE1,K)+ FELD2(IE2,K))
          ENDIF
        ENDDO
      ENDDO

! ----------------------------------------------------------------------

!*    4. DELETE LANDPOINTS WHICH ARE INTEPOLATED.                       
!        ----------------------------------------                       

      DO K=1,16
        FELD2(IL(K),KL(K)) = 0.
      ENDDO

! ----------------------------------------------------------------------

!*    5. COPY TO OUTPUT ARRAY.
!        ---------------------

      DO K=1,NGY
        DO I=1,NGX
          FELD1(I,K) = FELD2(I,K)
        ENDDO
      ENDDO
                                                                        
      RETURN 
      END SUBROUTINE CORMARS


