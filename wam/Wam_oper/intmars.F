      SUBROUTINE INTMARS (IU06, ITEST, FILNM, NKSEK0, NKSEK1, NKSEK2,
     &                   NKSEK3, NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                   NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
     &                   AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                   CDATE, IFORP, IPARAM, KZLEV, FIELD )

! -------------------------------------------------------------------

!***  *INTMARS* - READS IN AND UNPACKS FIELDS FROM MARS       
!                 AND PUTS IT ON THE WAVE MODEL GRID USING WAVE GRID
!                 MODIFIED BILINEAR INTERPOLATION IF NECESSARY.

!     J. BIDLOT  ECMWF  FEBRUARY 2000. 

!     PURPOSE.                                                          
!     --------                                                          


!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL INTMARS* (IU06, ITEST, FILNM, NKSEK0, NKSEK1, NKSEK2,
!    &                 NKSEK3, NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
!    &                 NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
!    &                 AMOWEP, AMOSOP, AMOEAP, AMONOP,
!    &                 CDATE, IFORP, IPARAM, KZLEV, FIELD )

!        *IU06*   _ OUTPUT UNIT.
!        *ITEST*  - TEST MESSAGE LEVEL.
!        *FILNM*  - DATA INPUT FILENAME. 
!        *NKSEK0* - SIZE OF KSEC0
!        *NKSEK1* - SIZE OF KSEC1
!        *NKSEK2* - SIZE OF KSEC2
!        *NKSEK3* - SIZE OF KSEC3
!        *NKSEK4* - SIZE OF KSEC4
!        *NPSEK2* - SIZE OF PSEC2
!        *NPSEK3* - SIZE OF PSEC3
!        *KSEK3*  - DEFAULT VALUES FOR KSEC3
!        *PSEK3*  - DEFAULT VALUES FOR PSEC3
!        WAVE MODEL GRID SPECIFICATION:
!        *NGX  *  - NUMBER OF COLUMNS IN ARRAY FIELD USED.
!        *NGY  *  - NUMBER OF ROWS    IN ARRAY FIELD USED.
!        *KRGG*   - GRID DEFINITION PARAMETER (O=REGULAR, 1=IRREGULAR)
!        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE
!        *XDELLA* - GRID POINT SPACING BETWEEN LATITUDES.
!        *ZDELL0* - GRID POINT SPACING PER LATITUDES.
!        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).
!        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).
!        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).
!        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).
!        OUTPUT:
!        *CDATE*  - DATE/TIME OF THE DATA READ.
!        *IFORP*  - FORCAST PERIOD IN SECONDS.
!        *IPARAM* - DATA CODE
!        *KZLEV*  - REFERENCE LEVEL IN full METER
!                   SHOULD BE 0 EXCEPT FOR 233 AND 245 WHERE IT
!                   MIGHT BE DIFFERENT THAN 0 PROVIDED IT WAS
!                   CODED IN KSEC1(7) AND KSEC1(8) SEE *GRIBPAC*
!        *FIELD*  - UNPACKED DATA ON WAVE MODEL GRID.


!     EXTERNALS.                                                        
!     ----------                                                        

!     *PBOPEN*         OPEN FILE TO BE READ BY PBGRIB
!     *PBGRIB*         READ GRIB FILE
!     *PBCLOSE*        CLOSE FILE

!     *GRB2WGRD*       UNPACKS GRIB and PUTS IT ON THE WAVE MODEL GRID.

! ----------------------------------------------------------------------

      PARAMETER (NBIT = 190000)

      INTEGER KLONRGG(NGY)                                          
      INTEGER KSEK3(NKSEK3)
      INTEGER, ALLOCATABLE ::  KGRIB(:)

      REAL FIELD(NGX,NGY)                                          
      REAL ZDELLO(NGY)
      REAL PSEK3(NPSEK3)

      CHARACTER CDATE*12, FILNM*24
                                                                        

!*    1.0 READ MARS DATA IN GRIB CODE FORMAT.                           
!         -----------------------------------                           

      CALL PBOPEN(IUNIT1,FILNM,'r',KRET)                            
      IF(KRET.LT.0) THEN
         IF(KRET.EQ.-1) THEN
           WRITE(IU06,*) '***********************************'
           WRITE(IU06,*) '*                                 *'
           WRITE(IU06,*) '*  FATAL ERROR IN SUB. INTMARS    *'
           WRITE(IU06,*) '*  ============================   *'
           WRITE(IU06,*) '*                                 *'
           WRITE(IU06,*) '*  COULD NOT OPEN '//FILNM
           WRITE(IU06,*) '*                                 *'
           WRITE(IU06,*) '***********************************'
           CALL ABORT1
         ENDIF
         IF(KRET.EQ.-2) THEN
           WRITE(IU06,*) '***********************************'
           WRITE(IU06,*) '*                                 *'
           WRITE(IU06,*) '*  FATAL ERROR IN SUB. INTMARS    *'
           WRITE(IU06,*) '*  ============================   *'
           WRITE(IU06,*) '*                                 *'
           WRITE(IU06,*) '*  INVALID FILENAME'//FILNM
           WRITE(IU06,*) '*                                 *'
           WRITE(IU06,*) '***********************************'
           CALL ABORT1
         ENDIF
      ENDIF

      ISIZE=NBIT
1001  KPLENG=ISIZE*KIND(ISIZE)
      IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))
      CALL PBGRIB(IUNIT1,KGRIB,KPLENG,KLEN,KRET)
      IF(KRET.EQ.-1) THEN 
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. INTMARS    *'
        WRITE(IU06,*) '*  ============================   *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  PBGRIB : END OF FILE ENCOUNTED' 
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF
      IF(KRET.EQ.-2) THEN 
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. INTMARS    *'
        WRITE(IU06,*) '*  ============================   *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  PBGRIB : FILE HANDLING ERROR'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF
      IF(KRET.EQ.-3) THEN 
        DEALLOCATE(KGRIB)
        CALL KGRIBSIZE(IU06, IUNIT1, KLEN, ISIZE, 'INTMARS')
        GOTO 1001
      ENDIF

      CALL PBCLOSE(IUNIT1,KRET)
      IF(KRET.EQ.-2) THEN 
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. INTMARS    *'
        WRITE(IU06,*) '*  ============================   *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  PBCLOSE :FILE HANDLING ERROR'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF


! ----------------------------------------------------------------------

!*    2.0 UNPACKS MARS FIELD AND PUTS IT ON THE WAVE MODEL GRID.
!         ------------------------------------------------------

      KTEST=ITEST+1
      ZDUM =0.
      CALL GRB2WGRD (IU06, KTEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &               NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &               KGRIB, ISIZE,
     &               NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
     &               AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &               ZDUM, ZDUM, 
     &               CDATE, IFORP, IPARAM, KZLEV, IDUM, IDUM, FIELD )


      DEALLOCATE(KGRIB)

      RETURN
      END SUBROUTINE INTMARS
