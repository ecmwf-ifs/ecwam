      SUBROUTINE FLDINTER (IU06, ITEST, NGPTOTG, NC, NR, NFIELDS,FIELDS,
     &                     NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
     &                     IXLG, KXLT, NIBLO, IGL, IJS, IJL,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                     RLONRGG, IKBLOCK, UWND, VWND, AIRD, PZIDL,  
     &                     LADEN, ROAIR, LGUST,LLINTERPOL,
     &                     DK1M, DII1M, DIIP1M, KKM, IIM, IIPM)
! ----------------------------------------------------------------------    

!***  *FLDINTER* - INTERPOLATION OF ATMOSPHERIC FIELDS OVER
!                  THE WAVE MODEL GRID USING BILINEAR
!                  INTERPOLATION.

!      S. ABDALLA   ECMWF   DECEMBER 2000. 

!     PURPOSE.                                                          
!     --------                                                          

!     IT INTERPOLATES ATMOSPHERIC DATA FIELDS TO PRODUCE A FIELD 
!     ON THE WAVE MODEL GRID USING BILINEAR INTERPOLATION.
!     IT MAKES USE OF THE INTERPOLATION COEFFICIENTS PRODUCED BY 
!     SUBROUTINE "INITIALINT" WHICH SHOULD BE CALLED BEFORE ANY
!     CALL TO "FLDINTER".
!     IT ONLY WORKS IF THE INPUT AND OUTPUT GRIDS ARE LATITUDE/LONGITUDE
!     (REGULAR OR IRREGULAR) OR GAUSSIAN GRIDS (FULL OR REDUCED) !!!!

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL* *FLDINTER* (IU06, ITEST, NGPTOTG, NC, NR, FIELDS,
!    &                    NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
!    &                    IXLG, KXLT, NIBLO, IGL, IJS, IJL, 
!    &                    AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
!    &                    RLONRGG, IKBLOCK, UWND, VWND, AIRD, PZIDL,  
!    &                    LADEN, ROAIR, LGUST,LLINTERPOL,
!    &                    DK1M, DII1M, DIIP1M, KKM, IIM, IIPM)
!
!        *IU06*   - OUTPUT UNIT.
!        *ITEST*  - TEST MESSAGE LEVEL.
!        ATMOSPHERIC MODEL GRID AND FIELD (INPUT):
!        *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!        *RLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE (INPUT FIELD)
!        *FIELDS* - ATMOSPHERIC FIELDS AS FOLLOWS:
!                   FIELDS(:,1) = U COMPONENT OF NEUTRAL WIND SPEED (U10)
!                   FIELDS(:,2) = V COMPONENT OF NEUTRAL WIND SPEED (V10)
!                   FIELDS(:,3) = AIR DENSITY
!                   FIELDS(:,4) = ZI/L USED FOR GUSTINESS
!                   FIELDS(:,5) = RESERVED FOR ICE BOUNDARY

!        WAVE MODEL GRID SPECIFICATION (INPUT):
!        *NGX*    - NUMBER OF COLUMNS IN ARRAY FIELD USED.              
!        *NGY*    - NUMBER OF ROWS    IN ARRAY FIELD USED.              
!        *KRGG*   - GRID DEFINITION PARAMETER (0=REGULAR, 1=IRREGULAR)
!        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE (OUTPUT FIELD)
!        *XDELLA* - GRID POINT SPACING BETWEEN LATITUDES.
!        *ZDELLO* - GRID POINT SPACING PER LATITUDES. 
!        *IXLG*   - INTEGER  LONG. GRID INDEX.
!        *KXLT*   - INTEGER  LAT. GRID INDEX.
!        *NIBLO*  - FIRST DIMENSION OF IXLG AND KXLT.
!        *IGL*    - SECOND DIMENSION OF IXLG AND KXLT.
!        *IJS*    - SMALLEST WAM GRID POINT INDEX USED BY THE PE.
!        *IJL*    - LARGEST WAM GRID POINT INDEX USED BY THE PE.
!        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).           
!        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).           
!        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).           
!        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).           
!        *IPERIODIC* - SPECIFIES IF WAM GRID IS PERIODIC (=1) OR NOT (=0)
!        INTERPOLATION COEFFICIENTS & INDECES (INPUT):
!        *IKBLOCK* - BLOCK INDEX FOR THE ATMOSPHERIC FIELDS.
!        *LLINTERPOL* - FLAG (TRUE=DO INTERPOLATION, FALSE=ASSIGN ONLY)
!        *DK1M*   - COEFFICIENT
!        *DII1M*  - COEFFICIENT
!        *DIIP1M* - COEFFICIENT
!        *KKM*    - INDEX
!        *IIM*    - INDEX
!        *IIPM*   - INDEX
!        OUTPUT:
!        *UWND*   - INTERPOLATED FIELD
!        *VWND*   - INTERPOLATED FIELD
!        *AIRD*   - INTERPOLATED FIELD
!        *PZIDL*  - INTERPOLATED FIELD


!     EXTERNALS.                                                        
!     ----------                                                        

! ----------------------------------------------------------------------

      INTEGER KLONRGG(NGY)

      INTEGER NGPTOTG, NC, NR, NIBLO, IGL, IJS, IJL 
      INTEGER, DIMENSION(NIBLO,IGL) :: IXLG, KXLT
      INTEGER RLONRGG(NR)
      INTEGER :: IKBLOCK(0:NC+1,NR)
      REAL :: ROAIR
      REAL FIELDS(NGPTOTG,NFIELDS)
      LOGICAL:: LADEN, LGUST, LLINTERPOL
      
      REAL, DIMENSION(NGX,NGY) :: UWND, VWND, AIRD, PZIDL

      REAL ZDELLO(NGY)

      INTEGER  KKM(NGY), IIM(NGX,NGY), IIPM(NGX,NGY)
      REAL DK1M(NGY), DII1M(NGX,NGY), DIIP1M(NGX,NGY)
! ----------------------------------------------------------------------

      IF(.NOT.LLINTERPOL) THEN

!       REARRANGE DATA FIELD.
!       --------------------

        DO IG=1,IGL
          DO IJ = IJS,IJL
            I = IXLG(IJ,IG)
            K = NGY-KXLT(IJ,IG) + 1

            UWND(I,K) = FIELDS(IKBLOCK(I,K),1) 
            VWND(I,K) = FIELDS(IKBLOCK(I,K),2) 
            IF (LADEN) THEN
              AIRD(I,K) = FIELDS(IKBLOCK(I,K),3) 
            ELSE
              AIRD(I,K) = ROAIR
            ENDIF
            IF (LGUST) THEN
              PZIDL(I,K) = FIELDS(IKBLOCK(I,K),4) 
            ELSE
              PZIDL(I,K) = 0.
            ENDIF

!!! nothing is done for the ice boundary !!!

          ENDDO
        ENDDO

      ELSE

!       INTERPOLATE TO WAVE MODEL GRID
!       ------------------------------

        DO IG=1,IGL
          DO IJ = IJS,IJL 
            I = IXLG(IJ,IG)
            K = NGY-KXLT(IJ,IG) + 1

            KK = KKM(K)
            KSN= NR-KK+1
            KK1= MIN(KK+1,NR)
            KSN1=NR-KK1+1
            DK1= DK1M(K) 
            DK2=1.-DK1

            II = IIM(I,K)
            II1 = MIN(II+1,RLONRGG(KSN)+IPERIODIC)
            DII1=DII1M(I,K) 
            DII2=1.-DII1

            IIP = IIPM(I,K) 
            IIP1 = MIN(IIP+1,RLONRGG(KSN1)+IPERIODIC)
            DIIP1= DIIP1M(I,K) 
            DIIP2=1.-DIIP1

            UWND(I,K)=DK2*( DII2*FIELDS(IKBLOCK(II,KK),1) +
     &                      DII1*FIELDS(IKBLOCK(II1,KK),1) ) +
     &                DK1*( DIIP2*FIELDS(IKBLOCK(IIP,KK1),1) +
     &                      DIIP1*FIELDS(IKBLOCK(IIP1,KK1),1) )

            VWND(I,K)=DK2*( DII2*FIELDS(IKBLOCK(II,KK),2) +
     &                      DII1*FIELDS(IKBLOCK(II1,KK),2) ) +
     &                DK1*( DIIP2*FIELDS(IKBLOCK(IIP,KK1),2) +
     &                      DIIP1*FIELDS(IKBLOCK(IIP1,KK1),2) )

            IF (LADEN) THEN
              AIRD(I,K)=DK2*( DII2*FIELDS(IKBLOCK(II,KK),3) +
     &                        DII1*FIELDS(IKBLOCK(II1,KK),3) ) +
     &                  DK1*( DIIP2*FIELDS(IKBLOCK(IIP,KK1),3) +
     &                        DIIP1*FIELDS(IKBLOCK(IIP1,KK1),3) )
            ELSE
              AIRD(I,K) = ROAIR
            ENDIF
            IF (LGUST) THEN
              PZIDL(I,K)=DK2*( DII2*FIELDS(IKBLOCK(II,KK),4) +
     &                         DII1*FIELDS(IKBLOCK(II1,KK),4) ) +
     &                   DK1*( DIIP2*FIELDS(IKBLOCK(IIP,KK1),4) +
     &                         DIIP1*FIELDS(IKBLOCK(IIP1,KK1),4) )
            ELSE
              PZIDL(I,K) = 0.
            ENDIF

!!! nothing is done for the ice boundary !!!


          ENDDO
        ENDDO

      ENDIF

      RETURN
      END SUBROUTINE FLDINTER

