
! ----------------------------------------------------------------------
      SUBROUTINE OUTWAMINPUT(FIELD, KSEC1, KKK, MMM, OFILENAME)

!**** *INTWAMINPUT* -

!     J. BIDLOT    ECMWF  AUGUST 2000 

!*    PURPOSE.
!     --------
!     TO ENCODE FIELD INTO GRIB AND WRITE THE RESULTS TO FILE 
!
!**   INTERFACE.
!     ----------
!     *CALL* *OUTWAMINPUT(FIELD, KSEC1, KKK, MMM, OFILENAME)

!     *FIELD*     FIELD TO BE ENCODED FIELD(NGX,NGY)
!     *KSEC1*     KSEC1 FOR FIELD (DIFFERENT THAN THE DEFAULT KSEK1
!     *KKK*       DIRECTION INDEX IF SPECTRAL FIELD (0 OTHERWISE)
!     *MMM*       FREQUENCY INDEX IF SPECTRAL FIELD (0 OTHERWISE)
!     *OFILENAME* OUTPUT FILENAME

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3   ,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3    ,
     &            PPMISS   ,PPEPS    ,PPREC    ,NTENCODE ,NGRBRESI ,
     &            IMDLGRBID_G,IMDLGRBID_M      ,NGRBRESS ,HOPERI   ,
     &            HOPERS 
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOSOP   ,AMONOP   ,XDELLA
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            CLDOMAIN
      USE YOWPCONS , ONLY : DEG
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IUOU, KRET
      INTEGER :: I, J, JSN, IC, IJ
      INTEGER :: KLENG, KLENP, KWORD
      INTEGER :: KKK, MMM, INDEX, KOUNT
      INTEGER :: KSEC0(NKSEK0), KSEC1(NKSEK1),KSEC2(NKSEK2),
     &           KSEC3(NKSEK3), KSEC4(NKSEK4)

      INTEGER, ALLOCATABLE :: KGRIB(:)

      REAL :: PSEC2(NPSEK2),PSEC3(NPSEK3)
      REAL :: FIELD(NGX,NGY)
      REAL, ALLOCATABLE :: PSEC4(:)

      CHARACTER* 1 HOPER
      CHARACTER*(*) OFILENAME

      LOGICAL LFRST

      DATA LFRST /.TRUE./
! ----------------------------------------------------------------------

      IF(LFRST) THEN 
        CALL PBOPEN(IUOU,OFILENAME,'w',KRET)
        IF(KRET.LT.0) THEN
          IF(KRET.EQ.-1) THEN
            WRITE(IU06,*) '***********************************'
            WRITE(IU06,*) '*  COULD NOT OPEN '//OFILENAME
            WRITE(IU06,*) '***********************************'
            CALL ABORT1
          ENDIF
          IF(KRET.EQ.-2) THEN
            WRITE(IU06,*) '***********************************'
            WRITE(IU06,*) '*  INVALID FILENAME'//OFILENAME
            WRITE(IU06,*) '***********************************'
            CALL ABORT1
          ENDIF
        ENDIF

        LFRST=.FALSE.
      ELSE
!       WE USE CLOSE AND REOPEN TO APPEND RESULTS TO MAKE SURE
!       THAT THE ORDER OF THE INPUT FILES IS MAINTAINED.
        CALL PBOPEN(IUOU,OFILENAME,'a',KRET)
        IF(KRET.LT.0) THEN
          IF(KRET.EQ.-1) THEN
            WRITE(IU06,*) '***********************************'
            WRITE(IU06,*) '*  FOR APPEND '
            WRITE(IU06,*) '*  COULD NOT OPEN '//OFILENAME
            WRITE(IU06,*) '***********************************'
            CALL ABORT1
          ENDIF
          IF(KRET.EQ.-2) THEN
            WRITE(IU06,*) '***********************************'
            WRITE(IU06,*) '*  INVALID FILENAME'//OFILENAME
            WRITE(IU06,*) '***********************************'
            CALL ABORT1
          ENDIF
        ENDIF
      ENDIF


      IF(IRGG.EQ.1 .OR. CLDOMAIN == 'm' ) THEN
         INDEX=1
         KLENP=NTENCODE
      ELSEIF(CLDOMAIN == 's' ) THEN
         INDEX=1
         KLENP=NTENCODE
      ELSEIF(IRGG.EQ.0) THEN
!        SEE COMMENT BELOW ON USING MODEL NORTH AND SOUTH BOUNDARIES.
         INDEX=1
         KLENP=NGY*NGX
      ELSE
         INDEX = (NINT((90. - AMONOP ) / XDELLA))*NGX + 1
         KLENP=NTENCODE
      ENDIF

      IF( CLDOMAIN == 'g' ) THEN
        KSEC1( 3) = IMDLGRBID_G   ! MODEL IDENTIFICATION.
      ELSE
        KSEC1( 3) = IMDLGRBID_M   ! MODEL IDENTIFICATION.
      ENDIF

      IF(KKK.NE.0 .AND. MMM.NE.0 ) THEN
!  ->    IF IT IS A SPECTRAL FIELD !!!!!

!        TRANSFFER FIELD TO ENCODING ARRAY
         ALLOCATE(PSEC4(KLENP))
         PSEC4=0.

         DO J=1,NGY
           JSN=NGY-J+1
           DO I=1,NLONRGG(JSN)
             IF(FIELD(I,J).NE.PSEK3(2)) THEN
               PSEC4(INDEX)=FIELD(I,J)
             ELSE
               PSEC4(INDEX)=0.
             ENDIF
             INDEX=INDEX+1
           ENDDO
         ENDDO

!        MODIFIED GRIB HEADERS

         KSEC1(44) = KKK
         KSEC1(45) = MMM
         KSEC1(46) = NANG
         KSEC1(47) = NFRE

         DO IC=1,NANG
           KSEC1(49+IC) = NINT(KSEC1(48)*TH(IC)*DEG)
         ENDDO
         DO IC=1,NFRE
           KSEC1(49+NANG+IC) = NINT(KSEC1(49)*FR(IC))
         ENDDO

         KSEC2=KSEK2
!        THE NORTH AND SOUTH BOUNDARY WILL BE SET TO THE MODEL ONE
!        WHEN THE GRID IS REGULAR IN ORDER TO REPRODUCE WHAT WAS DONE
!        WHEN THE DATA WERE INTERPOLATED WITH MARS ON THE MODEL DOMAIN.
         IF(IRGG.EQ.0) THEN
           KSEC2( 3) =  NGY
           KSEC2( 4) =  NINT(AMONOP*1000)
           KSEC2( 7) =  NINT(AMOSOP*1000)
         ENDIF

         KSEC3 = KSEK3
         PSEC3 = PSEK3

         KSEC4 = KSEK4
         KSEC4( 1) = KLENP
         KSEC4( 2) = NGRBRESS 
         KSEC4( 6) = 16
         KSEC4(50) =  1
         KSEC4(51) =  1
         KSEC4(53) =  NANG
         KSEC4(55) =  NFRE
         KSEC4(56) =  1
         KSEC4(57) =  2

         DO IC = 1, NANG
           KSEC4(60+IC-1) = TRANSFER (TH(IC), 1)
         ENDDO
         DO IC = 1, NFRE
           KSEC4(60+NANG+IC-1) = TRANSFER (FR(IC), 1) 
         ENDDO

         HOPER=HOPERS

!        LOG SCALING

         DO IJ=1,KLENP
           PSEC4(IJ) = LOG10(PSEC4(IJ)+PPEPS)+ABS(PPREC)
           IF ( PSEC4(IJ) .LE. PPMISS ) PSEC4(IJ)=PSEC3(2)
         ENDDO

      ELSE
!     ALL THE OTHER WAVE PARAMETERS

!        TRANSFER FIELD TO ENCODING ARRAY
         ALLOCATE(PSEC4(KLENP))
         PSEC4=PSEK3(2)

         DO J=1,NGY
           JSN=NGY-J+1
           DO I=1,NLONRGG(JSN)
             PSEC4(INDEX)=FIELD(I,J)
             INDEX=INDEX+1
           ENDDO
         ENDDO

!        MODIFIED GRIB HEADERS

         KSEC2=KSEK2
!        THE NORTH AND SOUTH BOUNDARY WILL BE SET TO THE MODEL ONE
!        WHEN THE GRID IS REGULAR IN ORDER TO REPRODUCE WHAT WAS DONE
!        WHEN THE DATA WERE INTERPOLATED WITH MARS ON THE MODEL DOMAIN.
         IF(IRGG.EQ.0) THEN
           KSEC2( 3) =  NGY
           KSEC2( 4) =  NINT(AMONOP*1000)
           KSEC2( 7) =  NINT(AMOSOP*1000)
         ENDIF

         KSEC3 = KSEK3
         PSEC3 = PSEK3

         KSEC4 = KSEK4
         KSEC4( 1) = KLENP
         KSEC4( 2) = NGRBRESI

         HOPER=HOPERI

       ENDIF

!      ENCODE AND OUTPUT 
!      -----------------

       KLENG = NINT(KLENP*(KSEC4(2)+1.)/(KIND(IU06)*8))+
     &         NKSEK0+NKSEK1+NKSEK2+NKSEK3+NKSEK4+NPSEK2+NPSEK3

       ALLOCATE(KGRIB(KLENG))

       CALL GRIBEX (KSEC0, KSEC1, KSEC2, PSEC2, KSEC3, PSEC3, KSEC4,
     &              PSEC4, KLENP, KGRIB, KLENG , KWORD, HOPER, KRET)
       IF(KRET.GT.0) THEN
         WRITE(IU06,*) '***************************************'
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '*  ENCODING PARAMETER ', KSEC1(6) 
         WRITE(IU06,*) '*  GRIBEX FAILED WITH ERROR CODE ',KRET
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
         WRITE(IU06,*) '***************************************'
         CALL ABORT1
       ENDIF

       DEALLOCATE(PSEC4)

!      WRITE OUT GRIB DATA

       KOUNT=KWORD*KIND(KWORD)
       CALL PBWRITE(IUOU,KGRIB,KOUNT,KRET)
       IF(KRET.LT.0) THEN
         WRITE (IU06,*) '**************************************'
         WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
         WRITE (IU06,*) '*  FILE ',OFILENAME
         WRITE (IU06,*) '*  KRET= ',KRET
         WRITE (IU06,*) '**************************************'
         CALL ABORT1
       ENDIF 

       CALL PBCLOSE(IUOU,KRET)
       IF(KRET.LT.0) THEN
         WRITE (IU06,*) '************************************'
         WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
         WRITE (IU06,*) '* FILE ',OFILENAME
         WRITE (IU06,*) '*  KRET= ',KRET
         WRITE (IU06,*) '************************************'
         CALL ABORT1
       ENDIF

       DEALLOCATE(KGRIB)

      END SUBROUTINE OUTWAMINPUT
