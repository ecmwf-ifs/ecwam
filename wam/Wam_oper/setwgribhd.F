
      SUBROUTINE SETWGRIBHD 

!----------------------------------------------------------------------

!**** *SETWGRIBHD* SETS DEFAULT VALUES FOR THE GRIB HEADERS

!     J. BIDLOT    ECMWF  OCTOBER 1998 

!*    PURPOSE.
!     --------

!     SET DEFAULT VALUES FOR THE GRIB HEADERS FOR INTEGER SECTION
!     1,2,3,4 AND REAL SECTION 2 AND 3.

!**   INTERFACE.
!     ----------

!     SUBROUTINE SETWGRIBHD  

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     *CLASS_WGRIB*


!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK0    ,KSEK1,
     &            KSEK2    ,KSEK3    ,KSEK4    ,PSEK2    ,PSEK3,
     &            NTENCODE ,IMDLGRBID_G,IMDLGRBID_M
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP,
     &            AMONOP   ,XDELLA   ,XDELLO
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY,
     &             CLDOMAIN
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : NENSFNB  ,NTOTENS  ,MARSTYPE ,YCLASS,
     &            YEXPVER  ,ISTREAM
      USE YOWTEST  , ONLY : IU06     ,ITEST 


      IMPLICIT NONE

      INTEGER JC, KST, JSN, CLASS_WGRIB, ICLASS

!-------------------------------------------------------------------

      NKSEK1 = 49+NANG+NFRE

      ALLOCATE(KSEK0(NKSEK0))
      ALLOCATE(KSEK1(NKSEK1))
      ALLOCATE(KSEK2(NKSEK2))
      ALLOCATE(KSEK3(NKSEK3))
      ALLOCATE(KSEK4(NKSEK4))

      ALLOCATE(PSEK2(NPSEK2))
      ALLOCATE(PSEK3(NPSEK3))
      
      KSEK0 = 0
      KSEK1 = 0
      KSEK2 = 0
      KSEK3 = 0
      KSEK4 = 0
      PSEK2 = 0.0
      PSEK3 = 0.0


!    *SECTION 1*  INTEGER  PRODUCT DEFINITION BLOCK.
!    -----------------------------------------------
!       WORD        VALUE   CONTENTS.
!       ----        -----   ---------
      KSEK1( 1) =   140     ! VERSION NUMBER OF CODE TABLE 2.
      KSEK1( 2) =    98     ! ORIGINATING CENTRE IDENTIFIER.
      IF( CLDOMAIN == 'g' ) THEN
        KSEK1( 3) = IMDLGRBID_G   ! MODEL IDENTIFICATION.
      ELSE
        KSEK1( 3) = IMDLGRBID_M   ! MODEL IDENTIFICATION.
      ENDIF
      KSEK1( 4) =   255     ! GRID DEFINITION.
      KSEK1( 5) =   192     ! FLAG (WMO CODE TABLE 1)
      KSEK1( 6) =     0     ! PARAMETER IDENTIFIER (WMO CODE TABLE 2)
      KSEK1( 7) =   102     ! TYPE OF LEVEL (WMO CODE TABLE 3)
      KSEK1( 8) =     0     ! TYPE OF LEVEL (WMO CODE TABLE 3)
      KSEK1( 9) =     0     ! VALUE(S) OF LEVEL (WMO CODE TABLE 3)
      KSEK1(10) =     0     ! YEAR OF DATA
      KSEK1(11) =     0     ! MONTH OF DATA
      KSEK1(12) =     0     ! DAY OF DATA
      KSEK1(13) =     0     ! HOUR OF DATA
      KSEK1(14) =     0     ! MINUTE OF DATA
      KSEK1(15) =     1     ! TIME UNIT (WMO CODE TABLE 4)
      KSEK1(16) =     0     ! TIME RANGE ONE (FCST TIME FROM BASE DTG)
      KSEK1(17) =     0     ! TIME RANGE TWO
      KSEK1(18) =    10     ! TIME RANGE FLAG (WMO CODE TABLE 5)
      KSEK1(19) =     0     ! NUMBER AVERAGED.
      KSEK1(20) =     0     ! NUMBER MISSING FROM AVERAGE.
      KSEK1(21) =     0     ! CENTURY OF REFERENCE TIME.
      KSEK1(22) =     0     ! 
      KSEK1(23) =     0     !
      KSEK1(24) =     1     ! INDICATES LOCAL USE OF SECTION 1.
!       25-36                                                  0
      KSEK1(37) =     1     ! INDICATES LABELLING WITH MARS PARAMETERS 
      ICLASS=CLASS_WGRIB(YCLASS,IU06)
      KSEK1(38) = ICLASS    ! INDICATES CLASS
      KSEK1(39) =     2     ! TYPE OF FIELD 1=FG 2=AN 9=FC 10=F 11=PF
      IF(ISTREAM.NE.0) THEN
        KSEK1(40)= ISTREAM  ! STREAM
      ELSE
        IF(MARSTYPE .EQ. 'cf' .OR. MARSTYPE .EQ. 'pf') THEN
          KSEK1(40)=1081      ! STREAM 1081=ENSEMBLE WAVE 
        ELSE
          KSEK1(40)=1045      ! STREAM 1045=WAVE
        ENDIF
      ENDIF
      KSEK1(41) =     0     ! EXPERIMENT VERSION
      DO JC = 1, 4
        KSEK1(41) = KSEK1(41)*256 + ICHAR(YEXPVER(JC:JC))
      ENDDO
      IF(MARSTYPE.EQ.'cf'.OR.MARSTYPE.EQ.'pf') THEN
        KSEK1(42) = NENSFNB ! ENSEMBLE FORECAST NUMBER
        KSEK1(43) = NTOTENS ! TOTAL ENSEMBLE FORECAST NUMBER
      ELSE
        KSEK1(42) = 0
        KSEK1(43) = 0
      ENDIF

!    *SECTION 2*  INTEGER  GRID DEFINITION BLOCK.

!       WORD        VALUE   CONTENTS.
!       ----        -----   ---------

      KSEK2( 1) =       0  ! DATA REPRESENTATION TYPE: WMO CODE TABLE 6
      KSEK2( 2) =     NGX  ! NUMBER OF POINTS ALONG A PARALLEL 
                           ! NUMBER OF POINTS ALONG A MERIDIAN
      IF( CLDOMAIN == 'g' ) THEN 
        KSEK2( 3) =  NINT(180./XDELLA) + 1
      ELSE
        KSEK2( 3) =     NGY
      ENDIF
                           ! LATITUDE OF THE FIRST GRID POINT
      IF( CLDOMAIN == 'g' ) THEN 
        KSEK2( 4) =  90000  
      ELSE
        KSEK2( 4) =  NINT(AMONOP*1000) 
      ENDIF
      KSEK2( 5) =  NINT(AMOWEP*1000) ! LONGITUDE OF ORIGIN (WEST -)
      KSEK2( 6) =     128  ! RESOLUTION FLAG (WMO CODE TABLE 7)
                           ! LATITUDE OF THE LAST GRID POINT
      IF( CLDOMAIN == 'g' ) THEN 
        KSEK2( 7) = -90000  
      ELSE
        KSEK2( 7) =  NINT(AMOSOP*1000) 
      ENDIF
      KSEK2( 8) =  NINT(AMOEAP*1000) ! LONGITUDE OF EXTREME POINT (WEST)
      KSEK2( 9) =  NINT(XDELLO*1000) ! LONGITUDE INCREMENT
      KSEK2(10) =  NINT(XDELLA*1000) ! LATITUDE INCREMENT
      KSEK2(11) =       0  ! SCANNING MODE FLAGS (WMO CODE TABLE 8)
!       12[17]-22                                               0
      KSEK2(17) =    IRGG  ! GRID TYPE.
      IF (IRGG .EQ. 1) THEN
        IF( CLDOMAIN == 'g' ) THEN
          KST = NINT((90. - AMONOP ) / XDELLA)
        ELSE
          KST = 0
        ENDIF
        DO JC = 1, NGY
          JSN=NGY-JC+1
          KSEK2(22+JC+KST) = NLONRGG(JSN)
        ENDDO
      ENDIF

!    *KSEK3*  INTEGER  BIT MAP SECTION.
!    ----------------------------------
      KSEK3(1) = 0          ! BITMAP IS INCLUDED

!    *PSEK3*  REAL BIT MAP SECTION.
!    ------------------------------
      PSEK3(2) = ZMISS      ! MISSING DATA VALUE FOR REAL

!    *SECTION 4*  INTEGER  DATA SECTION.
!    -----------------------------------
                           ! Number of bits used for each encoded value
      KSEK4( 2) = 16 ! 9 or 16 or 32 

!    *TOTAL NUMBER OF GRID POINTS FOR ENCODING.
!    ------------------------------------------

      IF(IRGG.EQ.1) THEN
        NTENCODE=0
        DO JC=1,NGY
          NTENCODE=NTENCODE+NLONRGG(JC)
        ENDDO
      ELSE
        IF (CLDOMAIN == 'g') THEN
          NTENCODE=(NGY+NINT((90. - AMONOP ) / XDELLA)
     &                 +NINT((AMOSOP + 90. ) / XDELLA) )*NGX
        ELSE
          NTENCODE=NGY*NGX
        ENDIF
      ENDIF

!-----------------------------------------------------------------------

      IF (ITEST.GT.3) THEN
        WRITE(IU06,'(" DEFAULT VALUES FOR GRIB HEADERS")') 
        WRITE(IU06,'(" I SECTION 0:", /,80(8x, 8I10,/))') KSEK0
        WRITE(IU06,'(" I SECTION 1:", /,80(8x, 8I10,/))') KSEK1
        WRITE(IU06,'(" I SECTION 2:", /,80(8x, 8I10,/))') KSEK2
        WRITE(IU06,'(" I SECTION 3:", /,80(8x, 8I10,/))') KSEK3
        WRITE(IU06,'(" I SECTION 4:", /,80(8x, 8I10,/))') KSEK4
        WRITE(IU06,'(" R SECTION 2:", /,80(8x, 8F10.5,/))') PSEK2
        WRITE(IU06,'(" R SECTION 3:", /,80(8x, 8F10.5,/))') PSEK3
        WRITE(IU06,'(" NTENCODE : ", 8I10)') NTENCODE 
        CALL FLUSH (IU06)
      ENDIF

      RETURN
      END SUBROUTINE SETWGRIBHD 
