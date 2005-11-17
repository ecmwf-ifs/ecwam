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
     &            NTENCODE ,IMDLGRBID_G,IMDLGRBID_M      ,NGRBRESI ,
     &            IFS_KSEC1,LGRHDIFS ,LNEWLVTP 
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,IQGAUSS  ,AMOWEP   ,AMOSOP   ,
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY,
     &             CLDOMAIN
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,
     &            MARSTYPE ,YCLASS   ,YEXPVER  ,ISTREAM  ,NLOCGRB  ,
     &            NCONSENSUS,NDWD    ,NMFR     ,NNCEP    ,NUKM     ,
     &            IREFDATE
      USE YOWTEST  , ONLY : IU06     ,ITEST 

      IMPLICIT NONE

      INTEGER IC, JC, KST, JSN, CLASS_WGRIB, ICLASS, ICENTRE, IFS_STREAM

      CHARACTER*2 MARSFCTYPE
      CHARACTER*4 CSTREAM
      CHARACTER*8 CLWORD

      LOGICAL :: LASTREAM

!-------------------------------------------------------------------

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
      IF(.NOT.LNEWLVTP) THEN
        KSEK1( 7) =   102     ! TYPE OF LEVEL (WMO CODE TABLE 3)
      ELSE
        KSEK1( 7) =   209   ! if this change check fld2fdb, grb2wgrd  
      ENDIF                 ! spec2fdb and intwaminput.
                            ! for spectra it will be reset to 212
      KSEK1( 8) =     0     ! TYPE OF LEVEL (WMO CODE TABLE 3)
      KSEK1( 9) =     0     ! VALUE(S) OF LEVEL (WMO CODE TABLE 3)
      KSEK1(10) =     0     ! YEAR OF DATA
      KSEK1(11) =     0     ! MONTH OF DATA
      KSEK1(12) =     0     ! DAY OF DATA
      KSEK1(13) =     0     ! HOUR OF DATA
      KSEK1(14) =     0     ! MINUTE OF DATA
      KSEK1(15) =     1     ! TIME UNIT (WMO CODE TABLE 4)
                            ! 1 = HOURS (if the number of hours > 2**16-1
                            ! MULTIPLE OF HOURS OR DAYS WILL BE USED.
      KSEK1(16) =     0     ! TIME RANGE ONE (FCST TIME FROM BASE DTG)
      KSEK1(17) =     0     ! TIME RANGE TWO
      KSEK1(18) =    10     ! TIME RANGE FLAG (WMO CODE TABLE 5)
      KSEK1(19) =     0     ! NUMBER AVERAGED.
      KSEK1(20) =     0     ! NUMBER MISSING FROM AVERAGE.
      KSEK1(21) =     0     ! CENTURY OF REFERENCE TIME.
      KSEK1(22) =     0     ! 
      KSEK1(23) =     0     !
      KSEK1(24) =     1     ! INDICATES LOCAL USE OF SECTION 1.
!         GET ISTREAM THAT CORRESPONDS TO IFS_STREAM
!       25-36                                                  0

!     DEFINE YOUR OWN LOCAL HEADER
      IF(.NOT. LGRHDIFS) THEN
      KSEK1(37) = NLOCGRB ! INDICATES LOCAL MARS TABLE USED.

        ICLASS=CLASS_WGRIB(YCLASS,IU06)
        KSEK1(38) = ICLASS    ! INDICATES CLASS
        KSEK1(39) =     2     ! TYPE OF FIELD 1=FG 2=AN 9=FC 10=F 11=PF

        IF(ISTREAM.GT.0) THEN
          KSEK1(40)= ISTREAM  ! STREAM
        ELSE
          WRITE(IU06,*) ''
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) ' ERROR IN SETWGRIBHD !!!!! '
          WRITE(IU06,*) ' ISTREAM MUST ALWAYS BE SPECIFIED > 0 !!!' 
          WRITE(IU06,*) ' SEE INPUT NAMELIST wam_input' 
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) ''
          CALL ABORT1
        ENDIF

        KSEK1(41) =     0     ! EXPERIMENT VERSION
        DO JC = 1, 4
          KSEK1(41) = KSEK1(41)*256 + ICHAR(YEXPVER(JC:JC))
        ENDDO

        IF(NLOCGRB == 18 ) THEN
          ! MULTI ANALYSIS ENSEMBLE RUNS
          KSEK1(42) = NENSFNB ! ENSEMBLE FORECAST NUMBER
          KSEK1(43) = NTOTENS ! TOTAL ENSEMBLE FORECAST NUMBER
          CLWORD='    ECMF'
          READ(CLWORD,'(A8)') KSEK1(45)
          IF(NCONSENSUS == 0) THEN
          !   Data from one centre (check also fld2fdb)
            KSEK1(46)=0
            IF(NDWD == 1) THEN
              KSEK1(44)=78
            ELSEIF(NMFR == 1) THEN
              KSEK1(44)=85
            ELSEIF(NNCEP == 1) THEN
              KSEK1(44)=7
            ELSEIF(NUKM == 1) THEN
              KSEK1(44)=74
            ENDIF
          ELSE
          ! Consensus analysis (always includes ECMWF)
            KSEK1(44)=255
            ICENTRE=1
            CLWORD='    ECMF'
            READ(CLWORD,'(A8)') KSEK1(46+ICENTRE)
            IF(NDWD == 1) THEN
              ICENTRE=ICENTRE+1
              CLWORD='    EDZW'
              READ(CLWORD,'(A8)') KSEK1(46+ICENTRE)
            ENDIF
            IF(NMFR == 1) THEN  
              ICENTRE=ICENTRE+1
              CLWORD='    LFPW'
              READ(CLWORD,'(A8)') KSEK1(46+ICENTRE)
            ENDIF
            IF(NNCEP == 1) THEN 
              ICENTRE=ICENTRE+1
              CLWORD='    KWBC'
              READ(CLWORD,'(A8)') KSEK1(46+ICENTRE)
            ENDIF     
            IF(NUKM == 1) THEN
              ICENTRE=ICENTRE+1
              CLWORD='    EGRR'
              READ(CLWORD,'(A8)') KSEK1(46+ICENTRE)
            ENDIF 
            KSEK1(46)=ICENTRE
          ENDIF
        ELSE IF(NLOCGRB == 15 .OR.
     &          NLOCGRB == 23 ) THEN
!         seasonal and monthly forecasts
          KSEK1(42) = NENSFNB ! ENSEMBLE FORECAST NUMBER
          KSEK1(43) = NTOTENS ! TOTAL ENSEMBLE FORECAST NUMBER
                              ! SHOULD BE 0 EXCEPT FOR ENSEMBLE RUNS ! 
          KSEK1(44) = NSYSNB ! SYSTEM NUMBER
          KSEK1(45) = NMETNB ! METHOD NUMBER

!         monthly forecast hindcasts need a reference date
          IF(KSEK1(40) == 1204 ) THEN
             KSEK1(49) = IREFDATE 
          ELSE
            KSEK1(49) = 0
          ENDIF
        ELSE IF(NLOCGRB == 26 ) THEN
!         EPS or DETERMINISTIC HINDCASTS
          KSEK1(42) = NENSFNB ! ENSEMBLE FORECAST NUMBER
          KSEK1(43) = NTOTENS ! TOTAL ENSEMBLE FORECAST NUMBER
                              ! SHOULD BE 0 EXCEPT FOR ENSEMBLE RUNS ! 
          KSEK1(44) = IREFDATE ! HINDCASTS NEED A REFERENCE DATE
          KSEK1(45) = 0
        ELSE
!         DETERMINISTIC OR EPS  
          KSEK1(42) = NENSFNB ! ENSEMBLE FORECAST NUMBER
          KSEK1(43) = NTOTENS ! TOTAL ENSEMBLE FORECAST NUMBER
                              ! SHOULD BE 0 EXCEPT FOR ENSEMBLE RUNS ! 
          KSEK1(44) = 0
          KSEK1(45) = 0
        ENDIF

      ELSE
!     USES GRIB HEADER INFORMATION FROM THE IFS COUPLING
        DO IC=24,NKSEK1
          KSEK1(IC) = IFS_KSEC1(IC)
        ENDDO
        IFS_STREAM=IFS_KSEC1(40)
        IF(.NOT.LNEWLVTP) THEN
!         GET ISTREAM THAT CORRESPONDS TO IFS_STREAM
          CALL WSTREAM_STRG(IFS_STREAM, CSTREAM, NENSFNB, NTOTENS,
     &                      MARSFCTYPE, ISTREAM, LASTREAM) 
          IF(CSTREAM.EQ.'****') THEN
            WRITE(IU06,*) '*****************************************'
            WRITE(IU06,*) ''
            WRITE(IU06,*) ' ERROR IN SETWGRIBHD !!!!'
            WRITE(IU06,*) ' IFS STREAM UNKNOWN '
            WRITE(IU06,*) ' INPUT ISTREAM = ', IFS_STREAM
            WRITE(IU06,*) ' BUT NOT DEFINED IN WSTREAM_STRG !!!!'
            WRITE(IU06,*) ''
            WRITE(IU06,*) '*****************************************'
            CALL ABORT1
          ENDIF
          KSEK1(40)= ISTREAM
        ELSE
          ISTREAM=IFS_STREAM
        ENDIF
      ENDIF
!       WORD        VALUE   CONTENTS.
!       ----        -----   ---------

      KSEK2( 1) =       0  ! DATA REPRESENTATION TYPE: WMO CODE TABLE 6
      KSEK2( 2) =     NGX  ! NUMBER OF POINTS ALONG A PARALLEL 
                           ! NUMBER OF POINTS ALONG A MERIDIAN
      IF( CLDOMAIN == 'g' .AND. IQGAUSS.NE.1 ) THEN 
        KSEK2( 3) =  NINT(180./XDELLA) + 1
      ELSE
        KSEK2( 3) =     NGY
      ENDIF
                           ! LATITUDE OF THE FIRST GRID POINT
      IF( CLDOMAIN == 'g' .AND. IQGAUSS.NE.1 ) THEN 
        KSEK2( 4) =  90000  
      ELSE
        KSEK2( 4) =  NINT(AMONOP*1000) 
      ENDIF
      KSEK2( 5) =  NINT(AMOWEP*1000) ! LONGITUDE OF ORIGIN (WEST -)
      IF (IRGG .EQ. 0) THEN
        KSEK2( 6) =     128  ! RESOLUTION FLAG (WMO CODE TABLE 7)
      ELSE
!!!!! until there is a fix in gribex, we will still code with
!!!!! 128 !!!!!
        KSEK2( 6) =     128 
      ENDIF
                           ! LATITUDE OF THE LAST GRID POINT
      IF( CLDOMAIN == 'g' .AND. IQGAUSS.NE.1 ) THEN 
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
        IF( CLDOMAIN == 'g' .AND. IQGAUSS.NE.1 ) THEN 
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
      KSEK4( 2) = NGRBRESI

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
