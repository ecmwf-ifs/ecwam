      SUBROUTINE PRESET_WGRIB_TEMPLATE(CT,IGRIB_HANDLE)

!----------------------------------------------------------------------

!**** *PRESET_WGRIB_TEMPLATE* SETS DEFAULT VALUES FOR GRIB TEMPLATES

!     J. BIDLOT    ECMWF JUNE 2009 

!*    PURPOSE.
!     --------

!     GET A GRIB TEMPLATE FROM FILE OR FROM THE IFS
!     AND MODIFY IT ACCORDINGLY.

!**   INTERFACE.
!     ----------

!     SUBROUTINE PRESET_WGRIB_TEMPLATE(CT,IGRIB_HANDLE)
!                INPUT:
!                CT           : "I" for INTEGRATED PARAMETERS AND
!                               "S" for SPECTRA
!                OUTPUT:
!                IGRIB_HANDLE : GRIB HANDLE THAT WILL BE CREATED. 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : 
     &            NTENCODE ,IMDLGRBID_G,IMDLGRBID_M      ,NGRBRESI ,
     &            NGRBRESS, LGRHDIFS ,LNEWLVTP 
      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_IFS
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,IQGAUSS  ,AMOWEP   ,AMOSOP   ,
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NANG     ,NFRE     ,
     &                      CLDOMAIN
      USE YOWPCONS , ONLY : ZMISS    ,DEG
      USE YOWSTAT  , ONLY : NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,
     &            MARSTYPE ,YCLASS   ,YEXPVER  ,ISTREAM  ,NLOCGRB  ,
     &            NCONSENSUS,NDWD    ,NMFR     ,NNCEP    ,NUKM     ,
     &            IREFDATE
      USE YOWTEST  , ONLY : IU06     ,ITEST

      USE GRIB_API_INTERFACE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

      IMPLICIT NONE

      INTEGER :: IC, JC, KST,JSN, KK, MM
      INTEGER :: ICLASS,ICENTRE,IFS_STREAM
      INTEGER :: IGRIB_HANDLE
      INTEGER :: IDIRSCALING, IFRESCALING
      INTEGER :: NY
      INTEGER :: KSYSNB, KMETNB, KREFDATE
      INTEGER :: IDUM, IRET 
      INTEGER :: ITHETA(NANG)
      INTEGER :: IFREQ(NFRE)
      INTEGER, DIMENSION(:), ALLOCATABLE :: PL

      REAL :: ZHOOK_HANDLE
! The following must NOT be changed from a 4 byte real
      REAL*4 REAL4
      REAL :: ZTHETA(NANG)
      REAL :: ZFREQ(NFRE)
      REAL, ALLOCATABLE :: SCFR(:), SCTH(:)

      CHARACTER*1 CT 
      CHARACTER*2 MARSFCTYPE
      CHARACTER*4 CSTREAM
      CHARACTER(LEN=96) :: CLWORD

      LOGICAL :: LASTREAM

!-------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PRESET_WGRIB_TEMPLATE',0,ZHOOK_HANDLE)

      IF(.NOT. LGRHDIFS) THEN
        IGRIB_HANDLE=-99
        CALL IGRIB_NEW_FROM_TEMPLATE(IGRIB_HANDLE,'gg_sfc')
      ELSE
         IGRIB_HANDLE=-99
         CALL IGRIB_CLONE(NGRIB_HANDLE_IFS,IGRIB_HANDLE)
      ENDIF
      IF(ITEST.GT.0) THEN
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'centre', ICENTRE)
        WRITE(IU06,*) ' PRESET_WGRIB_TEMPLATE: GRIB HANDLE CREATED: ',
     &               IGRIB_HANDLE, ' CENTRE = ',ICENTRE 
      ENDIF

!     PRODUCT DEFINITION.
!     -------------------

!     MODEL IDENTIFICATION.
      IF( CLDOMAIN == 'g' ) THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'generatingProcessIdentifier',
     &                       IMDLGRBID_G)
      ELSE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'generatingProcessIdentifier',
     &                       IMDLGRBID_M)
      ENDIF

      IF(.NOT.LNEWLVTP) THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'levtype',102)
      ELSE
        ! if this change check wgribout and intwaminput.
        IF(CT.EQ."S") THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'levtype',212)
        ELSE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'levtype',209)
        ENDIF
      ENDIF
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'level',0)

!     DEFINE YOUR OWN LOCAL HEADER
!     -----------------------------
      IF(.NOT. LGRHDIFS) THEN
        ! LOCAL MARS TABLE USED.

        IF(CT.EQ."S") THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'localDefinitionNumber',
     &                         13)
        ELSE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'localDefinitionNumber',
     &                         NLOCGRB)
        ENDIF

        ! EXPERIMENT VERSION
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'expver',YEXPVER)
        ! CLASS
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'class',YCLASS)
        ! TYPE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'type',2)
        ! STREAM
        IF(ISTREAM.GT.0) THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stream',ISTREAM)
        ELSE
          WRITE(IU06,*) ''
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) ' ERROR IN PRESET_WGRIB_TEMPLATE !!!!! '
          WRITE(IU06,*) ' ISTREAM MUST ALWAYS BE SPECIFIED > 0 !!!' 
          WRITE(IU06,*) ' SEE INPUT NAMELIST wam_input' 
          WRITE(IU06,*) '*******************************************'
          WRITE(IU06,*) ''
          CALL ABORT1
        ENDIF

        ! ENSEMBLE FORECAST NUMBER
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'perturbationNumber',NENSFNB)
        ! TOTAL ENSEMBLE FORECAST NUMBER
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'numberOfForecastsInEnsemble',
     &                       NTOTENS)

        IF(CT.EQ."S") THEN
!         special case for spectra....
!         (for seasonal or monthly forecast only)
          IF ( ISTREAM .EQ. 1082 .OR.
     &         ISTREAM .EQ. 1095 .OR.
     &         ISTREAM .EQ. 1203 .OR.
     &         ISTREAM .EQ. 1204 ) THEN
            KSYSNB=NSYSNB
            KMETNB=NMETNB
          ELSE
            KSYSNB=65535
            KMETNB=65535
          ENDIF
          ! SYSTEM NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'systemNumber',KSYSNB)
          ! METHOD NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'methodNumber',KMETNB)

          IF(ISTREAM == 1204 .OR.
     &       ISTREAM == 1078 .OR.
     &       ISTREAM == 1079 .OR.
     &       ISTREAM == 1084 .OR.
     &       ISTREAM == 1085
     &       ) THEN
            KREFDATE=IREFDATE
          ELSE
            KREFDATE=0 
          ENDIF
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'referenceDate',KREFDATE)

          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'climateDateFrom',0)
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'climateDateTo',0)
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'legBaseDate',0)
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'legBaseTime',0)
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'legNumber',0)
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'oceanAtmosphereCoupling',0)


        ELSE IF(NLOCGRB == 18 ) THEN
          ! MULTI ANALYSIS ENSEMBLE RUNS
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'modelIdentifier','ECMF')
          IF(NCONSENSUS == 0) THEN
            ! Data from one centre
            IF(NDWD == 1) THEN
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',78)
            ELSEIF(NMFR == 1) THEN
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',85)
            ELSEIF(NNCEP == 1) THEN
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',7)
            ELSEIF(NUKM == 1) THEN
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',74)
            ENDIF
          ELSE
            ! Consensus analysis (always includes ECMWF)
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',255)
            ICENTRE=1
            CLWORD='ECMF'
            IF(NDWD == 1) THEN
              ICENTRE=ICENTRE+1
              CLWORD=CLWORD//'EDZW'
            ENDIF
            IF(NMFR == 1) THEN
              ICENTRE=ICENTRE+1
              CLWORD=CLWORD//'LFPW'
            ENDIF
            IF(NNCEP == 1) THEN
             ICENTRE=ICENTRE+1
             CLWORD=CLWORD//'KWBC'
            ENDIF
            IF(NUKM == 1) THEN
             ICENTRE=ICENTRE+1
             CLWORD=CLWORD//'EGRR'
            ENDIF
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'consensusCount',ICENTRE)
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'ccccIdentifiers',
     &                           CLWORD(1:4*ICENTRE))
          ENDIF
        ELSE IF(NLOCGRB == 15 ) THEN
!         SEASONAL FORECASTS
          ! SYSTEM NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'system',NSYSNB)
          ! METHOD NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'method',NMETNB)
        ELSE IF(NLOCGRB == 23 ) THEN
!         Coupled atmospheric, wave and ocean means (with hindcast support)
          ! SYSTEM NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'system',NSYSNB)
          ! METHOD NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'method',NMETNB)
          ! REFERENCE DATE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'refdate',IREFDATE)
        ELSE IF(NLOCGRB == 26 ) THEN
!         EPS or DETERMINISTIC HINDCASTS
          ! REFERENCE DATE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'referenceDate',IREFDATE)
        ENDIF

      ELSE

!     USES GRIB HEADER INFORMATION FROM THE IFS COUPLING
!     --------------------------------------------------

        IF(CT.EQ."S") THEN
!         SPECTRA USE THEIR OWN GRIB TABLE !!!
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'localDefinitionNumber',13)
          CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'systemNumber',
     &                         IDUM,KRET=IRET)
          IF(IRET.NE.0) THEN
             KSYSNB=65535
             CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'systemNumber',KSYSNB)
          ENDIF
          CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'methodNumber',
     &                         IDUM,KRET=IRET)
          IF(IRET.NE.0) THEN
             KMETNB=65535
             CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'methodNumber',KMETNB)
          ENDIF
          CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'referenceDate',
     &                         IDUM,KRET=IRET)
          IF(IRET.NE.0) THEN
             KREFDATE=0 
             CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'referenceDate',KREFDATE)
          ENDIF

          CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'climateDateFrom',
     &                         IDUM,KRET=IRET)
          IF(IRET.NE.0)
     &      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'climateDateFrom',0)

          CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'climateDateTo',
     &                         IDUM,KRET=IRET)
          IF(IRET.NE.0)
     &      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'climateDateTo',0)

          CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'legBaseDate',
     &                         IDUM,KRET=IRET) 
          IF(IRET.NE.0)
     &      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'legBaseDate',0)

          CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'legBaseTime',
     &                         IDUM,KRET=IRET)
          IF(IRET.NE.0)
     &      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'legBaseTime',0)

          CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'legNumber', 
     &                         IDUM,KRET=IRET)
          IF(IRET.NE.0)
     &      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'legNumber',0)

          CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,
     &                         'oceanAtmosphereCoupling',
     &                         IDUM,KRET=IRET)
          IF(IRET.NE.0)
     &      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                           'oceanAtmosphereCoupling',0)
        ENDIF

!       RESET STREAM IF NEEDED
        CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'stream',IFS_STREAM)
        IF(.NOT.LNEWLVTP) THEN
!         GET ISTREAM THAT CORRESPONDS TO IFS_STREAM
          CALL WSTREAM_STRG(IFS_STREAM, CSTREAM, NENSFNB, NTOTENS,
     &                      MARSFCTYPE, ISTREAM, LASTREAM) 
          IF(CSTREAM.EQ.'****') THEN
            WRITE(IU06,*) '*****************************************'
            WRITE(IU06,*) ''
            WRITE(IU06,*) ' ERROR IN PRESET_WGRIB_TEMPLATE !!!!'
            WRITE(IU06,*) ' IFS STREAM UNKNOWN '
            WRITE(IU06,*) ' INPUT ISTREAM = ', IFS_STREAM
            WRITE(IU06,*) ' BUT NOT DEFINED IN WSTREAM_STRG !!!!'
            WRITE(IU06,*) ''
            WRITE(IU06,*) '*****************************************'
            CALL ABORT1
          ENDIF
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stream',ISTREAM)
        ENDIF

      ENDIF


!     SPECIFIC ENTRIES FOR SPECTRAL DATA
      IF(CT.EQ."S") THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'numberOfDirections',NANG)
        IDIRSCALING = 1000
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'directionScalingFactor',
     &                       IDIRSCALING)
        ALLOCATE(SCTH(NANG))
        DO KK=1,NANG
           SCTH(KK)=NINT(TH(KK)*IDIRSCALING*DEG)
        ENDDO
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'scaledDirections',SCTH)
        DEALLOCATE(SCTH)


        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'numberOfFrequencies',NFRE)
        IFRESCALING = 1000000
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'frequencyScalingFactor',
     &                       IFRESCALING)
        ALLOCATE(SCFR(NFRE))
        DO MM=1,NFRE
          SCFR(MM)=NINT(FR(MM)*IFRESCALING)
        ENDDO
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'scaledFrequencies',SCFR)
        DEALLOCATE(SCFR)

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'bitsPerValue',NGRBRESS)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'additionalFlagPresent',1)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'typeOfPacking',
     &                       'grid_simple_matrix')
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'NR',1)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'NC',1)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'NC1',NANG)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'NC2',NFRE)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'physicalFlag1',1)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'physicalFlag2',2)
        DO IC=1,NANG
          REAL4 = TH(IC)*DEG
!!!          ZTHETA(IC)=TRANSFER (REAL4, 1)
          ZTHETA(IC)=REAL4
        ENDDO
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'coefsFirst',ZTHETA)
        DO IC=1,NFRE
          REAL4 = FR(IC)
!!!          ZFREQ(IC)=TRANSFER (REAL4, 1)
          ZFREQ(IC)=REAL4
        ENDDO
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'coefsSecond',ZFREQ)

      ENDIF

!     GEOGRAPHY

      IF (IRGG .EQ. 0) THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'gridType','regular_ll')
      ELSE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'gridType','reduced_ll')
      ENDIF


      ! NUMBER OF POINTS ALONG A MERIDIAN
      IF( CLDOMAIN == 'g' .AND. IQGAUSS.NE.1 ) THEN 
        NY = NINT(180./XDELLA) + 1
      ELSE
        NY = NGY
      ENDIF
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'Nj',NY)

      ! NUMBER OF POINTS PER LATITUDE
      IF (IRGG .EQ. 0) THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'Ni',NGX)
      ELSE
        ALLOCATE(PL(NY))
        PL=0
        IF( CLDOMAIN == 'g' .AND. IQGAUSS.NE.1 ) THEN 
          KST = NINT((90. - AMONOP ) / XDELLA)
        ELSE
          KST = 0
        ENDIF
        DO JC = 1, NGY
          JSN=NGY-JC+1
          PL(JC+KST) = NLONRGG(JSN)
        ENDDO
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'pl',PL)
        DEALLOCATE(PL)
      ENDIF

!     RESOLUTION AND COMPONENT FLAGS
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                     'resolutionAndComponentFlags',128)

      ! LATITUDE OF THE FIRST GRID POINT
      IF( CLDOMAIN == 'g' .AND. IQGAUSS.NE.1 ) THEN 
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                       'latitudeOfFirstGridPointInDegrees',90.)
      ELSE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                       'latitudeOfFirstGridPointInDegrees',AMONOP)
      ENDIF

      ! LONGITUDE OF ORIGIN (WEST -)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                     'longitudeOfFirstGridPointInDegrees',AMOWEP)

      ! LATITUDE OF THE LAST GRID POINT
      IF( CLDOMAIN == 'g' .AND. IQGAUSS.NE.1 ) THEN 
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                       'latitudeOfLastGridPointInDegrees',-90.)
      ELSE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                       'latitudeOfLastGridPointInDegrees',AMOSOP)
      ENDIF

      ! LONGITUDE OF EXTREME POINT (WEST)
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                    'longitudeOfLastGridPointInDegrees',AMOEAP)

      ! LONGITUDE INCREMENT
      IF (IRGG .EQ. 0) THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                     'iDirectionIncrementInDegrees',XDELLO)
      ENDIF

      ! LATITUDE INCREMENT
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,
     &                     'jDirectionIncrementInDegrees',XDELLA)


      ! BITMAP PRESENT:
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'bitmapPresent',1)

      ! MISSING DATA VALUE FOR REAL
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'missingValue',ZMISS)

      ! Number of bits used for each encoded value
      IF(CT.EQ."S") THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'bitsPerValue',NGRBRESS)
      ELSE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'bitsPerValue',NGRBRESI)
      ENDIF

      ! TOTAL NUMBER OF GRID POINTS FOR ENCODING.
      IF(IRGG.EQ.1) THEN
        NTENCODE=0
        DO JC=1,NGY
          NTENCODE=NTENCODE+NLONRGG(JC)
        ENDDO
      ELSE
        NTENCODE=NY*NGX
      ENDIF

      IF (LHOOK) CALL DR_HOOK('PRESET_WGRIB_TEMPLATE',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE PRESET_WGRIB_TEMPLATE 
