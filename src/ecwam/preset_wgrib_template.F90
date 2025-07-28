! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PRESET_WGRIB_TEMPLATE(CT, IGRIB_HANDLE, LLCREATE, NBITSPERVALUE)

!----------------------------------------------------------------------

!**** *PRESET_WGRIB_TEMPLATE* SETS DEFAULT VALUES FOR GRIB TEMPLATES

!     J. BIDLOT    ECMWF JUNE 2009

!*    PURPOSE.
!     --------

!     GET A GRIB TEMPLATE FROM FILE OR FROM THE IFS
!     AND MODIFY IT ACCORDINGLY.

!**   INTERFACE.
!     ----------

!     SUBROUTINE PRESET_WGRIB_TEMPLATE(CT, IGRIB_HANDLE, LLCREATE, NBITSPERVALUE)
!                INPUT:
!                CT           : "I" for INTEGRATED PARAMETERS AND
!                               "S" for SPECTRA
!                OUTPUT:
!                IGRIB_HANDLE : GRIB HANDLE THAT WILL BE CREATED. 

!                OPTIONAL INPUT:
!                LLCREATE       IF TRUE, FORCE CREATION OF TEMPLATE FROM FILE
!                NBITSPERVALUE  NUMBER OF BITS FOR CODING. 
!                               IF PRESENT, IT WILL OVERRULE NGRBRESI AND NGRBRESS

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUP  , ONLY : LWCOUSAMEGRID
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIBHD, ONLY : LL_GRID_SIMPLE_MATRIX,                      &
     &            NTENCODE ,IMDLGRBID_G,IMDLGRBID_M      ,NGRBRESI ,    &
     &            NGRBRESS, LGRHDIFS ,LNEWLVTP,                         &
     &            NSPEC2TAB, NSPEC2TMPD, NSPEC2TMPP 
      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_IFS
      USE YOWMAP   , ONLY : IRGG     ,IQGAUSS  ,DAMOWEP   ,DAMOSOP   ,  &
     &            DAMOEAP  ,DAMONOP  ,DXDELLA  ,DXDELLO   ,NLONRGG   ,  &
     &            NGX      ,NGY      ,CLDOMAIN
      USE YOWPARAM , ONLY : NANG     ,NFRE_RED
      USE YOWPCONS , ONLY : ZMISS    ,DEG
      USE YOWSTAT  , ONLY : NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,    &
     &            MARSTYPE ,YCLASS   ,YEXPVER  ,ISTREAM  ,NLOCGRB  ,    &
     &            NCONSENSUS,NDWD    ,NMFR     ,NNCEP    ,NUKM     ,    &
     &            IREFDATE
      USE YOWTEST  , ONLY : IU06

      USE YOWGRIB  , ONLY : IGRIB_NEW_FROM_SAMPLES, &
                          & IGRIB_CLONE, &
                          & IGRIB_SET_VALUE, &
                          & IGRIB_GET_VALUE
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "wstream_strg.intfb.h"

      CHARACTER(LEN=1), INTENT(IN) :: CT 
      INTEGER(KIND=JWIM), INTENT(OUT) :: IGRIB_HANDLE
      LOGICAL, INTENT(IN), OPTIONAL :: LLCREATE
      INTEGER(KIND=JWIM) , INTENT(IN), OPTIONAL :: NBITSPERVALUE


      INTEGER(KIND=JWIM) :: IC, JC, KST,JSN, KK, MM
      INTEGER(KIND=JWIM) :: ICLASS,ICENTRE,IFS_STREAM
      INTEGER(KIND=JWIM) :: IREPR, IRESFLAGS
      INTEGER(KIND=JWIM) :: IBITSPERVALUE
      INTEGER(KIND=JWIM) :: IDIRSCALING, IFRESCALING
      INTEGER(KIND=JWIM) :: NJ
      INTEGER(KIND=JWIM) :: KSYSNB, KMETNB, KREFDATE
      INTEGER(KIND=JWIM) :: IDUM, IRET 
      INTEGER(KIND=JWIM) :: IGRIB_HANDLE_IFS
      INTEGER(KIND=JWIM) :: ITHETA(NANG)
      INTEGER(KIND=JWIM) :: IFREQ(NFRE_RED)
      INTEGER(KIND=JWIM), DIMENSION(:), ALLOCATABLE :: PL

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRU) :: RMOEAP
      REAL(KIND=JWRB) :: ZTHETA(NANG)
      REAL(KIND=JWRB) :: ZFREQ(NFRE_RED)
      REAL(KIND=JWRB), ALLOCATABLE :: SCFR(:), SCTH(:)

! The following must NOT be changed from a 4 byte real
      REAL(KIND=4) :: REAL4

      CHARACTER(LEN=2) :: MARSFCTYPE
      CHARACTER(LEN=4) :: CSTREAM
      CHARACTER(LEN=96) :: CLWORD

      LOGICAL :: LASTREAM
      LOGICAL :: LLCRT 

!-------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PRESET_WGRIB_TEMPLATE',0,ZHOOK_HANDLE)

      IF( PRESENT(LLCREATE) ) THEN
        LLCRT = LLCREATE
      ELSE
        LLCRT = .FALSE. 
      ENDIF


      IF( PRESENT(NBITSPERVALUE) ) THEN
        IBITSPERVALUE = NBITSPERVALUE
      ELSEIF(CT.EQ."S") THEN
        IBITSPERVALUE = NGRBRESS
      ELSE
        IBITSPERVALUE = NGRBRESI
      ENDIF

      IF (.NOT. LGRHDIFS .OR. LLCRT) THEN
        CALL IGRIB_NEW_FROM_SAMPLES(IGRIB_HANDLE,'gg_sfc_grib2')
      ELSE
        IGRIB_HANDLE_IFS = NGRIB_HANDLE_IFS
        CALL IGRIB_CLONE(IGRIB_HANDLE_IFS,IGRIB_HANDLE)
      ENDIF

!     PRODUCT DEFINITION.
!     -------------------

      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'discipline',10)


!     MODEL IDENTIFICATION.
      IF ( CLDOMAIN == 'g' ) THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'generatingProcessIdentifier', IMDLGRBID_G)
      ELSE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'generatingProcessIdentifier', IMDLGRBID_M)
      ENDIF

      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'level',0)

!     DEFINE YOUR OWN LOCAL HEADER
!     -----------------------------
      IF (.NOT. LGRHDIFS .OR. LLCRT) THEN
        ! LOCAL MARS TABLE USED.

        IF (CT == "S") THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'tablesVersion', NSPEC2TAB)
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'setLocalDefinition', 1)
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'localDefinitionNumber', 1)
          IF ( NTOTENS > 0 ) THEN
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'productDefinitionTemplateNumber', NSPEC2TMPP)
          ELSE
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'productDefinitionTemplateNumber', NSPEC2TMPD)
          ENDIF
        ELSE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'localDefinitionNumber', NLOCGRB)
        ENDIF

        ! EXPERIMENT VERSION
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'expver',YEXPVER)
        ! CLASS
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'class',YCLASS)
        ! TYPE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'type',2)
        ! STREAM
        IF (ISTREAM > 0) THEN
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

        IF ( NTOTENS > 0 ) THEN
          ! ENSEMBLE FORECAST NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'perturbationNumber',NENSFNB)
          ! TOTAL ENSEMBLE FORECAST NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'numberOfForecastsInEnsemble', NTOTENS)
        ENDIF

        IF (CT == "S") THEN
!         special case for spectra....

            IF ( ISTREAM .EQ. 1082 .OR.                                   &
     &           ISTREAM .EQ. 1095 .OR.                                   &
     &           ISTREAM .EQ. 1203 .OR.                                   &
     &           ISTREAM .EQ. 1204 ) THEN

              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'grib2LocalSectionNumber',30)

              KSYSNB=NSYSNB
              KMETNB=NMETNB
              ! SYSTEM NUMBER
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'systemNumber',KSYSNB)
              ! METHOD NUMBER
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'methodNumber',KMETNB)
            ENDIF

            IF (ISTREAM == 1204 .OR.                                       &
     &          ISTREAM == 1078 .OR.                                       &
     &          ISTREAM == 1079 .OR.                                       &
     &          ISTREAM == 1084 .OR.                                       &
     &          ISTREAM == 1085                                            &
     &          ) THEN

              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'grib2LocalSectionNumber',30)

              KREFDATE=IREFDATE
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'referenceDate',KREFDATE)

              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'climateDateFrom',0)
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'climateDateTo',0)
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'legBaseDate',0)
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'legBaseTime',0)
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'legNumber',0)
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'oceanAtmosphereCoupling',0)
            ENDIF

        ELSEIF (NLOCGRB == 18 ) THEN
          ! MULTI ANALYSIS ENSEMBLE RUNS
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'modelIdentifier','ECMF')
          IF (NCONSENSUS == 0) THEN
            ! Data from one centre
            IF (NDWD == 1) THEN
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',78)
            ELSEIF (NMFR == 1) THEN
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',85)
            ELSEIF (NNCEP == 1) THEN
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',7)
            ELSEIF (NUKM == 1) THEN
              CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',74)
            ENDIF
          ELSE
            ! Consensus analysis (always includes ECMWF)
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'origin',255)
            ICENTRE=1
            CLWORD='ECMF'
            IF (NDWD == 1) THEN
              ICENTRE=ICENTRE+1
              CLWORD=CLWORD//'EDZW'
            ENDIF
            IF (NMFR == 1) THEN
              ICENTRE=ICENTRE+1
              CLWORD=CLWORD//'LFPW'
            ENDIF
            IF (NNCEP == 1) THEN
             ICENTRE=ICENTRE+1
             CLWORD=CLWORD//'KWBC'
            ENDIF
            IF (NUKM == 1) THEN
             ICENTRE=ICENTRE+1
             CLWORD=CLWORD//'EGRR'
            ENDIF
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'consensusCount',ICENTRE)
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'ccccIdentifiers',CLWORD(1:4*ICENTRE))
          ENDIF
        ELSEIF (NLOCGRB == 15 ) THEN
!         SEASONAL FORECASTS
          ! SYSTEM NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'system',NSYSNB)
          ! METHOD NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'method',NMETNB)
        ELSEIF (NLOCGRB == 23 ) THEN
!         Coupled atmospheric, wave and ocean means (with hindcast support)
          ! SYSTEM NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'system',NSYSNB)
          ! METHOD NUMBER
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'method',NMETNB)
          ! REFERENCE DATE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'refdate',IREFDATE)
        ELSEIF (NLOCGRB == 26 ) THEN
!         EPS or DETERMINISTIC HINDCASTS
          ! REFERENCE DATE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'referenceDate',IREFDATE)
        ENDIF

      ELSE

!     USES GRIB HEADER INFORMATION FROM THE IFS COUPLING
!     --------------------------------------------------

        IF (CT == "S") THEN
!         SPECTRA USE THEIR OWN GRIB TABLE !!!
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'tablesVersion', NSPEC2TAB)
          IF ( NTOTENS > 0 ) THEN
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'productDefinitionTemplateNumber', NSPEC2TMPP)
          ELSE
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'productDefinitionTemplateNumber', NSPEC2TMPD)
          ENDIF
        ENDIF

!       RESET STREAM IF NEEDED
        CALL IGRIB_GET_VALUE(IGRIB_HANDLE_IFS,'stream',IFS_STREAM)
        IF (.NOT.LNEWLVTP) THEN
!         GET ISTREAM THAT CORRESPONDS TO IFS_STREAM
          CALL WSTREAM_STRG(IFS_STREAM, CSTREAM, NENSFNB, NTOTENS,       &
     &                      MARSFCTYPE, ISTREAM, LASTREAM) 
          IF (CSTREAM == '****') THEN
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
      IF (CT == "S") THEN

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'numberOfWaveDirections',NANG)
        IDIRSCALING = 2
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'scaleFactorOfWaveDirections',IDIRSCALING)
        ALLOCATE(SCTH(NANG))
        DO KK=1,NANG
          SCTH(KK)=NINT(TH(KK)*10**IDIRSCALING*DEG)
        ENDDO
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'scaledValuesOfWaveDirections',SCTH)
        DEALLOCATE(SCTH)

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'numberOfWaveFrequencies',NFRE_RED)
        IFRESCALING = 6
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'scaleFactorOfWaveFrequencies',IFRESCALING)
        ALLOCATE(SCFR(NFRE_RED))
        DO MM=1,NFRE_RED
          SCFR(MM)=NINT(FR(MM)*10**IFRESCALING)
        ENDDO
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'scaledValuesOfWaveFrequencies',SCFR)
        DEALLOCATE(SCFR)
        
      ENDIF

!     GEOGRAPHY

      IF (.NOT. LGRHDIFS .OR. LLCRT .OR. .NOT. LWCOUSAMEGRID ) THEN

        IF ( IQGAUSS == 1 ) THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'gridType','reduced_gg')
        ELSE
          IF (IRGG == 0) THEN
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'gridType','regular_ll')
          ELSE
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'gridType','reduced_ll')
          ENDIF
        ENDIF


        ! NUMBER OF POINTS ALONG A MERIDIAN
        IF ( CLDOMAIN == 'g' .AND. IQGAUSS /= 1 ) THEN 
          NJ = NINT(180.0_JWRU/DXDELLA) + 1
        ELSE
          NJ = NGY
        ENDIF
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'Nj',NJ)
        IF ( IQGAUSS == 1 ) THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'N',NJ/2)
        ENDIF

        ! NUMBER OF POINTS PER LATITUDE
        IF (IRGG == 0) THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'Ni',NGX)
        ELSE
          ALLOCATE(PL(NJ))
          PL(:)=0
          IF ( CLDOMAIN == 'g' .AND. IQGAUSS /= 1 ) THEN 
            KST = NINT((90.0_JWRU - DAMONOP ) / DXDELLA)
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

!       RESOLUTION AND COMPONENT FLAGS
        IF ( IQGAUSS == 1 ) THEN
          IRESFLAGS=0
        ELSE
          IRESFLAGS=128
        ENDIF
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'resolutionAndComponentFlags',IRESFLAGS)

        ! LATITUDE OF THE FIRST GRID POINT
        IF ( CLDOMAIN == 'g' .AND. IQGAUSS /= 1 ) THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'latitudeOfFirstGridPointInDegrees',90.)
        ELSE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'latitudeOfFirstGridPointInDegrees',DAMONOP)
        ENDIF

        ! LONGITUDE OF ORIGIN (WEST -)
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'longitudeOfFirstGridPointInDegrees',DAMOWEP)

        ! LATITUDE OF THE LAST GRID POINT
        IF ( CLDOMAIN == 'g' .AND. IQGAUSS /= 1 ) THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'latitudeOfLastGridPointInDegrees',-90.)
        ELSE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'latitudeOfLastGridPointInDegrees',DAMOSOP)
        ENDIF

        ! LONGITUDE OF EXTREME POINT (EAST)
        RMOEAP = DAMOEAP
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'longitudeOfLastGridPointInDegrees',RMOEAP)

        ! LONGITUDE INCREMENT
        IF (IRGG == 0) THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'iDirectionIncrementInDegrees',DXDELLO)
        ENDIF

        ! LATITUDE INCREMENT
        IF ( IQGAUSS /= 1 ) THEN 
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'jDirectionIncrementInDegrees',DXDELLA)
        ENDIF

      ENDIF


      ! BITMAP PRESENT:
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'bitmapPresent',1)

      ! MISSING DATA VALUE FOR REAL
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'missingValue',ZMISS)

      ! Number of bits used for each encoded value
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'bitsPerValue', IBITSPERVALUE)

      ! TOTAL NUMBER OF GRID POINTS FOR ENCODING.
      IF (IRGG == 1) THEN
        NTENCODE=0
        DO JC=1,NGY
          NTENCODE=NTENCODE+NLONRGG(JC)
        ENDDO
      ELSE
        NTENCODE=NJ*NGX
      ENDIF

IF (LHOOK) CALL DR_HOOK('PRESET_WGRIB_TEMPLATE',1,ZHOOK_HANDLE)

END SUBROUTINE PRESET_WGRIB_TEMPLATE 
