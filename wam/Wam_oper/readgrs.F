!------------------------------------------------------------------

      SUBROUTINE READGRS(IU06,FILENAME,KSPLIT,CDATE)

!----------------------------------------------------------------------

!**** *READGRS*  READS FROM GRIBBED WAVE MODEL SPECTRA 

!     J. BIDLOT   ECMWF   OCTOBER 1997 

!*    PURPOSE.
!     --------

!       INPUT FROM GRIBBED WAVE SPECTRA 

!**   INTERFACE.
!     ----------

!       *CALL* *READGRS*(IU06,FILENAME,KSPLIT,CDATE)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
!      *FILENAME*  CHARACTER FILENAME OF GRIBBED SPECTRA
!      *KSPLIT*    INTEGER   NUMBER OF SUBDOMAINS USED TO SAVE THE FULL
!                            SPECTRA IN GRIB (1 or 2)
!                  !!!       IF THERE ARE 2 DOMAINS, THE FIRST REQUEST
!                            MUST BE THE NORTHERN HEMISPHERE INCLUDING
!                            THE EQUATOR, AND THE SECOND DOMAIN THE
!                            SOUTHERN HEMISPHERE WIHTOUT THE EQUATOR.
!                            THE WHOLE REQUEST BEING SAVED IN FILENAME.
!      *CDATE*     CHARACTER DATE OF THE REQUESTED FIELD 

!     METHOD.
!     -------

!     USE INMARSS TO READ TO GRIBBED SPECTRA. TRANSFORM THE LOGARITHMIC
!     SPECTRAL COMPONENTS TO NORMAL ONES. PUT THE SPECTRUM INTO BLOCKS
!     USING MAKEBLOS.

!     EXTERNALS.
!     ----------

!      *ABORT1*
!      *INMARSS*
!      *ISPOFF*
!      *MAKEBLOS*
!      *PBCLOSE*
!      *PBOPEN*

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO
      USE YOWPCONS , ONLY : R        ,ZMISS
      USE YOWGRID  , ONLY : NLONRGG  ,IJS      ,IJL
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IRGG     ,AMOWEP   ,
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,ZDELLO
      USE YOWSPEC, ONLY : FL1
 
! -------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.




!-----------------------------------------------------------------------

      INTEGER,ALLOCATABLE :: KDOMRGG(:)
      REAL, ALLOCATABLE :: ZSEC4(:),ZTH(:),ZFR(:)

      CHARACTER*1 MODE
      CHARACTER*12 CZDATE,CDATE 
      CHARACTER*40 CLMSG, TITL*100
      CHARACTER*80 FILENAME
      LOGICAL LLEXIST
!-----------------------------------------------------------------------

!*    0. INITIALISTION 
!     -----------------

      WRITE(IU06,*)'            '
      WRITE(IU06,*)'SUB. READGRS: GRIBBED SPECTRA WILL BE READ FROM'
      WRITE(IU06,*) FILENAME
      WRITE(IU06,*)'            '

      IPARAM=250

      ALLOCATE(ZTH(NANG),ZFR(NFRE))
      ALLOCATE(KDOMRGG(NGY))

      IF (IRGG.EQ.1) THEN
        IDL=0
        DO JGY = 1, NGY 
          IDL = NLONRGG(JGY) + IDL 
        ENDDO
        IDL= IDL * NANG * NFRE 
      ELSE
        IDL = NGY * NGX * NFRE * NANG 
      ENDIF

      ALLOCATE (ZSEC4(IDL))
      ZSEC4=ZMISS


!*    1. INPUT OF GRIB DATA 
!     ----------------------


      LFILE=0
      LLEXIST=.FALSE.
      IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
      INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
      IF(LLEXIST) THEN
        MODE='r'
        CALL PBOPEN(IUGRS,FILENAME(1:LFILE),MODE,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN     *'
          WRITE (IU06,*) '*   IN READGRS                         *'
          IF(KRET.EQ.-1)
     &      WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAME
          IF(KRET.EQ.-2)
     &      WRITE (IU06,*) 'INVALID FILENAME ',FILENAME
          IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          CALL ABORT1
        ENDIF
      ELSE
        WRITE(IU06,*)'*********************************************'
        WRITE(IU06,*)'*                                           *'
        WRITE(IU06,*)'*SUB. READGRS: GRIBBED SPECTRA NOT FOUND IN *'
        WRITE(IU06,*)  FILENAME
        WRITE(IU06,*)'*PROGRAM WILL ABORT                         *'
        WRITE(IU06,*)'*                                           *'
        WRITE(IU06,*)'*********************************************'
        CALL ABORT1
      ENDIF


      DO ISPLIT=1,KSPLIT

        ZLATITUDE=0.0

        IF(ISPLIT.EQ.1) THEN
          ISTRT=1
          IF(KSPLIT.EQ.1) THEN
            SMOSOP=AMOSOP
            SMONOP=AMONOP
            KNGY=NGY
            KST=1
            KND=NGY
          ELSE
            SMOSOP=ZLATITUDE
            SMONOP=AMONOP
            KNGY=NGY/2+1
            KST=1
            KND=NGY/2+1
          ENDIF

        ELSE
          ZLONGITUDE=AMOEAP
          IY = NINT((AMONOP-ZLATITUDE)/XDELLA) + 1
          ZX = MOD(ZLONGITUDE-AMOWEP+720.,360.)
          ZX = ZX / ZDELLO(NGY-IY+1)
          IX = NINT(ZX) +1
          ITSP = ISPOFF(NLONRGG, NGY, NANG, NFRE, IY-1, IX) +
     &           (NANG*NFRE)
          ISTRT = 1+ITSP

          SMOSOP=AMOSOP
          SMONOP=ZLATITUDE-XDELLA
          KNGY=NGY/2
          KST=NGY/2+2
          KND=NGY
        ENDIF


        CALL INMARSS (IU06, IUGRS, NGY, IDL-ISTRT+1, NANG, NFRE,
     &               ZMISS, CZDATE, IFORP, KPARAM, KRGG, KDOMRGG,
     &               INGX, INGY, ZMOWEP, ZMOSOP, ZMOEAP, ZMONOP, 
     &               INANG, INFRE, ZSEC4(ISTRT), ZTH, ZFR, KRET, CLMSG)
        IF ( KRET .LE. 0 ) THEN
          WRITE(IU06,*)'            '
          WRITE(IU06,*)'SUB. READGRS: DATE=',CDATE,' IPARAM= ',IPARAM
          WRITE(IU06,*)'INMARSS PART ',ISPLIT,' DONE'
          WRITE(IU06,*)'            '
        ELSE
          WRITE(IU06,*)'SUB. READGRS: KRET=',KRET,' IPARAM= ',IPARAM
          WRITE(IU06,*)'INMARSS PART ',ISPLIT
          WRITE(IU06,*) CLMSG
          WRITE(IU06,*)'            '
        ENDIF


!*    1.1 ERROR MESSAGES
!     ------------------


        IF (KPARAM.NE.IPARAM) THEN
          WRITE(IU06,*)'********************************'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB READGRS   *'
          WRITE(IU06,*)'* ===========================  *'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* GRIB PARAMETER  ',KPARAM
          WRITE(IU06,*)'* WAS READ INSTEAD OF ',IPARAM
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'********************************'
          CALL ABORT1
        ENDIF
        IF (CZDATE.NE.CDATE) THEN
          WRITE(IU06,*)'**********************************'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB READGRS     *'
          WRITE(IU06,*)'* ===========================    *'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'* REQUESTED DATE IS NOT EQUAL TO *' 
          WRITE(IU06,*)'* RETRIEVED DATE.                *' 
          WRITE(IU06,*)'* CZDATE = ',CZDATE
          WRITE(IU06,*)'* CDATE = ',CDATE
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'**********************************'
          CALL ABORT1
        ENDIF
        IF (AMOWEP.NE.ZMOWEP .OR. AMOEAP.NE.ZMOEAP .OR.
     &      SMONOP.NE.ZMONOP .OR. SMOSOP.NE.ZMOSOP) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. READGRS    *'
          WRITE(IU06,*) '*  ============================   *'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  THIS IS NOT THE SAME GRID AREA *'
          WRITE(IU06,*) '*  AS YOU HAD DELARED IN THE      *'
          WRITE(IU06,*) '*  COURSE MODEL GRID              *'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  AMOSOP: ', SMOSOP, 'ZMOSOP : ',ZMOSOP
          WRITE(IU06,*) '*  AMONOP: ', SMONOP, 'ZMONOP : ',ZMONOP
          WRITE(IU06,*) '*  AMOWEP: ', AMOWEP, 'ZMOWEP : ',ZMOWEP
          WRITE(IU06,*) '*  AMOEAP: ', AMOEAP, 'ZMOEAP : ',ZMOEAP
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*     THE PROGRAM ABORTS          *'
          WRITE(IU06,*) '***********************************'
          CALL ABORT1
        ENDIF
        IF (INGX .NE. NGX) THEN
          WRITE(IU06,*)'********************************'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB READGRS   *'
          WRITE(IU06,*)'* ===========================  *'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* INGX SHOULD BE EQUAL TO NGX  *'
          WRITE(IU06,*)'* NGX  = ',NGX
          WRITE(IU06,*)'* INGX  = ',INGX
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'********************************'
          CALL ABORT1
        ENDIF
        IF (INGY .NE. KNGY) THEN
          WRITE(IU06,*)'********************************'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB READGRS   *'
          WRITE(IU06,*)'* ===========================  *'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* INGY SHOULD BE EQUAL TO KNGY *'
          WRITE(IU06,*)'* KNGY  = ',KNGY
          WRITE(IU06,*)'* INGY  = ',INGY
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'********************************'
          CALL ABORT1
        ENDIF
        IF (KRGG.NE.IRGG) THEN
          WRITE(IU06,*)'********************************'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* ERROR IN SUB READGRS         *'
          WRITE(IU06,*)'* =====================        *'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* KRGG IS NOT EQUAL TO IRGG    *'
          WRITE(IU06,*)'* KRGG  = ',KRGG
          WRITE(IU06,*)'* IRGG  = ',IRGG
          WRITE(IU06,*)'* PROGRAM CANNOT INTERPOLATE   *'
          WRITE(IU06,*)'* BETWEEN THOSE 2 GRIDS        *'
          WRITE(IU06,*)'* IT WILL ABORT                *'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'********************************'
          CALL ABORT1
        ENDIF
        IF (KRGG.EQ.1) THEN
          IC=0
          DO J=KST,KND
            JSN=NGY-J+1
            IC=IC+1
            IF(KDOMRGG(INGY+1-IC).NE.NLONRGG(JSN)) THEN
              WRITE(IU06,*)'***********************************'
              WRITE(IU06,*)'*                                 *'
              WRITE(IU06,*)'* ERROR IN SUB READGRS            *'
              WRITE(IU06,*)'* =====================           *'
              WRITE(IU06,*)'*                                 *'
              WRITE(IU06,*)'* THE IRREGULAR GRID DEFINITIONS  *'
              WRITE(IU06,*)'* OF THE INPUT FILE DO NOT MATCH  *'
              WRITE(IU06,*)'* WITH THE MODEL ONE.             *'
              WRITE(IU06,*)'* KDOMRGG IS NOT EQUAL TO NLONRGG *'
              IIC=0
              DO JJ=KST,KND
                JSN=NGY-JJ+1
                IIC=IIC+1
                WRITE(IU06,*)'* JSN = ',JSN
                WRITE(IU06,*)'* KDOMRGG  = ',KDOMRGG(INGY+1-IIC)
                WRITE(IU06,*)'* NLONRGG  = ',NLONRGG(JSN)
              ENDDO
              WRITE(IU06,*)'*                                 *'
              WRITE(IU06,*)'***********************************'
              CALL ABORT1
            ENDIF
          ENDDO
        ENDIF
        IF (INANG .NE. NANG) THEN
          WRITE(IU06,*)'**********************************'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB READGRS     *'
          WRITE(IU06,*)'* ===========================    *'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'* INANG SHOULD BE EQUAL TO NANG  *'
          WRITE(IU06,*)'* NANG  = ',NANG
          WRITE(IU06,*)'* INANG  = ',INANG
          WRITE(IU06,*)'* NO INTERPOLATION ON SPECTRAL   *'
          WRITE(IU06,*)'* DIRECTIONS IS YET IN PLACE !   *'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'**********************************'
          CALL ABORT1
        ENDIF
        IF (INFRE .NE. NFRE) THEN
          WRITE(IU06,*)'**********************************'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB READGRS     *'
          WRITE(IU06,*)'* ===========================    *'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'* INFRE SHOULD BE EQUAL TO NFRE  *'
          WRITE(IU06,*)'* NFRE  = ',NFRE
          WRITE(IU06,*)'* INFRE  = ',INFRE
          WRITE(IU06,*)'* NO INTERPOLATION ON SPECTRAL   *'
          WRITE(IU06,*)'* FREQUECIES IS YET IN PLACE !   *'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'**********************************'
          CALL ABORT1
        ENDIF
!
!!!!!  not yet in place
!
!       IF (INT(100*ZTH(1)).NE.INT(100*TH(1))) THEN
!          WRITE(IU06,*)'****************************************'
!          WRITE(IU06,*)'*                                      *'
!          WRITE(IU06,*)'* FATAL ERROR IN SUB READGRS           *'
!          WRITE(IU06,*)'* ===========================          *'
!          WRITE(IU06,*)'*                                      *'
!          WRITE(IU06,*)'* THE INITIAL DIRECTIONS ARE DIFFERENT *'
!          WRITE(IU06,*)'* TH(1)  = ',TH(1)
!          WRITE(IU06,*)'* ZTH(1)  = ',ZTH(1)
!          WRITE(IU06,*)'* NO INTERPOLATION ON SPECTRAL         *'
!          WRITE(IU06,*)'* DIRECTION IS YET IN PLACE !          *'
!          WRITE(IU06,*)'*                                      *'
!          WRITE(IU06,*)'****************************************'
!          CALL ABORT1
!       ENDIF
!       IF (INT(1000*ZFR(1)).NE.INT(1000*FR(1))) THEN
!          WRITE(IU06,*)'****************************************'
!          WRITE(IU06,*)'*                                      *'
!          WRITE(IU06,*)'* FATAL ERROR IN SUB READGRS           *'
!          WRITE(IU06,*)'* ===========================          *'
!          WRITE(IU06,*)'*                                      *'
!          WRITE(IU06,*)'* THE INITIAL FREQUENCY ARE DIFFERENT  *'
!          WRITE(IU06,*)'* FR(1)  = ',FR(1)
!          WRITE(IU06,*)'* ZFR(1)  = ',ZFR(1)
!          WRITE(IU06,*)'* NO INTERPOLATION ON SPECTRAL         *'
!          WRITE(IU06,*)'* FREQUECIES IS YET IN PLACE !         *'
!          WRITE(IU06,*)'*                                      *'
!          WRITE(IU06,*)'****************************************'
!          CALL ABORT1
!       ENDIF

      ENDDO

!     2.0 TRANSFORM SPECTRAL COMPONENT FROM THE LOGARITHMIC REPRESENTATION
!         TO ITS NORMAL ONE.
!         ----------------------------------------------------------------

      pprec=0.0
      xeps=1.0e-10

      DO JX=1,IDL
        IF (ZSEC4(JX) .NE. ZMISS) THEN
          ZSEC4(JX) = 10.**(ZSEC4(JX)-abs(pprec))- xeps
        ELSE
          ZSEC4(JX) = 0.0
        ENDIF
      ENDDO


!      3.0 TRANSFORM GRID DATA TO BLOCK DATA
!          ---------------------------------

      IF(.NOT.ALLOCATED(FL1)) ALLOCATE (FL1(0:NIBLO,NANG,NFRE))
      FL1=0.
      DO IG=1, NBLO
        CALL MAKEBLOS (FL1, ZSEC4, IDL,
     &   IJS(IG), IJL(IG), IG,
     &   NLONRGG, NGX, NGY, NANG, NFRE, IXLG, KXLT,
     &   NBLO, NIBLO)
      ENDDO

!      4.0 DEALLOCATION
!          ------------

      IF(ALLOCATED(ZFR))DEALLOCATE(ZFR)
      IF(ALLOCATED(ZTH))DEALLOCATE(ZTH)
      IF(ALLOCATED(KDOMRGG))DEALLOCATE(KDOMRGG)
      IF(ALLOCATED(ZSEC4))DEALLOCATE(ZSEC4)

      CALL PBCLOSE(IUGRS,KRET)
      IF(KRET.LT.0) THEN
        WRITE (IU06,*) '************************************'
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (IU06,*) '* IN READGRS                       *'
        WRITE (IU06,*) '* FILE ',FILENAME
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '************************************'
        CALL ABORT1
      ENDIF

      RETURN
      END SUBROUTINE READGRS
