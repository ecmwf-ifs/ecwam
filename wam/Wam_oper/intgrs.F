!----------------------------------------------------------------------

      SUBROUTINE INTGRS(IU06,INFILEN,KSPLIT,CDATE,LINTGRS,OUTFILEN)

!----------------------------------------------------------------------

!**** *INTGRS*  PERFORMS INTERPOLATION ON GRIBBED WAVE MODEL SPECTRA 

!     J. BIDLOT   ECMWF   JANUARY 1998 

!*    PURPOSE.
!     --------

!       PERFORMS INTERPOLATION ON GRIBBED WAVE SPECTRA IF NECESSARY 

!**   INTERFACE.
!     ----------

!       *CALL* *INTGRS*(IU06,INFILEN,KSPLIT,CDATE,LINTGRS,OUTFILEN)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
!      *INFILEN*   CHARACTER FILENAME OF INPUT GRIBBED SPECTRA
!      *KSPLIT*    INTEGER   NUMBER OF SUBDOMAINS USED TO SAVE THE FULL
!                            SPECTRA IN GRIB (1 or 2)
!                  !!!       IF THERE ARE 2 DOMAINS, THE FIRST REQUEST
!                            MUST BE THE NORTHERN HEMISPHERE INCLUDING
!                            THE EQUATOR, AND THE SECOND DOMAIN THE
!                            SOUTHERN HEMISPHERE WIHTOUT THE EQUATOR.
!                            THE WHOLE REQUEST BEING SAVED IN FILENAME.
!      *CDATE*     CHARACTER DATE OF THE REQUESTED FIELD
!      *LINTGRS*   LOGICAL   TRUE IF INTERPOLATION WAS NECESSARY 
!      *INFILEN*   CHARACTER FILENAME OF OUTPUT GRIBBED SPECTRA (WIll ONLY
!                            BE USED IF LINTGRS IS TRUE).

!     METHOD.
!     -------
!     READS THE GRIB HEADER FROM THE INPUT SPECTRA AND FIGURES OUT WHETHER
!     THERE IS A NEED FOR INTERPOLATION. IN THAT CASE, THE INTERPOLATION IS
!     CARRIED OUT USING CALLS TO INTOUT TO DEFINE THE OUTPUT GRID AND 
!     A CALL TO INTF TO PERFORM THE INTERPOLATION.


!     EXTERNALS.
!     ----------

!      *ABORT1*
!      *INCDATE*
!      *INMARSSI*
!      *INTF*
!      *INTOUT*
!      *ISPOFF*
!      *PBCLOSE*
!      *PBOPEN*
!      *PBSIZE*
!      *PBWRITE*

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWMPP   , ONLY : NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY
      USE YOWPCONS , ONLY : R        ,ZMISS

!-----------------------------------------------------------------------

      PARAMETER(NGYMAX=721)
      INTEGER,ALLOCATABLE :: KDOMRGG(:), KPTS(:)
      INTEGER INTV(4), IDUMMY, OUTLEN, OUTGRIBSIZE, OUTLENG
      INTEGER, ALLOCATABLE :: INGRIB(:),OUTGRIB(:)
      REAL, ALLOCATABLE :: ZTH(:),ZFR(:),PLATS(:)
      REAL REALV(4)

      CHARACTER*1 MODE
      CHARACTER*12 CZDATE,CDATE 
      CHARACTER*20 CHARV(4) 
      CHARACTER*40 CLMSG
      CHARACTER*80 INFILEN,OUTFILEN
      LOGICAL LINTGRS, LNEWAREA, LNINT, LLEXIST
!-----------------------------------------------------------------------

!*    0. INITIALISTION 
!     -----------------

      WRITE(IU06,*)'            '
      WRITE(IU06,*)'SUB. INTGRS: GRIBBED SPECTRA WILL BE READ FROM'
      WRITE(IU06,*) INFILEN
      WRITE(IU06,*)'            '

      IPARAM=250

      ALLOCATE(ZTH(NANG),ZFR(NFRE))
      ALLOCATE(KDOMRGG(NGYMAX),KPTS(NGY),PLATS(NGY))

      LINTGRS=.FALSE.
      LNINT=.FALSE.
      LNEWAREA=.FALSE.

!*    1. INPUT OF GRIB DATA 
!     ----------------------

      LFILE=0
      LLEXIST=.FALSE.
      IF (INFILEN.NE. ' ') LFILE=LEN_TRIM(INFILEN)
      INQUIRE(FILE=INFILEN(1:LFILE),EXIST=LLEXIST)
      IF(LLEXIST) THEN
        MODE='r'
        CALL PBOPEN(IUGRS,INFILEN(1:LFILE),MODE,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN     *'
          WRITE (IU06,*) '*   IN INTGRS                          *'
          IF(KRET.EQ.-1)
     &      WRITE (IU06,*) 'COULD NOT OPEN FILE ',INFILEN
          IF(KRET.EQ.-2)
     &      WRITE (IU06,*) 'INVALID FILENAME ',INFILEN
          IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          CALL ABORT1
        ENDIF
      ELSE
        WRITE(IU06,*)'*********************************************'
        WRITE(IU06,*)'*                                           *'
        WRITE(IU06,*)'*SUB. INTGRS : GRIBBED SPECTRA NOT FOUND IN *'
        WRITE(IU06,*)  INFILEN
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

!       GET SIZE OF INGRIB

        CALL PBSIZE (IUGRS, INLENG)
        ISIZE=(INLENG+NPRECI-1)/NPRECI
        ALLOCATE(INGRIB(ISIZE))

!       GET GRIB HEADERS

        CALL INMARSSI (IU06, IUGRS, NGYMAX, NANG, NFRE, ZMISS,
     &                CZDATE, IFORP, KPARAM, KRGG, KDOMRGG, KK, KM,
     &                INGX, INGY, ZMOWEP, ZMOSOP, ZMOEAP, ZMONOP, 
     &                INANG, INFRE, ZTH, ZFR, INLENG, ISIZE, INGRIB,
     &                KRET,CLMSG)
        IF ( KRET .LE. 0 ) THEN
          WRITE(IU06,*)'            '
          WRITE(IU06,*)'SUB. INTGRS : DATE=',CDATE,' IPARAM= ',IPARAM
          WRITE(IU06,*)'INMARSSI PART ',ISPLIT,' DONE'
          WRITE(IU06,*)'            '
        ELSE
          WRITE(IU06,*)'SUB. INTGRS : KRET=',KRET,' IPARAM= ',IPARAM
          WRITE(IU06,*)'INMARSSI PART ',ISPLIT
          WRITE(IU06,*) CLMSG
          WRITE(IU06,*)'            '
          CALL ABORT1
        ENDIF

!*    1.1 ERROR MESSAGES
!     ------------------

        IF (KPARAM.NE.IPARAM) THEN
          WRITE(IU06,*)'********************************'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB INTGRS    *'
          WRITE(IU06,*)'* ===========================  *'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* GRIB PARAMETER  ',KPARAM
          WRITE(IU06,*)'* WAS READ INSTEAD OF ',IPARAM
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'********************************'
          CALL ABORT1
        ENDIF
        CALL INCDATE (CZDATE,IFORP)
        IF (CZDATE.NE.CDATE) THEN
          WRITE(IU06,*)'**********************************'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB INTGRS      *'
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
        IF (INANG .NE. NANG) THEN
          WRITE(IU06,*)'**********************************'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB INTGRS      *'
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
          WRITE(IU06,*)'* FATAL ERROR IN SUB INTGRS      *'
          WRITE(IU06,*)'* ===========================    *'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'* INFRE SHOULD BE EQUAL TO NFRE  *'
          WRITE(IU06,*)'* NFRE  = ',NFRE
          WRITE(IU06,*)'* INFRE  = ',INFRE
          WRITE(IU06,*)'* NO INTERPOLATION ON SPECTRAL   *'
          WRITE(IU06,*)'* FREQUENCIES IS YET IN PLACE !  *'
          WRITE(IU06,*)'*                                *'
          WRITE(IU06,*)'**********************************'
          CALL ABORT1
        ENDIF
        IF (ABS(ZTH(1)-TH(1)).GT.0.01) THEN
          WRITE(IU06,*)'****************************************'
          WRITE(IU06,*)'*                                      *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB INTGRS            *'
          WRITE(IU06,*)'* ===========================          *'
          WRITE(IU06,*)'*                                      *'
          WRITE(IU06,*)'* THE INITIAL DIRECTIONS ARE DIFFERENT *'
          WRITE(IU06,*)'* TH(1)  = ',TH(1)
          WRITE(IU06,*)'* ZTH(1)  = ',ZTH(1)
          WRITE(IU06,*)'* NO INTERPOLATION ON SPECTRAL         *'
          WRITE(IU06,*)'* DIRECTION IS YET IN PLACE !          *'
          WRITE(IU06,*)'*                                      *'
          WRITE(IU06,*)'****************************************'
          CALL ABORT1
        ENDIF
        IF (ABS(ZFR(1)-FR(1)).GT.0.01*FR(1)) THEN
          WRITE(IU06,*)'****************************************'
          WRITE(IU06,*)'*                                      *'
          WRITE(IU06,*)'* FATAL ERROR IN SUB INTGRS            *'
          WRITE(IU06,*)'* ===========================          *'
          WRITE(IU06,*)'*                                      *'
          WRITE(IU06,*)'* THE INITIAL FREQUENCY ARE DIFFERENT  *'
          WRITE(IU06,*)'* FR(1)  = ',FR(1)
          WRITE(IU06,*)'* ZFR(1)  = ',ZFR(1)
          WRITE(IU06,*)'* NO INTERPOLATION ON SPECTRAL         *'
          WRITE(IU06,*)'* FREQUENCIES IS YET IN PLACE !        *'
          WRITE(IU06,*)'*                                      *'
          WRITE(IU06,*)'****************************************'
          CALL ABORT1
        ENDIF

!*    1.2 CHECK IF INTERPOLATION IS NECESSARY 
!     ----------------------------------------

        IF (AMOWEP.NE.ZMOWEP .OR. AMOEAP.NE.ZMOEAP .OR.
     &      SMONOP.NE.ZMONOP .OR. SMOSOP.NE.ZMOSOP) THEN
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  INTERPOLATION IS NECESSARY     *' 
          WRITE(IU06,*) '*  ==========================     *'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  THIS IS NOT THE SAME GRID AREA *'
          WRITE(IU06,*) '*  AS DELARED                     *'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  AMOSOP: ', SMOSOP, 'ZMOSOP : ',ZMOSOP
          WRITE(IU06,*) '*  AMONOP: ', SMONOP, 'ZMONOP : ',ZMONOP
          WRITE(IU06,*) '*  AMOWEP: ', AMOWEP, 'ZMOWEP : ',ZMOWEP
          WRITE(IU06,*) '*  AMOEAP: ', AMOEAP, 'ZMOEAP : ',ZMOEAP
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '***********************************'

          LINTGRS=.TRUE.
          LNEWAREA=.TRUE.
        ENDIF
        IF (INGX .NE. NGX) THEN
          WRITE(IU06,*)'********************************'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'*  INTERPOLATION IS NECESSARY  *' 
          WRITE(IU06,*)'* ===========================  *'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* INGX IS NOT EQUAL TO NGX     *'
          WRITE(IU06,*)'* NGX  = ',NGX
          WRITE(IU06,*)'* INGX  = ',INGX
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'********************************'

          LINTGRS=.TRUE.
        ENDIF
        IF (INGY .NE. KNGY) THEN
          WRITE(IU06,*)'********************************'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'*  INTERPOLATION IS NECESSARY  *' 
          WRITE(IU06,*)'* ===========================  *'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* INGY IS NOT EQUAL TO KNGY    *'
          WRITE(IU06,*)'* KNGY  = ',KNGY
          WRITE(IU06,*)'* INGY  = ',INGY
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'********************************'

          LINTGRS=.TRUE.
        ENDIF
        IF (KRGG.NE.IRGG) THEN
          WRITE(IU06,*)'********************************'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* INTERPOLATION IS NECESSARY   *' 
          WRITE(IU06,*)'* ==========================   *'
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'* KRGG IS NOT EQUAL TO IRGG    *'
          WRITE(IU06,*)'* KRGG  = ',KRGG
          WRITE(IU06,*)'* IRGG  = ',IRGG
          WRITE(IU06,*)'*                              *'
          WRITE(IU06,*)'********************************'

          LINTGRS=.TRUE.
        ENDIF
        IF (KRGG.EQ.1.AND.IRGG.EQ.1) THEN
          IC=0
          DO J=KST,KND
            JSN=NGY-J+1
            IC=IC+1
            IF(KDOMRGG(INGY+1-IC).NE.NLONRGG(JSN).AND..NOT.LINTGRS) THEN
              WRITE(IU06,*)'***********************************'
              WRITE(IU06,*)'*                                 *'
              WRITE(IU06,*)'* INTERPOLATION IS NECESSARY      *' 
              WRITE(IU06,*)'* ==========================      *'
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

              LINTGRS=.TRUE.
            ENDIF
          ENDDO
        ENDIF


!      2.0 INTERPOLATION 
!          -------------
        IF(LINTGRS) THEN

          IF (IRGG.EQ.0) THEN

!         DEFINE THE INTERPOLATION PARAMETERS FOR A REGULAR LAT LON GRID
!         --------------------------------------------------------------

            REALV(1)=XDELLO
            REALV(2)=XDELLA

            IRET = INTOUT('grid', INTV, REALV, CHARV)
            IF ( IRET .NE. 0 ) THEN
              WRITE(IU06,*)'*******************************************'
              WRITE(IU06,*)'*                                         *'
              WRITE(IU06,*)'*INTOUT regular grid CALL FAILED IN INTGRS*'
              WRITE(IU06,*)'*                                         *'
              WRITE(IU06,*)'*******************************************'
              CALL ABORT1
            ENDIF

          ELSE

!         DEFINE THE INTERPOLATION PARAMETERS FOR A REDUCED LAT LON GRID
!         --------------------------------------------------------------

            REALV(1)=XDELLO
            REALV(2)=XDELLA
            IRET = INTOUT('grid', INTV, REALV, CHARV)
            IF ( IRET .NE. 0 ) THEN
              WRITE(IU06,*)'*******************************************'
              WRITE(IU06,*)'*                                         *'
              WRITE(IU06,*)'*INTOUT reduced grid CALL FAILED IN INTGRS*'
              WRITE(IU06,*)'*                                         *'
              WRITE(IU06,*)'*******************************************'
              CALL ABORT1
            ENDIF

            INTV(1) = KNGY
            IRET = INTOUT('user_reduced_latlon', INTV, REALV, CHARV)
            IF ( IRET .NE. 0 ) THEN
              WRITE(IU06,*)'*******************************************'
              WRITE(IU06,*)'*                                         *'
              WRITE(IU06,*)'*INTOUT user_reduced CALL FAILED IN INTGRS*'
              WRITE(IU06,*)'*                                         *'
              WRITE(IU06,*)'*******************************************'
              CALL ABORT1
            ENDIF

            IC=0
            DO J=KST,KND
              JSN=NGY-J+1
              IC=IC+1
              KPTS(IC) = NLONRGG(JSN)
            ENDDO
            IRET = INTOUT('ll_pnts', KPTS, REALV, CHARV)
            IF ( IRET .NE. 0 ) THEN
              WRITE(IU06,*)'****************************************'
              WRITE(IU06,*)'*                                      *'
              WRITE(IU06,*)'* INTOUT ll_pnts CALL FAILED IN INTGRS *'
              WRITE(IU06,*)'*                                      *'
              WRITE(IU06,*)'****************************************'
              CALL ABORT1
            ENDIF

            IC=0
            DO J=KST,KND
              IC=IC+1
              PLATS(IC) = SMONOP-(IC-1)*XDELLA 
            ENDDO
            IRET = INTOUT('ll_lats', INTV, PLATS, CHARV)
            IF ( IRET .NE. 0 ) THEN
              WRITE(IU06,*)'****************************************'
              WRITE(IU06,*)'*                                      *'
              WRITE(IU06,*)'* INTOUT ll_lats CALL FAILED IN INTGRS *'
              WRITE(IU06,*)'*                                      *'
              WRITE(IU06,*)'****************************************'
              CALL ABORT1
            ENDIF

          ENDIF


!         DEFINE AREA RESOLUTION 
!         ----------------------
          IF(LNEWAREA) THEN
            REALV(1) = SMONOP 
            REALV(2) = AMOWEP 
            REALV(3) = SMOSOP 
            REALV(4) = AMOEAP 
          ELSE
            REALV(1) = 0. 
            REALV(2) = 0. 
            REALV(3) = 0. 
            REALV(4) = 0. 
          ENDIF
          IRET = INTOUT('area', INTV, REALV, CHARV)
          IF ( IRET .NE. 0 ) THEN
            WRITE(IU06,*)'*************************************'
            WRITE(IU06,*)'*                                   *'
            WRITE(IU06,*)'* INTOUT area CALL FAILED IN INTGRS *'
            WRITE(IU06,*)'*                                   *'
            WRITE(IU06,*)'*************************************'
            CALL ABORT1
          ENDIF

!         INTERPOLATE
!         -----------

          INLEN=ISIZE

          SCALE=(REAL(NGX)/REAL(INGX))*(REAL(KNGY)/REAL(INGY))
          IF(KRGG.EQ.1.AND.IRGG.EQ.0) SCALE=1.5*SCALE
          OUTGRIBSIZE=NINT(REAL(INLEN)*SCALE)+1
          ALLOCATE(OUTGRIB(OUTGRIBSIZE))
          OUTLEN=OUTGRIBSIZE
 
          IRET = INTF(INGRIB,INLEN,ZDUMMY,OUTGRIB,OUTLEN,ZDUMMY)
          IF ( IRET .NE. 0 ) THEN
            WRITE(IU06,*)'******************************'
            WRITE(IU06,*)'*                            *'
            WRITE(IU06,*)'* INTF CALL FAILED IN INTGRS *'
            WRITE(IU06,*)'* RETURN CODE = ',IRET 
            WRITE(IU06,*)'* INLEN = ', INLEN 
            WRITE(IU06,*)'* OUTGRIBSIZE = ', OUTGRIBSIZE 
            WRITE(IU06,*)'* OUTLEN = ', OUTLEN 
            WRITE(IU06,*)'*                            *'
            WRITE(IU06,*)'*******************************'
            CALL ABORT1
          ENDIF

          IF(IRET.EQ.0.AND.OUTLEN.EQ.0) THEN
            WRITE(IU06,*)'********************************************'
            WRITE(IU06,*)'*                                          *'
            WRITE(IU06,*)'* NO INTERPOLATION WAS NECESSARY           *'
            WRITE(IU06,*)'* GRIB DATA CAN BE FOUND IN THE INPUT FILE *'
            WRITE(IU06,*)'* ',INFILEN
            WRITE(IU06,*)'*                                          *'
            WRITE(IU06,*)'********************************************'
            LNINT=.TRUE.

          ELSE

!           OUTPUT INTERPOLATED DATA FILE
          
            WRITE(IU06,*)'**********************************'
            WRITE(IU06,*)'*                                *'
            WRITE(IU06,*)'* INTERPOLATION WAS CARRIED OUT  *'
            WRITE(IU06,*)'* GRIB DATA CAN BE FOUND IN FILE *'
            WRITE(IU06,*)'* ', OUTFILEN
            WRITE(IU06,*)'*                                *'
            WRITE(IU06,*)'**********************************'

            IF(ISPLIT.EQ.1) THEN
              LFILE=LEN_TRIM(OUTFILEN)
              MODE='w'
              CALL PBOPEN(IUOUT,OUTFILEN(1:LFILE),MODE,KRET)
              IF(KRET.LT.0) THEN
                WRITE (IU06,*) '*************************************'
                WRITE (IU06,*) '*                                   *'
                WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN  *'
                WRITE (IU06,*) '*   IN INTGRS                       *'
                IF(KRET.EQ.-1)
     &            WRITE (IU06,*) ' COULD NOT OPEN FILE ',OUTFILEN
                IF(KRET.EQ.-2)
     &            WRITE (IU06,*) ' INVALID FILENAME ',OUTFILEN
                IF(KRET.EQ.-3) 
     &            WRITE (IU06,*)' INVALID OPEN MODE SPECIFIED'
                WRITE (IU06,*) '*                                   *'
                WRITE (IU06,*) '*************************************'
                CALL ABORT1
              ENDIF
            ENDIF

            OUTLENG=OUTLEN*NPRECI
            CALL PBWRITE(IUOUT,OUTGRIB,OUTLENG,KRET)
            IF(KRET.LT.0) THEN
              WRITE (IU06,*) '**************************************'
              WRITE (IU06,*) '*                                    *'
              WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBWRITE   *'
              WRITE (IU06,*) '*  IN INTGRS                         *'
              WRITE (IU06,*) '*  FILE ',OUTFILEN
              WRITE (IU06,*) '*                                    *'
              WRITE (IU06,*) '**************************************'
              CALL ABORT1
            ENDIF

          ENDIF

        ELSE
          WRITE(IU06,*)'********************************************'
          WRITE(IU06,*)'*                                          *'
          WRITE(IU06,*)'* NO INTERPOLATION WAS NECESSARY           *'
          WRITE(IU06,*)'* GRIB DATA CAN BE FOUND IN THE INPUT FILE *'
          WRITE(IU06,*)'* ',INFILEN
          WRITE(IU06,*)'*                                          *'
          WRITE(IU06,*)'********************************************'
          LNINT=.TRUE.
        ENDIF

        IF(ALLOCATED(INGRIB)) DEALLOCATE(INGRIB)
        IF(ALLOCATED(OUTGRIB)) DEALLOCATE(OUTGRIB)

      ENDDO

      IF(LNINT) LINTGRS=.FALSE.

!      4.0 DEALLOCATION
!          ------------
 
      IF(ALLOCATED(ZFR))DEALLOCATE(ZFR)
      IF(ALLOCATED(ZTH))DEALLOCATE(ZTH)
      IF(ALLOCATED(KDOMRGG))DEALLOCATE(KDOMRGG)
      IF(ALLOCATED(KPTS))DEALLOCATE(KPTS)
      IF(ALLOCATED(PLATS))DEALLOCATE(PLATS)

      CALL PBCLOSE(IUGRS,KRET)
      IF(KRET.LT.0) THEN
        WRITE (IU06,*) '************************************'
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
        WRITE (IU06,*) '* IN INTGRS                        *'
        WRITE (IU06,*) '* FILE ',INFILEN
        WRITE (IU06,*) '*                                  *'
        WRITE (IU06,*) '************************************'
        CALL ABORT1
      ENDIF

      IF(LINTGRS) THEN
        CALL PBCLOSE(IUOUT,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '************************************'
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
          WRITE (IU06,*) '* IN INTGRS                        *'
          WRITE (IU06,*) '* FILE ',OUTFILEN
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '************************************'
          CALL ABORT1
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE INTGRS
