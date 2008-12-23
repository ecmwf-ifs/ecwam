      SUBROUTINE BUILDSTRESS(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,
     &                       RADOLD, RZIDLOLD, KICEMASK,
     &                       KEND)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF   APRIL 1998 

!     J. BIDLOT    ECMWF   FEBRUARY 1999 TAUT --> SQRT(TAUT)

!     S. ABDALLA   ECMWF   OCTOBER 1999 MODIFICATION THE CALL TO GETWND

!*    PURPOSE.
!     --------
!     CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND CD.

!**   INTERFACE.
!     ----------
!     CALL *BUILDSTRESS*(RU10OLD,RTHWOLD,RUSOLD,RTAUW,RZ0OLD,RADOLD,
!                                                      RZIDLOLD,KEND)
!     *RU10OLD*   WIND SPEED.
!     *RTHWOLD*   WIND DIRECTION (RADIANS).
!     *RUSOLD*    FRICTION VELOCITY.
!     *RTAUW*     WAVE STRESS.
!     *RZ0OLD*    ROUGHNESS LENGTH IN M.
!     *RAD0OLD*   AIR DENSITY IN KG/M3.
!     *RZIDL0OLD* Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!     *KEND*      USED UPPER LIMIT ON FIRST DIMENSION OF INPUT ARRAYS

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     *ABORT1*
!     *AIRSEA*
!     *GETWND*
!     *PBOPEN*
!     *PBREAD*
!     *PBCLOSE*
!     *READWGRIB*

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU    ,ALPHA    ,XKAPPA   ,XNLEV
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWPARAM , ONLY : NBLO     ,NIBLO
      USE YOWPCONS , ONLY : G        ,ROAIR    ,EPSUS    ,EPSU10
      USE YOWGRIBHD, ONLY : NKSEK1 
      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO   ,NWAM_BLKS
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,JPLEVT   ,TAUT
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL
      USE YOWMESPAS, ONLY : LMESSPASS,LNOCDIN  ,LWAVEWIND

! ----------------------------------------------------------------------
      LOGICAL LWNDFILE, LCLOSEWND
      REAL, ALLOCATABLE :: CD(:,:), TEMPTAUT(:,:,:)
      CHARACTER*24 FILNM
      REAL,DIMENSION(NIBLO,NBLO) :: RU10OLD,RTHWOLD,RUSOLD,RZ0OLD,RTAUW
     &                             ,RADOLD, RZIDLOLD
      INTEGER,DIMENSION(NIBLO,NBLO) :: KICEMASK 
      LOGICAL :: LLONLYPOS
! ----------------------------------------------------------------------
! "Dummy" arrays and parameters for compatability purposes:
      INTEGER NFIELDS, NGPTOTG, NCA, NRA
      PARAMETER (NFIELDS=1, NGPTOTG=1, NCA=1, NRA=1)
      INTEGER ISEC1(NKSEK1), ISEC2(22+NRA)
      REAL FIELDS(1)
      INTEGER JKGLO,KIJS,KIJL,NPROMA
      INTEGER MASK_IN(NGPTOTG)
! ----------------------------------------------------------------------

      CDATEWO = ' '
      CDAWIFL = ' '
      CDATEFL = ' '
      CDTPRO = CDATEA

!     GET ATMOSPHERIC MODEL 10m WINDS

      ILEN=1
      LWNDFILE=.TRUE.
      LCLOSEWND=.TRUE.

      DO IG=1,NBLO
        NPROMA=NIBLO/NWAM_BLKS+1
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=1,NIBLO,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,NIBLO)
          DO IJ=KIJS,KIJL
            RU10OLD(IJ,IG)=0.
            RTHWOLD(IJ,IG)=0.
            RUSOLD(IJ,IG)=0.
            RTAUW(IJ,IG)=0.
            RZ0OLD(IJ,IG)=0.
            RADOLD(IJ,IG)=ROAIR
            RZIDLOLD(IJ,IG)=0.
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
      ENDDO

      FCRANGE = 0.
      CALL GETWND (RU10OLD, RTHWOLD, RADOLD, RZIDLOLD, KICEMASK,
     &             1, NIBLO, CDTPRO,
!------------------ dummy parameters & arrays
     &             NFIELDS, NGPTOTG, NCA, NRA, FCRANGE,
     &             ISEC1, ISEC2, FIELDS, MASK_IN,
!---------------------------------------------
     &             LWNDFILE, LCLOSEWND)
      IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. GETWND DONE'

!     USE DATA FROM A FILE CONTAINING WIND SPEED MODIFIED BY
!     A PREVIOUS WAVE MODEL RUN ON WHICH THIS RESTART IS BASED
!     IF IT EXISTS, OTHERWISE PRESET WILL USE THE INFORMATION FROM
!     THE ATMOSPHERIC MODEL WIND SPEED ARCHIVE.

      LWAVEWIND=.FALSE.

      FILNM='uwavein'
      LIU = LEN_TRIM(FILNM)
      FILNM=FILNM(1:LIU)
      INQUIRE(FILE=FILNM,EXIST=LWAVEWIND)

      IF(LWAVEWIND) THEN

        IPARAM=245
        LLONLYPOS=.FALSE.
        CALL READWGRIB(IU06,FILNM,IPARAM,CDTPRO,RU10OLD,1,NIBLO,KEND,
     &                 KZLEVUWAVE, LLONLYPOS)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. READWGRIB DONE FOR ',FILNM

        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' A DATA FILE CONTAINING WIND SPEED INFORMATION'
        WRITE(IU06,*) ' AS PROVIDED BY A PREVIOUS WAVE MODEL RUN'
        WRITE(IU06,*) ' WAS USED TO UPDATE THE INPUT ATMOSPHERIC WINDS'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE INPUT WINDS AND DRAG COEFFICIENT ARE FOUND'
        WRITE(IU06,*) ' TO HAVE BEEN DETERMINED FOR HEIGHT AT '
        WRITE(IU06,*) ' ',KZLEVUWAVE,' m'
        IF (ITEST.GT.0) CALL FLUSH(IU06) 

        ILEV=1


      ELSE

        ILEV=1
        KZLEVUWAVE=10
        WRITE(IU06,*) ' '
        WRITE(IU06,*) '          !!!! NOTE !!!!'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' NO INFORMATION ON WIND SPEEDS FROM THE'
        WRITE(IU06,*) ' WAVE MODEL WAS PROVIDED'
        WRITE(IU06,*) ' NO UPDATE OF THE INPUT ATMOSPHERIC WINDS'
        WRITE(IU06,*) ' WAS POSSIBLE'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE INPUT WIND AND DRAG COEFFICIENT ARE ASSUMED'
        WRITE(IU06,*) ' TO HAVE BEEN DETERMINED FOR HEIGHT AT '
        WRITE(IU06,*) ' ',KZLEVUWAVE,' m'

        IF (ITEST.GT.0) CALL FLUSH(IU06) 
      ENDIF
 
!     TEST WHETHER THE HEIGHT OF THE INPUT WINDS IS THE SAME AS DEFINED
!     FOR THE REST OF THE RUN (SEE CALL TO STRESS IN PREPROC). IF NOT
!     RECOMPUTE THE TABLE. THE ORIGINAL TABLE WILL BE SWAP BACK AT THE
!     END OF THE ROUTINE.

      IF (KZLEVUWAVE.NE.NINT(XNLEV(1))) THEN
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE REFERENCE HEIGHT TO BE USED IN WAMODEL'
        WRITE(IU06,*) ' ',XNLEV(1) 
        WRITE(IU06,*) ' IS DIFFERENT THAN THE INPUT FIELDS HEIGHT'
        WRITE(IU06,*) ' ',FLOAT(KZLEVUWAVE)
        WRITE(IU06,*) ' THE NECESSARY ADJUSTMENTS WILL BE MADE'
        WRITE(IU06,*) ' TO DETERMINE THE INITIAL FIELDS.'
        WRITE(IU06,*) ' '
        IF (ITEST.GT.0) CALL FLUSH(IU06) 

        ALLOCATE(TEMPTAUT(0:ITAUMAX,0:JUMAX,JPLEVT))
        TEMPTAUT=TAUT
        TEMPXNLEV=XNLEV(1)
        XNLEV(1)=KZLEVUWAVE

        CALL STRESS(IU06,ITEST)
        IF (ITEST.GE.2)
     1    WRITE(IU06,*) ' SUB. BUILDSTRESS: CALL STRESS DONE'
        IF (ITEST.GT.0) CALL FLUSH(IU06) 

      ENDIF


      IF(.NOT.ALLOCATED(CD)) ALLOCATE(CD(NIBLO,NBLO))

!     INITIALISE CD USING THE FRICTION VELOCITY FOR TAUW=0.

      DO IG=1,IGL
        IF(LMESSPASS) THEN
          CALL AIRSEA (RU10OLD(1,IG), RTAUW(1,IG),RUSOLD(1,IG),
     &     RZ0OLD(1,IG), 1, KEND, ILEV)
          IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. AIRSEA DONE AT 1'

          DO IJ=1,KEND
!! the numerical relation between usold and u10old should always be as
!! in outgrid
            CDINV = MAX(RU10OLD(IJ,IG)**2,EPSU10)/
     &              MAX(RUSOLD(IJ,IG)**2,EPSUS)
            CDINV = MIN(CDINV,10000.) 
            CD(IJ,IG) = 1./CDINV
          ENDDO
        ELSE
          CALL AIRSEA (RU10OLD(1,IG), RTAUW(1,IG),RUSOLD(1,IG),
     &     RZ0OLD(1,IG), IJS(IG), IJL(IG), ILEV)
          IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. AIRSEA DONE AT 1'

          DO IJ=IJS(IG),IJL(IG)
!! the numerical relation between usold and u10old should always be as
!! in outgrid
            CDINV = MAX(RU10OLD(IJ,IG)**2,EPSU10)/
     &              MAX(RUSOLD(IJ,IG)**2,EPSUS)
            CDINV = MIN(CDINV,10000.) 
            CD(IJ,IG) = 1./CDINV
          ENDDO
        ENDIF
      ENDDO

      IF(.NOT.LNOCDIN) THEN
        IPARAM=233
        LLONLYPOS=.TRUE.
        FILNM='cdwavein'
        CALL READWGRIB(IU06,FILNM,IPARAM,CDTPRO,CD,1,NIBLO,KEND,
     &                 KZLEVCD, LLONLYPOS)
        IF (ITEST.GT.0) WRITE (IU06,*) ' SUB. READWGRIB DONE FOR ',FILNM

!       TEST REFERENCE LEVEL FOR UWAVE AND CD

        IF(KZLEVUWAVE.NE.0.AND.KZLEVCD.NE.0) THEN
          IF(KZLEVUWAVE.NE.KZLEVCD) THEN
            WRITE(IU06,*)'************************************'
            WRITE(IU06,*)'*                                  *'
            WRITE(IU06,*)'* FATAL ERROR IN SUB BUILDSTRESS   *'
            WRITE(IU06,*)'* ==============================   *'
            WRITE(IU06,*)'* REFERENCE LEVEL FOR CD AND UWAVE *'
            WRITE(IU06,*)'* DIFFER. PROGRAM WILL ABORT       *'
            WRITE(IU06,*)'* REF. LEVEL CD = ',KZLEVCD
            WRITE(IU06,*)'* REF. LEVEL UWAVE = ',KZLEVUWAVE
            WRITE(IU06,*)'*                                  *'
            WRITE(IU06,*)'************************************'
            CALL ABORT1
          ENDIF
        ENDIF

!       COMPUTE TAUW,USOLD AND Z0OLD

        ILEV=1

        DO IG=1,IGL
           IF(LMESSPASS) THEN
             IJSTART=1
             IJSTOP=KEND
           ELSE
             IJSTART=IJS(IG)
             IJSTOP=IJL(IG)
           ENDIF
           DO IJ=IJSTART,IJSTOP
            CDSQRTINV = MIN(1./SQRT(CD(IJ,IG)),100.0)
            RZ0OLD(IJ,IG) = XNLEV(ILEV)*EXP(-XKAPPA*CDSQRTINV)
!! RUSOLD will first contain its square
!! the numerical relation between usold and u10old should always be as
!! in outgrid
            RUSOLD(IJ,IG) = CD(IJ,IG)*MAX(RU10OLD(IJ,IG)**2,EPSU10)
            RUSOLD(IJ,IG) = MAX(RUSOLD(IJ,IG),EPSUS)
            CHARNOCK = G*RZ0OLD(IJ,IG)/RUSOLD(IJ,IG)
            CHARNOCK = MAX(CHARNOCK,ALPHA)
            RTAUW(IJ,IG) = RUSOLD(IJ,IG)*(1.-(ALPHA/CHARNOCK)**2)
            RUSOLD(IJ,IG) = SQRT(RUSOLD(IJ,IG))
          ENDDO
        ENDDO
      ENDIF

      IF(ALLOCATED(CD)) DEALLOCATE(CD)

      IF(ALLOCATED(TEMPTAUT)) THEN
        XNLEV(1)=TEMPXNLEV
        TAUT=TEMPTAUT
        DEALLOCATE(TEMPTAUT)
      ENDIF

      IF (ITEST.GE.2)
     &    WRITE(IU06,*) ' SUB. BUILDSTRESS: INPUT OF RESTART FILES DONE'

      RETURN
      END SUBROUTINE BUILDSTRESS
