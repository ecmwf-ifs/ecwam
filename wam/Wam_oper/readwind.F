
      SUBROUTINE READWIND (CDTWIR, NFIELDS, NGPTOTG, NCA, NRA, 
     &                     ISEC1, ISEC2, FIELDS, LWNDFILE,LCLOSEWND)

!***  *READWIND* - PROGRAM TO REFORMAT WIND FIELDS FROM MARS 
!                 EITHER FROM A FILE OR 
!                 FROM THE ARRAY "FIELDS" PASSED FROM ATM. MODEL

!    PETER JANSSEN      ECMWF        FEBRUARY 1987                


!    MODIFIED BY:                                                 
!    ------------                                                       

!    LIANA ZAMBRESKY     GKSS/ECMWF   JULY 1988                  

!        1. PASSING OF ARRAYS AS FORMAL PARMETERS TO SAVE MEMORY
!        2. PROCESSING WINDS FOR ONLY ONE WIND TIME STEP       

!    PEDRO VITERBO       ECMWF        OCTOBER 1988            

!        DEALING WITH MARS PACKED DATA.                      

!    H. GUNTHER      ECMWF/GKSS       OCTOBER 1991          

!        GRIB EDITION 1.                                   

!     J. DOYLE    NRL/ECMWF     OCT. 1996
!                 DECODES GRID WIND FIELDS WHICH ARE PASSED AS ARRAYS
!                 WHEN MODEL IS COUPLED TO ATMOPHERIC MODEL
!                 SEE LWNDFILE

!    B. HANSEN    ECMWF 1997
!                 RENAME SUBROUTINE TO MATCH THE NAME OF THE FILE.
!                 USE INCLUDE FILES FOR ALL COMMON BLOCKS.

!    J. BIDLOT    ECMWF 1998 RENAMED READWND_MARS AS READWIND

!    S. ABDALLA   ECMWF 2001
!                 GETS ATMOSPHERIC PARAMETERS FOR GUSTINESS
!                 AND AIR DENSITY COMPUTATIONS

!     PURPOSE                                                       
!     -------                                                      

!     *READWIND*
!         READS TARGET FILES FROM THE MARS ARCHIVE AND INFERS FROM THIS 
!             A) THE DATE OF THE WIND FIELD                             
!             B) THE DEFINITION OF THE MARS GRID                        
!             C) THE U AND V COMPONENT OF THE WIND FIELD                
!             D) THE TYPE OF WIND FIELD AND                             
!             E) THE UNPACKED WINDS                                     

!     INTERFACE                                                         
!     ---------                                                         

!     *CALL* *READWIND (CDTWIR, NFIELDS, NGPTOTG, NCA, NRA, 
!    &                  ISEC1, ISEC2, FIELDS, LWNDFILE,LCLOSEWND)*
!
!        *CDTWIR* - DATE/TIME OF THE DATA READ                          
!        *NFIELDS*- NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!        *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!        *ISEC1*  - GRIB ISEC1 (TO EXTRACT DATE/TIME OF ATM. FIELDS)
!        *ISEC2*  - GRIB ISEC2 (TO EXTRACT ATM. GRID CHARACTERISTICS)
!        *FIELDS* - ATMOSPHERIC FIELDS AS FOLLOWS:
!                   FIELDS(:,1) = U component of wind speed (U10)
!                   FIELDS(:,2) = V component of wind speed (V10)
!                   FIELDS(:,3) = surface specific humidity (at 2m)
!                   FIELDS(:,4) = meas surface atm. pressure (MSL)
!                   FIELDS(:,5) = air temperature at surface (SKT)
!                   FIELDS(:,6) = Reserved for ice boundary
!                   FIELDS(:,7) = Zi/L used for gustiness
!        *LWNDFILE*  IF TRUE DATA ARE ACCESSED FROM THE INPUT FILE
!                    sfcwind
!                    ELSE THE DATA WILL BE SUPPLIED VIA  "FIELDS"
!        *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND THE UNIT 
!                    RESET

!     EXTERNALS                                                         
!     ---------                                                         

!     *ABORT1*         TERMINATES PROCESSING
!     *INCDATE*        INCREMENTS DATE 
!     *PBGRIB*
!     *PBOPEN*
!     *PBCLOSE*
!     *KGRIBSIZE*
!     *FLUSH*
!     *GRB2WGRD*       TRANSFORM GRIB FIELDS TO WAM GRID FIELDS
!     *INITIALINT*     COMPUTES INTERPOLATION COEFFICIENTS (FOR FLDINTER)
!     *FLDINTER*       INTERPOLATES ATMOSPHERIC FIELDS TO WAM GRID

! --------------------------------------------------------------------- 

      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWPARAM , ONLY : NGX      ,NGY      ,SWAMPWIND,CLDOMAIN ,
     &            SWAMPWIND2,DTNEWWIND,LTURN90
      USE YOWSTAT  , ONLY : CDATEA   ,IDELWI   ,LADEN    ,LGUST
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWNDG  , ONLY : DLAM     ,DPHI     ,RLATS    ,RLATN    ,
     &            RLONL    ,RLONR    ,KCOL     ,KROW     ,ICODE    ,
     &            IWPER    ,ICOORD
      USE YOWWIND  , ONLY : NC       ,NR       ,UWND     ,VWND     ,
     &            AIRD     ,PZIDL

      USE YOWPCONS , ONLY : ROAIR, Rconstant, EpsWaterVapor

! --------------------------------------------------------------------- 

      INTEGER, ALLOCATABLE :: KGRIB(:)

      INTEGER NFIELDS, NGPTOTG, NCA, NRA
      INTEGER ISEC1(255), ISEC2(22+NRA)
      INTEGER, ALLOCATABLE, SAVE :: RLONRGG(:)
      REAL FIELDS(NGPTOTG,NFIELDS)
      LOGICAL, SAVE :: LLINTERPOL

      INTEGER, ALLOCATABLE, SAVE :: KK(:), II(:,:), IIP(:,:)
      REAL, ALLOCATABLE, SAVE :: DK1(:), DII1(:,:), DIIP1(:,:)

      REAL, DIMENSION(NFIELDS) :: FAVG,FMIN,FMAX

      CHARACTER FILNM*24
      CHARACTER*12 CDTWIR, CDTTURN
      LOGICAL FRSTIME, LWNDFILE, LCLOSEWND

      LOGICAL FRSTMATM
      REAL, ALLOCATABLE :: AirDenLocal (:)

      DATA FRSTIME / .TRUE. /
      DATA FRSTMATM / .TRUE. /
      SAVE IUNIT1, IPERIODIC
      SAVE FRSTIME, FRSTMATM

! --------------------------------------------------------------------  

!     for the swamp case create your own wind fields
      IF(CLDOMAIN.EQ.'s') goto 8000

      NC=NGX
      NR=NGY

 
      IF(LWNDFILE) THEN

!       WIND FIELDS ARE SUPPLIED FROM INPUT FILE sfcwindin
!       AS IS THE CASE IN THE STAND ALONE VERSION OR AT INITIAL INPUT
!       TIME IN COUPLED MODEL.

!       1.0 INITIALISE VALUES. 
!           ------------------ 

        NBIT=175000

        ICOORD = 1

        IF (FRSTIME) THEN
          FILNM='sfcwindin'
          CALL PBOPEN(IUNIT1,FILNM,'r',KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '****************************************'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBOPEN      *'
            WRITE (IU06,*) '*  IN READWIND                         *'
            IF(KRET.EQ.-1)
     &       WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILNM
            IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILNM
            IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '****************************************'
            CALL ABORT1
          ENDIF

!         TRYING TO FIGURE OUT THE SIZE OF EACH FIELD

          FRSTIME = .FALSE.
        ENDIF  


! --------------------------------------------------------------------  

!*      2.0 READ MARS U and V WIND COMPONENTS  IN GRIB CODE FORMAT.
!           -------------------------------------------------------

        ISIZE=NBIT

        IWPER = 1


        IF(.NOT.ALLOCATED(UWND)) ALLOCATE(UWND(NGX,NGY)) 
        IF(.NOT.ALLOCATED(VWND)) ALLOCATE(VWND(NGX,NGY))
        IF(.NOT.ALLOCATED(AIRD)) ALLOCATE(AIRD(NGX,NGY))
        IF(.NOT.ALLOCATED(PZIDL)) ALLOCATE(PZIDL(NGX,NGY))

        AIRD=ROAIR ! Initializes AirDensity field by the default value.
                ! Actual values can be computed here as well but needs
                ! to get the corresponding fields in the file sfcwindin.
        PZIDL = 0.

        WND: DO IVAR=1,2

2001      KPLENG=ISIZE*KIND(ISIZE)
          IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))

          CALL PBGRIB(IUNIT1,KGRIB,KPLENG,KLEN,KRET)

          IF (KRET.EQ.-1) THEN
            WRITE(IU06,*) ' *******************************************'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *       FATAL ERROR IN SUB. READWIND      *'
            WRITE(IU06,*) ' *       ===========================       *'
            WRITE(IU06,*) ' * PBGRIB EOF ON FILE ',FILNM 
            WRITE(IU06,*) ' *       NO MORE WINDS                     *'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS      *'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *******************************************'
            CALL ABORT1
          ELSE IF (KRET.EQ.-2) THEN
            WRITE(IU06,*) ' *******************************************'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *       FATAL ERROR IN SUB. READWIND      *'
            WRITE(IU06,*) ' *       ===========================       *'
            WRITE(IU06,*) ' * PBGRIB ERROR ON FILE ',FILNM
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS      *'
            WRITE(IU06,*) ' *                                         *'
            WRITE(IU06,*) ' *******************************************'
            CALL ABORT1
          ELSE IF (KRET.EQ.-3) THEN
            DEALLOCATE(KGRIB)
            CALL KGRIBSIZE(IU06, IUNIT1, KLEN, ISIZE, 'READWIND')
            GOTO 2001
          ENDIF

          IF (ITEST.GT.1) THEN
            WRITE(IU06,*) ' SUB. READWIND - READ FROM ',FILNM
          ENDIF

! ----------------------------------------------------------------------

!*        3.0 UNPACK MARS FIELDS.
!             -------------------
          KLENG = KLEN 
          ZDUM=0.
          IF(IVAR.EQ.1) THEN
            CALL GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                     NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                     KGRIB, KLENG ,
     &                     NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                     ZDUM, ZDUM,
     &                     CDTWIR, IFORP, IPARAM, KZLEV, IDM, IDM, UWND)
          ELSE
            CALL GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                     NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                     KGRIB, KLENG ,
     &                     NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                     ZDUM, ZDUM,
     &                     CDTWIR, IFORP, IPARAM, KZLEV, IDM, IDM, VWND)
          ENDIF

!*      3.1.2 DETERMINE CODE FOR WIND FIELD TYPE:
!             165 = U10 FIELD , 180 = STRESS FIELD.
!             ------------------------------------ 

          IF (IPARAM.EQ.165 .OR. IPARAM.EQ.166) THEN
            ICODE = 3
          ELSEIF (IPARAM.EQ.180 .OR. IPARAM.EQ.181) THEN
            ICODE = 2
          ELSE
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +      WARNING ERROR IN SUB. READWIND     +'
            WRITE(IU06,*) ' +      ==============================     +'
            WRITE(IU06,*) ' + SUSPICIOUS WIND FIELD TYPE IN MARS BLOCK+'
            WRITE(IU06,*) ' + TYPE IS = ', IPARAM
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +    PROGRAM CONTINUES                    +'
            WRITE(IU06,*) ' +                                         +'
            WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
            CALL ABORT1
          ENDIF


!*      3.2.1 DETERMINE GRID PARAMETERS.                
!             -------------------------- 

          KCOL = NGX
          KROW = NGY
          RLATN = AMONOP
          RLATS = AMOSOP
          RLONL = AMOWEP
          RLONR = AMOEAP
          IF(KCOL.NE.1) THEN
            DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
          ELSE
            DLAM  = 0.
          ENDIF
          IF(KROW.NE.1) THEN
            DPHI  = (RLATN-RLATS)/REAL(KROW-1)
          ELSE
            DPHI  = 0.
          ENDIF

        ENDDO WND

! --------------------------------------------------------------------  

!*      5. WRITE TEST OUTPUT.
!          ------------------

        IF (ITEST.GT.2) THEN
          WRITE(IU06,*) ' READWIND -  WIND FIELD FOR THE CDTWIR = ',
     &     CDTWIR
          WRITE(IU06,66) ((UWND(I,J),I=1,24),J=1,5)
          WRITE(IU06,66) ((VWND(I,J),I=1,24),J=1,5)
 66       FORMAT(1X,24F5.2)
          CALL FLUSH(IU06)
        ENDIF

! --------------------------------------------------------------------  
                                                                        
        DEALLOCATE(KGRIB)

        IF (LCLOSEWND) THEN
          CALL PBCLOSE(IUNIT1,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '****************************************'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBCLOSE     *'
            WRITE (IU06,*) '*  IN READWIND                         *'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '****************************************'
            CALL FLUSH(IU06)
            CALL ABORT1
          ENDIF
          FRSTIME = .TRUE.
          IF (ITEST.GT.1)
     &     WRITE(IU06,*) ' SUB. READWIND - PBCLOSE ON ',FILNM
        ENDIF  

      ELSE


!*      6.0 READ ATMOS. MODEL FIELDS FOR COUPLED MODEL
!           ------------------------------------------

        IF (FRSTMATM) THEN
!         6.1 First time ONLY, compute interpolation coefficients:
!             ----------------------------------------------------
          IF (.NOT.ALLOCATED(RLONRGG)) ALLOCATE (RLONRGG(NRA))
          IF (.NOT.ALLOCATED(KK)) ALLOCATE (KK(NGY))
          IF (.NOT.ALLOCATED(II)) ALLOCATE (II(NGX,NGY))
          IF (.NOT.ALLOCATED(IIP)) ALLOCATE (IIP(NGX,NGY))
          IF (.NOT.ALLOCATED(DK1)) ALLOCATE (DK1(NGY))
          IF (.NOT.ALLOCATED(DII1)) ALLOCATE (DII1(NGX,NGY))
          IF (.NOT.ALLOCATED(DIIP1)) ALLOCATE (DIIP1(NGX,NGY))
          CALL INITIALINT(IU06, ITEST, NCA, NRA, ISEC2,
     &                    NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                    AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                    RLONRGG, 
     &                    LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)
          FRSTMATM = .FALSE.
          IF (ITEST.GT.1)
     &     WRITE(IU06,*) ' WAM_READWIND - INTERPOLATION COEFFICIENTS'
     &                     //' COMPUTED'
          IF (LADEN.AND.NFIELDS.LT.5) THEN 
            LADEN = .FALSE.
            WRITE(IU06,*) ' WAM_READWIND - WARNING: '
            WRITE(IU06,*) '      AIR DENSITY RUN WAS REQUESTED BUT '
     &                 //'NOT ENOUGH FIELDS PASSED.'
            WRITE(IU06,*) '      NFIELDS = ', NFIELDS
            WRITE(IU06,*) '      NFIELDS MUST BE AT LEAST 5'
            WRITE(IU06,*) '      RESET  LADEN  TO  .FALSE.'
          ENDIF
          IF (LADEN) THEN 
            WRITE(IU06,*) ' WAM_READWIND - VARIABLE AIR DENSITY RUN'
          ELSE
            WRITE(IU06,*) ' WAM_READWIND - AIR DENSITY IS '
     &                     //' CONSTANT = ', ROAIR
          ENDIF
          IF (LGUST.AND.NFIELDS.LT.7) THEN 
            LGUST = .FALSE.
            WRITE(IU06,*) ' SUB. READWIND - WARNING: '
            WRITE(IU06,*) '      GUSTINESS RUN WAS REQUESTED BUT '
     &                 //'NOT ENOUGH FIELDS PASSED.'
            WRITE(IU06,*) '      NFIELDS = ', NFIELDS
            WRITE(IU06,*) '      NFIELDS MUST BE AT LEAST 7'
            WRITE(IU06,*) '      RESET  LGUST  TO  .FALSE.'
          ENDIF
          IF (LGUST) THEN 
            WRITE(IU06,*) ' WAM_READWIND - GUSTINESS RUN'
          ELSE
            WRITE(IU06,*) ' WAM_READWIND - NO GUSTINESS EFFECT.'
          ENDIF
        ENDIF
     
        IWPER = 1
        ICOORD = 1

        IF(.NOT.ALLOCATED(UWND)) ALLOCATE(UWND(NGX,NGY)) 
        IF(.NOT.ALLOCATED(VWND)) ALLOCATE(VWND(NGX,NGY))
        IF(.NOT.ALLOCATED(AIRD)) ALLOCATE(AIRD(NGX,NGY))
        IF(.NOT.ALLOCATED(PZIDL)) ALLOCATE(PZIDL(NGX,NGY))

        CALL FLDINTER (IU06, ITEST, NGPTOTG, NCA, NRA, FIELDS(:,1), 
     &                 NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                 AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                 RLONRGG, UWND, 
     &                 LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)

        CALL FLDINTER (IU06, ITEST, NGPTOTG, NCA, NRA, FIELDS(:,2),
     &                 NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                 AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                 RLONRGG, VWND,
     &                 LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)

!       computation of the norms of input fields
        WRITE(IU06,*) ' NORM OF INPUT ATMOSPHERIC FIELDS : ' 
        DO IFLD=1,NFIELDS
          FAVG(IFLD)=0.
          FMIN(IFLD)=999999999999.
          FMAX(IFLD)=-99999999999.
          DO I=1,NGPTOTG 
            FAVG(IFLD)=FAVG(IFLD)+FIELDS(I,IFLD)
            FMIN(IFLD)=MIN(FMIN(IFLD),FIELDS(I,IFLD))
            FMAX(IFLD)=MAX(FMAX(IFLD),FIELDS(I,IFLD))
          ENDDO
          FAVG(IFLD)=FAVG(IFLD)/NGPTOTG
          WRITE(IU06,*) ' ',FAVG(IFLD),FMIN(IFLD),FMAX(IFLD),NGPTOTG
        ENDDO
        WRITE(IU06,*) ' '
        IF(ITEST.GT.0) CALL FLUSH(IU06)

! Computations of the Air Density:
        IF (LADEN) THEN 
          ALLOCATE (AirDenLocal(NGPTOTG))
          Dummy = Rconstant * EpsWaterVapor

!         The density computation below assumes:
!           FIELDS(:,3) = surface specific humidity (at 2m)
!           FIELDS(:,4) = meas surface atmospheric pressure (MSL)
!           FIELDS(:,5) = air temperature at surface (skin temp., SKT)

          AirDenLocal = ROAIR
          DO I=1,NGPTOTG 
            IF (FIELDS(I,4).GT.0. )
     &        AirDenLocal(I) = FIELDS(I,4)/
     &              (Rconstant*FIELDS(I,5)+Dummy*FIELDS(I,3))
          ENDDO
          CALL FLDINTER (IU06, ITEST, NGPTOTG, NCA, NRA, AirDenLocal,
     &                   NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                   AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                   RLONRGG, AIRD,
     &                   LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)

          DEALLOCATE (AirDenLocal)
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. READWIND - AIR DENSITY COMPUTED'
     &                     //' AND INTERPOLATED'
        ELSE
          AIRD = ROAIR
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. READWIND - AIR DENSITY USED AS '
     &                    //' CONSTANT = ', ROAIR
        ENDIF

! The following is for ice boundary interpolation 
! un-comment when things are ready 
! (A flag may be needed for this option)
!       IF (NFIELDS.GE.6) THEN
!         CALL FLDINTER (IU06, ITEST, NGPTOTG, NCA, NRA, FIELDS(:,6),
!    &                    NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
!    &                    AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
!    &                    RLONRGG, ICEBND,
!    &                    LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)
!       ELSE
!         ICEBND = 0. !!!!
!       ENDIF


        IF (LGUST) THEN 
          CALL FLDINTER (IU06, ITEST, NGPTOTG, NCA, NRA, FIELDS(:,7),
     &                   NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                   AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                   RLONRGG, PZIDL,
     &                   LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. READWIND - Zi/L INTERPOLATED'

        ELSE
          PZIDL = 0.
        ENDIF


! Decode the date of the fields:
!       START DATE.

        IYYYY=(ISEC1(21)-1)*100+ISEC1(10)
        WRITE(CDTWIR,'(I4.4,4I2.2)') IYYYY,ISEC1(11),ISEC1(12),
     &                             ISEC1(13),ISEC1(14)

!*      DETERMINE TIME UNIT.

        IF (ISEC1(15).EQ.0) THEN
          JCONS = 60
        ELSEIF (ISEC1(15).EQ.1) THEN
          JCONS = 3600
        ELSEIF (ISEC1(15).EQ.2) THEN
          JCONS = 86400
        ELSE
          WRITE(IU06,*) '***********************************'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*    PROBLEM IN SUB. READWIND     *'
          WRITE(IU06,*) '*    =========================    *'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '*  UNABLE TO DECODE THE TIME UNIT *'
          WRITE(IU06,*) '*                                 *'
          WRITE(IU06,*) '***********************************'
          CALL FLUSH(IU06)
        ENDIF

!       DETERMINE FORECAST PERIOD.

        IFORP = JCONS*ISEC1(16)

!       GET DATE OF THE FORECAST INSTEAD OF STARTING DATE

        CALL INCDATE (CDTWIR,IFORP)


        IPARAMU=165
        IPARAMV=166

!*        6.1.2 DETERMINE CODE FOR WIND FIELD TYPE:
!               165 = U10 FIELD , 180 = STRESS FIELD.
!               ------------------------------------

        IF (IPARAMU.EQ.165 .AND. IPARAMV.EQ.166) THEN
          ICODE = 3
        ELSEIF (IPARAMU.EQ.180 .AND. IPARAMV.EQ.181) THEN
          ICODE = 2
        ELSE
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
          WRITE(IU06,*) ' +                                         +'
          WRITE(IU06,*) ' +      WARNING ERROR IN SUB. READWIND     +'
          WRITE(IU06,*) ' +      ==============================     +'
          WRITE(IU06,*) ' + SUSPICIOUS WIND FIELD TYPE IN MARS BLOCK+'
          WRITE(IU06,*) ' + U TYPE IS = ', IPARAMU
          WRITE(IU06,*) ' + V TYPE IS = ', IPARAMV
          WRITE(IU06,*) ' +                                         +'
          WRITE(IU06,*) ' +    PROGRAM CONTINUES                    +'
          WRITE(IU06,*) ' +                                         +'
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
          CALL ABORT1
        ENDIF

!*      6.2 GRID DEFINITION BLOCK.
!           ---------------------

!*        6.2.1 DETERMINE GRID PARAMETERS.
!               --------------------------

        KCOL = NGX
        KROW = NGY
        RLATN = AMONOP
        RLATS = AMOSOP
        RLONL = AMOWEP
        RLONR = AMOEAP
        IF(KCOL.NE.1) THEN
          DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
        ELSE
          DLAM  = 0.
        ENDIF
        IF(KROW.NE.1) THEN
          DPHI  = (RLATN-RLATS)/REAL(KROW-1)
        ELSE
          DPHI  = 0.
        ENDIF


!*      6.3. WRITE TEST OUTPUT.
!            ------------------

        IF (ITEST.GT.2) THEN
          WRITE(IU06,*) ' READWIND -  WIND FIELD FOR THE CDTWIR = ',
     &     CDTWIR
          WRITE(IU06,63) ((UWND(I,J),I=1,24),J=51,55)
          WRITE(IU06,63) ((VWND(I,J),I=1,24),J=51,55)
          WRITE(IU06,*)' READWIND - AIR DENSITY FIELD FOR THE CDTWIR= ',
     &     CDTWIR
          WRITE(IU06,63) ((AIRD(I,J),I=1,24),J=51,55)
          WRITE(IU06,*) ' READWIND - ZI/L FIELD FOR THE CDTWIR = ',
     &     CDTWIR
          WRITE(IU06,63) ((PZIDL(I,J),I=1,24),J=51,55)
          CALL FLUSH(IU06)
 63       FORMAT(1X,24F5.2)
        ENDIF

      ENDIF

      RETURN


!     SWAMP CASE : CREATE YOUR OWN WIND FIELD

8000  CONTINUE

!     THE DIMENSIONS ARE HARD CODED !!
      NC=NGX
      NR=NGY

      IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
      ALLOCATE(UWND(NC,NR))                                
      IF(ALLOCATED(VWND)) DEALLOCATE(VWND)
      ALLOCATE(VWND(NC,NR))                                
      IF(ALLOCATED(AIRD)) DEALLOCATE(AIRD)
      ALLOCATE(AIRD(NC,NR))                                
      IF(ALLOCATED(PZIDL)) DEALLOCATE(PZIDL)
      ALLOCATE(PZIDL(NC,NR))                                

      IF (FRSTIME) THEN
         RLATS = AMOSOP
         RLATN = AMONOP 
         RLONL = AMOWEP 
         RLONR = AMOEAP 
         ICODE = 3
         CDTWIR = CDATEA
         CALL INCDATE(CDTWIR,-IDELWI)
         KCOL = NC
         KROW = NR
         IWPER = 0
         ICOORD = 1
         IF(KCOL.NE.1) THEN
           DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
         ELSE
           DLAM  = 0.
         ENDIF
         IF(KROW.NE.1) THEN
           DPHI  = (RLATN-RLATS)/REAL(KROW-1)
         ELSE
           DPHI  = 0.
         ENDIF
         FRSTIME = .FALSE.
      ENDIF

      CALL INCDATE(CDTWIR,IDELWI)

      CDTTURN = CDATEA
      IDTTURN=DTNEWWIND*3600
      CALL INCDATE(CDTTURN,IDTTURN)

      IF(SWAMPWIND2.LE.0. .OR. CDTWIR.LE.CDTTURN) THEN
        WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',
     &                  SWAMPWIND 
        DO J=1,NR
          DO I=1,NC
            UWND(I,J)=0.0
          ENDDO
        ENDDO
        DO J=1,NR
          DO I=1,NC
           VWND(I,J)=SWAMPWIND
          ENDDO
        ENDDO
      ELSE
        WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',
     &                  SWAMPWIND2 
        IF(LTURN90) THEN
          DO J=1,NR
            DO I=1,NC
              UWND(I,J)=SWAMPWIND2
            ENDDO
          ENDDO
          DO J=1,NR
            DO I=1,NC
             VWND(I,J)=0.
            ENDDO
          ENDDO
        ELSE
          DO J=1,NR
            DO I=1,NC
              UWND(I,J)=0.0
            ENDDO
          ENDDO
          DO J=1,NR
            DO I=1,NC
             VWND(I,J)=SWAMPWIND2
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      DO J=1,NR
        DO I=1,NC
          AIRD(I,J)=ROAIR
          PZIDL(I,J)=0.
        ENDDO
      ENDDO

!*      WRITE TEST OUTPUT.
!       ------------------

        IF (ITEST.GT.2) THEN
          WRITE(IU06,*) ' READWIND -  WIND FIELD FOR THE CDTWIR = ',
     &     CDTWIR
          WRITE(IU06,67) ((UWND(I,J),I=1,24),J=1,5)
          WRITE(IU06,67) ((VWND(I,J),I=1,24),J=1,5)
 67       FORMAT(1X,24F5.2)
        ENDIF

      RETURN


      END SUBROUTINE READWIND
