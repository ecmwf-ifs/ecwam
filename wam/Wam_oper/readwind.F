      SUBROUTINE READWIND (CDTWIR, NFIELDS, NGPTOTG, NCA, NRA, FCRANGE,
     &                     ISEC1, ISEC2, FIELDS, MASK_IN,
     &                     LWNDFILE,LCLOSEWND,
     &                     KINF, KSUP,
     &                     FILNM, LLNOTOPENED )

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
!
!    P.A.E.M. JANSSEN   ECMWF 2003
!                 USE NEUTRAL WINDS RATHER THAN REAL WINDS
!                 TRANSFER AIR DENSITY CALCULATION TO IFS

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

!     *CALL* *READWIND (CDTWIR, NFIELDS, NGPTOTG, NCA, NRA, FCRANGE,
!    &                  ISEC1, ISEC2, FIELDS, LWNDFILE,LCLOSEWND,
!    &                  KINF, KSUP,
!    &                  FILNM, LLNOTOPENED )*
!
!        *CDTWIR* - DATE/TIME OF THE DATA READ                          
!        *NFIELDS*- NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!        *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!        *FCRANGE*- FORECAST RANGE IN SECONDS
!        *ISEC1*  - GRIB ISEC1 (TO EXTRACT DATE/TIME OF ATM. FIELDS)
!        *ISEC2*  - GRIB ISEC2 (TO EXTRACT ATM. GRID CHARACTERISTICS)
!        *FIELDS* - ATMOSPHERIC FIELDS AS FOLLOWS:
!                   FIELDS(:,1) = U COMPONENT OF WIND SPEED (U10)
!                   FIELDS(:,2) = V COMPONENT OF WIND SPEED (V10)
!                   FIELDS(:,3) = AIR DENSITY
!                   FIELDS(:,4) = ZI/L USED FOR GUSTINESS
!                   FIELDS(:,5) = RESERVED FOR ICE BOUNDARY
!                                 NOT YET USED !
!        *LWNDFILE*  IF TRUE DATA ARE ACCESSED FROM THE INPUT FILE
!                    sfcwind
!                    ELSE THE DATA WILL BE SUPPLIED VIA  "FIELDS"
!        *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND THE UNIT 
!                    RESET
!        *KINF*   - SMALLEST WAM GRID POINT INDEX USED BY THE PE 
!        *KSUP*   - LARGEST WAM GRID POINT INDEX USED BY THE PE 
!        *FILNM*    FILENAME OF INPUT FILE
!        *LLNOTOPENED*  TRUE IF THE INPUT FILE HAS TO BE OPENED

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

      USE YOWCOUP  , ONLY : LWCOUNORMS
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3
      USE YOWGRID  , ONLY : IGL      ,NLONRGG
      USE YOWICE   , ONLY : LICERUN  ,IPARAMCI ,CI
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,IXLG     , 
     &            KXLT
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO    ,CLDOMAIN ,
     &            SWAMPWIND,SWAMPWIND2,DTNEWWIND,LTURN90 ,LWDINTS
      USE YOWSTAT  , ONLY : CDATEA   ,IDELWI   ,LADEN    ,LGUST    ,
     &            NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWNDG  , ONLY : DLAM     ,DPHI     ,RLATS    ,RLATN    ,
     &            RLONL    ,RLONR    ,KCOL     ,KROW     ,ICODE    ,
     &            IWPER    ,ICOORD
      USE YOWWIND  , ONLY : NC       ,NR       ,UWND     ,VWND     ,
     &            AIRD     ,PZIDL    ,IUNITW   ,CWDFILE
      USE YOWPCONS , ONLY : RAD      ,ZMISS    ,
     &                ROAIR,Rconstant,EpsWaterVapor

! --------------------------------------------------------------------- 

      INTEGER, ALLOCATABLE :: KGRIB(:)

      INTEGER NFIELDS, NGPTOTG, NCA, NRA
      INTEGER ISEC1(NKSEK1), ISEC2(22+NRA)
      INTEGER, ALLOCATABLE, SAVE :: RLONRGG(:), IKBLOCK(:,:)
! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA,I,JF
! End Mod for OPENMP

      REAL FIELDS(NGPTOTG,NFIELDS)
      INTEGER MASK_IN(NGPTOTG)
      LOGICAL, SAVE :: LLINTERPOL

      INTEGER, ALLOCATABLE, SAVE :: KK(:), II(:,:), IIP(:,:)
      REAL, ALLOCATABLE, SAVE :: DK1(:), DII1(:,:), DIIP1(:,:)
      REAL, DIMENSION(NFIELDS) :: FAVG,FMIN,FMAX

      CHARACTER FILNM*24
      CHARACTER*14 CDTWIR, CDTTURN
      CHARACTER*14, SAVE :: CWDDATE
      CHARACTER*16 :: AFLABEL(NFIELDS)
      LOGICAL LLNOTOPENED, LWNDFILE, LCLOSEWND

      LOGICAL FRSTMATM
      REAL, ALLOCATABLE :: AirDenLocal (:)

      DATA FRSTMATM / .TRUE. /
      SAVE IPERIODIC
      SAVE FRSTMATM
! --------------------------------------------------------------------  

!     for the swamp case create your own wind fields
!     or if the file windforcing_time_series is present then
!     the wind forcing will be extracted from the wind speed
!     and direction time series specified herewith.
      IF(CLDOMAIN.EQ.'s' .OR. LWDINTS) goto 8000

      NC=NGX
      NR=NGY

 
      IF(LWNDFILE) THEN

!       WIND FIELDS ARE SUPPLIED FROM INPUT FILE sfcwindin
!       AS IS THE CASE IN THE STAND ALONE VERSION OR AT INITIAL INPUT
!       TIME IN COUPLED MODEL.

!       1.0 INITIALISE VALUES. 
!           ------------------ 

        NBIT=520000

        ICOORD = 1

        FILNM='sfcwindin'
        IF (LLNOTOPENED) THEN
          CALL PBOPEN(IUNITW,FILNM,'r',KRET)
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

          LLNOTOPENED = .FALSE.
        ENDIF  


! --------------------------------------------------------------------  

!*      2.0 READ MARS U and V WIND COMPONENTS  IN GRIB CODE FORMAT.
!           AND SEA ICE FACTION (OR SST) IF LICERUN=TRUE
!           -------------------------------------------------------

        ISIZE=NBIT

        IWPER = 1

        AIRD=ROAIR ! Initializes AirDensity field by the default value.
                ! Actual values can be computed here as well but needs
                ! to get the corresponding fields in the file sfcwindin.
        PZIDL = 0.

        IF(LICERUN) THEN
          NFLD=3
        ELSE
          CI = 0.
          NFLD=2
          IPARAMCI=31
        ENDIF

        WND: DO IVAR=1,NFLD

2001      KPLENG=ISIZE*KIND(ISIZE)
          IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))

          CALL PBGRIB(IUNITW,KGRIB,KPLENG,KLEN,KRET)

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
            CALL KGRIBSIZE(IU06, IUNITW, KLEN, ISIZE, 'READWIND')
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

            IF (IPARAM.EQ.165) THEN
              ICODE = 3
            ELSEIF (IPARAM.EQ.33) THEN
              ICODE = 3
            ELSEIF (IPARAM.EQ.131) THEN
              ICODE = 3
            ELSEIF (IPARAM.EQ.180) THEN
              ICODE = 2
            ELSE
              WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
              WRITE(IU06,*) ' +                                       +'
              WRITE(IU06,*) ' +    WARNING ERROR IN SUB. READWIND     +'
              WRITE(IU06,*) ' +    ==============================     +'
              WRITE(IU06,*) ' + SUSPICIOUS U_WIND FIELD TYPE IN MARS  +'
              WRITE(IU06,*) ' + TYPE IS = ', IPARAM
              WRITE(IU06,*) ' +                                       +'
              WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
              CALL ABORT1
            ENDIF

          ELSE IF (IVAR.EQ.2) THEN
            CALL GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                     NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                     KGRIB, KLENG ,
     &                     NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                     ZDUM, ZDUM,
     &                     CDTWIR, IFORP, IPARAM, KZLEV, IDM, IDM, VWND)

            IF (IPARAM.EQ.166) THEN
              ICODE = 3
            ELSEIF (IPARAM.EQ.34) THEN
              ICODE = 3
            ELSEIF (IPARAM.EQ.132) THEN
              ICODE = 3
            ELSEIF (IPARAM.EQ.181) THEN
              ICODE = 2
            ELSE
              WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
              WRITE(IU06,*) ' +                                       +'
              WRITE(IU06,*) ' +    WARNING ERROR IN SUB. READWIND     +'
              WRITE(IU06,*) ' +    ==============================     +'
              WRITE(IU06,*) ' + SUSPICIOUS V_WIND FIELD TYPE IN MARS  +'
              WRITE(IU06,*) ' + TYPE IS = ', IPARAM
              WRITE(IU06,*) ' +                                       +'
              WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
              CALL ABORT1
            ENDIF

          ELSE
            CALL GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                     NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                     KGRIB, KLENG ,
     &                     NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                     ZDUM, ZDUM,
     &                     CDTWIR, IFORP, IPARAMCI, KZLEV, IDM, IDM, CI)

            IF (IPARAMCI.NE.31 .AND. IPARAMCI.NE.139) THEN
              WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
              WRITE(IU06,*) ' +                                       +'
              WRITE(IU06,*) ' +    WARNING ERROR IN SUB. READWIND     +'
              WRITE(IU06,*) ' +    ==============================     +'
              WRITE(IU06,*) ' + SUSPICIOUS SEA ICE FIELD TYPE IN MARS +'
              WRITE(IU06,*) ' + TYPE IS = ', IPARAMCI
              WRITE(IU06,*) ' +                                       +'
              WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
              CALL ABORT1
            ENDIF

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
          WRITE(IU06,66) ((CI(I,J),I=1,24),J=1,5)
 66       FORMAT(1X,24F5.2)
          CALL FLUSH(IU06)
        ENDIF

! --------------------------------------------------------------------  
                                                                        
        DEALLOCATE(KGRIB)

      ELSE


!*      6.0 READ ATMOS. MODEL FIELDS FOR COUPLED MODEL
!           ------------------------------------------

        IF (FRSTMATM) THEN
!        
!         6.1 FIRST TIME ONLY, COMPUTE INTERPOLATION COEFFICIENTS:
!             ----------------------------------------------------
!
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

          IF (.NOT.ALLOCATED(IKBLOCK)) ALLOCATE(IKBLOCK(0:NCA+1,NRA))

          L = 0
          DO K=1,NRA
            JSN=NRA-K+1
            DO I=1,RLONRGG(JSN)
              L = L+1
              IKBLOCK(I,K) = L
            ENDDO
          ENDDO
          IF(NGPTOTG.NE.L)THEN
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' ####  WARNING IN  READWIND WARNING  ###'
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' ##  NGPTOTG  IS  NOT  EQUAL TO  L  ###'
            WRITE(IU06,*)' ##  NGPTOTG = ', NGPTOTG
            WRITE(IU06,*)' ##        L = ', L
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' ##  NCA = ', NCA, '   NRA = ', NRA
            WRITE(IU06,*)' ##  SUM RLONRGG = ',SUM(RLONRGG(1:NRA))
            WRITE(IU06,*)' ##  RLONRGG = ', RLONRGG(1:NRA)
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' #######################################'
            CALL FLUSH(IU06)
          ENDIF

          IF(IPERIODIC.EQ.1) THEN
            DO K=1,NRA
            JSN=NRA-K+1
            IKBLOCK(0,K)= IKBLOCK(RLONRGG(JSN),K)
            IKBLOCK(RLONRGG(JSN)+1,K)= IKBLOCK(1,K)
            ENDDO
          ENDIF
     
          FRSTMATM = .FALSE.
          
          IF (ITEST.GT.1)
     &     WRITE(IU06,*) ' WAM_READWIND - INTERPOLATION COEFFICIENTS'
     &                     //' COMPUTED'
     
          IF (LADEN.AND.NFIELDS.LT.3) THEN 
            LADEN = .FALSE.
            WRITE(IU06,*) ' WAM_READWIND - WARNING: '
            WRITE(IU06,*) '      AIR DENSITY RUN WAS REQUESTED BUT '
     &                 //'NOT ENOUGH FIELDS PASSED.'
            WRITE(IU06,*) '      NFIELDS = ', NFIELDS
            WRITE(IU06,*) '      NFIELDS MUST BE AT LEAST 3'
            WRITE(IU06,*) '      RESET  LADEN  TO  .FALSE.'
          ENDIF
          
          IF (LADEN) THEN 
            WRITE(IU06,*) ' WAM_READWIND - VARIABLE AIR DENSITY RUN'
          ELSE
            WRITE(IU06,*) ' WAM_READWIND - AIR DENSITY IS '
     &                     //' CONSTANT = ', ROAIR
          ENDIF
          
          IF (LGUST.AND.NFIELDS.LT.4) THEN 
            LGUST = .FALSE.
            WRITE(IU06,*) ' SUB. READWIND - WARNING: '
            WRITE(IU06,*) '      GUSTINESS RUN WAS REQUESTED BUT '
     &                 //'NOT ENOUGH FIELDS PASSED.'
            WRITE(IU06,*) '      NFIELDS = ', NFIELDS
            WRITE(IU06,*) '      NFIELDS MUST BE AT LEAST 4'
            WRITE(IU06,*) '      RESET  LGUST  TO  .FALSE.'
          ENDIF
          
          IF (LGUST) THEN 
            WRITE(IU06,*) ' WAM_READWIND - GUSTINESS RUN'
          ELSE
            WRITE(IU06,*) ' WAM_READWIND - NO GUSTINESS EFFECT.'
          ENDIF
        ENDIF
!        
!       COMPUTATION OF THE NORMS OF INPUT FIELDS
!
        AFLABEL(1)='neutral 10m U   '
        AFLABEL(2)='neutral 10m V   '
        IF (NFIELDS.GE.3) AFLABEL(3)='air density     '
        IF (NFIELDS.GE.4) AFLABEL(4)='Zi/L            '
        IF (NFIELDS.GE.5) AFLABEL(5)='sea ice fraction'

!       COMPUTATION OF THE NORMS OF INPUT FIELDS
        IF( LWCOUNORMS )THEN
          WRITE(IU06,*) ' NORM OF INPUT ATMOSPHERIC FIELDS : ' 
          CALL GSTATS(1436,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JF,I)
          DO IFLD=1,NFIELDS
            JF=IFLD
            FAVG(JF)=FIELDS(1,JF)
            FMIN(JF)=FIELDS(1,JF)
            FMAX(JF)=FIELDS(1,JF)
            DO I=2,NGPTOTG 
              FAVG(JF)=FAVG(JF)+FIELDS(I,JF)
              FMIN(JF)=MIN(FMIN(JF),FIELDS(I,JF))
              FMAX(JF)=MAX(FMAX(JF),FIELDS(I,JF))
            ENDDO
            FAVG(JF)=FAVG(JF)/NGPTOTG
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1436,1)
          DO IFLD=1,NFIELDS
            WRITE(IU06,*) AFLABEL(IFLD),
     &                    FAVG(IFLD),FMIN(IFLD),FMAX(IFLD),NGPTOTG
          ENDDO
          WRITE(IU06,*) ' '
          IF(ITEST.GT.0) CALL FLUSH(IU06)
        ENDIF

        IWPER = 1
        ICOORD = 1

! Mod for OPENMP
        NPROMA=(KSUP-KINF+1)/NWAM_BLKS+1
        CALL GSTATS(1437,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=KINF, KSUP,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,KSUP)

          CALL FLDINTER (IU06, ITEST, NGPTOTG, NCA, NRA, NFIELDS,FIELDS,
     &                   NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                   IXLG, KXLT, NIBLO, IGL, KIJS, KIJL,
     &                   AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                   RLONRGG, IKBLOCK, UWND, VWND, AIRD, PZIDL,
     &                   CI, ZMISS,
     &                   LADEN, ROAIR, LGUST,LLINTERPOL,
     &                   DK1, DII1, DIIP1, KK, II, IIP, MASK_IN)
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1437,1)

        IF (LADEN) THEN 
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. READWIND - AIR DENSITY COMPUTED'
     &                     //' AND INTERPOLATED'
        ELSE
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. READWIND - AIR DENSITY USED AS '
     &                    //' CONSTANT = ', ROAIR
        ENDIF

        IF (LGUST) THEN 
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. READWIND - Zi/L INTERPOLATED'
        ENDIF


! Decode the date of the fields:
!       START DATE.

        IYYYY=(ISEC1(21)-1)*100+ISEC1(10)
        WRITE(CDTWIR(1:12),'(I4.4,4I2.2)') IYYYY,ISEC1(11),ISEC1(12),
     &                             ISEC1(13),ISEC1(14)
        CDTWIR(13:14)='00'


!       DETERMINE FORECAST PERIOD

        IFORP = NINT(FCRANGE)

!       GET DATE OF THE FORECAST INSTEAD OF STARTING DATE

        CALL INCDATE (CDTWIR,IFORP)


        IPARAMU=165
        IPARAMV=166
        IPARAMCI=31

!*        6.1.2 DETERMINE CODE FOR WIND FIELD TYPE:
!               165 = U10 FIELD , 180 = STRESS FIELD.
!               ------------------------------------

        IF (IPARAMU.EQ.165 .AND. IPARAMV.EQ.166) THEN
          ICODE = 3
        ELSEIF (IPARAMU.EQ.33 .AND. IPARAMV.EQ.34) THEN
          ICODE = 3
        ELSEIF (IPARAMU.EQ.131 .AND. IPARAMV.EQ.132) THEN
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
          WRITE(IU06,*) ' READWIND -  FIELDS FOR THE CDTWIR = ',
     &     CDTWIR
          DO IG=1,IGL
            DO IJ = KINF,KSUP
              I = IXLG(IJ,IG)
              K = NGY-KXLT(IJ,IG) + 1
              WRITE(IU06,63) I,K,
     &                       UWND(I,K),VWND(I,K),AIRD(I,K),PZIDL(I,K),
     &                       CI(I,K)
            ENDDO
          ENDDO
          CALL FLUSH(IU06)
 63       FORMAT(1X,2I4,1X,5F5.2)
        ENDIF

      ENDIF

      RETURN


!     SWAMP CASE : CREATE YOUR OWN WIND FIELD

8000  CONTINUE

!     THE DIMENSIONS ARE HARD CODED !!
      NC=NGX
      NR=NGY

      IF (LLNOTOPENED) THEN
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

         IF(LWDINTS) THEN
           IUNITW=I_GET_UNIT(IU06,CWDFILE, 'r', 'f', 0)
           OPEN(IUNITW,FILE=CWDFILE,FORM='FORMATTED')
           CWDDATE=CDATEA
         ELSE
           IUNITW=-1
         ENDIF

         LLNOTOPENED = .FALSE.

      ENDIF

      CALL INCDATE(CDTWIR,IDELWI)

      IF(LWDINTS) THEN

         READ(IUNITW,100,END=110,ERR=120) IWTIME,WSPEED,WTHETA
100      FORMAT(i8,1x,f6.2,1x,f6.1)

         CALL INCDATE(CWDDATE,IWTIME)
         WRITE(IU06,'(a28,a14,1x,f6.2,1x,f6.1)')
     &   '  WIND INPUT TIME SERIES AT ', CWDDATE,WSPEED,WTHETA
  
         WTHETA=RAD*WTHETA
         UWIND=-WSPEED*SIN(WTHETA)
         VWIND=-WSPEED*COS(WTHETA)
      ELSE
!       THE TRADITIONAL ONE GRID POINT FORCING
        CDTTURN = CDATEA
        IDTTURN=DTNEWWIND*3600
        CALL INCDATE(CDTTURN,IDTTURN)
        IF(SWAMPWIND2.LE.0. .OR. CDTWIR.LE.CDTTURN) THEN
          WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',
     &                  SWAMPWIND 
          UWIND=0.0
          VWIND=SWAMPWIND
        ELSE
          WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',
     &                  SWAMPWIND2 
          IF(LTURN90) THEN
            UWIND=SWAMPWIND2
            VWIND=0.0
          ELSE
            UWIND=0.0
            VWIND=SWAMPWIND2
          ENDIF
        ENDIF
      ENDIF

      IPARAMCI=31
      DO J=1,NR
        DO I=1,NC
          UWND(I,J)=UWIND
          VWND(I,J)=VWIND
        ENDDO
      ENDDO


      DO J=1,NR
        DO I=1,NC
          AIRD(I,J)=ROAIR
          PZIDL(I,J)=0.
          CI(I,J)=0.
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

110   CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*  ERROR IN READWIND:                  *'
      WRITE(IU06,*) '*  END OF FILE REACHED FOR ',CWDFILE
      WRITE(IU06,*) '*  THE WIND TIME SERIES SHOULD AT LEAST*'
      WRITE(IU06,*) '*  BE AS LONG AS THE INTENDED RUN !!!! *'
      WRITE(IU06,*) '*  THE RUN HAS ABORTED                 *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1

120   CONTINUE
      WRITE(IU06,*) '****************************************'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '*  ERROR IN READWIND:                  *'
      WRITE(IU06,*) '*  ERROR READING FILE ',CWDFILE
      WRITE(IU06,*) '*  THE RUN HAS ABORTED                 *'
      WRITE(IU06,*) '*                                      *'
      WRITE(IU06,*) '****************************************'
      CALL ABORT1


      END SUBROUTINE READWIND
