      SUBROUTINE READWIND(CDTWIR,IUGB,IVGB,ILEN,lWNDFILE,LCLOSEWND)

!***  *READWIND* - PROGRAM TO REFORMAT WIND FIELDS FROM MARS 
!                 EITHER FROM A FILE OR FROM GRIB DATA IN IUGB, IVGB

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

!     *CALL* *READWIND (CDTWIR, IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)*
!        *CDTWIR* - DATE/TIME OF THE DATA READ                          
!        *IUGB*   - ARRAY CONTAINING THE GRIB U COMPONENT
!        *IVGB*   - ARRAY CONTAINING THE GRIB V COMPONENT
!        *ILEN*   - SIZE OF IUGB AND IVGB
!        *LWNDFILE*  IF TRUE DATA ARE ACCESSED FROM THE INPUT FILE
!                    sfcwind
!                    ELSE THE DATA WILL BE SUPPLIED IN IUGB AND IVGB
!        *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND THE UNIT 
!                    RESET

!     EXTERNALS                                                         
!     ---------                                                         

!     *ABORT1*         TERMINATES PROCESSING
!     *ADJUST*         ADJUSTS EAST AND WEST BOUNDARIES
!     *GRIBEX*         UNPACKS MARS WINDS                               
!     *GRPRS1*         PRINTS PRODUCT DEFINITION BLOCK                  
!     *GRPRS2*         PRINTS GRID  DEFINITION BLOCK                    
!     *INCDATE*        INCREMENTS DATE                                  
!     *PBCLOSE*
!     *PBGRIB*
!     *PBOPEN*

! --------------------------------------------------------------------- 

      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWPARAM , ONLY : NGX      ,NGY      ,SWAMPWIND,CLDOMAIN
      USE YOWSTAT  , ONLY : CDATEA   ,IDELWI 
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWNDG  , ONLY : DLAM     ,DPHI     ,RLATS    ,RLATN    ,
     &            RLONL    ,RLONR    ,KCOL     ,KROW     ,ICODE    ,
     &            IWPER    ,ICOORD
      USE YOWWIND  , ONLY : NC       ,NR       ,UWND     ,VWND  

! --------------------------------------------------------------------- 

      INTEGER, ALLOCATABLE :: KGRIB(:)
      INTEGER IUGB(ILEN), IVGB(ILEN)
      CHARACTER FILNM*24
      CHARACTER*12 CDTWIR
      LOGICAL FRSTIME, LWNDFILE, LCLOSEWND

      DATA FRSTIME / .TRUE. /
      SAVE IUNIT1      

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
            CALL ABORT1
          ENDIF
          FRSTIME = .TRUE.
          IF (ITEST.GT.1)
     &     WRITE(IU06,*) ' SUB. READWIND - PBCLOSE ON ',FILNM
        ENDIF  


      ELSE


!*      6.0 READ ATMOS. MODEL GRIB FIELDS FOR COUPLED MODEL
!           -----------------------------------------------

        IWPER = 1
        ICOORD = 1

        IF(.NOT.ALLOCATED(UWND)) ALLOCATE(UWND(NGX,NGY)) 
        IF(.NOT.ALLOCATED(VWND)) ALLOCATE(VWND(NGX,NGY))

        CALL GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                 NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                 IUGB, ILEN,
     &                 NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                 AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                 ZDUM, ZDUM,
     &                 CDTWIR, IFORP, IPARAMU, KZLEV, IDM, IDM, UWND)
        CALL GRB2WGRD (IU06, ITEST, NKSEK0, NKSEK1, NKSEK2, NKSEK3,
     &                 NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                 IVGB, ILEN,
     &                 NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                 AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                 ZDUM, ZDUM,
     &                 CDTWIR, IFORP, IPARAMV, KZLEV, IDM, IDM, VWND)

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


      ENDIF

      RETURN


!     SWAMP CASE : CREATE YOUR OWN WIND FIELD

8000  CONTINUE
      WRITE(IU06,*) ' READWIND - SWAMP CASE WITH WIND SPEED = ',
     &                SWAMPWIND 

!     THE DIMENSIONS ARE HARD CODED !!
      NC=NGX
      NR=NGY

      IF(ALLOCATED(UWND)) DEALLOCATE(UWND)
      ALLOCATE(UWND(NC,NR))                                
      IF(ALLOCATED(VWND)) DEALLOCATE(VWND)
      ALLOCATE(VWND(NC,NR))                                

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
