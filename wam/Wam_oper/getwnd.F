      SUBROUTINE GETWND (US, DS, ADS, ZIDL, KINF, KSUP, KEND, CDTWIS,
     &                   NFIELDS,NGPTOTG, NC,NR, ISEC1,ISEC2, FIELDS,
     &                   LWNDFILE, LCLOSEWND)

! ----------------------------------------------------------------------

!**** *GETWND* - ROUTINE TO READ AND PROCESS ONE WINDFIELD.

!*    PURPOSE.
!     --------

!        READ A WINDFIELD FROM THE WINDFILE (SEARCH FOR IT)
!        AND CALCULATES THE WIND VELOCITY  AND DIRECTION
!        FOR ALL WAM BLOCKS.

!**   INTERFACE.
!     ----------

!       *CALL* *GETWND (US, DS, ADS, ZIDL, KINF, KSUP, KEND, CDTWIS,
!                       NFIELDS,NGPTOTG, NC,NR, ISEC1,ISEC2, FIELDS,
!                       LWNDFILE, LCLOSEWND)*
!         *US*     - MAGNITUDE OF USTAR IN EACH POINT AND BLOCK.
!         *DS*     - DIRECTION OF USTAR IN EACH POINT AND BLOCK.
!         *ADS*    - AIR DENSITY IN EACH POINT AND BLOCK.
!         *ZIDL*   - Zi/L  IN EACH POINT AND BLOCK.
!         *KINF*   - LOWER BOUND ON FIRST DIEMENSION OF FIELD
!         *KSUP*   - UPPER BOUND ON FIRST DIEMENSION OF FIELD
!         *KEND*   - USED UPPER LIMIT ON FIRST DIMENSION OF FIELD
!         *CDTWIS* - DATE OF WIND FIELD TO BE LOOKED FOR.
!         *NFIELDS*- NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!         *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!         *NC*     - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!         *NR*     - NUMBER OF ATM. ROWS OF LATITUDES
!         *ISEC1*  - GRIB ISEC1 (TO EXTRACT DATE/TIME OF ATM. FIELDS)
!         *ISEC2*  - GRIB ISEC2 (TO EXTRACT ATM. GRID CHARACTERISTICS)
!         *FIELDS* - ATMOSPHERIC FIELDS (U10, V10, Q2, MSL, SKT, ice, Zi/L)
!         *LWNDFILE - FLAG USED TO DETERMINE WHETHER WINDS ARE READ FROM
!                     FILE OR PASSED AS ARGUMENTS IN IUGB AND IVGB
!         *LCLOSEWND* IF TRUE THE INPUT FILE WILL BE CLOSED AND
!                     THE UNIT RESET

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *READWIND*   - READING WINDS.
!       *WAMWND*    - CALCULATE WIND IN WAM POINTS.

!     REFERENCE.
!     ----------

!       NONE.


!    MODIFIED BY:
!    ------------
!    B. HANSEN    ECMWF 1997
!                 RESTRUCTURE CALL TO READWIND.
!
!    S. ABDALLA   ECMWF OCTOBER 2001
!                 MODIFICATION THE CALL TO READWIND WAMWND; AND 
!                 INCLUSION OF AIR DENSITY AND Zi/L.
! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NBLO     ,CLDOMAIN
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------
      INTEGER NFIELDS, NGPTOTG, NC, NR
      INTEGER ISEC1(255), ISEC2(22+NR)
      REAL FIELDS(NGPTOTG,NFIELDS)
! ----------------------------------------------------------------------

      REAL, DIMENSION (KINF:KSUP,NBLO) :: US, DS, ADS, ZIDL
      INTEGER :: IUNIT1
      LOGICAL FRSTIME, LWNDFILE, LCLOSEWND

      CHARACTER*12 CDTWIS, CDTWIR
      CHARACTER*24 :: FILNM

      DATA FRSTIME / .TRUE. /

      SAVE CDTWIR
      SAVE IUNIT1

! ----------------------------------------------------------------------

!*    1. WIND DATA ARE READ
!        ------------------

 1000 CONTINUE

      CALL READWIND (CDTWIR, NFIELDS, NGPTOTG, NC, NR, 
     &               ISEC1, ISEC2, FIELDS,LWNDFILE, LCLOSEWND,
     &               KINF, KSUP,
     &               IUNIT1, FILNM, FRSTIME )

!*    2. CHECK WIND FIELD DATE
!        ----------------------

      IF (CDTWIR.LT.CDTWIS) THEN

!     2.1 DATE OF INPUT FIELD IS BEFORE REQUESTED DATE
!         TRY AGAIN
!         --------------------------------------------
        IF(LWNDFILE) THEN
          GOTO 1000
        ELSE
          WRITE (IU06,*) ' ******************************************'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND            *'
          WRITE (IU06,*) ' *      =======================            *'
          WRITE (IU06,*) ' * WIND DATE IS EARLIER THAN EXPECTED      *'
          WRITE (IU06,*) ' * DECODED DATE IS  CDTWIR = ', CDTWIR
          WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS        *'
          WRITE (IU06,*) ' *                                         *'
          WRITE (IU06,*) ' ******************************************'
          CALL ABORT1
        ENDIF
      ELSEIF (CDTWIR.GT.CDTWIS) THEN

!     2.2 DATE OF INPUT FIELD IS LATER THAN REQUESTED DATE
!         ------------------------------------------------
        WRITE (IU06,*) ' ******************************************'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *      FATAL ERROR SUB. GETWND           *'
        WRITE (IU06,*) ' *      =======================           *'
        WRITE (IU06,*) ' * WIND DATE IS LATER THAN EXPECTED       *'
        IF(LWNDFILE) THEN
          WRITE (IU06,*) ' * DATE READ IS     CDTWIR = ', CDTWIR
        ELSE
          WRITE (IU06,*) ' * DECODED DATE IS  CDTWIR = ', CDTWIR
        ENDIF
        WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
        WRITE (IU06,*) ' *                                        *'
        WRITE (IU06,*) ' ******************************************'
        CALL ABORT1
      ENDIF


      IF (LCLOSEWND .AND. LWNDFILE .AND. (CLDOMAIN.NE.'s') ) THEN
        CALL PBCLOSE(IUNIT1,KRET)
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '****************************************'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBCLOSE     *'
          WRITE (IU06,*) '*  IN GETWND                           *'
          WRITE (IU06,*) '*                                      *'
          WRITE (IU06,*) '****************************************'
          CALL ABORT1
        ENDIF
        FRSTIME = .TRUE.
        IF (ITEST.GT.1)
     &    WRITE(IU06,*) ' SUB. GETWND - PBCLOSE ON ',FILNM
      ENDIF


! ----------------------------------------------------------------------

!*    3. INTERPOLATE AND BLOCK WINDFIELD
!        -------------------------------

      CALL WAMWND (US, DS, ADS, ZIDL, KINF, KSUP, KEND)

      RETURN
      END SUBROUTINE GETWND
