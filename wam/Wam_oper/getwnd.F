      SUBROUTINE GETWND (US, DS, KINF, KSUP, KEND, CDTWIS,
     &                   IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)

! ----------------------------------------------------------------------

!**** *GETWND* - ROUTINE TO READ AND PROCESS ONE WINDFIELD.

!*    PURPOSE.
!     --------

!        READ A WINDFIELD FROM THE WINDFILE (SEARCH FOR IT)
!        AND CALCULATES THE WIND VELOCITY  AND DIRECTION
!        FOR ALL WAM BLOCKS.

!**   INTERFACE.
!     ----------

!       *CALL* *GETWND (US, DS, KINF, KSUP, KEND, CDTWIS,
!                       IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)*
!         *US*     - MAGNITUDE OF USTAR IN EACH POINT AND BLOCK.
!         *DS*     - DIRECTION OF USTAR IN EACH POINT AND BLOCK.
!         *KINF*   - LOWER BOUND ON FIRST DIEMENSION OF FIELD
!         *KSUP*   - UPPER BOUND ON FIRST DIEMENSION OF FIELD
!         *KEND*   - USED UPPER LIMIT ON FIRST DIMENSION OF FIELD
!         *CDTWIS* - DATE OF WIND FIELD TO BE LOOKED FOR.
!         *IUGB*   - INPUT GRIB WIND FIELD U COMPONENT.
!         *IVGB*   - INPUT GRIB WIND FIELD V COMPONENT.
!         *ILEN*   - LENGTH OF IUGB AND IVGB
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

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NBLO
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------
      INTEGER IUGB(ILEN), IVGB(ILEN)
! ----------------------------------------------------------------------

      REAL, DIMENSION (KINF:KSUP,NBLO) :: US, DS
      LOGICAL LWNDFILE, LCLOSEWND

      CHARACTER*12 CDTWIS, CDTWIR

      SAVE CDTWIR

! ----------------------------------------------------------------------

!*    1. WIND DATA ARE READ
!        ------------------

 1000 CONTINUE

      CALL READWIND (CDTWIR, IUGB, IVGB, ILEN, LWNDFILE, LCLOSEWND)

!*    2. CHECK WIND FIELD DATE
!        ----------------------

      IF (CDTWIR.LT.CDTWIS) THEN

!     2.1 DATE OF INPUT FIELD IS BEFORE REQUESTED DATE
!         TRY AGAIN
!         --------------------------------------------
        IF(LWNDFILE) THEN

          IF(LCLOSEWND) THEN
            WRITE (IU06,*) ' ***************************************'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *    FATAL ERROR SUB. GETWND          *'
            WRITE (IU06,*) ' *    =======================          *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' *  THE INPUT FILE HAS BEEN CLOSED     *'
            WRITE (IU06,*) ' *  IT CANNOT BE ACCESSED AGAIN !      *'
            WRITE (IU06,*) ' *  PROGRAM ABORTS  PROGRAM ABORTS     *'
            WRITE (IU06,*) ' *                                     *'
            WRITE (IU06,*) ' ***************************************'
            CALL ABORT1
          ENDIF

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

! ----------------------------------------------------------------------

!*    3. INTERPOLATE AND BLOCK WINDFIELD
!        -------------------------------

      CALL WAMWND (US, DS, KINF, KSUP, KEND)

      RETURN
      END SUBROUTINE GETWND
