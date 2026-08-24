! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

! ======================================================================

MODULE YOWDATE_UTILS

! ======================================================================

!****     *YOWDATE_UTILS* - MODULE FOR DATE TIME GROUP UTILITIES

!     J. BIDLOT   FEB 2007    RE-INRODUCING THE OLD WAY WITHOUT
!                             THE NEED FOR ECLIB. ADDING SECONDS
!     REFACTORED  MAR 2026    COMBINES DATE INCREMENT AND DIFFERENCE
!                             UTILITIES WITH GENERIC 32/64 BIT SUPPORT

!**   PURPOSE.
!     --------
!       PROVIDES GENERIC INTERFACES FOR UPDATING DATE TIME GROUPS
!       AND COMPUTING DATE DIFFERENCES WITH BOTH 32-BIT AND 64-BIT
!       INTEGER ARGUMENTS.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWIB, JWRB, JWRU
      USE YOWABORT, ONLY : WAM_ABORT
      USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: INCDATE, DIFDATE

      INTERFACE INCDATE
        MODULE PROCEDURE INCDATE_JWIM, INCDATE_JWIB
      END INTERFACE INCDATE

      INTERFACE DIFDATE
        MODULE PROCEDURE DIFDATE_JWIM, DIFDATE_JWIB
      END INTERFACE DIFDATE

      CONTAINS

! ======================================================================

      SUBROUTINE INCDATE_JWIM (CDATE, ISHIFT)

! ======================================================================

!****     *INCDATE_JWIM* - UPDATE DATE TIME GROUP (32-BIT VERSION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(INOUT) :: CDATE
      INTEGER(KIND=JWIM), INTENT(IN) :: ISHIFT

      INTEGER(KIND=JWIM) :: IL, IRET
      INTEGER(KIND=JWIM) :: IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
      INTEGER(KIND=JWIM) :: MON(12)

      CHARACTER(LEN=80) :: CLFMT

      LOGICAL :: LLND

      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/

! ----------------------------------------------------------------------

#define ISHIFT_KIND JWIM
#include "incdate_template.inc"
#undef ISHIFT_KIND

      RETURN

      END SUBROUTINE INCDATE_JWIM

! ======================================================================

      SUBROUTINE INCDATE_JWIB (CDATE, ISHIFT)

! ======================================================================

!****     *INCDATE_JWIB* - UPDATE DATE TIME GROUP (64-BIT VERSION)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(INOUT) :: CDATE
      INTEGER(KIND=JWIB), INTENT(IN) :: ISHIFT

      INTEGER(KIND=JWIM) :: IL, IRET
      INTEGER(KIND=JWIM) :: IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
      INTEGER(KIND=JWIM) :: MON(12)

      CHARACTER(LEN=80) :: CLFMT

      LOGICAL :: LLND

      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/

! ----------------------------------------------------------------------

#define ISHIFT_KIND JWIB
#include "incdate_template.inc"
#undef ISHIFT_KIND

      RETURN

      END SUBROUTINE INCDATE_JWIB

! ======================================================================

      SUBROUTINE DIFDATE_JWIM (CDATE1, CDATE2, KSHIFT)

! ======================================================================

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: CDATE1, CDATE2
      INTEGER(KIND=JWIM), INTENT(OUT) :: KSHIFT

      INTEGER(KIND=JWIB) :: KSHIFT64
      INTEGER(KIND=JWIB) :: IMAXSHIFT

      CALL DIFDATE_CORE(CDATE1, CDATE2, KSHIFT64)

      IMAXSHIFT = INT(HUGE(KSHIFT), KIND=JWIB)
      IF (KSHIFT64 < -IMAXSHIFT .OR. KSHIFT64 > IMAXSHIFT) THEN
        CALL WAM_ABORT('DIFDATE: DATE DIFFERENCE EXCEEDS HUGE(JWIM)',   &
     &                 __FILENAME__, __LINE__)
      ENDIF

      KSHIFT = INT(KSHIFT64, KIND=JWIM)

      END SUBROUTINE DIFDATE_JWIM

! ======================================================================

      SUBROUTINE DIFDATE_JWIB (CDATE1, CDATE2, KSHIFT)

! ======================================================================

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: CDATE1, CDATE2
      INTEGER(KIND=JWIB), INTENT(OUT) :: KSHIFT

      CALL DIFDATE_CORE(CDATE1, CDATE2, KSHIFT)

      END SUBROUTINE DIFDATE_JWIB

! ======================================================================

      SUBROUTINE DIFDATE_CORE (CDATE1, CDATE2, KSHIFT)

! ======================================================================

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: CDATE1, CDATE2
      INTEGER(KIND=JWIB), INTENT(OUT) :: KSHIFT

      INTEGER(KIND=JWIM) :: M, ISI, IDUM, IRET
      INTEGER(KIND=JWIM) :: IYEAR1, IMON1, IDAY1, IHOUR1, IMIN1, ISEC1
      INTEGER(KIND=JWIM) :: IYEAR2, IMON2, IDAY2, IHOUR2, IMIN2, ISEC2
      INTEGER(KIND=JWIM) :: MON(12)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=14) :: CDT1, CDT2

      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('DIFDATE',0,ZHOOK_HANDLE)

      CALL PARSE_DIFDATE_INPUT(CDATE1, CDT1, IYEAR1, IMON1, IDAY1,      &
     &                        IHOUR1, IMIN1, ISEC1, 'DATE1')
      CALL PARSE_DIFDATE_INPUT(CDATE2, CDT2, IYEAR2, IMON2, IDAY2,      &
     &                        IHOUR2, IMIN2, ISEC2, 'DATE2')

      IRET=0
      MON(2)=MFEB_LENGTH(IYEAR1)
      IF (IMON1 < 1 .OR. IMON1 > 12) IRET=-1
      IF (IDAY1 < 1 .OR. IDAY1 > MON(MIN(MAX(IMON1,1),12))) IRET=-2
      IF (IHOUR1 < 0 .OR. IHOUR1 > 23) IRET=-3
      IF (IMIN1 < 0 .OR. IMIN1 > 59) IRET=-4
      IF (ISEC1 < 0 .OR. ISEC1 > 59) IRET=-5

      MON(2)=MFEB_LENGTH(IYEAR2)
      IF (IMON2 < 1 .OR. IMON2 > 12) IRET=-1
      IF (IDAY2 < 1 .OR. IDAY2 > MON(MIN(MAX(IMON2,1),12))) IRET=-2
      IF (IHOUR2 < 0 .OR. IHOUR2 > 23) IRET=-3
      IF (IMIN2 < 0 .OR. IMIN2 > 59) IRET=-4
      IF (ISEC2 < 0 .OR. ISEC2 > 59) IRET=-5

      IF (IRET /= 0) THEN
        WRITE(6, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '("           INPUT DATE(S) INCORRECT,")')
        WRITE(6, '("           IYEAR1 ", I4 )') IYEAR1
        WRITE(6, '("           IMON1 ", I2 )') IMON1
        WRITE(6, '("           IDAY1 ", I2 )') IDAY1
        WRITE(6, '("           IHOUR1 ", I2 )') IHOUR1
        WRITE(6, '("           IMIN1 ", I2 )') IMIN1
        WRITE(6, '("           ISEC1 ", I2 )') ISEC1
        WRITE(6, '("           IYEAR2 ", I4 )') IYEAR2
        WRITE(6, '("           IMON2 ", I2 )') IMON2
        WRITE(6, '("           IDAY2 ", I2 )') IDAY2
        WRITE(6, '("           IHOUR2 ", I2 )') IHOUR2
        WRITE(6, '("           IMIN2 ", I2 )') IMIN2
        WRITE(6, '("           ISEC2 ", I2 )') ISEC2
        WRITE(*, '(" DIFDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '("           INPUT DATE(S) INCORRECT,")')
        WRITE(*, '("           IYEAR1 ", I4 )') IYEAR1
        WRITE(*, '("           IMON1 ", I2 )') IMON1
        WRITE(*, '("           IDAY1 ", I2 )') IDAY1
        WRITE(*, '("           IHOUR1 ", I2 )') IHOUR1
        WRITE(*, '("           IMIN1 ", I2 )') IMIN1
        WRITE(*, '("           ISEC1 ", I2 )') ISEC1
        WRITE(*, '("           IYEAR2 ", I4 )') IYEAR2
        WRITE(*, '("           IMON2 ", I2 )') IMON2
        WRITE(*, '("           IDAY2 ", I2 )') IDAY2
        WRITE(*, '("           IHOUR2 ", I2 )') IHOUR2
        WRITE(*, '("           IMIN2 ", I2 )') IMIN2
        WRITE(*, '("           ISEC2 ", I2 )') ISEC2
        CALL WAM_ABORT('DIFDATE: INPUT DATE(S) INCORRECT',              &
     &                 __FILENAME__, __LINE__)
      ENDIF

!*    CHANGE DATE TIME GROUPS TO ENSURE THAT THE SECOND IS LARGER.
!     ------------------------------------------------------------

      ISI = 1
      IF (CDT1 > CDT2) THEN
         ISI  = -1
         IDUM=IYEAR1
         IYEAR1=IYEAR2
         IYEAR2=IDUM
         IDUM=IMON1
         IMON1=IMON2
         IMON2=IDUM
         IDUM=IDAY1
         IDAY1=IDAY2
         IDAY2=IDUM
         IDUM=IHOUR1
         IHOUR1=IHOUR2
         IHOUR2=IDUM
         IDUM=IMIN1
         IMIN1=IMIN2
         IMIN2=IDUM
         IDUM=ISEC1
         ISEC1=ISEC2
         ISEC2=IDUM
      ENDIF

! ----------------------------------------------------------------------

!*    2.0 CORRECT DAYS IN FEBRUARY OF FIRST YEAR FOR LEAP-YEAR.
!         -----------------------------------------------------

      MON(2)=MFEB_LENGTH(IYEAR1)

! ----------------------------------------------------------------------

!*    3.0 COMPUTE TIME DIFFERENCE IN SECONDS.
!         -----------------------------------

!     3.1 DIFFERENCE BETWEEN DAY, HOUR ,MINUTE, SECOND.

      KSHIFT = (((INT(IDAY2-IDAY1, JWIB)*24_JWIB +                      &
     &            INT(IHOUR2-IHOUR1, JWIB))*60_JWIB +                   &
     &            INT(IMIN2-IMIN1, JWIB))*60_JWIB) +                    &
     &         INT(ISEC2-ISEC1, JWIB)

!     3.2 ADD DIFFERENCE FROM MONTH.

      IF (IYEAR2 > IYEAR1) THEN

!     3.2.1 START AND END MONTHS ARE IN DIFFERENT YEARS.

         DO WHILE (IYEAR2 > IYEAR1)
            DO M=IMON1,12
               KSHIFT = KSHIFT + INT(MON(M), JWIB)*86400_JWIB
            ENDDO
            IMON1 = 1
            IYEAR1 = IYEAR1 + 1
            MON(2)=MFEB_LENGTH(IYEAR1)
         ENDDO
         IF (IMON2 > 1) THEN
            MON(2)=MFEB_LENGTH(IYEAR2)
            DO M=1,IMON2-1
               KSHIFT = KSHIFT + INT(MON(M), JWIB)*86400_JWIB
            ENDDO
         ENDIF
      ELSE

!     3.2.2 START AND END MONTHS ARE IN THE SAME YEAR.

         IF (IMON2 > IMON1) THEN
            DO M=IMON1,IMON2-1
               KSHIFT = KSHIFT + INT(MON(M), JWIB)*86400_JWIB
            ENDDO
         ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    4.0 CHANGE SIGN OF DIFFERENCE IF DATES HAVE BEEN EXCHANGED.
!         -------------------------------------------------------

      KSHIFT = KSHIFT*INT(ISI, JWIB)

      IF (LHOOK) CALL DR_HOOK('DIFDATE',1,ZHOOK_HANDLE)

      END SUBROUTINE DIFDATE_CORE

! ======================================================================

      SUBROUTINE PARSE_DIFDATE_INPUT (CDATE, CDTNORM, IYEAR, IMON, IDAY,&
     &                                IHOUR, IMIN, ISEC, CLABEL)

! ======================================================================

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: CDATE, CLABEL
      CHARACTER(LEN=14), INTENT(OUT) :: CDTNORM
      INTEGER(KIND=JWIM), INTENT(OUT) :: IYEAR, IMON, IDAY
      INTEGER(KIND=JWIM), INTENT(OUT) :: IHOUR, IMIN, ISEC

      INTEGER(KIND=JWIM) :: IL, ISHIFT
      INTEGER(KIND=JWIM) :: IYEAR_TMP, IMON_TMP, IDAY_TMP
      INTEGER(KIND=JWIM) :: IHOUR_TMP, IMIN_TMP, ISEC_TMP

      CHARACTER(LEN=80) :: CLFMT

      IL = LEN_TRIM(CDATE)
      IF (IL == 10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)')                                  &
     &     '(" DIFDATE : NON-Y2K COMPLIANT ", A, " IN USE ", A',        &
     &      IL, ' )'
        WRITE(6, CLFMT ) TRIM(CLABEL), CDATE
        WRITE(*, CLFMT ) TRIM(CLABEL), CDATE
        CALL WAM_ABORT('DIFDATE: NON-Y2K COMPLIANT INPUT DATE',         &
     &                 __FILENAME__, __LINE__)
      ELSEIF (IL == 12) THEN
        READ(CDATE,'(I4, 5I2)')IYEAR, IMON, IDAY, IHOUR, IMIN
        ISEC = 0
        CDTNORM = CDATE(1:12)//'00'
      ELSEIF (IL == 14) THEN
        READ(CDATE,'(I4, 5I2)')IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
!       with satellite data it might happen that the seconds were
!       rounded to 60. One needs to increment the date properly.
        IF (ISEC == 60) THEN
          IYEAR_TMP = IYEAR
          IMON_TMP = IMON
          IDAY_TMP = IDAY
          IHOUR_TMP = IHOUR
          IMIN_TMP = IMIN
          ISEC_TMP = ISEC-1
          ISHIFT = 1
          WRITE(CDTNORM,'(I4.4, 5I2.2)')IYEAR_TMP, IMON_TMP, IDAY_TMP,  &
     &                                  IHOUR_TMP, IMIN_TMP, ISEC_TMP
          CALL INCDATE(CDTNORM, ISHIFT)
          READ(CDTNORM,'(I4, 5I2)')IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
        ELSE
          CDTNORM = CDATE
        ENDIF
      ELSEIF (IL == 0) THEN
        WRITE(6, *) ' '
        WRITE(6, *) ' DIFDATE:  FATAL@! ', TRIM(CLABEL),                &
     &              ' IS 0 CHARACTER LONG'
        WRITE(6, *) ' '
        CALL WAM_ABORT('DIFDATE: EMPTY INPUT DATE',                     &
     &                 __FILENAME__, __LINE__)
      ELSE
        WRITE(6, *) ' '
        WRITE(CLFMT,'(A48, A, A13, I2.2, A2)')                          &
     & '(" DIFDATE :  FATAL@! ", A, " IS ",I2," CHARACTERS LONG!!= ",A',&
     &  IL, ' )'
        WRITE(6, CLFMT) TRIM(CLABEL), IL, CDATE
        WRITE(6, *) ' '
        CALL WAM_ABORT('DIFDATE: INVALID INPUT DATE LENGTH',            &
     &                 __FILENAME__, __LINE__)
      ENDIF

      END SUBROUTINE PARSE_DIFDATE_INPUT

! ======================================================================

      INTEGER(KIND=JWIM) FUNCTION MFEB_LENGTH(IYEAR)

! ======================================================================

!     LENGTH OF FEBRUARY IN DAYS FOR YEAR IYEAR (YYYY)

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IYEAR

      IF (MOD(IYEAR,400) == 0) THEN
        MFEB_LENGTH=29
      ELSEIF (MOD(IYEAR,100) == 0) THEN
        MFEB_LENGTH=28
      ELSEIF (MOD(IYEAR,4) == 0) THEN
        MFEB_LENGTH=29
      ELSE
        MFEB_LENGTH=28
      ENDIF

      END FUNCTION MFEB_LENGTH

END MODULE YOWDATE_UTILS