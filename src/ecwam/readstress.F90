! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE READSTRESS(IJINF, IJSUP, NREAL, RFIELD, FILENAME, LRSTPARAL, &
     &                      NNEMO, RNEMO, LDNEMO)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     READS THE RESTART WIND AND STRESS FIELDS

!**   INTERFACE.
!     ----------
!     CALL *READSTRESS*(NREAL, RFIELD, FILENAME)  
!     *NREAL*    NUMBER OF FIELDS IN RFIELD
!     *RFIELD*     REAL RESTART FIELDS
!     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE
!     *NNEMO*    (OPTIONAL) NUMBER OF NEMO RECORDS THAT MAY FOLLOW THE NREAL ONES
!     *RNEMO*    (OPTIONAL) NEMO FIELDS (KIND JWRO), PER-TASK FILES ONLY
!     *LDNEMO*   (OPTIONAL) TRUE IF THE FILE HAD THE NEMO RECORDS

!     METHOD.
!     -------
!     READS ARRAYS FROM FILE. FORTRAN READS

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D     ,LLUNSTR
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL
      USE YOWABORT,  ONLY : WAM_ABORT
#ifdef WAM_HAVE_UNWAM
      USE YOWUNBLKRORD, ONLY : UNBLKRORD
#endif

      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "iwam_get_unit.intfb.h"
      INTEGER(KIND=JWIM), INTENT(IN) :: IJINF, IJSUP, NREAL
      REAL(KIND=JWRB),DIMENSION(IJINF:IJSUP,NREAL),INTENT(OUT) :: RFIELD
      CHARACTER(LEN=296), INTENT(IN) :: FILENAME
      LOGICAL, INTENT(IN) :: LRSTPARAL
      INTEGER(KIND=JWIM), INTENT(IN), OPTIONAL :: NNEMO
      REAL(KIND=JWRO),DIMENSION(IJINF:IJSUP,*),INTENT(OUT), OPTIONAL :: RNEMO
      LOGICAL, INTENT(OUT), OPTIONAL :: LDNEMO


      INTEGER(KIND=JWIM) :: LFILE, IUNIT
      INTEGER(KIND=JWIM) :: IFLD, IJ, IOS

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(IJINF:IJSUP) :: RFIELD_G 

      LOGICAL :: LLEXIST

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('READSTRESS',0,ZHOOK_HANDLE)

      LFILE=0
      IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
      LLEXIST=.FALSE.
      INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
      IF (.NOT. LLEXIST) THEN
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (IU06,*) '*  IN READSTRESS :                  *'
        WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILENAME
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        WRITE (*,*) '*************************************'
        WRITE (*,*) '*                                   *'
        WRITE (*,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (*,*) '*  IN READSTRESS :                  *'
        WRITE (*,*) '*  COULD NOT FIND FILE ',FILENAME
        WRITE (*,*) '*                                   *'
        WRITE (*,*) '*************************************'
        CALL ABORT1
      ENDIF
      IUNIT=IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0, 'READWRITE')

      READ(IUNIT) CDTPRO,CDATEWO,CDAWIFL,CDATEFL


      IF (LLUNSTR .AND. .NOT.LRSTPARAL) THEN
#ifdef WAM_HAVE_UNWAM
        DO IFLD=1,NREAL
          READ(IUNIT)(RFIELD_G(IJ),IJ=IJINF,IJSUP)
          CALL UNBLKRORD(-1,IJINF,IJSUP,1,1,1,1,RFIELD(1,IFLD),RFIELD_G(1))
        ENDDO
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSEIF (LRSTPARAL .OR. LL1D .OR. NPROC == 1) THEN
!     ALL PE'S READ OR 1D DECOMPOSITION 
        DO IFLD=1,NREAL
          READ(IUNIT)(RFIELD(IJ,IFLD),IJ=IJINF,IJSUP)
        ENDDO
      ELSE
!     WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE
!     RE-LABELLED BUT THE BINARY INPUT FILES ARE IN THE OLD MAPPING
        DO IFLD=1,NREAL
          READ(IUNIT)(RFIELD(IJ2NEWIJ(IJ),IFLD),IJ=IJINF,IJSUP)
        ENDDO
      ENDIF

!     OPTIONAL NEMO RECORDS (SEE RECVNEMOFIELDS): ALL OR NONE. A FILE WITHOUT
!     THEM (END OF FILE AT THE FIRST ONE) IS A RESTART FROM BEFORE THEY EXISTED.
      IF (PRESENT(LDNEMO)) THEN
        LDNEMO = .FALSE.
        IF (PRESENT(NNEMO) .AND. PRESENT(RNEMO)) THEN
          DO IFLD=1,NNEMO
            READ(IUNIT, IOSTAT=IOS)(RNEMO(IJ,IFLD),IJ=IJINF,IJSUP)
            IF (IFLD == 1 .AND. IS_IOSTAT_END(IOS)) EXIT
            IF (IOS /= 0) THEN
              WRITE (IU06,*) ' READSTRESS: CANNOT READ NEMO RECORD ', IFLD, ' OF ', NNEMO, &
     &                       ' (IOSTAT=', IOS, ') FROM ', FILENAME(1:LFILE)
              CALL WAM_ABORT("Incomplete or wrong-kind NEMO records in LAW restart file", &
     &                       __FILENAME__, __LINE__)
            ENDIF
            IF (IFLD == NNEMO) LDNEMO = .TRUE.
          ENDDO
        ENDIF
      ENDIF

      CLOSE(IUNIT)

      IF (LHOOK) CALL DR_HOOK('READSTRESS',1,ZHOOK_HANDLE)

      END SUBROUTINE READSTRESS
