! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WGRIBENCODE_VALUES ( I1, I2, &
&                               FIELD, &
&                               ITABPAR, LLSPECNOT251, &
&                               PPMISS, PPEPS, PPREC, PPRESOL, PPMIN_RESET, NTENCODE, &
&                               NGRBRESS, LPADPOLES, &
&                               NLONRGG_SIZE, NLONRGG, IRGG, &
&                               AMONOP, AMOSOP, XDELLA, CLDOMAIN, &
&                               ZMISS, &
&                               VALUES )


! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIBHD, ONLY : NTRG2TMPD, NTRG2TMPP, LLRSTGRIBPARAM

      USE YOWGRIB  , ONLY : IGRIB_GET_VALUE, IGRIB_SET_VALUE, JPGRIB_SUCCESS
      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR
      USE OML_MOD  , ONLY : OML_GET_MAX_THREADS

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: I1, I2
      INTEGER(KIND=JWIM), INTENT(IN) :: ITABPAR
      LOGICAL, INTENT(IN) :: LLSPECNOT251

      REAL(KIND=JWRB), INTENT(INOUT) :: FIELD(I1,I2)

      ! From yowgribhd
      REAL(KIND=JWRB), INTENT(IN)      :: PPMISS  ! every spectral values less or equal ppmiss are replaced by the missing data indicator
      REAL(KIND=JWRB), INTENT(IN)      :: PPEPS   ! Small number used in spectral packing of 140251
      REAL(KIND=JWRB), INTENT(IN)      :: PPREC   ! Reference value for spectral packing of 140251
      REAL(KIND=JWRB), INTENT(IN)      :: PPRESOL ! Maximun resolution possible when encoding spectra (parameter 140251).
      REAL(KIND=JWRB), INTENT(IN)      :: PPMIN_RESET      ! Can be used to set the minimum of ppmin in wgribout to a lower value.
      INTEGER(KIND=JWIM), INTENT(IN)   :: NTENCODE         ! Total number of grid points for encoding
      INTEGER(KIND=JWIM), INTENT(IN)   :: NGRBRESS         ! Number of bits used to encode spectra
      LOGICAL, INTENT(IN)   :: LPADPOLES        ! True if poles are padded when savind to grib.

      ! From yowgrid
      INTEGER(KIND=JWIM), INTENT(IN)   :: NLONRGG_SIZE
      INTEGER(KIND=JWIM), INTENT(IN)   :: NLONRGG(NLONRGG_SIZE)

      ! From yowmap
      INTEGER(KIND=JWIM), INTENT(IN) :: IRGG                ! Grid code: 0 = regular, 1 = irregular.
      REAL(KIND=JWRB),    INTENT(IN) :: AMONOP              ! Most Northern latitude in grid (degree).
      REAL(KIND=JWRB),    INTENT(IN) :: AMOSOP              ! Most Southern latitude in grid (degree).
      REAL(KIND=JWRB),    INTENT(IN) :: XDELLA              ! Grid increment for latitude (degree).

      ! From yowparam
      CHARACTER(LEN=1), INTENT(IN) :: CLDOMAIN   ! Defines the domain of the model (for the fdb and for selection of some variables)

      ! From yowpcons
      REAL(KIND=JWRB), INTENT(IN)        :: ZMISS           ! Missing data indicator (set in chief or via the ifs).

      REAL(KIND=JWRB), DIMENSION(:), INTENT(INOUT) :: VALUES


! ----------------------------------------------------------------------
      INTEGER(KIND=JWIM) :: ICOUNT, NN, I, J, JSN
      INTEGER(KIND=JWIM) :: NPROMA, MTHREADS, JC, JCS, JCL, JJ, ITHRS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS

      REAL(KIND=JWRB) :: TEMP
      REAL(KIND=JWRB) :: ZMINSPEC, PPMAX, PPMIN, DELTAPP, ABSPPREC



      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE :: VALM(:)

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('WGRIBENCODE_VALUES',0,ZHOOK_HANDLE)

      MTHREADS=OMP_GET_MAX_THREADS()
      NPROMA=NTENCODE/MTHREADS + 1


!*    0. PUT FIELD INTO GLOBAL MATRIX VALUES.
!        -----------------------------------

      IF (IRGG == 1 .OR. CLDOMAIN == 'm' ) THEN
        ICOUNT=1
      ELSEIF (CLDOMAIN == 's' ) THEN
        ICOUNT=1
      ELSE
        ICOUNT = (NINT((90._JWRB - AMONOP ) / XDELLA))*I1 + 1
        VALUES(:)=ZMISS
      ENDIF

!     PAD THE POLES IF INCLUDED IN THE GRID
      IF (LPADPOLES) THEN
        IF ((NINT((90._JWRB - AMONOP ) / XDELLA)) == 0) THEN
          TEMP=0._JWRB
          NN=0
          J=2
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            IF (FIELD(I,J) /= ZMISS) THEN
              TEMP=TEMP+FIELD(I,J)
              NN=NN+1
            ENDIF
          ENDDO
          IF (NN > 0) THEN
            TEMP=TEMP/NN
          ELSE
            TEMP=ZMISS
          ENDIF
          J=1
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            FIELD(I,J)=TEMP
          ENDDO
        ENDIF
        IF ((NINT((-90._JWRB - AMOSOP ) / XDELLA)) == 0) THEN
          TEMP=0._JWRB
          NN=0
          J=I2-1
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            IF (FIELD(I,J) /= ZMISS) THEN
              TEMP=TEMP+FIELD(I,J)
              NN=NN+1
            ENDIF
          ENDDO
          IF (NN > 0) THEN
            TEMP=TEMP/NN
          ELSE
            TEMP=ZMISS
          ENDIF
          J=I2
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            FIELD(I,J)=TEMP
          ENDDO
        ENDIF
      ENDIF

!     FILL THE ENCODING ARRAY
!     -----------------------
      DO J = 1, I2
        JSN = I2-J+1
        VALUES(ICOUNT:ICOUNT+NLONRGG(JSN)-1)=FIELD(1:NLONRGG(JSN), J)
        ICOUNT=ICOUNT+NLONRGG(JSN)
      ENDDO

      IF (ITABPAR == 140251 .AND. .NOT. LLSPECNOT251 ) THEN

!       spectra are encoded on a log10 scale and small values below a certain threshold are set to missing
        DELTAPP=(2**NGRBRESS-1)*PPRESOL
        ABSPPREC=ABS(PPREC)
        ZMINSPEC = LOG10(PPEPS)+ABSPPREC

        ALLOCATE(VALM(MTHREADS))
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JC, JCS, JCL, JJ, ITHRS)
        DO JC= 1, NTENCODE, NPROMA
          JCS = JC
          JCL = MIN(JCS+NPROMA-1, NTENCODE)
          ITHRS = JC/NPROMA + 1
          VALM(ITHRS) = ZMINSPEC
          DO JJ = JCS, JCL
            IF(VALUES(JJ) /= ZMISS ) THEN
              VALUES(JJ) = LOG10(VALUES(JJ)+PPEPS)+ABSPPREC
              VALM(ITHRS) = MAX(VALM(ITHRS), VALUES(JJ))
            ELSE
              VALUES(JJ) = ZMINSPEC
            ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        PPMAX=MAXVAL(VALM(:))
        DEALLOCATE(VALM)

        PPMIN=MIN(PPMISS, PPMAX-DELTAPP)
        PPMIN=MIN(PPMIN, PPMIN_RESET)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JC, JCS, JCL, JJ)
        DO JC= 1, NTENCODE, NPROMA
          JCS = JC
          JCL = MIN(JCS+NPROMA-1, NTENCODE)
          DO JJ = JCS, JCL
            IF ( VALUES(JJ) <= PPMIN ) VALUES(JJ)=ZMISS
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

      ENDIF

      IF (LHOOK) CALL DR_HOOK('WGRIBENCODE_VALUES',1,ZHOOK_HANDLE)

END SUBROUTINE WGRIBENCODE_VALUES
