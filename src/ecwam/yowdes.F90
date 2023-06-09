! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      MODULE YOWDES

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

!     MODULE FOR SAR INVERSION SOFTWARE

      INTEGER(KIND=JWIM) :: JPMAX1

      INTEGER(KIND=JWIM) :: JPWAM
      INTEGER(KIND=JWIM) :: JPKSAR
      INTEGER(KIND=JWIM) :: JPTSAR
      INTEGER(KIND=JWIM) :: JPISAR
      INTEGER(KIND=JWIM) :: JPNSAR
      INTEGER(KIND=JWIM) :: JPDATA

      INTEGER(KIND=JWIM) :: KSX
      INTEGER(KIND=JWIM) :: KSY
      INTEGER(KIND=JWIM) :: KSXT
      INTEGER(KIND=JWIM) :: KSYT
      INTEGER(KIND=JWIM) :: LOT
      INTEGER(KIND=JWIM) :: KSXLOT
      INTEGER(KIND=JWIM) :: KSXP1
      INTEGER(KIND=JWIM) :: KSXM1
      INTEGER(KIND=JWIM) :: KSYP1
      INTEGER(KIND=JWIM) :: KSYM1
      INTEGER(KIND=JWIM) :: NTABLE
      INTEGER(KIND=JWIM) :: NWORK
      INTEGER(KIND=JWIM) :: NOR

      INTEGER(KIND=JWIM) :: JPSPEC
      INTEGER(KIND=JWIM) :: JPIMAF
      INTEGER(KIND=JWIM) :: JPIM2
      INTEGER(KIND=JWIM) :: JPHALF

      INTEGER(KIND=JWIM) :: IDTMAX

      INTEGER(KIND=JWIM)              :: NXSUP
      INTEGER(KIND=JWIM)              :: NXELE
      INTEGER(KIND=JWIM), ALLOCATABLE :: NXSUPL(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: NTDLST(:)

      REAL(KIND=JWRB),    ALLOCATABLE :: XSDATA(:)
      REAL(KIND=JWRB),    ALLOCATABLE :: XODATA(:)

      REAL(KIND=JWRB), PARAMETER :: XEPS=1.0E-12_JWRB


!*    VARIABLE.   TYPE.     PURPOSE.
!     ---------   -------   --------
!     JPMAX1      INTEGER   NUMBER OF MAXIMA SELECTED
!     JPWAM       INTEGER   MAXIMUM NUMBER OF WAM WAVE SPECTRA PER BLOCK
!     JPKSAR      INTEGER   WAVE NUMBER GRID POINTS OF OBSERVED SAR SPECTRUM
!     JPTSAR      INTEGER   DIRECTION GRID POINTS OF OBSERVED SAR SPECTRUM
!     JPISAR      INTEGER   NUMBER OF GRID POINTS OF A SAR SPECTRUM
!                           JPISAR = JPKSAR * JPTSAR
!     JPNSAR      INTEGER   MAXIMUM NUMBER OF OBSERVED SAR SPECTRA
!     JPDATA      INTEGER   DIMENSION OF THE ARRAY CONTAINING ALL DECODED
!                           SAR DATA RELEVANT FOR THE INVERSION. 
!                           JPDATA = 24 + JPKSAR + JPTSAR + JPISAR + 1
!     KSX         INTEGER   NUMBER OF WAVE NUMBERS IN SAR SPECTRUM(X-D)
!     KSY         INTEGER   NUMBER OF WAVE NUMBERS IN SAR SPECTRUM(Y-D)
!     KSXT        INTEGER   2*KSX
!     KSYT        INTEGER   2*KSY
!     LOT         INTEGER   KSYT
!     KSXLOT      INTEGER   4*KSXT*LOT
!     KSXP1       INTEGER   KSX+1
!     KSXM1       INTEGER   KSX-1
!     KSYP1       INTEGER   KSY+1
!     KSYM1       INTEGER   KSY-1
!     NTABLE      INTEGER   100+2*(KSXT*2)
!     NWORK       INTEGER   4*KSXT*KSYT
!     NOR         INTEGER   DIMENSION OF WORK ARRAYS IN FIMAGEP
!     JPSPEC      INTEGER   SPECTRAL LENGTH (JPSPEC=2*KSX+1)
!     JPIMAF      INTEGER   LENGTH OF THE REAL IMAGE FIELDS (JPSPEC - 1) 
!     JPHALF      INTEGER   JPIMAF/2 + 1
!     IDTMAX      INTEGER   MAXIMUN TIME SEPARATION IN SECONDS BETWEEN
!                           MODEL FIRST GUESS AND SAR SPECTRA.
!     NXSUP       INTEGER   DIMENSION OF NXSUPL
!     NXELE       INTEGER   DIMENSION OF NTDLST
!     NXSUPL      INTEGER   ARRAY USED TO MANIPULATE DECODED SAR DATA
!     NTDLST      INTEGER   ARRAY USED TO MANIPULATE DECODED SAR DATA
!     XSDATA      REAL      ARRAY USED TO MANIPULATE DECODED SAR DATA
!     XODATA      REAL      ARRAY USED TO MANIPULATE DECODED SAR DATA
!     XEPS        REAL      SMALL NUMBER

      END MODULE YOWDES
