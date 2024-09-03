! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

PROGRAM combine_bathy_laked

! ----------------------------------------------------------------------

!**** *COMBINE_BATHY_LAKED* 

!*    PURPOSE.
!     --------

!     TO COMBINE BATHY WITH IFS LAND SEA MASK and LAKE COVER USING LAKE DEPTH DATA

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWABORT , ONLY : WAM_ABORT
USE YOWGRIBINFO, ONLY : WVGETGRIDINFO
USE YOWMAP   , ONLY : NGX, NGY, IPER, IRGG, IQGAUSS,                         &
                   &  DAMOWEP, DAMOSOP, DAMOEAP, DAMONOP, DXDELLA, DXDELLO,  &
                   &  NLONRGG
USE YOWPCONS , ONLY : ZMISS
USE YOWGRIB  , ONLY : IGRIB_GET_VALUE, IGRIB_CLOSE_FILE, IGRIB_RELEASE,       &
                    & IGRIB_SET_VALUE, IGRIB_OPEN_FILE, IGRIB_NEW_FROM_FILE,  &
                    & IGRIB_CLONE, IGRIB_GET_MESSAGE_SIZE, IGRIB_GET_MESSAGE, &
                    & IGRIB_WRITE_BYTES, &
                    & JPGRIB_END_OF_FILE, JPKSIZE_T

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "wvopenbathy.intfb.h"

INTEGER(KIND=JWIM), PARAMETER :: NPARAM=3
INTEGER(KIND=JWIM) :: IC, K, IP, LFILE, IRET, IVAL
INTEGER(KIND=JWIM) :: IU06, IU07, KGRIB_HANDLE_BATHY, KGRIB_HANDLE_NEW_BATHY
INTEGER(KIND=JWIM) :: NUMBEROFVALUES, ISIZE, NPRECI
INTEGER(KIND=JWIM) :: NGX_LAKE, NGY_LAKE, IPER_LAKE, IRGG_LAKE, IQGAUSS_LAKE
INTEGER(KIND=JWIM) :: I4(2)
INTEGER(KIND=JWIM) :: NPROMA, MTHREADS, JC, JCS, JCL, IPR
!$ INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS

INTEGER(KIND=JWIM), DIMENSION(NPARAM) :: IULAKE, KGRIB_HANDLE_LAKE, IPARAMID
INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: INEWLAKE, IAVG
INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: NLONRGG_LAKE
INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: KGRIB_BUFR
INTEGER(KIND=JPKSIZE_T) :: KBYTES

REAL(KIND=JWRU), PARAMETER :: BATHYMAX = 999.0_JWRU !! ecWAM maximum depth
REAL(KIND=JWRU), PARAMETER :: THRSLSM = 0.5_JWRU    !! points with LSM below THRSLSM are assumed to be sea/ocean
REAL(KIND=JWRU), PARAMETER :: THRSLAKE = 0.65_JWRU  !! points with lake cover above THRSLAKE are assumed to be resolved lake

REAL(KIND=JWRU) :: DAMOWEP_LAKE, DAMOSOP_LAKE, DAMOEAP_LAKE, DAMONOP_LAKE, DXDELLA_LAKE, DXDELLO_LAKE
REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:) :: VALUES_BATHY
REAL(KIND=JWRU), ALLOCATABLE, DIMENSION(:,:) ::  VALUES_LAKE

CHARACTER(LEN=80) :: FILENAME, OUTFILENAME
CHARACTER(LEN=80), DIMENSION(NPARAM) :: INFILENAME

LOGICAL :: LLEXIST

LOGICAL :: LLSCANNS, LLSCANNS_LAKE, LLSAMEGRID

! ----------------------------------------------------------------------

! The new bathymetry will be writen out to
OUTFILENAME='new_bathy'


NPRECI = KIND(I4)

IU06 = 6

IU07 = -1
KGRIB_HANDLE_BATHY = -99

IULAKE(:) = -1
KGRIB_HANDLE_LAKE(:) = -99

LLSAMEGRID = .TRUE.


WRITE(IU06,*) ''
WRITE(IU06,*) ' Combining model bathymetry with lsm and lake info'
WRITE(IU06,*) ''

!  INPUT FILE:
!  -----------

! LAND SEA MASK
INFILENAME(1)='lsm'
IPARAMID(1)=172
! LAKE COVER
INFILENAME(2)='clake'
IPARAMID(2)=26
! LAKE DEPTH
INFILENAME(3)='lakedl'
IPARAMID(3)=228007

! BATHY (wam_grid_tables):
CALL WVOPENBATHY (IU06, IU07, KGRIB_HANDLE_BATHY)

DO IP = 1, NPARAM
  FILENAME = INFILENAME(IP)
  LFILE=0
  LLEXIST=.FALSE.
  IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
  INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)

  IF (LLEXIST) THEN
    CALL IGRIB_OPEN_FILE(IULAKE(IP),FILENAME(1:LFILE),'r')
  ELSE
    WRITE(IU06,*) '*****************************************************************'
    WRITE(IU06,*) '*                                                               *'
    WRITE(IU06,*) '*  FATAL ERROR IN                               *'
    WRITE(IU06,*) '*  =============================                                *'
    WRITE(IU06,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE), ' IS MISSING !!!!'
    WRITE(IU06,*) '*                                                               *'
    WRITE(IU06,*) '*****************************************************************'
    CALL WAM_ABORT("INPUT FILE '"//FILENAME(1:LFILE)//"' IS MISSING !!!",__FILENAME__,__LINE__)
  ENDIF

ENDDO


IF ( KGRIB_HANDLE_BATHY > 0 ) THEN

! GRID INFO:
! ---------
  CALL WVGETGRIDINFO(IU06, KGRIB_HANDLE_BATHY, &
 &                   NGX, NGY, IPER, IRGG, IQGAUSS, NLONRGG, LLSCANNS, &
 &                   DAMOWEP, DAMOSOP, DAMOEAP, DAMONOP, DXDELLA, DXDELLO )

! LAKE GRID INFO:
! ---------------
  DO IP = 1, NPARAM

    CALL IGRIB_NEW_FROM_FILE(IULAKE(IP), KGRIB_HANDLE_LAKE(IP), IRET)

    IF (IRET /= JPGRIB_END_OF_FILE) THEN

      CALL WVGETGRIDINFO(IU06, KGRIB_HANDLE_LAKE(IP), &
 &                       NGX_LAKE, NGY_LAKE, IPER_LAKE, IRGG_LAKE, IQGAUSS_LAKE, NLONRGG_LAKE, LLSCANNS_LAKE, &
 &                       DAMOWEP_LAKE, DAMOSOP_LAKE, DAMOEAP_LAKE, DAMONOP_LAKE, DXDELLA_LAKE, DXDELLO_LAKE )

      !! Check that it is the same as input BATHY
      IF ( NGX_LAKE /= NGX .OR. NGY_LAKE /= NGY .OR. IPER_LAKE .NE. IPER .OR. IQGAUSS_LAKE .NE. IQGAUSS .OR. &
 &         DAMOWEP_LAKE /= DAMOWEP .OR. DAMOSOP_LAKE /= DAMOSOP .OR.                                         &
 &         DAMONOP_LAKE /= DAMONOP .OR. DAMOEAP_LAKE /= DAMOEAP .OR.                                         &
 &         DXDELLA_LAKE /= DXDELLA .OR. DXDELLO_LAKE /= DXDELLO .OR. LLSCANNS_LAKE .NEQV. LLSCANNS ) THEN
           
         LLSAMEGRID = .FALSE.
      ELSE
         DO K = 1, NGY
           IF (NLONRGG_LAKE(MIN(K,NGY_LAKE)) /= NLONRGG(K) ) THEN
              LLSAMEGRID = .FALSE.
             EXIT
           ENDIF
         ENDDO
      ENDIF

      IF( .NOT.  LLSAMEGRID ) THEN
         LFILE=0
         IF (INFILENAME(IP) /= ' ') LFILE=LEN_TRIM(INFILENAME(IP))
         WRITE(IU06,*) "Checking that file ", INFILENAME(IP)(1:LFILE)
         WRITE(IU06,*) "has the same grid as input bathymetry"
         WRITE(IU06,*) "But it is not the case !!!"
         WRITE(IU06,*)  NGX, NGY, IPER, IRGG, IQGAUSS, LLSCANNS
         WRITE(IU06,*)  NGX_LAKE, NGY_LAKE, IPER_LAKE, IRGG_LAKE, IQGAUSS_LAKE, LLSCANNS_LAKE
         DO K = 1, NGY
           WRITE(IU06,*)  K, NLONRGG(K), NLONRGG_LAKE(MIN(K,NGY_LAKE))
         ENDDO
         WRITE(IU06,*) DAMOWEP, DAMOSOP, DAMOEAP, DAMONOP, DXDELLA, DXDELLO
         WRITE(IU06,*) DAMOWEP_LAKE, DAMOSOP_LAKE, DAMOEAP_LAKE, DAMONOP_LAKE, DXDELLA_LAKE, DXDELLO_LAKE
         CALL FLUSH(IU06)
         CALL WAM_ABORT("Not same grid !",__FILENAME__,__LINE__)
       ENDIF


       ! Check that it is the input as expected
       CALL IGRIB_GET_VALUE(KGRIB_HANDLE_LAKE(IP),'paramId', IVAL)
       IF ( IPARAMID(IP) /= IVAL ) THEN
         LFILE=0
         IF (INFILENAME(IP) /= ' ') LFILE=LEN_TRIM(INFILENAME(IP))
         WRITE(IU06,*) "Checking that file ", INFILENAME(IP)(1:LFILE)
         WRITE(IU06,*) "contains the correct paramId"
         WRITE(IU06,*) "But it is not the case !!!"
         WRITE(IU06,*) 'IPARAMID = ',IPARAMID(IP), 'paramId = ',IVAL
         CALL WAM_ABORT("paramId not matching !",__FILENAME__,__LINE__)
       ENDIF

    ELSE
      WRITE(IU06,*) "Trying to read from file ", INFILENAME(IP)
      WRITE(IU06,*) "but end of file reached !"
      CALL WAM_ABORT("end of file reached !",__FILENAME__,__LINE__)
    ENDIF

  ENDDO

! READING DATA
! ------------

  ! BATHY
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE_BATHY,'missingValue',ZMISS)
  CALL IGRIB_GET_VALUE(KGRIB_HANDLE_BATHY,'getNumberOfValues',NUMBEROFVALUES)
  ALLOCATE(VALUES_BATHY(NUMBEROFVALUES))
  CALL IGRIB_GET_VALUE(KGRIB_HANDLE_BATHY,'values',VALUES_BATHY)

  ! LAKE VARIABLES
  ALLOCATE(VALUES_LAKE(NUMBEROFVALUES,NPARAM))
  DO IP = 1, NPARAM
    CALL IGRIB_SET_VALUE(KGRIB_HANDLE_LAKE(IP),'missingValue',ZMISS)
    CALL IGRIB_GET_VALUE(KGRIB_HANDLE_LAKE(IP),'getNumberOfValues',IVAL)
    IF ( NUMBEROFVALUES /= IVAL ) THEN
      LFILE=0
      IF (INFILENAME(IP) /= ' ') LFILE=LEN_TRIM(INFILENAME(IP))
      WRITE(IU06,*) "Checking that file ", INFILENAME(IP)(1:LFILE)
      WRITE(IU06,*) "contains the correct number of values"
      WRITE(IU06,*) "But it is not the case !!!"
      WRITE(IU06,*) 'NUMBEROFVALUES = ',NUMBEROFVALUES, 'paramId = ',IVAL
      CALL WAM_ABORT("number of values not matching !",__FILENAME__,__LINE__)
    ENDIF
    CALL IGRIB_GET_VALUE(KGRIB_HANDLE_LAKE(IP),'values',VALUES_LAKE(:,IP))
  ENDDO


! COMBINING DATA
! ---------------

  MTHREADS=1
!$ MTHREADS=OMP_GET_MAX_THREADS()
  NPROMA=NUMBEROFVALUES/MTHREADS + 1

  ALLOCATE(INEWLAKE(MTHREADS))
  ALLOCATE(IAVG(MTHREADS))
  INEWLAKE(:) = 0
  IAVG(:) = 0

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JC, JCS, JCL, IC, IPR)
  DO JC = 1, NUMBEROFVALUES, NPROMA
    IPR = (JC+NPROMA-1)/NPROMA
    JCS = JC
    JCL = MIN(JCS+NPROMA-1, NUMBEROFVALUES)

    DO IC = JCS, JCL
      !! If land sea mask <= THRSLSM or lake cover > THRSLAKE then add that point to BATHY
      !  -----------------------------------------------------------------------
      IF ( VALUES_LAKE(IC,2) >= THRSLAKE ) THEN
      !! Lake point with cover above and equal to THRSLAKE, take the lake value
        IF ( VALUES_BATHY(IC) /= ZMISS ) THEN
!!!! for now only update lakes that were already in BATHY !!!!!!!!!!!!
          VALUES_BATHY(IC) = MIN(VALUES_LAKE(IC,3), BATHYMAX)
          INEWLAKE(IPR) = INEWLAKE(IPR) + 1
        ENDIF
      ELSEIF ( VALUES_LAKE(IC,1) <= THRSLSM .AND. VALUES_LAKE(IC,2) <= 0.01_JWRU ) THEN
      !! Not a lake point with a land sea mask below and equal THRSLSM (i.e. assumed to be ocean)
        IF ( VALUES_BATHY(IC) /= ZMISS ) THEN
          !! Take the average between BATHY and lake depth
          IF ( VALUES_BATHY(IC) <= 3.0_JWRU ) THEN
            VALUES_BATHY(IC) = 0.25_JWRU * VALUES_BATHY(IC) + 0.75_JWRU * MIN(VALUES_LAKE(IC,3), BATHYMAX)
          ELSEIF (VALUES_LAKE(IC,3) <= 3.0_JWRU ) THEN
            VALUES_BATHY(IC) = 0.75_JWRU * VALUES_BATHY(IC) + 0.25_JWRU * MIN(VALUES_LAKE(IC,3), BATHYMAX)
          ELSE
            VALUES_BATHY(IC) = 0.5_JWRU * ( VALUES_BATHY(IC) + MIN(VALUES_LAKE(IC,3), BATHYMAX) )
          ENDIF
          IAVG(IPR) = IAVG(IPR) + 1
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!$OMP END PARALLEL DO



! CLEANING UP
! -----------
  CALL IGRIB_CLOSE_FILE(IU07)
  DO IP = 1, NPARAM
    CALL IGRIB_CLOSE_FILE(IULAKE(IP))
    CALL IGRIB_RELEASE(KGRIB_HANDLE_LAKE(IP))
  ENDDO

  DEALLOCATE(VALUES_LAKE)

! SAVING NEW BATHY
! ----------------

  ! Reencode BATHY
  CALL IGRIB_CLONE(KGRIB_HANDLE_BATHY,KGRIB_HANDLE_NEW_BATHY)
  CALL IGRIB_RELEASE(KGRIB_HANDLE_BATHY)

  CALL IGRIB_SET_VALUE(KGRIB_HANDLE_NEW_BATHY,'missingValue',ZMISS)
  CALL IGRIB_SET_VALUE(KGRIB_HANDLE_NEW_BATHY,'values',VALUES_BATHY)

  DEALLOCATE(VALUES_BATHY)

  ! Save to file
  LFILE=0
  IF (OUTFILENAME /= ' ') LFILE=LEN_TRIM(OUTFILENAME)
  CALL IGRIB_OPEN_FILE(IU07,OUTFILENAME(1:LFILE),'w')
  WRITE(IU06,*) ''
  WRITE(IU06,*) ' New bathymetry file can be found in ',OUTFILENAME(1:LFILE)
  WRITE(IU06,*) ''
  WRITE(IU06,*) '  Number of new lake points      = ', SUM(INEWLAKE(:))
  WRITE(IU06,*) '  Number of average ocean points = ', SUM(IAVG(:))
  WRITE(IU06,*) ''

  CALL IGRIB_GET_MESSAGE_SIZE(KGRIB_HANDLE_NEW_BATHY,KBYTES)
  ISIZE=(KBYTES+NPRECI-1)/NPRECI
  ALLOCATE(KGRIB_BUFR(ISIZE))
  CALL IGRIB_GET_MESSAGE(KGRIB_HANDLE_NEW_BATHY,KGRIB_BUFR)

  CALL IGRIB_WRITE_BYTES(IU07,KGRIB_BUFR,KBYTES)
  CALL IGRIB_CLOSE_FILE(IU07)

ELSE

  CALL WAM_ABORT("Adding Lake information to binary bathymetry not available",__FILENAME__,__LINE__)

ENDIF ! GRIB OR BINARY INPUT


END PROGRAM 
