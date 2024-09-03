! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE OUTCOM (IU07, BATHY, LLGRIB_BATHY_OUT)

! ----------------------------------------------------------------------

!**** *OUTCOM* - ROUTINE TO WRITE MODULES TO DISK

!     H.GUNTHER            ECMWF       04/04/1990

!     J. BIDLOT            ECMWF       10/1998.
!!!!!!                     COMMON BLOCKS HAVE BEEN CONVERTED TO MODULES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT            ECMWF       11/2003
!                          IF YOUR ARE RUNNING AT ECMWF:
!                          BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                          STRUCTURE OF THE OUTPUT FILE YOU WILL HAVE TO
!                          MAKE SURE THAT IT IS CREATED FOR YOUR RUN, 
!                          OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                          THAT IS ALREADY ON DISK.
!                          YOU ALSO HAVE TO CHANGE READPRE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*    PURPOSE.
!     -------

!       TO WRITE OUT THE COMPUTED MODULES
!       (MODULE UBUF IS WRITTEN IN MUBUF)

!**   INTERFACE.
!     ----------

!       *CALL* *OUTCOM (IU07, BATHY, LLGRIB_BATHY_OUT)*
!          *IU07*   - LOGICAL UNIT FOR  UNFORMATED WRITE.

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : LFDB, FFLAG, GFLAG, NFLAG, IRBATHY, INFOBOUT
      USE YOWGRIBHD, ONLY : CDATECLIM,IMDLGRBID_G
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWMAP   , ONLY : NGX      ,NGY      , IPER     ,IRGG    ,    &
     &                      AMOWEP   ,AMOSOP   ,AMOEAP    ,AMONOP  ,    &
     &                      XDELLA   ,XDELLO   ,NLONRGG
      USE YOWSTAT  , ONLY : MARSTYPE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWABORT, ONLY : WAM_ABORT

#ifdef WAM_HAVE_UNWAM
      USE UNWAM , ONLY : UNWAM_OUT
      USE YOWUNPOOL ,ONLY : LPREPROC
#endif

      USE YOWGRIB  , ONLY : IGRIB_CLOSE_FILE

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "wgribenout.intfb.h"
#include "outnam.intfb.h"
#include "mpcrtbl.intfb.h"
 
      INTEGER(KIND=JWIM), INTENT(IN) :: IU07
      REAL(KIND=JWRB), INTENT(INOUT) :: BATHY(NGX, NGY)
      LOGICAL, INTENT(IN) :: LLGRIB_BATHY_OUT

      INTEGER(KIND=JWIM) :: IDUM, K, M, L, JY
      INTEGER(KIND=JWIM) :: ITABLE, IPARAM, IZLEV, IFCST, ITMIN, ITMAX
      INTEGER(KIND=JWIM) :: NKIND !Precision used when writing

      REAL(KIND=JWRB), ALLOCATABLE :: ZDUM(:)

      CHARACTER(LEN=14) :: CDATE

! ----------------------------------------------------------------------

!*    GRID INFORMATION
!     ----------------

      IF ( LLGRIB_BATHY_OUT ) THEN
        ! Grib output

        !!! check MPCRTBL
        FFLAG(:)=.FALSE.
        NFLAG(:)=.FALSE.
        GFLAG(:)=.TRUE.

        CALL MPCRTBL

        ITABLE=INFOBOUT(IRBATHY,1)
        IPARAM=INFOBOUT(IRBATHY,2)
        IZLEV=INFOBOUT(IRBATHY,3)
        ITMIN=INFOBOUT(IRBATHY,4)
        ITMAX=INFOBOUT(IRBATHY,5)

        CDATE=CDATECLIM
        IFCST=0

        !!! For grib output the bathymetry file has to be re-order from north to south
        !!! ecWAM internally has BATHY defined from south to north !!!
        ALLOCATE(ZDUM(NGX))
        DO JY=1,NGY/2
          ZDUM(1:NGX) = BATHY(1:NGX,JY)
          BATHY(1:NGX,JY) = BATHY(1:NGX,NGY-JY+1)
          BATHY(1:NGX,NGY-JY+1) = ZDUM(1:NGX)
        ENDDO
        DEALLOCATE(ZDUM)

        CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, BATHY,               &
     &                  ITABLE, IPARAM, IZLEV, ITMIN, ITMAX, 0, 0,  &
     &                  CDATE, IFCST, MARSTYPE, LFDB, IU07)

      ELSE
        ! Binary output

!       WRITE IDENTIFIERS (MAKE SURE TO UPDATE THE VALUES IF YOU CHANGE
!       ANYTHING TO THE MODEL). 
        NKIND = KIND(AMOSOP)
        WRITE(IU07) NKIND, IMDLGRBID_G

        WRITE(IU07) NGX, NGY

        WRITE(IU07) NLONRGG

        WRITE(IU07) IPER, IRGG, AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO

        WRITE(IU07) BATHY

      ENDIF

! ----------------------------------------------------------------------

!*    WRITE NAMELIST PARWAM.
!     ----------------------
      
      CALL OUTNAM  (NGX, NGY)

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        IF (LPREPROC) THEN
          CALL UNWAM_OUT(IU07)
        ENDIF
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ENDIF

      IF ( LLGRIB_BATHY_OUT ) THEN
        CALL IGRIB_CLOSE_FILE(IU07)
      ELSE
        CLOSE(IU07)
      ENDIF

      END SUBROUTINE OUTCOM
