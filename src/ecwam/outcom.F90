! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE OUTCOM (IU07, IU17, IFORM, BATHY)

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

!       *CALL* *OUTCOM (IU07, IU17, IFORM)*
!          *IU07*   - LOGICAL UNIT FOR  UNFORMATED WRITE.
!          *IU17*   - LOGICAL UNIT FOR    FORMATED WRITE.
!          *IFORM*   - FORMAT OPTION  = 1  UNFORMATED WRITE.
!                                     = 2  FORMATED WRITE.
!                                     OTHERWISE BOTH.

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

      USE YOWGRIBHD, ONLY : IMDLGRBID_G,IMDLGRBID_M 
      USE YOWPARAM , ONLY : NGX      ,NGY      ,CLDOMAIN ,IMDLGRDID,LLUNSTR
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,    &
     &            IJS      ,IJL
      USE YOWMAP   , ONLY : NX       ,NY       ,    &
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,NLONRGG
      USE YOWABORT, ONLY : WAM_ABORT

#ifdef WAM_HAVE_UNWAM
      USE UNWAM , ONLY : UNWAM_OUT
      USE YOWUNPOOL ,ONLY : LPREPROC
#endif

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "outnam.intfb.h"
 
      INTEGER(KIND=JWIM), INTENT(IN) :: IU07, IU17, IFORM
      REAL(KIND=JWRB), INTENT(IN) :: BATHY(NGX, NGY)

      INTEGER(KIND=JWIM) :: IDUM, K, M, L
      INTEGER(KIND=JWIM) :: NKIND !Precision used when writing

! ----------------------------------------------------------------------

  995 FORMAT(5I8)
  996 FORMAT(I8,1X,2E16.7)
  997 FORMAT(14I8,1X,A1)
  998 FORMAT(10I8)
  999 FORMAT(5E16.7)

! ----------------------------------------------------------------------
!        WRITE IDENTIFIERS (MAKE SURE TO UPDATE THE VALUES IF YOU CHANGE
!        ANYTHING TO THE MODEL). 
      
      NKIND = KIND(DELPHI)

      IF (IFORM /= 2) THEN
        WRITE(IU07) NKIND, IMDLGRDID, IMDLGRBID_G, IMDLGRBID_M
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE(IU17,998) NKIND, IMDLGRDID, IMDLGRBID_G, IMDLGRBID_M
      ENDIF
!*    0. WRITE YOWPARAM (BLOCK SIZES).
!        ----------------------------

!     IDUM REPLACES IREFRA WHICH IS NO LONGER USED IN PREPROC.
      IDUM=0

      IF (IFORM /= 2) THEN
        WRITE(IU07) NGX, NGY, CLDOMAIN
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE(IU17,997) NGX, NGY, CLDOMAIN 
      ENDIF


! ----------------------------------------------------------------------

!*    2. WRITE MODULE YOWGRID.
!        ---------------------

      IF (IFORM /= 2) THEN
        WRITE (IU07) DELPHI, (DELLAM(L),L=1,NY), (NLONRGG(L),L=1,NY),   &
     &   (SINPH(L),L=1,NY), (COSPH(L),L=1,NY),                          &
     &   IJS, IJL
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,999) DELPHI,(DELLAM(L),L=1,NY),(REAL(NLONRGG(L),JWRB),L=1,NY), &
     &   (SINPH(L),L=1,NY), (COSPH(L),L=1,NY)
        WRITE (IU17,998) IJS, IJL
      ENDIF

! ----------------------------------------------------------------------

!*    3. WRITE MODULE YOWMAP.
!        --------------------

      IF (IFORM /= 2) THEN
        WRITE (IU07) NX, NY, IPER,          &
     &   AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,                &
     &   ZDELLO, IRGG
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,998) NX, NY, IPER
        WRITE (IU17,999) AMOWEP, AMOSOP, AMOEAP, AMONOP,                &
     &   XDELLA, XDELLO,ZDELLO
      ENDIF

! ----------------------------------------------------------------------

!*    8. WRITE MODULE YOWSHAL.
!        --------------------

      IF (IFORM /= 2) THEN
        WRITE (IU07) BATHY
      ENDIF
      IF (IFORM /= 1) THEN
        WRITE (IU17,999) BATHY
      ENDIF

! ----------------------------------------------------------------------

!*   11. WRITE NAMELIST PARWAM.
!        ----------------------
      
      CALL OUTNAM  (NGX, NGY)

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        IF (LPREPROC) THEN
          CALL UNWAM_OUT(IU07)
        ENDIF
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      END IF

      END SUBROUTINE OUTCOM
