! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MBOUNF (IU03, IU10, IU20, IFORM, IINPC)

! ----------------------------------------------------------------------

!**** *MBOUNF* - MAKE FINE GRID BOUNDARY.

!     R. PORTZ     MPI         15/01/1991
!     Oyvind Breivik (oyvind.breivik@ecmwf.int) 2012-11-23: 
!       Added tolerance to grid comparison
!                   

!*    PURPOSE.
!     -------

!       COMPUTE ALL INFORMATION FOR FINE GRID BOUNDARY VALUE
!       INPUT (COMMON FBOUND).

!**   INTERFACE.
!     ----------

!       *CALL* *MBOUNF (IU03, IU10, IU20, IFORM, IINPC)
!          *IU03*  -  INPUT UNIT OF COARSE GRID BOUNDARY INFORMATION
!                     COMMON CBOUND GENERATED BY A COARSE GRID PREPROC
!                     (UNFORMATED).
!          *IU10*   - LOGICAL OUTPUT UNIT FOR  UNFORMATED WRITE OF
!                     COMMON FBOUND.
!          *IU20*   - LOGICAL OUTPUT UNIT FOR    FORMATED WRITE OF
!                     COMMON FBOUND.
!          *IFORM*  - FORMAT OPTION  = 1  UNFORMATED WRITE/READ.
!                                    = 2  FORMATED WRITE/READ
!                                    OTHERWISE WRITE BOTH
!                                              READ UNFORMATED.
!          *IINPC*  - NUMBER OF INPUT POINTS FROM COARSE GRID.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.
!       *FINDB*     - FIND BLOCK AND SEA POINT NUMBERS.
!       *MBOXB*     - FIND LAT AND LONG OF BOUNDARY POINTS.
!       *MINTF*     - MAKES INTERPOLATION TABLE FOR FINE TO COARSE GRID.
!       *PACKI*     - PACKS AN INTEGER ARRAY.
!       *PACKC*     - PACKS A  REAL    ARRAY.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCPBO  , ONLY : NBOUNC   ,IJARC    ,IGARC    ,DLAMAC   ,    &
     &            DPHIAC   ,BLATC    ,BLNGC    ,                        &
     &            GBOUNC   ,IPOGBO
      USE YOWFPBO  , ONLY : NBOUNF   ,IJARF    ,IGARF    ,IBFL     ,    &
     &            IBFR     ,BFW
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,    &
     &            XDELLA   ,XDELLO   ,NGX      ,NGY
      USE YOWPRPROC,ONLY : NBMAX
      USE YOWTEST  ,ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "findb.intfb.h"
#include "mboxb.intfb.h"
#include "mintf.intfb.h"
#include "packi.intfb.h"
#include "packr.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU03, IU10, IU20, IFORM
      INTEGER(KIND=JWIM), INTENT(OUT) :: IINPC

      INTEGER(KIND=JWIM) :: I, II, IRATIO, IO
      INTEGER(KIND=JWIM) :: NBOUNC_LOC
      INTEGER(KIND=JWIM) :: NBOUNEW
      INTEGER(KIND=JWIM) :: GBOUNC_LOC
      INTEGER(KIND=JWIM) :: ListSTART(1), ListEND(1)

      REAL(KIND=JWRB), PARAMETER :: GRDTOL = 0.000001_JWRB
      REAL(KIND=JWRB) :: GRDDIFF
      REAL(KIND=JWRB) :: AMOSOA, AMONOA, AMOEAA, AMOWEA
      REAL(KIND=JWRB), ALLOCATABLE :: BLATF(:), BLNGF(:)

! ----------------------------------------------------------------------

  998 FORMAT(10I8)
  999 FORMAT(5E16.7)

! ----------------------------------------------------------------------

!*    1. INITIAL.
!     -----------

      NBMAX = (NGX+NGY)*2-4

      ALLOCATE(IJARF(NBMAX))
      ALLOCATE(IGARF(NBMAX))
      ALLOCATE(IBFL(NBMAX))
      ALLOCATE(IBFR(NBMAX))
      ALLOCATE(BFW(NBMAX))
      ALLOCATE(BLATF(NBMAX))
      ALLOCATE(BLNGF(NBMAX))

      DO I = 1,NBMAX
        IJARF(I) = 0
        ! IGARF obsolete
        IGARF(I) = 1
        IBFL (I) = 0
        IBFR (I) = 0
        BFW  (I) = 0.0_JWRB
        BLATF(I) = 0.0_JWRB
        BLNGF(I) = 0.0_JWRB
      ENDDO

! ----------------------------------------------------------------------

!*    2. READ COMMON BOUNDC
!     ---------------------

      IF (IFORM.NE.2) THEN
        READ (IU03) GBOUNC_LOC
        ALLOCATE(IPOGBO(0:GBOUNC_LOC))
        READ(IU03)(IPOGBO(I),I=1,GBOUNC_LOC)
        IPOGBO(0)=0
        DO II=1,GBOUNC_LOC
!         FIND THE RIGHT DOMAIN
          READ (IU03) NBOUNC_LOC
          ALLOCATE(IGARC(NBOUNC_LOC))
          ALLOCATE(IJARC(NBOUNC_LOC))
          ALLOCATE(BLNGC(NBOUNC_LOC))
          ALLOCATE(BLATC(NBOUNC_LOC))
          READ (IU03) (IGARC(I),I=1,NBOUNC_LOC)
          READ (IU03) (IJARC(I),I=1,NBOUNC_LOC)
          READ (IU03) DLAMAC, DPHIAC, AMOSOA, AMONOA, AMOEAA, AMOWEA,   &
     &     (BLNGC(I),I=1,NBOUNC_LOC), (BLATC(I),I=1,NBOUNC_LOC)

          GRDDIFF = ABS(AMOWEP-AMOWEA)+ABS(AMOEAP-AMOEAA)+              &
     &              ABS(AMONOP-AMONOA)+ABS(AMOSOP-AMOSOA)
          IF (GRDDIFF < GRDTOL) THEN
!           FOUND THE RIGHT ONE
            EXIT
          ELSE
            DEALLOCATE (IJARC,IGARC,BLATC,BLNGC)
          ENDIF
        ENDDO
      ELSE
        READ (IU03,998) GBOUNC_LOC
        ALLOCATE(IPOGBO(0:GBOUNC_LOC))
        READ(IU03,998)(IPOGBO(I),I=1,GBOUNC_LOC)
        IPOGBO(0)=0
        DO II=1,GBOUNC
!         FIND THE RIGHT DOMAIN
          READ (IU03,998) NBOUNC_LOC
          ALLOCATE(IGARC(NBOUNC_LOC))
          ALLOCATE(IJARC(NBOUNC_LOC))
          ALLOCATE(BLNGC(NBOUNC_LOC))
          ALLOCATE(BLATC(NBOUNC_LOC))
          READ (IU03,998) (IGARC(I),I=1,NBOUNC_LOC)
          READ (IU03,998) (IJARC(I),I=1,NBOUNC_LOC)
          READ (IU03,999) DLAMAC, DPHIAC, AMOSOA,                       &
     &      AMONOA, AMOEAA, AMOWEA,                                     &
     &     (BLNGC(I),I=1,NBOUNC_LOC), (BLATC(I),I=1,NBOUNC_LOC)
          GRDDIFF = ABS(AMOWEP-AMOWEA)+ABS(AMOEAP-AMOEAA)+              &
     &              ABS(AMONOP-AMONOA)+ABS(AMOSOP-AMOSOA)
          IF (GRDDIFF < GRDTOL) THEN
!           FOUND THE RIGHT ONE
            EXIT
          ELSE
            DEALLOCATE (IJARC,IGARC,BLATC,BLNGC)
          ENDIF
        ENDDO
      ENDIF
      NBOUNC = NBOUNC_LOC
      IINPC = NBOUNC

      DEALLOCATE(IPOGBO)

! ----------------------------------------------------------------------

!*    3. CHECK THE INPUT
!     ------------------

!     IS THE FINE GRID THE SAME AS IN THE COARSE GRID ?

      IF (GRDDIFF >= GRDTOL) THEN
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. MBOUNF     *'
        WRITE(IU06,*) '*  ==========================     *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  THIS IS NOT ANY SAME GRID      *'
        WRITE(IU06,*) '*  AS YOU HAD DELARED IN THE      *'
        WRITE(IU06,*) '*  COARSE GRID                    *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  AMOSOP: ', AMOSOP, ' AMOSOA: ', AMOSOA
        WRITE(IU06,*) '*  AMONOP: ', AMONOP, ' AMONOA: ', AMONOA
        WRITE(IU06,*) '*  AMOWEP: ', AMOWEP, ' AMOWEA: ', AMOWEA
        WRITE(IU06,*) '*  AMOEAP: ', AMOEAP, ' AMOEAA: ', AMOEAA
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS          *'
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF

!     IS THE STEP OF LAT. AND LONG. OF THE FINE GRID A MULTIPLE
!              OF THE COARSE GRID ?
      IRATIO=NINT(DPHIAC/XDELLA)
      IF (ABS(IRATIO*XDELLA-DPHIAC)/DPHIAC.GT.1.E-6_JWRB) THEN
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. MBOUNF     *'
        WRITE(IU06,*) '*  --------------------------     *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  THE STEP OF THE LAT. OF THE    *'
        WRITE(IU06,*) '*  FINE GRID IS NOT A MULTIPLE    *'
        WRITE(IU06,*) '*  OF THE COARSE GRID             *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  FINE GRID XDELLA: ', XDELLA
        WRITE(IU06,*) '*  COARSE GRID DPHIAC: ', DPHIAC
        WRITE(IU06,*) '*  IRATIO: ',IRATIO
        WRITE(IU06,*) '*  DIFF = ',ABS(IRATIO*XDELLA-DPHIAC)/DPHIAC
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS          *'
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF

      IRATIO=NINT(DLAMAC/XDELLO)
      IF (ABS(IRATIO*XDELLO-DLAMAC)/DLAMAC.GT.1.E-6_JWRB) THEN
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. MBOUNF     *'
        WRITE(IU06,*) '*  ==========================     *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  THE STEP OF THE LON. OF THE    *'
        WRITE(IU06,*) '*  FINE GRID IS NOT A MULTIPLE    *'
        WRITE(IU06,*) '*  OF THE COARSE GRID             *'
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*  FINE GRID   XDELLO: ', XDELLO
        WRITE(IU06,*) '*  COARSE GRID DLAMAC: ', DLAMAC
        WRITE(IU06,*) '*  IRATIO: ',IRATIO
        WRITE(IU06,*) '*  DIFF = ',ABS(IRATIO*XDELLO-DLAMAC)
        WRITE(IU06,*) '*                                 *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS          *'
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF

! ----------------------------------------------------------------------

!*    4. COMPUTED THE SQUARE BOX
!     --------------------------

      CALL MBOXB (NBOUNF, AMOWEP, AMOSOP, AMOEAP, AMONOP, BLATF, BLNGF)

! ----------------------------------------------------------------------

!*    5. SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
!        -----------------------------------------
      ListSTART(1)=0
      ListEND(1)=0
      CALL FINDB (NBMAX, NBOUNF, BLATF, BLNGF, IJARF,            &
     &            1, ListSTART, ListEND, 1, 1)

! ----------------------------------------------------------------------

!*    6. MAKE INTERPOLATED ARRAYS
!     ----------------------------

      CALL MINTF

! ----------------------------------------------------------------------

!*    7. PACKED ALL ARRAYS
!     --------------------

      CALL PACKR (NBOUNF, NBOUNEW, NBMAX, IGARF, BLNGF)
      CALL PACKR (NBOUNF, NBOUNEW, NBMAX, IGARF, BLATF)
      CALL PACKR (NBOUNF, NBOUNEW, NBMAX, IGARF, BFW)
      CALL PACKI (NBOUNF, NBOUNEW, NBMAX, IGARF, IBFL)
      CALL PACKI (NBOUNF, NBOUNEW, NBMAX, IGARF, IBFR)
      CALL PACKI (NBOUNF, NBOUNEW, NBMAX, IGARF, IJARF)
      CALL PACKI (NBOUNF, NBOUNEW, NBMAX, IGARF, IGARF)
      NBOUNF = NBOUNEW

! ----------------------------------------------------------------------

!*    8. PRINTER PROTOCOL.
!        -----------------

      WRITE (IU06,'(''1BOUNDARY INPUT POINTS FOR THIS GRID:'')')
      WRITE (IU06,*) '===================================='
      WRITE (IU06,*) '  NUMBER OF BOUNDARY POINTS IS NBOUNF = ', NBOUNF
      WRITE (IU06,*) '        |-----FINE GRID INPUT POINTS-----|',      &
     &              '-RELATED COARSE GRID INDICES--|'
      WRITE (IU06,*) '      NO.    LAT.   LONG.  BLOCK.  POINT. ',      &
     &              '   LEFT   RIGHT   WEIGHT '

      DO IO = 1, NBOUNF
        WRITE (IU06,'(4X,I5,2F8.2,4I8,F10.4)')                          &
     &   IO, BLATF(IO), BLNGF(IO), IGARF(IO), IJARF(IO),                &
     &   IBFL(IO), IBFR(IO), BFW(IO)
      ENDDO

      IF (ALLOCATED(BLATF)) DEALLOCATE(BLATF)
      IF (ALLOCATED(BLNGF)) DEALLOCATE(BLNGF)

      IF (ALLOCATED(IJARC)) DEALLOCATE(IJARC)
      IF (ALLOCATED(IGARC)) DEALLOCATE(IGARC)
      IF (ALLOCATED(BLATC)) DEALLOCATE(BLATC)
      IF (ALLOCATED(BLNGC)) DEALLOCATE(BLNGC)

! ----------------------------------------------------------------------
!*    9. WRITE COMMON FBOUND

      IF (IFORM.NE.2) THEN
        WRITE(IU10) NBOUNF
        WRITE(IU10) (IGARF(I),I=1,NBOUNF), (IJARF(I),I=1,NBOUNF),       &
     &   (IBFL(I),I=1,NBOUNF), (IBFR(I),I=1,NBOUNF),                    &
     &   (BFW(I),I=1,NBOUNF)
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE(IU20,998) NBOUNF
        WRITE(IU20,998) (IGARF(I),I=1,NBOUNF), (IJARF(I),I=1,NBOUNF),   &
     &   (IBFL(I),I=1,NBOUNF), (IBFR(I),I=1,NBOUNF)
        WRITE(IU20,999) (BFW(I),I=1,NBOUNF)
      ENDIF

      END SUBROUTINE MBOUNF
