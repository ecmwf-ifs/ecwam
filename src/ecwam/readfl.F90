! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE READFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,       &
     &                  FILENAME, IUNIT, LOUNIT, LCUNIT, LRSTPARAL)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997

!*    PURPOSE.
!     --------
!     READS ARRAY FL

!**   INTERFACE.
!     ----------
!     CALL *READFL*(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,
!                   FILENAME, IUNIT, LOUNIT, LCUNI, LRSTPARAL)
!     *FL*       ARRAY TO BE WRITTEN TO FILE
!     *IJINF*    FIRST LOWER DIMENSION OF FL
!     *IJSUP*    FIRST UPPER DIMENSION OF FL
!     *KINF*     SECOND LOWER DIMENSION BOUND OF FL
!     *KSUP*     SECOND UPPER DIMENSION BOUND OF FL
!     *MINF*     THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*     THIRD UPPER DIMENSION BOUND OF FL
!     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE
!     *IUNIT*    FILE UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LOUNIT*  LOGICAL, TRUE IF FREE UNIT HAS TO BE FOUND
!     *LCUNIT* LOGICAL, TRUE IF UNIT HAs TO BE CLOSED 
!     *LRSTPARAL* LOGICAL, TRUE THEN READ IN PARALLEL


!     METHOD.
!     -------
!     READS ARRAY FL FROM FILE (FORTRAN READ)

!     IF PBIO IS USED, THEN PBREAD IS CALLED TO READ (MSUP-MINF+1) 
!     FREQUENCY CONTRIBUTION TO FL.


!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE ISO_C_BINDING, ONLY : C_INT, C_LOC, C_PTR

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D     ,LLUNSTR
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06
#ifdef WAM_HAVE_UNWAM
      USE YOWUNBLKRORD, ONLY : UNBLKRORD
#endif

      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
      USE YOWABORT, ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "iwam_get_unit.intfb.h"
      INTEGER(KIND=JWIM), INTENT(IN) :: IJINF, IJSUP, KINF, KSUP, MINF, MSUP
      INTEGER(KIND=JWIM), INTENT(INOUT) :: IUNIT

      REAL(KIND=JWRB), DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: FL

      CHARACTER(LEN=296), INTENT(IN) :: FILENAME

      LOGICAL, INTENT(IN) :: LOUNIT, LCUNIT, LRSTPARAL

      INTEGER(KIND=JWIM) :: LFILE, IJ, J2, J3

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP) :: FL_G

      LOGICAL :: LLEXIST

!     The restart is big-endian Fortran sequential: each record is a
!     leading 4-byte length marker, the payload, then the same marker
!     again. Opening with CONVERT='BIG_ENDIAN' makes the runtime swap
!     element by element, which runs ~28x slower than reading the bytes
!     raw and swapping them in one pass. LLFAST takes the raw path; it
!     needs STREAM access so the markers can be stepped over by hand.
      INTEGER(KIND=JWIM), PARAMETER :: IELEMBYTES = STORAGE_SIZE(1.0_JWRB)/8
#if defined(__amdflang__)
      LOGICAL, PARAMETER :: LLFAST = (IELEMBYTES == 4 .OR. IELEMBYTES == 8)
#else
      LOGICAL, PARAMETER :: LLFAST = .FALSE.
#endif
!     Record markers stay 4-byte whatever the real kind, so with 8-byte
!     reals markers and payload no longer swap at the same width and the
!     staging layout has to differ. See READ_RECORD_WIDE.
      LOGICAL, PARAMETER :: LLWIDE = (IELEMBYTES == 8)
!     Ceiling on the staging buffer used to batch whole planes into one
!     read. Large enough that the read length stops mattering, small
!     enough that it is irrelevant beside the spectra themselves.
      INTEGER(KIND=8), PARAMETER :: MAXRAW = 1073741824_8

!     Staging buffer for the batched reads, shared by both record readers
!     and held across calls: every band is the same size, and handing
!     several hundred MB back to the allocator each time costs more in
!     page faults than the batching saves.
      REAL(KIND=JWRB), ALLOCATABLE, TARGET, SAVE :: RAW(:)

!     Passed as C_PTR rather than a typed buffer so one interface serves
!     both real kinds.
      INTERFACE
        SUBROUTINE BSWAP32_ARRAY(BUF, N) BIND(C, NAME='bswap32_array')
          USE ISO_C_BINDING, ONLY : C_PTR, C_INT
          TYPE(C_PTR), VALUE :: BUF
          INTEGER(KIND=C_INT), VALUE, INTENT(IN) :: N
        END SUBROUTINE BSWAP32_ARRAY
        SUBROUTINE BSWAP64_ARRAY(BUF, N) BIND(C, NAME='bswap64_array')
          USE ISO_C_BINDING, ONLY : C_PTR, C_INT
          TYPE(C_PTR), VALUE :: BUF
          INTEGER(KIND=C_INT), VALUE, INTENT(IN) :: N
        END SUBROUTINE BSWAP64_ARRAY
      END INTERFACE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('READFL',0,ZHOOK_HANDLE)

      LFILE=0
      IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)

      IF (LOUNIT) THEN
        LLEXIST=.FALSE.
        INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
        IF (.NOT. LLEXIST) THEN
          WRITE (IU06,*) '*************************************'
          WRITE (IU06,*) '*                                   *'
          WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
          WRITE (IU06,*) '*  IN READFL :                      *'
          WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILENAME
          WRITE (IU06,*) '*                                   *'
          WRITE (IU06,*) '*************************************'
          WRITE (*,*) '*************************************'
          WRITE (*,*) '*                                   *'
          WRITE (*,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
          WRITE (*,*) '*  IN READFL :                      *'
          WRITE (*,*) '*  COULD NOT FIND FILE ',FILENAME
          WRITE (*,*) '*                                   *'
          WRITE (*,*) '*************************************'
          CALL ABORT1
        ENDIF
        IF (LLFAST) THEN
          OPEN(NEWUNIT=IUNIT, FILE=FILENAME(1:LFILE),                   &
     &         FORM='UNFORMATTED', ACCESS='STREAM',                     &
     &         STATUS='OLD', ACTION='READ')
        ELSE
          IUNIT=IWAM_GET_UNIT(IU06, FILENAME(1:LFILE), 'r', 'u',0,'READWRITE')
        ENDIF
      ENDIF

      IF (LLUNSTR .AND. .NOT.LRSTPARAL) THEN
#ifdef WAM_HAVE_UNWAM
        IF (LLFAST) THEN
          CALL READ_RECORD(FL_G)
        ELSE
          CALL READ_PLANES(FL_G)
        ENDIF

        CALL UNBLKRORD(-1,IJINF,IJSUP,KINF,KSUP,MINF,MSUP,              &
     &                 FL(IJINF:IJSUP,KINF:KSUP,MINF:MSUP),             &
     &               FL_G(IJINF:IJSUP,KINF:KSUP,MINF:MSUP))
#else
      CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSEIF (LRSTPARAL .OR. LL1D .OR. NPROC == 1) THEN
        IF (LLFAST) THEN
          CALL READ_RECORD(FL)
        ELSE
          CALL READ_PLANES(FL)
        ENDIF
      ELSE
!       WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!       BUT THE BINARY INPUT FILES ARE IN THE OLD MAPPING

        IF (LLFAST) THEN
!         Permute straight out of the staging buffer. The separate FL_G
!         copy and the 640 MB automatic array it needs both disappear:
!         the file is read, swapped, then permuted, three passes over the
!         spectra instead of four.
          CALL READ_RECORD_REORDER(FL)
        ELSE
          CALL READ_PLANES(FL_G)

!         RE-ORDER
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(J3, J2, IJ)
          DO J3 = MINF, MSUP
            DO J2 = KINF, KSUP
              DO IJ = IJINF, IJSUP
                FL(IJ2NEWIJ(IJ),J2,J3) = FL_G(IJ,J2,J3)
              ENDDO
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
        ENDIF


      ENDIF
      IF (LCUNIT) CLOSE(IUNIT)

      IF (LHOOK) CALL DR_HOOK('READFL',1,ZHOOK_HANDLE)

      CONTAINS

      INTEGER(KIND=JWIM) FUNCTION ISWAP32(I)
!     Byte-reverse a single 32-bit record marker.
      INTEGER(KIND=4), INTENT(IN) :: I
      ISWAP32 = IOR(IOR(ISHFT(IBITS(I, 0,8),24), ISHFT(IBITS(I, 8,8),16)),  &
     &              IOR(ISHFT(IBITS(I,16,8), 8),        IBITS(I,24,8)))
      END FUNCTION ISWAP32

      SUBROUTINE SWAP_ELEMS(P, N)
!     Byte-reverse N reals at P, at whatever width the build uses.
      TYPE(C_PTR), INTENT(IN) :: P
      INTEGER(KIND=JWIM), INTENT(IN) :: N
      IF (LLWIDE) THEN
        CALL BSWAP64_ARRAY(P, INT(N, C_INT))
      ELSE
        CALL BSWAP32_ARRAY(P, INT(N, C_INT))
      ENDIF
      END SUBROUTINE SWAP_ELEMS

      SUBROUTINE BAD_MARKER(ILEAD, ITRAIL, NEXPECT)
      INTEGER(KIND=JWIM), INTENT(IN) :: ILEAD, ITRAIL, NEXPECT
      WRITE (IU06,*) '* READFL: BAD RECORD MARKER IN ',FILENAME(1:LFILE)
      WRITE (IU06,*) '*   LEADING  =',ILEAD
      WRITE (IU06,*) '*   TRAILING =',ITRAIL
      WRITE (IU06,*) '*   EXPECTED =',NEXPECT
      WRITE (*,*)    '* READFL: BAD RECORD MARKER IN ',FILENAME(1:LFILE)
      CALL ABORT1
      END SUBROUTINE BAD_MARKER

      SUBROUTINE READ_PLANES(BUF)
!     One record per (direction, frequency) plane, as the file holds and
!     as the fast path emits. Callers ask for a whole band, so a single
!     READ over the span would expect its planes merged into one record.
      REAL(KIND=JWRB), INTENT(OUT) :: BUF(IJINF:IJSUP,KINF:KSUP,MINF:MSUP)
      INTEGER(KIND=JWIM) :: J2, J3

      DO J3 = MINF, MSUP
        DO J2 = KINF, KSUP
          READ(IUNIT) BUF(:,J2,J3)
        ENDDO
      ENDDO
      END SUBROUTINE READ_PLANES

      SUBROUTINE READ_RECORD(BUF)
!     The restart holds one big-endian sequential record per (direction,
!     frequency) plane, so a request spanning several planes is several
!     records laid end to end. Rather than issue a read per plane with two
!     4-byte marker reads in between - which drops the device back to idle
!     between every plane - pull as many whole records as the staging
!     budget allows in a single read and unpick them in memory. Byte order
!     is fixed for the span in one BSWAP32_ARRAY pass, markers included,
!     which leaves each marker sitting in native order ready to verify.
      REAL(KIND=JWRB), TARGET, INTENT(OUT) :: BUF(IJINF:IJSUP,KINF:KSUP,MINF:MSUP)
      INTEGER(KIND=4) :: IMARK1, IMARK2
      INTEGER(KIND=JWIM) :: NELEM, NEXPECT, NPLANE, NPERR, NREAD, NWORD
      INTEGER(KIND=JWIM) :: IPL, IP0, IB, JP, J2, J3, NK

      NELEM   = IJSUP-IJINF+1
      NEXPECT = NELEM*IELEMBYTES
      NK      = KSUP-KINF+1
      NPLANE  = NK*(MSUP-MINF+1)

!     One plane: read straight into BUF, no staging copy needed.
      IF (NPLANE == 1) THEN
        READ(IUNIT) IMARK1
        READ(IUNIT) BUF
        READ(IUNIT) IMARK2
        IF (ISWAP32(IMARK1) /= NEXPECT .OR. ISWAP32(IMARK2) /= NEXPECT) THEN
          CALL BAD_MARKER(ISWAP32(IMARK1), ISWAP32(IMARK2), NEXPECT)
        ENDIF
        CALL SWAP_ELEMS(C_LOC(BUF(IJINF,KINF,MINF)), NELEM)
        RETURN
      ENDIF

      IF (LLWIDE) THEN
        CALL READ_RECORD_WIDE(BUF, NELEM, NEXPECT, NPLANE, NK, .FALSE.)
        RETURN
      ENDIF

      NWORD = NELEM+2
      NPERR = MAX(1, INT(MAXRAW/(INT(NWORD,8)*IELEMBYTES), JWIM))
      NPERR = MIN(NPERR, NPLANE)
      IF (ALLOCATED(RAW)) THEN
        IF (SIZE(RAW) < NPERR*NWORD) DEALLOCATE(RAW)
      ENDIF
      IF (.NOT.ALLOCATED(RAW)) ALLOCATE(RAW(NPERR*NWORD))

      IPL = 0
      DO WHILE (IPL < NPLANE)
        NREAD = MIN(NPERR, NPLANE-IPL)

        READ(IUNIT) RAW(1:NREAD*NWORD)
        CALL BSWAP32_ARRAY(C_LOC(RAW(1)), INT(NREAD*NWORD, C_INT))

        DO JP = 1, NREAD
          IB = (JP-1)*NWORD
          IF (TRANSFER(RAW(IB+1), 0_4) /= NEXPECT .OR.                  &
     &        TRANSFER(RAW(IB+NWORD), 0_4) /= NEXPECT) THEN
            CALL BAD_MARKER(TRANSFER(RAW(IB+1), 0_4),                   &
     &                      TRANSFER(RAW(IB+NWORD), 0_4), NEXPECT)
          ENDIF
        ENDDO

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JP, IB, IP0, J2, J3)
        DO JP = 1, NREAD
          IP0 = IPL + JP - 1
          J3  = MINF + IP0/NK
          J2  = KINF + MOD(IP0, NK)
          IB  = (JP-1)*NWORD
          BUF(:,J2,J3) = RAW(IB+2:IB+NELEM+1)
        ENDDO
!$OMP   END PARALLEL DO

        IPL = IPL + NREAD
      ENDDO
      END SUBROUTINE READ_RECORD

      SUBROUTINE READ_RECORD_WIDE(FLOUT, NELEM, NEXPECT, NPLANE, NK, LLREORDER)
!     Batched reader for 8-byte reals, mirroring WRITE_RECORD_WIDE. The
!     marker between two planes is read as one 8-byte slot holding the
!     trailing and leading markers back to back, so a single
!     BSWAP64_ARRAY pass brings payload and markers alike into native
!     order; the pair slot's halves come back exchanged, which for the
!     equal-length records of one span is no change. The batch's first
!     leading and last trailing marker sit outside the pairing and are
!     read separately.
      REAL(KIND=JWRB), INTENT(OUT) :: FLOUT(IJINF:IJSUP,KINF:KSUP,MINF:MSUP)
      INTEGER(KIND=JWIM), INTENT(IN) :: NELEM, NEXPECT, NPLANE, NK
      LOGICAL, INTENT(IN) :: LLREORDER
      INTEGER(KIND=JWIM) :: NPERR, NREAD, NSLOT, NCAP, NSTRIDE
      INTEGER(KIND=JWIM) :: IPL, IP0, IB, JP, J2, J3, IJ
      INTEGER(KIND=4) :: IMARK1, IMARK2, IPAIR(2)

      NSTRIDE = NELEM+1
      NPERR = MAX(1, INT(MAXRAW/(INT(NSTRIDE,8)*IELEMBYTES), JWIM))
      NPERR = MIN(NPERR, NPLANE)
      NCAP  = NPERR*NSTRIDE
      IF (ALLOCATED(RAW)) THEN
        IF (SIZE(RAW) < NCAP) DEALLOCATE(RAW)
      ENDIF
      IF (.NOT.ALLOCATED(RAW)) ALLOCATE(RAW(NCAP))

      IPL = 0
      DO WHILE (IPL < NPLANE)
        NREAD = MIN(NPERR, NPLANE-IPL)
        NSLOT = NREAD*NSTRIDE - 1

        READ(IUNIT) IMARK1
        READ(IUNIT) RAW(1:NSLOT)
        READ(IUNIT) IMARK2
        CALL BSWAP64_ARRAY(C_LOC(RAW(1)), INT(NSLOT, C_INT))

        IF (ISWAP32(IMARK1) /= NEXPECT .OR. ISWAP32(IMARK2) /= NEXPECT) THEN
          CALL BAD_MARKER(ISWAP32(IMARK1), ISWAP32(IMARK2), NEXPECT)
        ENDIF
        DO JP = 1, NREAD-1
          IPAIR = TRANSFER(RAW(JP*NSTRIDE), IPAIR, 2)
          IF (IPAIR(1) /= NEXPECT .OR. IPAIR(2) /= NEXPECT) THEN
            CALL BAD_MARKER(IPAIR(1), IPAIR(2), NEXPECT)
          ENDIF
        ENDDO

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JP, IB, IP0, J2, J3, IJ)
        DO JP = 1, NREAD
          IP0 = IPL + JP - 1
          J3  = MINF + IP0/NK
          J2  = KINF + MOD(IP0, NK)
          IB  = (JP-1)*NSTRIDE
          IF (LLREORDER) THEN
            DO IJ = IJINF, IJSUP
              FLOUT(IJ2NEWIJ(IJ),J2,J3) = RAW(IB+1+IJ-IJINF)
            ENDDO
          ELSE
            FLOUT(:,J2,J3) = RAW(IB+1:IB+NELEM)
          ENDIF
        ENDDO
!$OMP   END PARALLEL DO

        IPL = IPL + NREAD
      ENDDO
      END SUBROUTINE READ_RECORD_WIDE

      SUBROUTINE READ_RECORD_REORDER(FLOUT)
!     As READ_RECORD, but applies the 2-D decomposition relabelling on the
!     way out of the staging buffer instead of afterwards. Saves a full
!     pass over the spectra and the FL_G automatic array that pass needed.
      REAL(KIND=JWRB), TARGET, INTENT(OUT) :: FLOUT(IJINF:IJSUP,KINF:KSUP,MINF:MSUP)
      INTEGER(KIND=JWIM) :: NELEM, NEXPECT, NPLANE, NPERR, NREAD, NWORD
      INTEGER(KIND=JWIM) :: IPL, IP0, IB, JP, J2, J3, NK, IJ

      NELEM   = IJSUP-IJINF+1
      NEXPECT = NELEM*IELEMBYTES
      NK      = KSUP-KINF+1
      NPLANE  = NK*(MSUP-MINF+1)

      IF (LLWIDE) THEN
        CALL READ_RECORD_WIDE(FLOUT, NELEM, NEXPECT, NPLANE, NK, .TRUE.)
        RETURN
      ENDIF

      NWORD   = NELEM+2

      NPERR = MAX(1, INT(MAXRAW/(INT(NWORD,8)*IELEMBYTES), JWIM))
      NPERR = MIN(NPERR, NPLANE)
      IF (ALLOCATED(RAW)) THEN
        IF (SIZE(RAW) < NPERR*NWORD) DEALLOCATE(RAW)
      ENDIF
      IF (.NOT.ALLOCATED(RAW)) ALLOCATE(RAW(NPERR*NWORD))

      IPL = 0
      DO WHILE (IPL < NPLANE)
        NREAD = MIN(NPERR, NPLANE-IPL)

        READ(IUNIT) RAW(1:NREAD*NWORD)
        CALL BSWAP32_ARRAY(C_LOC(RAW(1)), INT(NREAD*NWORD, C_INT))

        DO JP = 1, NREAD
          IB = (JP-1)*NWORD
          IF (TRANSFER(RAW(IB+1), 0_4) /= NEXPECT .OR.                  &
     &        TRANSFER(RAW(IB+NWORD), 0_4) /= NEXPECT) THEN
            CALL BAD_MARKER(TRANSFER(RAW(IB+1), 0_4),                   &
     &                      TRANSFER(RAW(IB+NWORD), 0_4), NEXPECT)
          ENDIF
        ENDDO

!       Each plane owns a distinct (J2,J3), so the scatter cannot collide
!       across threads. Kept as a scatter rather than a gather through the
!       inverse permutation: RAW is cache-warm from the swap just above, so
!       streaming it sequentially beats making those reads random.
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JP, IB, IP0, J2, J3, IJ)
        DO JP = 1, NREAD
          IP0 = IPL + JP - 1
          J3  = MINF + IP0/NK
          J2  = KINF + MOD(IP0, NK)
          IB  = (JP-1)*NWORD + 2 - IJINF
          DO IJ = IJINF, IJSUP
            FLOUT(IJ2NEWIJ(IJ),J2,J3) = RAW(IB+IJ)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

        IPL = IPL + NREAD
      ENDDO
      END SUBROUTINE READ_RECORD_REORDER

      END SUBROUTINE READFL
