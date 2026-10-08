! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WRITEFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,      &
     &                   FILENAME, IUNIT, LOUNIT, LCUNIT, LRSTPARAL)

! ---------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     WRITES ARRAY FL TO FILE.

!**   INTERFACE.
!     ----------
!     CALL *WRITEFL(FL, IJINF, IJSUP, KINF, KSUP, MINF, MSUP,
!    &              FILENAME,IUNIT,LOUNIT,LCUNIT,LRSTPARAL)
!     *FL*        ARRAY TO BE WRITTEN TO FILE
!     *IJINF*     FIRST LOWER DIMENSION OF FL
!     *IJSUP*     FIRST UPPER DIMENSION OF FL
!     *KINF*      SECOND LOWER DIMENSION BOUND OF FL
!     *KSUP*      SECOND UPPER DIMENSION BOUND OF FL
!     *MINF*      THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*      THIRD UPPER DIMENSION BOUND OF FL
!     *FILENAME*  FILENAME (INCLUDING PATH) OF TARGET FILE
!     *IUNIT*     PBIO UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LOUNIT*    LOGICAL, TRUE IF FREE UNIT HAS TO BE FOUND
!     *LCUNIT*    LOGICAL, TRUE ON THE LAST CALL FOR THIS FILE
!     *LRSTPARAL* LOGICAL, TRUE THEN WRITE IN PARALLEL

!     METHOD.
!     -------
!     WRITES ARRAY FL TO FILE (A FORTRAN WRITE)

!     FL IS WRITTEN AS UNFORMATTED BINARY TO FILENAME CONNECTED
!     TO UNIT IUNIT.
!     IT IS IN A FORM THAT IS INDEPENDENT OF THE MODEL DECOMPOSITION.



!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE ISO_C_BINDING, ONLY : C_INT, C_LOC, C_PTR

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : LL1D, LLUNSTR
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06
      USE YOWABORT , ONLY : WAM_ABORT
#ifdef WAM_HAVE_UNWAM
      USE YOWUNBLKRORD, ONLY : UNBLKRORD
#endif

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "iwam_get_unit.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJINF, IJSUP, KINF, KSUP, MINF, MSUP
      INTEGER(KIND=JWIM), INTENT(INOUT) :: IUNIT

      REAL(KIND=JWRB), DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP), INTENT(INOUT) :: FL

      CHARACTER(LEN=296), INTENT(IN) :: FILENAME

      LOGICAL, INTENT(IN) :: LOUNIT, LCUNIT, LRSTPARAL

      INTEGER(KIND=JWIM) :: IJ, K, M, J2, J3
      INTEGER(KIND=JWIM) :: KRET, KOUNT, LFILE

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJINF:IJSUP,KINF:KSUP,MINF:MSUP) :: FL_G
      REAL(KIND=JWRB), ALLOCATABLE :: RFL(:,:,:)

!     Mirror of the READFL fast path: writing through CONVERT='BIG_ENDIAN'
!     swaps element by element, so instead swap the whole record in one
!     pass and emit the record markers by hand over STREAM access.
      INTEGER(KIND=JWIM), PARAMETER :: IELEMBYTES = STORAGE_SIZE(1.0_JWRB)/8
#if defined(__amdflang__)
      LOGICAL, PARAMETER :: LLFAST = (IELEMBYTES == 4 .OR. IELEMBYTES == 8)
#else
      LOGICAL, PARAMETER :: LLFAST = .FALSE.
#endif
!     Record markers stay 4-byte whatever the real kind, so with 8-byte
!     reals markers and payload no longer swap at the same width and the
!     staging layout has to differ. See WRITE_RECORD.
      LOGICAL, PARAMETER :: LLWIDE = (IELEMBYTES == 8)
!     Matches READFL: ceiling on the staging buffer used to batch whole
!     planes into one write.
      INTEGER(KIND=8), PARAMETER :: MAXRAW = 1073741824_8

!     Held across calls, as in READFL: re-faulting the staging buffer every
!     band costs more than the batching saves.
      REAL(KIND=JWRB), ALLOCATABLE, TARGET, SAVE :: RAW(:)

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

      IF (LHOOK) CALL DR_HOOK('WRITEFL',0,ZHOOK_HANDLE)


      LFILE=0
      IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
      IF (LLFAST) THEN
!       Opened once and held to the end of the file. Reopening per plane
!       costs an append seek every call for no benefit.
        IF (LOUNIT) THEN
          OPEN(NEWUNIT=IUNIT, FILE=FILENAME(1:LFILE),                   &
     &         FORM='UNFORMATTED', ACCESS='STREAM',                     &
     &         STATUS='REPLACE', ACTION='WRITE')
        ENDIF
      ELSEIF (LOUNIT) THEN
        IUNIT=IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'w', 'u', 0, 'READWRITE')
      ELSE
        IUNIT=IWAM_GET_UNIT(-1, FILENAME(1:LFILE) , 'a', 'u', 0, 'READWRITE')
      ENDIF

      IF (LLUNSTR .AND. .NOT.LRSTPARAL) THEN
#ifdef WAM_HAVE_UNWAM
        FL_G(0,:,:)=0.0_JWRB
        CALL UNBLKRORD(1,IJINF,IJSUP,KINF,KSUP,MINF,MSUP,               &
     &                 FL(IJINF:IJSUP,KINF:KSUP,MINF:MSUP),             &
     &               FL_G(IJINF:IJSUP,KINF:KSUP,MINF:MSUP))

        IF (LLFAST) THEN
          CALL WRITE_RECORD(FL_G)
        ELSE
          CALL WRITE_PLANES(FL_G)
        ENDIF
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSEIF (LRSTPARAL .OR. LL1D .OR. NPROC == 1) THEN
        IF (LLFAST) THEN
!         Swap a copy so the caller's FL is left untouched.
          FL_G = FL
          CALL WRITE_RECORD(FL_G)
        ELSE
          CALL WRITE_PLANES(FL)
        ENDIF
      ELSE
!     WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!     BUT THE SINGLE BINARY INPUT FILES SHOULD BE IN THE OLD MAPPING

        IF (LLFAST) THEN
!         Relabel on the way into the staging buffer. Saves a full pass
!         over the spectra and the FL_G automatic array it needed.
          CALL WRITE_RECORD_REORDER(FL)
        ELSE
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(J3, J2, IJ)
          DO J3 = MINF, MSUP
            DO J2 = KINF, KSUP
              DO IJ = IJINF, IJSUP
                FL_G(IJ,J2,J3) = FL(IJ2NEWIJ(IJ),J2,J3)
              ENDDO
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO

          CALL WRITE_PLANES(FL_G)
        ENDIF
      ENDIF

      IF (LCUNIT .OR. .NOT.LLFAST) CLOSE(IUNIT)

      IF (LHOOK) CALL DR_HOOK('WRITEFL',1,ZHOOK_HANDLE)

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

      SUBROUTINE WRITE_PLANES(BUF)
!     One record per (direction, frequency) plane, as the fast path emits
!     and as the restart format requires. Callers hand over a whole band,
!     so a single WRITE over the span would merge its planes into one
!     record and change the file format.
      REAL(KIND=JWRB), INTENT(IN) :: BUF(IJINF:IJSUP,KINF:KSUP,MINF:MSUP)
      INTEGER(KIND=JWIM) :: J2, J3

      DO J3 = MINF, MSUP
        DO J2 = KINF, KSUP
          WRITE(IUNIT) BUF(:,J2,J3)
        ENDDO
      ENDDO
      END SUBROUTINE WRITE_PLANES

      SUBROUTINE WRITE_RECORD(BUF)
!     Emit one big-endian sequential record per (direction, frequency)
!     plane on the STREAM-opened unit. Mirror of READFL: when the caller
!     hands over several planes, assemble the whole span - markers
!     interleaved with payload - in native order, flip it in one
!     BSWAP32_ARRAY pass, and push it out in a single write rather than
!     five writes per plane. In the single-plane case BUF is swapped in
!     place, so it must be scratch the caller does not need afterwards;
!     the batched path stages through RAW and leaves BUF untouched.
      REAL(KIND=JWRB), TARGET, INTENT(INOUT) :: BUF(IJINF:IJSUP,KINF:KSUP,MINF:MSUP)
      INTEGER(KIND=JWIM) :: NELEM, NBYTES, NPLANE, NPERW, NWRITE, NWORD
      INTEGER(KIND=JWIM) :: IPL, IP0, IB, JP, J2, J3, NK, NSLOT
      REAL(KIND=JWRB) :: ZMARK

      NELEM  = IJSUP-IJINF+1
      NBYTES = NELEM*IELEMBYTES
      NK     = KSUP-KINF+1
      NPLANE = NK*(MSUP-MINF+1)

      IF (NPLANE == 1) THEN
        CALL SWAP_ELEMS(C_LOC(BUF(IJINF,KINF,MINF)), NELEM)
        WRITE(IUNIT) INT(ISWAP32(INT(NBYTES,4)),4)
        WRITE(IUNIT) BUF
        WRITE(IUNIT) INT(ISWAP32(INT(NBYTES,4)),4)
        RETURN
      ENDIF

      IF (LLWIDE) THEN
        CALL WRITE_RECORD_WIDE(BUF, NELEM, NBYTES, NPLANE, NK, .FALSE.)
        RETURN
      ENDIF

!     Marker carried in native order; the bulk swap below puts it on disk
!     big-endian along with the payload.
      ZMARK = TRANSFER(INT(NBYTES,4), ZMARK)

      NWORD = NELEM+2
      NPERW = MAX(1, INT(MAXRAW/(INT(NWORD,8)*IELEMBYTES), JWIM))
      NPERW = MIN(NPERW, NPLANE)
      IF (ALLOCATED(RAW)) THEN
        IF (SIZE(RAW) < NPERW*NWORD) DEALLOCATE(RAW)
      ENDIF
      IF (.NOT.ALLOCATED(RAW)) ALLOCATE(RAW(NPERW*NWORD))

      IPL = 0
      DO WHILE (IPL < NPLANE)
        NWRITE = MIN(NPERW, NPLANE-IPL)

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JP, IB, IP0, J2, J3)
        DO JP = 1, NWRITE
          IP0 = IPL + JP - 1
          J3  = MINF + IP0/NK
          J2  = KINF + MOD(IP0, NK)
          IB  = (JP-1)*NWORD
          RAW(IB+1)            = ZMARK
          RAW(IB+2:IB+NELEM+1) = BUF(:,J2,J3)
          RAW(IB+NWORD)        = ZMARK
        ENDDO
!$OMP   END PARALLEL DO

        CALL BSWAP32_ARRAY(C_LOC(RAW(1)), INT(NWRITE*NWORD, C_INT))
        WRITE(IUNIT) RAW(1:NWRITE*NWORD)

        IPL = IPL + NWRITE
      ENDDO
      END SUBROUTINE WRITE_RECORD

      SUBROUTINE WRITE_RECORD_WIDE(FLIN, NELEM, NBYTES, NPLANE, NK, LLREORDER)
!     Batched writer for 8-byte reals. Markers stay 4-byte, so markers and
!     payload cannot share a swap width. The trailing marker of a plane and
!     the leading marker of the next are adjacent in the stream and, every
!     record in a span being the same length, identical - so the pair fills
!     exactly one 8-byte slot. Staging them that way keeps every slot
!     aligned and lets one BSWAP64_ARRAY pass cover markers and payload
!     together: the pair slot comes back with its halves exchanged, which
!     for equal halves is no change. Only the first leading and last
!     trailing marker of a batch fall outside the pairing, and those two go
!     out on their own, so the payload still leaves in a single write.
      REAL(KIND=JWRB), INTENT(IN) :: FLIN(IJINF:IJSUP,KINF:KSUP,MINF:MSUP)
      INTEGER(KIND=JWIM), INTENT(IN) :: NELEM, NBYTES, NPLANE, NK
      LOGICAL, INTENT(IN) :: LLREORDER
      INTEGER(KIND=JWIM) :: NPERW, NWRITE, NSLOT, NCAP, NSTRIDE
      INTEGER(KIND=JWIM) :: IPL, IP0, IB, JP, J2, J3, IJ
      REAL(KIND=JWRB) :: ZPAIR
      INTEGER(KIND=4) :: IMARKBE

      ZPAIR   = TRANSFER([INT(NBYTES,4), INT(NBYTES,4)], ZPAIR)
      IMARKBE = ISWAP32(INT(NBYTES,4))
      NSTRIDE = NELEM+1

      NPERW = MAX(1, INT(MAXRAW/(INT(NSTRIDE,8)*IELEMBYTES), JWIM))
      NPERW = MIN(NPERW, NPLANE)
      NCAP  = NPERW*NSTRIDE
      IF (ALLOCATED(RAW)) THEN
        IF (SIZE(RAW) < NCAP) DEALLOCATE(RAW)
      ENDIF
      IF (.NOT.ALLOCATED(RAW)) ALLOCATE(RAW(NCAP))

      IPL = 0
      DO WHILE (IPL < NPLANE)
        NWRITE = MIN(NPERW, NPLANE-IPL)
        NSLOT  = NWRITE*NSTRIDE - 1

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JP, IB, IP0, J2, J3, IJ)
        DO JP = 1, NWRITE
          IP0 = IPL + JP - 1
          J3  = MINF + IP0/NK
          J2  = KINF + MOD(IP0, NK)
          IB  = (JP-1)*NSTRIDE
          IF (LLREORDER) THEN
            DO IJ = IJINF, IJSUP
              RAW(IB+1+IJ-IJINF) = FLIN(IJ2NEWIJ(IJ),J2,J3)
            ENDDO
          ELSE
            RAW(IB+1:IB+NELEM) = FLIN(:,J2,J3)
          ENDIF
          IF (JP < NWRITE) RAW(IB+NSTRIDE) = ZPAIR
        ENDDO
!$OMP   END PARALLEL DO

        CALL BSWAP64_ARRAY(C_LOC(RAW(1)), INT(NSLOT, C_INT))

        WRITE(IUNIT) IMARKBE
        WRITE(IUNIT) RAW(1:NSLOT)
        WRITE(IUNIT) IMARKBE

        IPL = IPL + NWRITE
      ENDDO
      END SUBROUTINE WRITE_RECORD_WIDE

      SUBROUTINE WRITE_RECORD_REORDER(FLIN)
!     As WRITE_RECORD, but applies the 2-D decomposition relabelling while
!     filling the staging buffer, so the permutation costs no extra pass
!     and FLIN is left untouched.
      REAL(KIND=JWRB), INTENT(IN) :: FLIN(IJINF:IJSUP,KINF:KSUP,MINF:MSUP)
      INTEGER(KIND=JWIM) :: NELEM, NBYTES, NPLANE, NPERW, NWRITE, NWORD
      INTEGER(KIND=JWIM) :: IPL, IP0, IB, JP, J2, J3, NK, IJ
      REAL(KIND=JWRB) :: ZMARK

      NELEM  = IJSUP-IJINF+1
      NBYTES = NELEM*IELEMBYTES
      NK     = KSUP-KINF+1
      NPLANE = NK*(MSUP-MINF+1)

      IF (LLWIDE) THEN
        CALL WRITE_RECORD_WIDE(FLIN, NELEM, NBYTES, NPLANE, NK, .TRUE.)
        RETURN
      ENDIF

      ZMARK  = TRANSFER(INT(NBYTES,4), ZMARK)

      NWORD = NELEM+2
      NPERW = MAX(1, INT(MAXRAW/(INT(NWORD,8)*IELEMBYTES), JWIM))
      NPERW = MIN(NPERW, NPLANE)
      IF (ALLOCATED(RAW)) THEN
        IF (SIZE(RAW) < NPERW*NWORD) DEALLOCATE(RAW)
      ENDIF
      IF (.NOT.ALLOCATED(RAW)) ALLOCATE(RAW(NPERW*NWORD))

      IPL = 0
      DO WHILE (IPL < NPLANE)
        NWRITE = MIN(NPERW, NPLANE-IPL)

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JP, IB, IP0, J2, J3, IJ)
        DO JP = 1, NWRITE
          IP0 = IPL + JP - 1
          J3  = MINF + IP0/NK
          J2  = KINF + MOD(IP0, NK)
          IB  = (JP-1)*NWORD
          RAW(IB+1)     = ZMARK
          RAW(IB+NWORD) = ZMARK
          DO IJ = IJINF, IJSUP
            RAW(IB+2+IJ-IJINF) = FLIN(IJ2NEWIJ(IJ),J2,J3)
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

        CALL BSWAP32_ARRAY(C_LOC(RAW(1)), INT(NWRITE*NWORD, C_INT))
        WRITE(IUNIT) RAW(1:NWRITE*NWORD)

        IPL = IPL + NWRITE
      ENDDO
      END SUBROUTINE WRITE_RECORD_REORDER

      END SUBROUTINE WRITEFL
