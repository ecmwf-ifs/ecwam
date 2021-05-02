      SUBROUTINE MSTART (IU12, IU14, IU15, IOPTI, FETCH, FRMAX,         &
     &                   IJS, IJL, FL1, U10OLD, THWOLD)
! ----------------------------------------------------------------------

!**** *MSTART* - MAKES START FIELDS FOR WAMODEL.

!      H. GUNTHER    ECMWF    MAY 1990
!      H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.
!      J. BIDLOT     ECMWF    FEBRUARY 96  MESSAGE PASSING
!      J. BIDLOT     ECMWF    MARCH 97     REMOVE ALL OUTPUT TO IU12
!                                          and IU15 

!*    PURPOSE.
!     --------

!       TO GENERATE WAMODEL START FIELDS.

!**   INTERFACE.
!     ----------

!   *CALL* *MSTART (IU12, IU14, IU15, IOPTI, FETCH, FRMAX,
!    1              FL1,U10OLD,THWOLD)*
!      *IU12*   INTEGER    OUTPUT UNIT BLOCKS OF SPECTRA.
!      *IU14*   INTEGER    OUTPUT UNIT SECOND LAT OF BLOCKS.
!      *IU15*   INTEGER    OUTPUT UNIT LAST WINDFIELDS.
!      *IOPTI*  INTEGER    START FIELD OPTION
!                          = 0 FROM PARAMETERS.
!                          = 1 FROM WINDS CALM ENERGY=0.
!                          = 2 FROM WINDS CALM FROM PARAMETERS.
!      *FETCH*  REAL       FETCH IN METERS.
!      *FRMAX*  REAL       MAXIMUM PEAK FREQUENCY IN HERTZ.
!      *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT 
!      *U10NEW*    NEW WIND SPEED IN M/S.
!      *U10OLD*    INTERMEDIATE STORAGE OF MODULUS OF WIND
!                  VELOCITY.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
!                  WIND VELOCITY.
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *ABORT*     - TERMINATES PROCESSING.
!       *PEAK*      - COMPUTE PARAMETERS FROM WIND FOR A BLOCK.
!       *SPECTRA*   - COMPUTES SPECTRA OF A BLOCK.

!    REFERENCE.
!    ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO    ,NBLO
      USE YOWPCONS , ONLY : DEG
      USE YOWCOUT  , ONLY : NGOUT    ,IGAR     ,IJAR
      USE YOWJONS  , ONLY : FP       ,ALPHJ    ,THES     ,FM       ,    &
     &            ALFA     ,GAMMA    ,SA       ,SB       ,THETAQ
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "peak.intfb.h"
#include "spectra.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU12, IU14, IU15, IOPTI
      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL

      REAL(KIND=JWRB), INTENT(IN) :: FETCH, FRMAX
      REAL(KIND=JWRB),DIMENSION(NINF:NSUP,NBLO), INTENT(IN) :: U10OLD, THWOLD
      REAL(KIND=JWRB),DIMENSION(IJS:IJL, NANG, NFRE), INTENT(INOUT) :: FL1

      INTEGER(KIND=JWIM) :: IG
      INTEGER(KIND=JWIM) :: M, K, IJ, NGOU

      CHARACTER(LEN=14), PARAMETER :: ZERO='              '

! ----------------------------------------------------------------------

!     0. ALLOCATE ARRAYS
!        ---------------

      IG = 1

      ALLOCATE(FP(NIBLO))
      ALLOCATE(ALPHJ(NIBLO))
      ALLOCATE(THES(NIBLO))

!*    1. DEFINE SPECTRUM FOR LAND POINTS AND WRITE OUTPUT.
!        -------------------------------------------------

      WRITE (IU06,'(1X,/,1X,''  PARAMETER AT OUTPUT SITES:'')')
      WRITE (IU06,'(1X,''  NGOU    IG    IJ     U10    UDIR'',          &
     &            ''      FP   ALPHAJONS   GAMMA      SA      SB'')')

! ----------------------------------------------------------------------

!*    2.1 COMPUTE PEAK FREQUENCIES AND ALPHAJONS PARAMETERS.
!         ----------------------------------------------


!*    2.1.1 INITIAL VALUES DUE TO OPTION.
!           -----------------------------

        IF (IOPTI.EQ.1) THEN
          DO IJ = IJS, IJL
            FP(IJ) = 0.0_JWRB
            ALPHJ(IJ) = 0.0_JWRB
            THES(IJ) = THWOLD(IJ,IG)
          ENDDO
        ELSE IF (IOPTI.EQ.0) THEN
          DO IJ = IJS, IJL
            FP(IJ) = FM
            ALPHJ(IJ) = ALFA
            THES(IJ) = THETAQ
          ENDDO
        ELSE
          DO IJ = IJS, IJL
            FP(IJ) = FM
            ALPHJ(IJ) = ALFA
            IF (U10OLD(IJ,IG) .GT. 0.1E-08_JWRB) THEN
              THES(IJ) = THWOLD(IJ,IG)
            ELSE
              THES(IJ) = 0.0_JWRB
            ENDIF
          ENDDO
        ENDIF

!*    2.1.2 PEAK FREQUENCY AND ALPHAJONS FROM FETCH LAW.
!           --------------------------------------------

        IF (IOPTI.NE.0) THEN
          CALL PEAK (IJS, IJL, IG, FETCH, FRMAX, U10OLD)
        ENDIF
        IF (ITEST.GT.1) THEN
          IF (IG.LE.ITESTB) WRITE (IU06,*) '    SUB. PEAK DONE'
        ENDIF

!*    2.1.3 PRINT PARAMETERS AT OUTPUT POINTS.
!           ----------------------------------

        DO NGOU = 1, NGOUT
          IF (IG.EQ.IGAR(NGOU)) THEN
            IJ = IJAR(NGOU)
            WRITE (IU06,'(1X,3I6,F8.2,F8.2,5F8.4)')  NGOU, IG, IJ,      &
     &       U10OLD(IJ,IG), THWOLD(IJ,IG)*DEG, FP(IJ), ALPHJ(IJ),       &
     &       GAMMA, SA, SB
          ENDIF
        ENDDO

!*    2.2 COMPUTE SPECTRA FROM PARAMETERS.
!         --------------------------------

        CALL SPECTRA (IJS, IJL, FL1)
        IF (ITEST.GT.1) THEN
          IF (IG.LE.ITESTB) WRITE (IU06,*) '    SUB. SPECTRA DONE'
        ENDIF


! ----------------------------------------------------------------------

!*    3. PREPARE OUTPUT OF WIND.
!        ----------------------

      CDTPRO  = ZERO
      CDATEWO = ZERO
      CDAWIFL = ZERO
      CDATEFL = ZERO

      END SUBROUTINE MSTART
