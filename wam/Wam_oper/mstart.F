      SUBROUTINE MSTART (IU12, IU14, IU15, IOPTI, FETCH, FRMAX,
     1                   FL1,U10OLD,THWOLD,USOLD,Z0OLD,TAUW)
C ----------------------------------------------------------------------
C
C**** *MSTART* - MAKES START FIELDS FOR WAMODEL.
C
C      H. GUNTHER    ECMWF    MAY 1990
C      H. GUNTHER    ECMWF    DECEMBER 90  MODIFIED FOR CYCLE_4.
C      J. BIDLOT     ECMWF    FEBRUARY 96  MESSAGE PASSING
C      J. BIDLOT     ECMWF    MARCH 97     REMOVE ALL OUTPUT TO IU12 and IU15 
C
C*    PURPOSE.
C     --------
C
C       TO GENERATE WAMODEL START FIELDS.
C
C**   INTERFACE.
C     ----------
C
C   *CALL* *MSTART (IU12, IU14, IU15, IOPTI, FETCH, FRMAX,
C    1              FL1,U10OLD,THWOLD,USOLD,Z0OLD,TAUW)*
C      *IU12*   INTEGER    OUTPUT UNIT BLOCKS OF SPECTRA.
C      *IU14*   INTEGER    OUTPUT UNIT SECOND LAT OF BLOCKS.
C      *IU15*   INTEGER    OUTPUT UNIT LAST WINDFIELDS.
C      *IOPTI*  INTEGER    START FIELD OPTION
C                          = 0 FROM PARAMTERS.
C                          = 1 FROM WINDS CALM ENERGY=0.
C                          = 2 FROM WINDS CALM FROM PARAMETERS.
C      *FETCH*  REAL       FETCH IN METERS.
C      *FRMAX*  REAL       MAXIMUM PEAK FREQUENCY IN HERTZ.
C      *FL1*      REAL      2-D SPECTRUM FOR EACH GRID POINT 
C      *U10NEW*    NEW WIND SPEED IN M/S.
C      *U10OLD*    INTERMEDIATE STORAGE OF MODULUS OF WIND
C                  VELOCITY.
C      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                  NOTATION (POINTING ANGLE OF WIND VECTOR,
C                  CLOCKWISE FROM NORTH).
C      *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
C                  WIND VELOCITY.
C      *USNEW*     NEW FRICTION VELOCITY IN M/S.
C      *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
C                  VELOCITY.
C      *Z0OLD*     INTERMEDIATE STORAGE OF ROUGHNESS LENGTH IN
C                  M.
C      *TAUW*      WAVE STRESS IN (M/S)**2
C
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT*     - TERMINATES PROCESSING.
C       *PEAK*      - COMPUTE PARAMETERS FROM WIND FOR A BLOCK.
C       *SPECTRA*   - COMPUTES SPECTRA OF A BLOCK.
C       *SPLITBL*   - STORING OVERLAPPING LATITUDES OF BLOCKS.
C
C    REFERENCE.
C    ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comcout.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "comjons.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL1 

C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      CHARACTER*12 ZERO
      CHARACTER FILENA*16,CYMPPATH*80,FILENAME*80,MODE*1
      PARAMETER(ZERO=' ')
C
C ----------------------------------------------------------------------
C
C*    1. DEFINE SPECTRUM FOR LAND POINTS AND WRITE OUTPUT.
C        -------------------------------------------------
C
 1000 CONTINUE
      DO M=1,NFRE
        DO K=1,NANG
          FL1(NINF-1,K,M) = 0.
        ENDDO
      ENDDO

      WRITE (IU06,'(1X,/,1X,''  PARAMETER AT OUTPUT SIDES:'')')
      WRITE (IU06,'(1X,''  NGOU    IG    IJ     U10    UDIR'',
     1            ''      FP   ALPHAJONS   GAMMA      SA      SB'')')
C
C ----------------------------------------------------------------------
C
C*    2. LOOP FOR BLOCKS.
C        ----------------
C
 2000 CONTINUE
      DO 2001 IG=1,IGL
C
C*    2.1 COMPUTE PEAK FREQUENCIES AND ALPHAJONS PARAMETERS.
C         ----------------------------------------------
C
C
C*    2.1.1 INITIAL VALUES DUE TO OPTION.
C           -----------------------------
C
 2110 CONTINUE
         IF (IOPTI.EQ.1) THEN
            DO 2111 IJ = IJS(IG), IJL(IG)
               FP(IJ) = 0.
               ALPHJ(IJ) = 0.
               THES(IJ) = THWOLD(IJ,IG)
 2111       CONTINUE
         ELSE
            DO 2112 IJ = IJS(IG), IJL(IG)
               FP(IJ) = FM
               ALPHJ(IJ) = ALFA
               IF (U10OLD(IJ,IG) .GT. 0.1E-08) THEN
                  THES(IJ) = THWOLD(IJ,IG)
               ELSE
                  THES(IJ) = THES(IJ)
               ENDIF
 2112       CONTINUE
         ENDIF
C
C*    2.1.2 PEAK FREQUENCY AND ALPHAJONS FROM FETCH LAW.
C           --------------------------------------------
C
         IF (IOPTI.NE.0) THEN
            CALL PEAK (IJS(IG), IJL(IG), IG, FETCH, FRMAX,U10OLD)
         ENDIF
         IF (ITEST.GT.1) THEN
            IF (IG.LE.ITESTB) WRITE (IU06,*) '    SUB. PEAK DONE'
         ENDIF
C
C*    2.1.3 PRINT PARAMETERS AT OUTPUT POINTS.
C           ----------------------------------
C
         DO 2131 NGOU = 1, NGOUT
            IF (IG.EQ.IGAR(NGOU)) THEN
               IJ = IJAR(NGOU)
               WRITE (IU06,'(1X,3I6,F8.2,F8.2,5F8.4)')  NGOU, IG, IJ,
     1            U10OLD(IJ,IG), THWOLD(IJ,IG)*DEG, FP(IJ), ALPHJ(IJ),
     2            GAMMA, SA, SB
            ENDIF
 2131    CONTINUE
C
C*    2.2 COMPUTE SPECTRA FROM PARAMETERS.
C         --------------------------------
C
         CALL SPECTRA (IJS(IG), IJL(IG), IG, FL1)
         IF (ITEST.GT.1) THEN
            IF (IG.LE.ITESTB) WRITE (IU06,*) '    SUB. SPECTRA DONE'
         ENDIF
C
C*    2.3 OUTPUT OF OVERLAPPING LATITUDES.
C         --------------------------------
C
         IF (IG.GT.1) CALL SPLITBL (FL1, IJS(IG), IJL2(IG), IU14)
         IF (ITEST.GT.1) THEN
            IF (IG.LE.ITESTB) WRITE (IU06,*) '    SUB. SPLITBL DONE'
         ENDIF
C
C
C*    BRANCHING BACK TO 2. FOR NEXT BLOCK.
C
 2001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. PREPARE OUTPUT OF WIND.
C        ----------------------
C
 3000 CONTINUE
      CDTPRO  = ZERO
      CDATEWO = ZERO
      CDAWIFL = ZERO
      CDATEFL = ZERO

      RETURN
      END
