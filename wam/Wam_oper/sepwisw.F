      SUBROUTINE SEPWISW (FL3, FL1, IJS, IJL, IG,THWOLD,USOLD,USNEW)

C ----------------------------------------------------------------------
C
C**** *SEPWISW* - COMPUTES THE SWELL ENERGY, THE MEAN SWELL DIRECTION,
C****             THE MEAN SWELL FREQUENCY AND THE MEAN WINDSEA DIR.
C
C     P.LIONELLO     FEBRUARY 87
C
C     L.ZAMBRESKY    NOVEMBER 87   GKSS/ECMWF   OPTIMIZED SUB.
!     J. BIDLOT      FEBRARY 1996       ECMWF   MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       TO SEPARATE THE SWELL FROM THE WIND INTERACTING SEE
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SEPWISW (FL3, FL1, IJS, IJL, IG,THWOLD,USOLD,USNEW)*
C          *FL3* - BLOCK OF SPECTRA
C          *FL1* - SWELL FLAG ARRAY
C                  = 1 IF COMPONENT IS SWELL
C                  = 0 IF COMPONENT IS SEA
C          *IJS* - INDEX OF FIRST GRIDPOINT
C          *IJL* - INDEX OF LAST GRIDPOINT
C          *IG * - BLOCK NUMBER
C          *THWOLD*    INTERMEDIATE STORAGE OF ANGLE (RADIANS) OF
C                      WIND VELOCITY.
C          *USNEW*     NEW FRICTION VELOCITY IN M/S.
C          *USOLD*     INTERMEDIATE STORAGE OF MODULUS OF FRICTION
C                      VELOCITY.
C
C     METHOD.
C     -------
C
C       THE WAVES WHICH DO NOT INTERACT WITH THE WIND ARE
C       CONSIDERED SWELL.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comfred.h"
C
#include "commean.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comswel.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE)  :: FL1,FL3

      REAL,DIMENSION(NINF:NSUP,NBLO) :: THWOLD,USOLD
      REAL,DIMENSION(NINF:NSUP) :: USNEW
C
C ----------------------------------------------------------------------
      PARAMETER (FRIC = 28.)
      REAL,DIMENSION(IJS:IJL) :: CIS, SIS, CIWS, SIWS, SF3SW, SF3WIS
C
#include "parcons.h"
C
C ----------------------------------------------------------------------
C
C*    1. THE SWELL DISTRIBUTION IS COMPUTED.
C        -----------------------------------
C
 1000 CONTINUE
      DO M=1,NFRE
        DO K=1,NANG
          DO IJ=IJS,IJL
            FL1(IJ,K,M) = 0.
          ENDDO
        ENDDO
      ENDDO

      DO 1002 K=1,NANG
         DO 1003 IJ=IJS,IJL
            SIS(IJ) = 1.2*USOLD(IJ,IG)*COS(TH(K)-THWOLD(IJ,IG))
 1003    CONTINUE
         DO M=1,NFRE
           CM = FRIC/C(M)
           DO IJ=IJS,IJL
             CHECKTA = CM*SIS(IJ)
             IF (CHECKTA.LT.1.) THEN
               FL1(IJ,K,M)=1.
             ENDIF
           ENDDO
         ENDDO
 1002 CONTINUE
C
C*    2. COMPUTATION OF SWELL TOTAL ENERGY,MEAN FREQUENCY,MEAN DIRECTION
C*       AND OF WIND-SEA MEAN DIRECTION.
C        ---------------------------------------------------------------
C
 2000 CONTINUE
      DO 2001 IJ=IJS,IJL
         ESWELL(IJ) = 0.
         FSWELL(IJ) = 0.
         SIS(IJ)    = 0.
         CIS(IJ)    = 0.
         ESEA(IJ)   = 0.
         FSEA(IJ)   = 0.
         CIWS(IJ)   = 0.
         SIWS(IJ)   = 0.
 2001 CONTINUE
C
      DO 2002 K=1,NANG
         CIK = COSTH(K)
         SIK = SINTH(K)
         DO 2003 IJ=IJS,IJL
            SF3SW(IJ)  = 0.
            SF3WIS(IJ) = 0.
 2003    CONTINUE
         DO 2004 M=1,NFRE
            DF = DFIM(M)
            FD = DF/FR(M)
            DO 2005 IJ=IJS,IJL
               F3SW  = FL3(IJ,K,M)*FL1(IJ,K,M)
               F3WIS = FL3(IJ,K,M)-F3SW
               SF3SW(IJ) = SF3SW(IJ)+F3SW*DF
               SF3WIS(IJ) = SF3WIS(IJ)+F3WIS*DF
               ESWELL(IJ) = ESWELL(IJ)+DF*F3SW
               FSWELL(IJ) = FSWELL(IJ)+FD*F3SW
               ESEA(IJ)   = ESEA(IJ)+DF*F3WIS
               FSEA(IJ)   = FSEA(IJ)+FD*F3WIS
 2005       CONTINUE
 2004    CONTINUE
         DO 2006 IJ=IJS,IJL
            SIS(IJ)  = SIS(IJ)+SF3SW(IJ)*SIK
            CIS(IJ)  = CIS(IJ)+SF3SW(IJ)*CIK
            SIWS(IJ) = SIWS(IJ)+SF3WIS(IJ)*SIK
            CIWS(IJ) = CIWS(IJ)+SF3WIS(IJ)*CIK
 2006    CONTINUE
 2002 CONTINUE

      DO 2007 IJ=IJS,IJL
         IF (CIWS(IJ).EQ.0.) CIWS(IJ) = EPSMIN
         IF (CIS(IJ).EQ.0.) CIS(IJ) = EPSMIN
         IF (FSWELL(IJ).EQ.0.) FSWELL(IJ) = EPSMIN
         IF (FSEA(IJ).EQ.0.) FSEA(IJ) = EPSMIN
         IF (USNEW(IJ).EQ.0.) ESWELL(IJ) = EMEAN(IJ)
 2007 CONTINUE

      DO 2008 IJ=IJS,IJL
         THWISEA(IJ) = ATAN2(SIWS(IJ),CIWS(IJ))
         THSWELL(IJ) = ATAN2(SIS(IJ),CIS(IJ))
         FSWELL(IJ)  = ESWELL(IJ)/(FSWELL(IJ))
         FSEA(IJ)    = ESEA(IJ)/(FSEA(IJ))
 2008 CONTINUE

      RETURN
      END
