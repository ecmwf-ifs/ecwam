      SUBROUTINE OUTGRID (FL3, IJS, IJL, IG,U10NEW,THWNEW,USNEW,TAUW,
     1                    NSTART,NEND)
C ----------------------------------------------------------------------
C
C**** *OUTGRID* - OUTPUT OF WAVE ENERGY,MEAN ANGLE, MEAN FREQUENCY
C****             FRICTION VELOCITY AND WIND DIRECTION IN GRIDDED
C****             FORMAT.
C
C     P.JANSSEN      KNMI
C
C
C     P.LIONELLO  FEB. 87
C                 OUTPUT OF SWELL ENERGY ,MEAN SWELL DIRECTION ,
C                 MEAN SWELL FREQUENCY AND MEAN WIND-SEA WAVE
C                 DIRECTION AT ALL ACTIVE GRID POINTS
C
C     H.GUNTHER   GKSS/ECMWF     NOVEMBER 1989
C     J. BIDLOT   ECMWF          APRIL 1996
C                 MESSAGE PASSING
C
C                 OUTPUT FORMAT IS NOW GRIDDED FIELDS.
C
C*    PURPOSE.
C     --------
C
C       OUTPUT OF WAVE AND WIND FIELDS INTO COMMON BLOCK.
C
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *OUTGRID (FL3, IJS, IJL, IG,U10NEW,THWNEW,USNEW,TAUW)*
C         *FL3*    - BLOCK OF SPECTRA.
C         *IJS*    - INDEX OF FIRST GRIDPOINT.
C         *IJL*    - INDEX OF LAST GRIDPOINT.
C         *IG*     - BLOCK NUMBER.
C         *U10NEW*    NEW WIND SPEED IN M/S.
C         *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                     NOTATION (POINTING ANGLE OF WIND VECTOR,
C                     CLOCKWISE FROM NORTH).
C         *USNEW*     NEW FRICTION VELOCITY IN M/S.
C         *TAUW*      WAVE STRESS IN (M/S)**2
C         *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C         *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C
C
C     EXTERNALS.
C     ----------
C
C       *MAKEGRID*  - GENERATES GRIDDED FIELDS.
C       *PEAKFR*    - COMPUTATION OF PEAK FREQUENCY.
C       *MPGATHERSCFLD* - GATHER SCALAR FIELD DISTRIBUTED ON DIFFERENT
C                         PROCESSES
C
C     METHOD.
C     -------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
C      
#include "param.h"
C
#include "parcons.h"
C
#include "comcout.h"
C
#include "comintp.h"
C
#include "comints.h"
C
#include "commean.h"
C
#include "comice.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "comstat.h"
C
#include "comswel.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL3
      REAL,DIMENSION(NINF:NSUP,NBLO) :: TAUW
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW

C ----------------------------------------------------------------------

      REAL,ALLOCATABLE,DIMENSION(:) :: PEAKF,CD,WVSTR,THW,US,U10 
C
C ----------------------------------------------------------------------
C

C*    1. INTEGRATED PARAMETERS OF TOTAL SEA.
C        -----------------------------------
C
      IF (CDTINTT.EQ.CDTPRO) THEN
C
C*    1.1 COMPUTE PEAK FREQUENCY.
C         -----------------------
C
         ALLOCATE (PEAKF(NIBLO))
         CALL PEAKFR (FL3, IJS, IJL, PEAKF)
C
C*    1.2 COMPUTE NORMALISED WAVESTRESS AND DRAG. 
C         ON ICEPOINTS GIVE MISSING VALUE INDICATOR.
C         ---------------------------------------
C
         ALLOCATE (WVSTR(NIBLO))
         ALLOCATE (CD(NIBLO))
         DO 1200 IJ = IJS,IJL
            TAU = USNEW(IJ)**2 +0.0001
            WVSTR(IJ) = TAUW(IJ,IG)/TAU
CC       if the numerical computation of TAU and CD change, a similar
CC       modification has to be put in buildstress where the friction
cc       velocity is determined from U10 and CD.
            CD(IJ) = TAU/(U10NEW(IJ)**2+0.01)
 1200    CONTINUE
C
         IF (LICERUN) THEN
            DO 1210 IJ = IJS,IJL
               IF (ICEMASK(IJ,IG).EQ.0) THEN
                   EMEAN(IJ)  = ZMISS
                   THQ(IJ)    = ZMISS
                   PEAKF(IJ)  = ZMISS
                   FMEAN(IJ)  = ZMISS
                   SMEAN(IJ)  = ZMISS
                   THWNEW(IJ) = ZMISS
                   USNEW(IJ)  = ZMISS
                   CD(IJ)     = ZMISS
                   WVSTR(IJ)  = ZMISS   
               ENDIF   
 1210      CONTINUE
         ENDIF
C
C*    1.3 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
C         -------------------------------------------------
C
         IF(LMESSPASS) THEN
           CALL MPGATHERSCFLD(IPFGTBL(1),1,NSTART,NEND,EMEAN)
           IF(IPFGTBL(1).EQ.IRANK) 
     1        CALL MAKEGRID (EMEAN , WHGTTG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(2),2,NSTART,NEND,THQ)
           IF(IPFGTBL(2).EQ.IRANK) 
     1        CALL MAKEGRID (THQ , WDIRTG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(3),3,NSTART,NEND,FMEAN)
           IF(IPFGTBL(3).EQ.IRANK) 
     1        CALL MAKEGRID (FMEAN , WMNFTG, 1, NEND(NPROC), IG)

           IF(IPFGTBL(4).NE.0) THEN
              ALLOCATE(US(NIBLO))
              US = 0.
              DO IJ = IJS,IJL
                US(IJ) = USNEW(IJ)
              ENDDO
           ENDIF
           CALL MPGATHERSCFLD(IPFGTBL(4),4,NSTART,NEND,US)
           IF(IPFGTBL(4).EQ.IRANK) 
     1        CALL MAKEGRID (US, USTARG, 1, NEND(NPROC), IG)
           IF(ALLOCATED(US)) DEALLOCATE(US)

           IF(IPFGTBL(5).NE.0) THEN 
              ALLOCATE(THW(NIBLO))
              THW = 0.
              DO IJ = IJS,IJL
                 THW(IJ) = THWNEW(IJ) 
              ENDDO
           ENDIF
           CALL MPGATHERSCFLD(IPFGTBL(5),5,NSTART,NEND,THW)
           IF(IPFGTBL(5).EQ.IRANK) 
     1        CALL MAKEGRID (THW, UDIRG, 1, NEND(NPROC), IG)
           IF(ALLOCATED(THW)) DEALLOCATE(THW)

           CALL MPGATHERSCFLD(IPFGTBL(6),6,NSTART,NEND,PEAKF)
           IF(IPFGTBL(6).EQ.IRANK)
     1        CALL MAKEGRID (PEAKF , WPKFTG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(7),7,NSTART,NEND,CD)
           IF(IPFGTBL(7).EQ.IRANK)
     1        CALL MAKEGRID (CD   ,  CDG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(8),8,NSTART,NEND,WVSTR)
           IF(IPFGTBL(8).EQ.IRANK)
     1        CALL MAKEGRID (WVSTR,  TAUWG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(9),9,NSTART,NEND,SMEAN)
           IF(IPFGTBL(9).EQ.IRANK)
     1        CALL MAKEGRID (SMEAN , SMEANG, 1, NEND(NPROC), IG)

           IF(IPFGTBL(10).NE.0) THEN
              ALLOCATE(U10(NIBLO))
              U10 = 0.
              DO IJ = IJS,IJL
                U10(IJ) = U10NEW(IJ)
                IF (LICERUN.AND.ICEMASK(IJ,IG).EQ.0) U10(IJ)=ZMISS
              ENDDO
           ENDIF
           CALL MPGATHERSCFLD(IPFGTBL(10),10,NSTART,NEND,U10)
           IF(IPFGTBL(10).EQ.IRANK) 
     1        CALL MAKEGRID (U10, U10G, 1, NEND(NPROC), IG)
           IF(ALLOCATED(U10)) DEALLOCATE(U10)

         ELSE
           CALL MAKEGRID (EMEAN , WHGTTG, IJS, IJL, IG)
           CALL MAKEGRID (THQ   , WDIRTG, IJS, IJL, IG)
           CALL MAKEGRID (PEAKF , WPKFTG, IJS, IJL, IG)
           CALL MAKEGRID (FMEAN , WMNFTG, IJS, IJL, IG)
           CALL MAKEGRID (SMEAN , SMEANG, IJS, IJL, IG)
           CALL MAKEGRID (THWNEW, UDIRG , IJS, IJL, IG)
           CALL MAKEGRID (USNEW,  USTARG, IJS, IJL, IG)
           CALL MAKEGRID (CD   ,  CDG   , IJS, IJL, IG)
           CALL MAKEGRID (WVSTR,  TAUWG , IJS, IJL, IG)
           CALL MAKEGRID (U10NEW, U10G  , IJS, IJL, IG)
         ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. INTEGRATED PARAMETERS OF SWELL.
C        -------------------------------
C
      IF (CDTINTS.EQ.CDTPRO) THEN
C
C*    2.1 STORE PARAMETERS FOR ACTUAL BLOCK IN GRID ARRAYS.
C         -------------------------------------------------
C
         IF (LICERUN) THEN
            DO 2100 IJ = IJS,IJL
               IF (ICEMASK(IJ,IG).EQ.0) THEN
                  ESEA(IJ)    = ZMISS
                  ESWELL(IJ)  = ZMISS
                  THWISEA(IJ) = ZMISS
                  THSWELL(IJ) = ZMISS
                  FSEA(IJ)    = ZMISS
                  FSWELL(IJ)  = ZMISS
               ENDIF
 2100      CONTINUE
         ENDIF
C         
         IF(LMESSPASS) THEN
           CALL MPGATHERSCFLD(IPFGTBL(11),11,NSTART,NEND,ESEA)
           IF(IPFGTBL(11).EQ.IRANK) 
     1        CALL MAKEGRID (ESEA,   WHGTWG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(12),12,NSTART,NEND,ESWELL)
           IF(IPFGTBL(12).EQ.IRANK) 
     1        CALL MAKEGRID (ESWELL, WHGTSG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(13),13,NSTART,NEND,THWISEA)
           IF(IPFGTBL(13).EQ.IRANK) 
     1        CALL MAKEGRID (THWISEA,WDIRWG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(14),14,NSTART,NEND,THSWELL)
           IF(IPFGTBL(14).EQ.IRANK) 
     1        CALL MAKEGRID (THSWELL,WDIRSG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(15),15,NSTART,NEND,FSEA)
           IF(IPFGTBL(15).EQ.IRANK) 
     1        CALL MAKEGRID (FSEA,   WMNFWG, 1, NEND(NPROC), IG)

           CALL MPGATHERSCFLD(IPFGTBL(16),16,NSTART,NEND,FSWELL)
           IF(IPFGTBL(16).EQ.IRANK) 
     1        CALL MAKEGRID (FSWELL, WMNFSG, 1, NEND(NPROC), IG)

         ELSE
           CALL MAKEGRID(ESWELL, WHGTSG, IJS, IJL, IG)
           CALL MAKEGRID(THSWELL,WDIRSG, IJS, IJL, IG)
           CALL MAKEGRID(FSWELL, WMNFSG, IJS, IJL, IG)
           CALL MAKEGRID(ESEA,   WHGTWG, IJS, IJL, IG)
           CALL MAKEGRID(THWISEA,WDIRWG, IJS, IJL, IG)
           CALL MAKEGRID(FSEA,   WMNFWG, IJS, IJL, IG)
         ENDIF
      ENDIF

C
C 3. DEALLOCATE LOCAL ARRAYS
C --------------------------
      IF(ALLOCATED(PEAKF)) DEALLOCATE (PEAKF)
      IF(ALLOCATED(CD)) DEALLOCATE(CD)
      IF(ALLOCATED(WVSTR)) DEALLOCATE(WVSTR)
C
      IF(LMESSPASS) THEN
        CALL MPE_BARRIER(IERR)
        IF(IERR.LT.0)CALL MPEI_ABORT('MPE_BARRIER ERROR in OUTGRID')
      ENDIF

      RETURN
      END
