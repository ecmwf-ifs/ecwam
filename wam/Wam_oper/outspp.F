      SUBROUTINE OUTSPP (FL3, FL1, IG, IU25, IU26,NSTART,NEND,
     1                   U10NEW,THWNEW,USNEW)

C ----------------------------------------------------------------------
C
C**** *OUTSPP* -  OUTPUT OF SPECTRA , WAVE HEIGHT , WAVE DIRECTION,
C****             USTAR AND WIND DIRECTION AT GIVEN OUTPUT POINTS.
C
C      S.D.HASSELMANN.
C      P.JANSSEN       KNMI         FEBRUARY 1986
C
C      P. LIONELLO                  FEBRUARY 1987
C                  OUTPUT OF SWELL 2-D DISTRIBUTION , SWELL WAVE HEIGHT
C                  MEAN SWELL FREQUENCY , MEAN SWELL DIRECTION ,WIND-SE
C                  WAVE HEIGTH AND WIND SEA WAVE DIRECTION.
C
C      L. ZAMBRESKY   GKSS/ECMWF    JULY 89
C                 VARIANCE ENERGY SUMMED OVER FREQUENCY AND DIRECTION.
C
C      H. GUNTHER     ECMWF         OCTOBER 1990
C                 EXTENDED HEADER AND OUTPUT OF 1D SPECTRA IN PRINT.
!      J. BIDLOT      ECMWF         APRIL 1996 MESSAGE PASSING
!      J. BIDLOT      ECMWF         JANURARY 1997 MODIFY DEFINITION OF TH0 
C
C**   INTERFACE.
C     ----------
C
C        *CALL  OUTSPP (FL3, FL1, IG, IU25, IU26,,NSTART,NEND,
C                       U10NEW,THWNEW,USNEW)*
C           *FL3*    -   BLOCK OF SPECTRA.
C           *FL1*    -   SWELL SPECTRUM.
C           *IG*     -   BLOCK NUMBER.
C           *IU25*   -   OUTPUT UNIT FOR SPECTRA.
C           *IU26*   -   OUTPUT UNIT FOR SEA AND SWELL SPECTRA.
C           *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C           *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C           *U10NEW*    NEW WIND SPEED IN M/S.
C           *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
C                       NOTATION (POINTING ANGLE OF WIND VECTOR,
C                       CLOCKWISE FROM NORTH).
C           *USNEW*     NEW FRICTION VELOCITY IN M/S.
C
C     EXTERNALS.
C     ----------
C
C       *PRSPP*     - PRINT A SPECTRUM.
C       *PRSPPS*    - PRINT A SWELL SPECTRUM.
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
#include "param.h"
C
#include "comfred.h"
C
#include "commap.h"
C
#include "commean.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "commesspass.h"
C
#include "comcout.h"
C
#include "comstat.h"
C
#include "comswel.h"
C
#include "comtest.h"
C
#include "comwind.h"
C
C ----------------------------------------------------------------------
C
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL3,FL1
      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW

C ----------------------------------------------------------------------



#include "parcons.h"
C
      PARAMETER ( CO = 1.1)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C        *CO*      REAL      FREQUENCY RATIO.
C
      CHARACTER NAME*20
      DIMENSION SPEC(NANG,NFRE)
C
C ----------------------------------------------------------------------
C  
C
C     WORK ARRAYS
      REAL,ALLOCATABLE :: FLPTS(:,:,:,:)
      REAL,ALLOCATABLE,DIMENSION(:) :: U10,THW,US

C     0. GATHER THE FIELDS ON PROCESS 1 
C        ------------------------------
C
      IF(NGOUT.GE.1) THEN 
        ALLOCATE(U10(NIBLO),THW(NIBLO),US(NIBLO))
        DO IJ=NSTART(IRANK),NEND(IRANK)
           U10(IJ)=U10NEW(IJ)
           THW(IJ)=THWNEW(IJ)
           US(IJ)=USNEW(IJ)
        ENDDO
      ENDIF 
      IF((LMESSPASS).AND.(NGOUT.GE.1)) THEN 
        IRECV=1
        NSPFLD=2 
        NSCFLD=11 
        ALLOCATE(FLPTS(NSPFLD,NGOUT,NANG,NFRE))
        CALL MPGATHERSPP(IRECV,25,NSTART,NEND,NSPFLD,NSCFLD,
     1                   FL3,FL1,FLPTS,
     2                   EMEAN,FMEAN,THQ,SMEAN,ESWELL,FSWELL,THSWELL,
     3                   THWISEA,U10,THW,US)

      ENDIF
C     
C*    1. LOOP OVER OUTPUT POINTS.
C        ------------------------
C
      IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
      DO 1001 NGOU=1,NGOUT
         IF (IG.EQ.IGAR(NGOU)) THEN
            IJ = IJAR(NGOU)
            IX = IXLG(IJ,IG)
            KX = KXLT(IJ,IG)
            XLON = AMOWEP + REAL(IX-1)*ZDELLO(KX)
            XLAT = AMOSOP + REAL(KX-1)*XDELLA
C
C*    1.1 TOTAL SPECTRUM.
C         ---------------
C
            IF (CDTSPT.EQ.CDTPRO) THEN
C
C*    1.1.1 FILE OUTPUT TO IU25.
C           --------------------
C
               IF (FFLAG(17)) THEN
                  XANG = REAL(NANG)
                  XFRE = REAL(NFRE)
                  TH0 = TH(1)*DEG
                  WRITE(IU25) XLON, XLAT, CDTPRO, XANG, XFRE,
     1                        TH0, FR(1), CO
                  WRITE(IU25) 4.*SQRT(EMEAN(IJ)), DEG*THQ(IJ),
     1                        FMEAN(IJ), US(IJ), DEG*THW(IJ),
     2                        U10(IJ)
                  IF(LMESSPASS) THEN
                    WRITE(IU25) ((FLPTS(1,NGOU,K,M),K=1,NANG),M=1,NFRE)
                  ELSE
                    WRITE(IU25) ((FL3(IJ,K,M),K=1,NANG),M=1,NFRE)
                  ENDIF
               ENDIF
C
C*    1.1.2 PRINT SPECTRUM.
C           ---------------
C
               IF (PFLAG(17)) THEN
                  DO M=1,NFRE
                    DO K=1,NANG
                      IF(LMESSPASS) THEN
                      SPEC(K,M) = FLPTS(1,NGOU,K,M)
                      ELSE
                      SPEC(K,M) = FL3(IJ,K,M)
                      ENDIF
                    ENDDO
                  ENDDO
                  WRITE (NAME,'(''BL= '',I3,'' POINT= '',I4)') IG,IJ
                  CALL PRSPP (IU06, CDTPRO, XLON, XLAT, NAME,
     1                   4.*SQRT(EMEAN(IJ)), DEG*THQ(IJ),
     2                   FMEAN(IJ), US(IJ), DEG*THW(IJ),
     3                   U10(IJ),
     4                   FR, TH, SPEC, NANG, NFRE, NANG, NFRE)
               ENDIF
            ENDIF
CPL****
C*    1.2 OUTPUT SWELL INFORMATION.
C         -------------------------
C
            IF (CDTSPS.EQ.CDTPRO) THEN
               EWSEA=EMEAN(IJ)-ESWELL(IJ)
               FWSEA=EWSEA/(EMEAN(IJ)/FMEAN(IJ)-ESWELL(IJ)/FSWELL(IJ))
               EWSEA = 4.*SQRT(EWSEA)
                  IF(LMESSPASS) THEN
                    DO M=1,NFRE
                      DO K=1,NANG
                        SPEC(K,M) = FLPTS(2,NGOU,K,M)*FLPTS(1,NGOU,K,M)
                      ENDDO
                    ENDDO
                  ELSE
                    DO M=1,NFRE
                      DO K=1,NANG
                        SPEC(K,M) = FL1(IJ,K,M)*FL3(IJ,K,M)
                      ENDDO
                    ENDDO
                  ENDIF
C
C*    1.2.1 FILE OUTPUT TO IU26.
C           --------------------
C
               IF (FFLAG(18)) THEN
                  XANG = REAL(NANG)
                  XFRE = REAL(NFRE)
                  TH0 = TH(1)*DEG
                  WRITE(IU26) XLON, XLAT, CDTPRO, XANG, XFRE,
     1                        TH0, FR(1), CO
                  WRITE(IU26) EWSEA, 4.*SQRT(ESWELL(IJ)),
     1                        DEG*THWISEA(IJ), DEG*THSWELL(IJ),
     2                        FWSEA, FSWELL(IJ)
                  WRITE(IU26) ((SPEC(K,M),K=1,NANG),M=1,NFRE)
               ENDIF
C
C*    1.2.2 PRINT OUTPUT TO IU26.
C           ---------------------
C
               IF (PFLAG(18)) THEN
                  WRITE (NAME,'(''BL= '',I3,'' POINT= '',I4)') IG,IJ
                  CALL PRSPPS (IU06, CDTPRO, XLON, XLAT, NAME,
     1               EWSEA, DEG*THWISEA(IJ), FWSEA,
     2               4.*SQRT(ESWELL(IJ)), DEG*THSWELL(IJ),
     3               FSWELL(IJ), FR, TH, SPEC, NANG, NFRE, NANG, NFRE)
               ENDIF
            ENDIF
         ENDIF
 1001 CONTINUE
      ENDIF
      IF((LMESSPASS).AND.(NGOUT.GE.1)) DEALLOCATE(FLPTS)
      IF(NGOUT.GE.1) DEALLOCATE(U10,THW,US)

      RETURN
      END
