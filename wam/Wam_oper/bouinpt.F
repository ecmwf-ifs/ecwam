      SUBROUTINE BOUINPT (FL3, IG, IU02)

C ----------------------------------------------------------------------
C
C**** *BOUINPT* - BOUNDARY VALUE INPUT INTO A BLOCK.
C
C     H. GUNTHER    GKSS/ECMWF         JANUARY 1991
C
C*    PURPOSE.
C     --------
C
C       INPUT AND SPACE INTERPOLATION OF BOUNDARY SPECTRA.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *BOUINPT (FL3, IG, IU02)*
C          *FL3*     REAL         BLOCK OF SPECTRA.
C          *IG*      INTEGER      BLOCK NUMBER.
C          *IU02*    INTEGER      UNIT FOR INPUT OF BOUNDARY VALUES.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT*     - TERMINATES PROCESSING.
C       *GSFILE*    - GETS OR SAVES A FILE.
C       *INTSPEC*   - INTERPOLATE A SPECTRUM.
C
C     METHOD.
C     -------
C
C       IN THE FIRST CALL OF THE SUB. THE FILE HEADER IS READ
C       AND THE CONSISTENCY IS CHECKED. THE SUB. READS A COMPLETE
C       SET OF BOUNDARY VALUES EACH PROPAGATION TIMESTEP WHEN
C       IT IS CALLED FOR THE FIRST BLOCK. THE SPECTRA REQUIRED FOR
C       THE ACTUAL BLOCK ARE INTERPOLATED AND STORED IN THE BLOCK.
C       INDICES AND WEIGHTS NECESSARY FOR THE INTERPOLATION AND
C       STORAGE ARE PRECOMPUTED IN PROG. PREPROC.
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
#include "comfbou.h"
C
#include "comfred.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION FL3(NINF-1:NSUP,NANG,NFRE)
C
C ----------------------------------------------------------------------
C
      DIMENSION F1(NANG,NFRE,0:NBINP), FMEAN1(0:NBINP), EMEAN1(0:NBINP),
     1          THQ1(0:NBINP), FL(NANG,NFRE)
      SAVE F1, FMEAN1, EMEAN1, THQ1, NBOINP
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *F1*        REAL      SPECTRA FROM COARSE GRID.
C      *FMEAN1*    REAL      MEAN FREQUENCIES FROM COARSE GRID.
C      *EMEAN1*    REAL      TOTAL ENERGIES FROM COARSE GRID.
C      *THQ1*      REAL      MEAN DIRECTIONS FROM COARSE GRID (RAD).
C      *NBOINP*    REAL      NUMBER OF INPUT SPECTRA.
C      *FL*        REAL      INTERPOLATED SPECTRUM.
C
      CHARACTER*12 CDTLAST, CDATE1, ZERO
      DATA IJAF, CDTLAST /0, '000000000000'/  
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *IJFA*      REAL      NUMBER OF FIRST INPUT POINT IN BLOCK.
C      *CDTLAST*   CHAR*12   DATE OF LAST FETCH FILE.
C
C
C ----------------------------------------------------------------------
C
C*    1. INITIALIZATION.
C        ---------------
C
 1000 CONTINUE
C
C*    1.1 FIRST CALL OF SUB?
C         ------------------
C
      ZERO = ' '
      IF (CDTLAST.EQ.ZERO) THEN
C
C*    1.1.1 INITIALISE FOR LANDPOINTS.
C           --------------------------
C
         DO 1111 M=1,NFRE
            DO 1112 K=1,NANG
               F1(K,M,0) = 0.
 1112       CONTINUE
 1111    CONTINUE
         FMEAN1(0) = 0.
         EMEAN1(0) = 0.
         THQ1(0) = 0.
      ENDIF
C
C*    1.2 HAS A NEW BOUNDARY VALUE INPUT FILE TO BE FETCH?
C         ------------------------------------------------
C
      IF (CDTLAST.LT.CDTRES) THEN
C
C*    1.2.1 FETCH INPUT FILE.
C           -----------------
C
         CDTLAST = CDTRES
         CALL GSFILE (IU06, IU02, 0, CDTRES, CDATEF, 'FBI', 'G')
C
C*    1.2.2 READ BOUNDARY FILE HEADER.
C           --------------------------
C
         READ (IU02, ERR=5000, END=5000)
     1        XANG, XFRE, TH0, FR1, CO, XNBO, XDELIN
         KL = NINT(XANG)
         ML = NINT(XFRE)
         NBOINP  = NINT(XNBO)
         IDELINP = NINT(XDELIN)

         IF (ITEST.GT.3) THEN
            WRITE (IU06,*) ' '
            WRITE (IU06,*) ' BOUNDARY VALUE INPUT FILE HEADER IS:'
            WRITE (IU06,*) ' NO. OF DIRECTIONS IS    KL     = ', KL
            WRITE (IU06,*) ' NO. OF FREQUENCIES IS   ML     = ', ML
            WRITE (IU06,*) ' FIRST DIRECTION IS      TH0    = ', TH0
            WRITE (IU06,*) ' FIRST FREQUENCY IS      FR1    = ', FR1
            WRITE (IU06,*) ' FREQUENCY RATIO IS      CO     = ', CO
            WRITE (IU06,*) ' NO. OF BOUNDRAY POINTS  NBOINP = ', NBOINP
            WRITE (IU06,*) ' TIME STEP OF DATA IS    IDELINP= ', IDELINP
         ENDIF
C
C*    1.2.3 CHECK CONSISTENCY.
C           ------------------
C
         IF (KL.NE.NANG .OR. ML.NE.NFRE .OR. NBOINP.GT.NBINP .OR.
     1      FR1.NE.FR(1) .OR. TH0.NE.TH(1) .OR. IDELINP.NE.IDELPRO) THEN
            WRITE (IU06,*) '****************************************'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '*    FATAL ERROR SUB. BOUINP           *'
            WRITE (IU06,*) '*    =======================           *'
            WRITE (IU06,*) '* VALUES IN BOUNDARY FILE HEADER ARE   *'
            WRITE (IU06,*) '* INCONSISTENT WITH MODEL SET-UP.      *'
            WRITE (IU06,*) '* MODEL VALUES ARE:                    *'
            WRITE (IU06,*) '* NO. OF DIRECTIONS     NANG   = ', NANG
            WRITE (IU06,*) '* NO. OF FREQUENCIES    NFRE   = ', NFRE
            WRITE (IU06,*) '* FIRST DIRECTION       TH0    = ', TH0
            WRITE (IU06,*) '* FIRST FREQUENCY       FR(1)  = ', FR(1)
            WRITE (IU06,*) '* DIMENSION             NBINP  = ', NBINP
            WRITE (IU06,*) '* PROPAGATION TIMESTEP  IDELPRO= ', IDELPRO
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.    *'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '****************************************'
            CALL ABORT1
         ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. CALL OF SUB FOR FIRST BLOCK?
C        ----------------------------
C
 2000 CONTINUE
      IF (IG.EQ.1) THEN
         IJAF = 1
C
C*    2.1 READ BOUNDARY SPECTRA.
C         ----------------------
C
 2100 CONTINUE
         DO 2101 IJ=1,NBOINP
            READ (IU02, ERR=5001, END=5001)
     1           DUM1, DUM2, CDATE1, EMEAN1(IJ),
     2           THQ1(IJ), FMEAN1(IJ), DUM3, DUM4
            READ (IU02, ERR=5002, END=5002)
     1           ((F1(K,M,IJ),K=1,NANG),M=1,NFRE)
 2101    CONTINUE
         IF (ITEST.GT.3) THEN
            WRITE (IU06,*) ' '
            WRITE (IU06,*) ' BOUNDARY VALUES READ CDATE1 = ', CDATE1
         ENDIF
C
C*    2.2 CHECK DATES.
C         ------------
C
         IF (CDATE1.LT.CDTPRO) THEN
            WRITE (IU06,*) '****************************************'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '*    WARNING ERROR SUB. BOUINP         *'
            WRITE (IU06,*) '*    =========================         *'
            WRITE (IU06,*) '* BOUNDARY VALUE INPUT DATE IS BEFORE  *'
            WRITE (IU06,*) '* MODEL DATE.                          *'
            WRITE (IU06,*) '* MODEL DATE IS          CDTPRO = ', CDTPRO
            WRITE (IU06,*) '* BOUNDARY VALUE DATE IS CDATE1 = ', CDATE1
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '* PROGRAM WILL TRY NEXT INPUT.         *'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '****************************************'
            GOTO 2100
         ELSE IF (CDATE1.GT.CDTPRO) THEN
            WRITE (IU06,*) '****************************************'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '*    FATAL ERROR SUB. BOUINP           *'
            WRITE (IU06,*) '*    =======================           *'
            WRITE (IU06,*) '* DATES DO NOT MATCH.                  *'
            WRITE (IU06,*) '* MODEL DATE IS          CDTPRO = ', CDTPRO
            WRITE (IU06,*) '* BOUNDARY VALUE DATE IS CDATE1 = ', CDATE1
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.    *'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '****************************************'
            CALL ABORT1
         ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. LOOP OVER BOUNDARY POINTS.
C        --------------------------
C
 3000 CONTINUE
      DEL12 = 1.
      DO 3001 IJ=IJAF,NBOUNF
C
C*    3.1 GET INDICES.
C         ------------
C
         IGF = IGARF(IJ)
         IF (IGF.EQ.IG) THEN
            IJF = IJARF(IJ)
            IBCL = IBFL(IJ)
            IBCR = IBFR(IJ)
            DEL1L = BFW(IJ)
C
C*    3.2 CHECK INTERPOLATION WEIGHT.
C         ---------------------------
C
            IF (DEL1L.GT.0.) THEN
C
C*    3.2.1. WEIGHT IS GT ZERO INTERPOLATE.
C            ------------------------------
C
               CALL INTSPEC (NFRE, NANG, NFRE, NANG, FR, DEL12, DEL1L,
     1            F1(1,1,IBCL), FMEAN1(IBCL), EMEAN1(IBCL), THQ1(IBCL),
     2            F1(1,1,IBCR), FMEAN1(IBCR), EMEAN1(IBCR), THQ1(IBCR),
     3            FL, FMEAN, EMEAN, THQ)
               DO 3201 M=1,NFRE
                  DO 3202 K=1,NANG
                     FL3(IJF,K,M) = FL(K,M)
 3202             CONTINUE
 3201          CONTINUE
            ELSE
C
C*    3.2.2. WEIGHT IS ZERO COPY LEFT POINT.
C            -------------------------------
C
               DO 3301 M=1,NFRE
                  DO 3302 K=1,NANG
                     FL3(IJF,K,M) = F1(K,M,IBCL)
 3302             CONTINUE
 3301          CONTINUE
            ENDIF
         ELSE IF (IGF.GT.IG) THEN
            IJAF = IJ
            GOTO 4000
         ENDIF
 3001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    4. RETURN TO CALLING PROG.
C        -----------------------
C
 4000 CONTINUE

      RETURN
C
C ----------------------------------------------------------------------
C
C*    5. ERROR MESSAGES.
C        ---------------
C
 5000 CONTINUE
         WRITE (IU06,*) '*******************************************'
         WRITE (IU06,*) '*                                         *'
         WRITE (IU06,*) '*      FATAL ERROR SUB. BOUNINP.          *'
         WRITE (IU06,*) '*      =========================          *'
         WRITE (IU06,*) '* END OF FILE OR READ ERROR.              *'
         WRITE (IU06,*) '* PROGRAM TRIES TO READ                   *'
         WRITE (IU06,*) '* HEADER OF BOUNDARY VALUES               *'
         WRITE (IU06,*) '* UNIT IS IU02 = ', IU02
         WRITE (IU06,*) '*                                         *'
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE (IU06,*) '*                                         *'
         WRITE (IU06,*) '*******************************************'
         CALL ABORT1
 5001 CONTINUE
         WRITE (IU06,*) '*******************************************'
         WRITE (IU06,*) '*                                         *'
         WRITE (IU06,*) '*      FATAL ERROR SUB. BOUNINP.          *'
         WRITE (IU06,*) '*      =========================          *'
         WRITE (IU06,*) '* END OF FILE OR READ ERROR.              *'
         WRITE (IU06,*) '* PROGRAM TRIES TO READ                   *'
         WRITE (IU06,*) '* A SPECTRUM HEADER                       *'
         WRITE (IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
         WRITE (IU06,*) '* UNIT IS          IU02 = ', IU02
         WRITE (IU06,*) '*                                         *'
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE (IU06,*) '*                                         *'
         WRITE (IU06,*) '*******************************************'
         CALL ABORT1
 5002 CONTINUE
         WRITE (IU06,*) '*******************************************'
         WRITE (IU06,*) '*                                         *'
         WRITE (IU06,*) '*      FATAL ERROR SUB. BOUNINP.          *'
         WRITE (IU06,*) '*      =========================          *'
         WRITE (IU06,*) '* END OF FILE OR READ ERROR.              *'
         WRITE (IU06,*) '* PROGRAM TRIES TO READ A SPECTRUM        *'
         WRITE (IU06,*) '* DATE IS        CDATE1 = ', CDATE1
         WRITE (IU06,*) '* SPECTRA COUNTER IS IJ = ', IJ
         WRITE (IU06,*) '* UNIT IS          IU02 = ', IU02
         WRITE (IU06,*) '*                                         *'
         WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.       *'
         WRITE (IU06,*) '*                                         *'
         WRITE (IU06,*) '*******************************************'
         CALL ABORT1

      END
