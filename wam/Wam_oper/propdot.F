      SUBROUTINE PROPDOT(SL)

C ----------------------------------------------------------------------
C
C**** *PROPDOT* - PROPAGATION DOT TERMS FROM DEPTH AND CURRENT GRADIENT.
C
C     H. GUNTHER   GKSS/ECMWF   17/02/91
C     J. BIDLOT    ECMWF        FEBRUARY 97 MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF COMMON REFDOT FOR PROPAGATION.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *PROPDOT*
C
C     METHOD.
C     -------
C
C       IN A LOOP OVER THE BLOCKS THE COMMON UBUF IS READ,
C       THE DEPTH AND CURRENT GRADIENTS ARE COMPUTED,
C       COMMON REFDOT (DEPTH AND CURRENT REFRACTION FOR THETA DOT)
C       IS COMPUTED AND WRITTEN TO MASS STORAGE (IU16).
C       IN CASE OF CURRENT REFRACTION THE COMPLETE SIGMA DOT TERM
C       IS COMPUTED AND WRITTEN TO IU16 ADDITIONALLY.
C       WRITE OPERATIONS ARE NOT DONE FOR COMMON UBUF AND REFDOT
C       IF THIS IS A ONE BLOCK MODEL.
C
C     EXTERNALS.
C     ----------
C
C       *GRADI*     - COMPUTES DEPTH AND CURRENT GRADIENTS.
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
#include "comcurr.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "comrefd.h"
C
#include "comshal.h"
C
#include "comstat.h"
C
#include "comubuf.h"
C
#include "comunit.h"
C
C ----------------------------------------------------------------------
C
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: SL
C
C ----------------------------------------------------------------------
C
      DIMENSION DDPHI(NIBLD), DDLAM(NIBLD), DUPHI(NIBLC), DULAM(NIBLC),
     1          DVPHI(NIBLC), DVLAM(NIBLC), DCO(NIBLD), OMDD(NIBLC)
C
C ----------------------------------------------------------------------
C
      CHARACTER*30 FILENAME

C*    1. IF CARTESIAN PROPAGATION SET COSINE OF LAT TO 1.
C         -----------------------------------------------
C
 1000 CONTINUE
      IF (ICASE.NE.1) THEN
         DO 1001 IJ = 1,NIBLD
            DCO(IJ) = 1.
 1001    CONTINUE
      ENDIF
C
C*    2. LOOP OVER BLOCKS.
C        -----------------
C
      IF (IGL.NE.1) REWIND IU08
      DO 2000 IG = 1,IGL
C
C*    2.1 IF MULTI BLOCK VERSION.
C         -----------------------
C
         IF (IGL.NE.1) THEN
C
C*    2.1.1 READ COMMON BLOCK UBUF.
C           -----------------------
C
            READ (IU08) KLAT, KLON
CSHALLOW
C
C*    2.1.2 COMPUTE SHALLOW WATER TABLE INDICES.
C           ------------------------------------
C
            IF (ISHALLO.NE.1) THEN
               DO 2121 IJ=NINF,IJLT(IG)
                  XD = LOG(DEPTH(IJ,IG)/DEPTHA)/LOG(DEPTHD)+1.
                  ID = NINT(XD)
                  ID = MAX(ID,1)
                  INDEP(IJ) = MIN(ID,NDEPTH)
 2121          CONTINUE
            ENDIF
CSHALLOW
         ENDIF
C
C*    2.2 DEPTH AND CURRENT GRADIENTS.
C         ----------------------------
C
         CALL GRADI (IG, IREFRA, DDPHI, DDLAM, DUPHI,
     1               DULAM, DVPHI, DVLAM)
C
C*    2.3 COSINE OF LATITUDES IF SPHERICAL PROPAGATION.
C         ---------------------------------------------
C
         IF (ICASE.EQ.1) THEN
            DO 2301 IJ = IJS(IG),IJL(IG)
               JH = KXLT(IJ,IG)
               DCO(IJ) = 1./COSPH(JH)
 2301       CONTINUE
         ENDIF
C
C*    2.4 DEPTH GRADIENT PART OF SIGMA DOT.
C         ---------------------------------
C
         IF (ISHALLO.NE.1 .AND. IREFRA.EQ.2) THEN
            DO 2401 IJ = IJS(IG),IJL(IG)
               OMDD(IJ) = V(IJ,IG)*DDPHI(IJ) +
     1                    U(IJ,IG)*DDLAM(IJ)*DCO(IJ)
 2401       CONTINUE
         ENDIF
C
C*    2.5. LOOP OVER DIRECTIONS.
C          ---------------------
C
         DO 2501 K=1,NANG
            SD = SINTH(K)
            CD = COSTH(K)
C
C*    2.5.1. DEPTH GRADIENT OF THETA DOT.
C            ----------------------------
C
            IF (ISHALLO.NE.1) THEN
               DO 2511 IJ = IJS(IG),IJL(IG)
                  THDD(IJ,K) = SD*DDPHI(IJ) - CD*DDLAM(IJ)*DCO(IJ)
 2511          CONTINUE
            ENDIF
C
C*    2.5.2 SIGMA DOT AND THETA DOT PART FROM CURRENT GRADIENT.
C           ---------------------------------------------------
C
            IF (IREFRA.EQ.2) THEN
               SS  = SD**2
               SC  = SD*CD
               CC  = CD**2
               DO 2521 IJ = IJS(IG),IJL(IG)
                  SL(IJ,K,NFRE) = -SC*DUPHI(IJ) - CC*DVPHI(IJ)
     1                          - (SS*DULAM(IJ) + SC*DVLAM(IJ))*DCO(IJ)
                  THDC(IJ,K) =  SS*DUPHI(IJ) + SC*DVPHI(IJ)
     1                       - (SC*DULAM(IJ) + CC*DVLAM(IJ))*DCO(IJ)
 2521          CONTINUE
C
C*    2.5.3 LOOP OVER FREQUENCIES.
C           ----------------------
C
               IF (ISHALLO.NE.1) THEN
               DO 2530 M=1,NFRE
                  DO 2531 IJ=IJS(IG),IJL(IG)
                     SL(IJ,K,M) = (SL(IJ,K,NFRE)*TCGOND(INDEP(IJ),M)
     1                          + OMDD(IJ)*TSIHKD(INDEP(IJ),M))
     2                          * TFAK(INDEP(IJ),M)
 2531             CONTINUE
C
C*    BRANCH BACK TO 2.5.3 FOR NEXT FREQUENCY.
C
 2530          CONTINUE
               ENDIF
            ENDIF
C
C*    BRANCH BACK TO 2.5 FOR NEXT DIRECTION.
C
 2501    CONTINUE
C
C*    2.6. WRITE COMPUTED DOT TERMS.
C           ------------------------
C
         IF (IGL.NE.1) THEN
            IF (ISHALLO.NE.1) THEN
               IF (IREFRA.EQ.1) THEN
                  WRITE(IU16) THDD
               ELSE
                  WRITE(IU16) THDD, THDC, SL
               ENDIF
            ELSE
               WRITE(IU16) THDC, ((SL(IJ,K,NFRE),IJ=1,NIBLC),K=1,NANG)
            ENDIF
         ELSE
            IF(LMESSPASS) THEN
              FILENAME='fort.16.%p'
              CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAME,1)
              OPEN(IU16,FILE=FILENAME,STATUS='UNKNOWN')

              WRITE(IU16) SL
              CLOSE(IU16)
            ELSE
              IF (ISHALLO.NE.1) THEN
                IF (IREFRA.EQ.2) THEN
                   WRITE(IU16) SL
                ENDIF
              ELSE
                WRITE(IU16) ((SL(IJ,K,NFRE),IJ=1,NIBLC),K=1,NANG)
              ENDIF
            ENDIF
         ENDIF
C
C*    BRANCH BACK TO 2. FOR NEXT BLOCK.
C
 2000 CONTINUE

      RETURN
      END
