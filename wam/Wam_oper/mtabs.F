      SUBROUTINE MTABS (ML, KL)

C ----------------------------------------------------------------------
C
C**** *MTABS* - ROUTINE TO COMPUTE TABLES USED FOR SHALLOW WATER.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       TO COMPUTE TABLES USED FOR SHALLOW WATER.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MTABS (ML, KL)*
C          *ML*      - NUMBER OF FREQUENCIES.
C          *KL*      - NUMBER OF DIRECTIONS.
C
C     METHOD.
C     -------
C
C       TABLES FOR GROUP VELOCITY, WAVE NUMBER AND OMEGA/SINH(2KD)
C       ARE COMPUTED AT ALL FREQUENCIES AND FOR A DEPTH TABLE
C       OF LENGTH NDEPTH, STARTING AT DEPTHA METERS AND INCREMENTED
C        BY DEPTHD METRES.
C
C     EXTERNALS.
C     ----------
C
C       *AKI*       - FUNCTION TO COMPUTE WAVE NUMBER.
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
#include "comshal.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
#include "parcons.h"
C
      DEPTHA = 5.
      DEPTHD = 1.1
C
C ----------------------------------------------------------------------
C
C*    1. GROUP VELOCITY AND WAVE NUMBER.
C        -------------------------------
C
 1000 CONTINUE
C
C*    1.1 LOOP OVER FREQUENCIES.
C         ----------------------
C
 1100 CONTINUE
      GH = G/(4.*PI)
      DO 1101 M=1,ML
         OM=ZPI*FR(M)
C
C*    1.1.1 LOOP OVER DEPTH.
C           -----------------
C
 1110 CONTINUE
         DO 1111 JD=1,NDEPTH
            AD=DEPTHA*DEPTHD**(JD-1)
            AK=AKI(OM,AD)
            TFAK(JD,M)=AK
            AKD=AK*AD
            IF(AKD.LE.10.0) THEN
               TCGOND(JD,M) = 0.5*SQRT(G*TANH(AKD)/AK)*
     1                       (1.0+2.0*AKD/SINH(2.*AKD))
               TSIHKD(JD,M) = OM/SINH(2.*AKD)
            ELSE
               TCGOND(JD,M) = GH/FR(M)
               TSIHKD(JD,M) = 0.
            ENDIF
 1111    CONTINUE
 1101 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. PRINT TABLES.
C        -------------
C
 2000 CONTINUE
      NAN  = 10
      NSTP = NDEPTH/NAN
      NSTP = MAX(NSTP,1)
      DEPTHE = DEPTHA*DEPTHD**(NDEPTH-1)
      WRITE (IU06,'(1H1, '' SHALLOW WATER TABLES:'',/)')
      WRITE (IU06,'(''  LOGARITHMIC DEPTH FROM: DEPTHA = '',F5.1,
     1  '' TO DEPTHE  = '',F5.1, ''IN STEPS OF DEPTHD = '',F5.1)')
     2    DEPTHA, DEPTHE, DEPTHD
      WRITE (IU06,'(''  PRINTED IN STEPS OF '',I3,'' ENTRIES'',/)') NSTP
      DO 2001 JD=1,NDEPTH,NSTP
         AD=DEPTHA*DEPTHD**(JD-1)
         WRITE (IU06,'(1X,''DEPTH = '',F7.1,'' METRES '')') AD
         WRITE (IU06,'(1X,''GROUP VELOCITY IN METRES/SECOND'')')
         WRITE (IU06,'(1x,13F10.5)') (TCGOND(JD,M),M=1,ML)
         WRITE (IU06,'(1X,''WAVE NUMBER IN 1./METRES'')')
         WRITE (IU06,'(1x,13F10.5)') (TFAK(JD,M),M=1,ML)
         WRITE (IU06,'(1X,''OMEGA/SINH(2KD) IN 1./SECOND'')')
         WRITE (IU06,'(1x,13F10.5)') (TSIHKD(JD,M),M=1,ML)
 2001 CONTINUE

      RETURN
      END
