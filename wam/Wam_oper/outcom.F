      SUBROUTINE OUTCOM (IU07, IU17, IFORM, IREFRA, ML, KL)

C ----------------------------------------------------------------------
C
C**** *OUTCOM* - ROUTINE TO WRITE COMMON TO DISK
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       TO WRITE OUT THE COMPUTED COMMON BLOCKS
C       (COMMON UBUF IS WRITTEN IN OUTUBUF)
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *OUTCOM (IU07, IU17, IFORM, IREFRA, ML, KL)*
C          *IU07*   - LOGICAL UNIT FOR  UNFORMATED WRITE.
C          *IU17*   - LOGICAL UNIT FOR    FORMATED WRITE.
C          *IFORM*   - FORMAT OPTION  = 1  UNFORMATED WRITE.
C                                     = 2  FORMATED WRITE.
C                                     OTHERWISE BOTH.
C          *IREFRA*  - REFRACTION OPTION.
C          *ML*      - NUMBER OF FREQUENCIES.
C          *KL*      - NUMBER OF DIRECTIONS.
C
C     METHOD.
C     -------
C
C       COMMON BLOCKS COUPLE, CURRENT, FREDIR, INDNL, GRIDPAR, MAP,
C       COUT, TABLE, AND SHALLOW ARE WRITTEN TO UNIT.
C       ALL FREQUENCY AND DIRECTION DEPENDENT ARRAYS
C       ARE WRITTEN FROM 1 TO THE USED NUMBER OF FREQUENCIES (ML),
C       AND THE USED NUMBER OF DIRECTIONS (KL). OTHER ARRAYS ARE
C       WRITTEN ACCORDING TO THEIR DIMENSIONS.
C
C     EXTERNALS.
C     ----------
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
#include "comcoup.h"
C
#include "comcurr.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "comindn.h"
C
#include "commap.h"
C
#include "comcout.h"
C
#include "comshal.h"
C
#include "comtabl.h"
C
C ----------------------------------------------------------------------
C
  998 FORMAT(10I8)
  999 FORMAT(5E16.7)
C
C ----------------------------------------------------------------------
C
C*    1. WRITE COMMON FREDIR.
C        --------------------
C
 1000 CONTINUE
      IF (IFORM.NE.2) THEN
         WRITE (IU07) (FR(M),M=1,ML), (DFIM(M),M=1,ML),
     1                (GOM(M),M=1,ML), (C(M),M=1,ML),
     2                DELTH, DELTR, (TH(K),K=1,KL),
     3                (COSTH(K),K=1,KL), (SINTH(K),K=1,KL)
      ENDIF
      IF (IFORM.NE.1) THEN
         WRITE (IU17,999) (FR(M),M=1,ML), (DFIM(M),M=1,ML),
     1                    (GOM(M),M=1,ML), (C(M),M=1,ML),
     2                    DELTH, DELTR, (TH(K),K=1,KL),
     3                    (COSTH(K),K=1,KL), (SINTH(K),K=1,KL)
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. WRITE COMMON GRIDPAR.
C        ---------------------
C
 2000 CONTINUE
      IF (IFORM.NE.2) THEN
         WRITE (IU07) DELPHI, (DELLAM(L),L=1,NY), (NLONRGG(L),L=1,NY),
     1                (SINPH(L),L=1,NY), (COSPH(L),L=1,NY),
     2                IGL, IJS, IJL2, IJLS, IJL, IJLT
      ENDIF
      IF (IFORM.NE.1) THEN
         WRITE (IU17,999) DELPHI,(DELLAM(L),L=1,NY),(NLONRGG(L),L=1,NY),
     1                    (SINPH(L),L=1,NY), (COSPH(L),L=1,NY)
         WRITE (IU17,998) IGL, IJS, IJL2, IJLS, IJL, IJLT
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. WRITE COMMON MAP.
C        -----------------
C
 3000 CONTINUE
      IF (IFORM.NE.2) THEN
         WRITE (IU07) IXLG, KXLT, NX, NY, IPER,
     1                AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,
     2                ZDELLO, IRGG
      ENDIF
      IF (IFORM.NE.1) THEN
         WRITE (IU17,998) IXLG, KXLT, NX, NY, IPER
         WRITE (IU17,999) AMOWEP, AMOSOP, AMOEAP, AMONOP,
     1                    XDELLA, XDELLO,ZDELLO
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    4. WRITE COMMON INDNL.
C        -------------------
C
 4000 CONTINUE
      IF (IFORM.NE.2) THEN
         WRITE (IU07) (IKP(M),M=1,ML+4), (IKP1(M),M=1,ML+4),
     1                (IKM(M),M=1,ML+4), (IKM1(M),M=1,ML+4),
     2                ((K1W(K,L),K=1,KL),L=1,2),
     3                ((K2W(K,L),K=1,KL),L=1,2),
     4                ((K11W(K,L),K=1,KL),L=1,2),
     5                ((K21W(K,L),K=1,KL),L=1,2),
     6                (AF11(M),M=1,ML+4), (FKLAP(M),M=1,ML+4),
     7                (FKLAP1(M),M=1,ML+4), (FKLAM(M),M=1,ML+4),
     8                (FKLAM1(M),M=1,ML+4),
     9                ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
      ENDIF
      IF (IFORM.NE.1) THEN
         WRITE (IU17,998) (IKP(M),M=1,ML+4), (IKP1(M),M=1,ML+4),
     1                    (IKM(M),M=1,ML+4), (IKM1(M),M=1,ML+4),
     2                    ((K1W(K,L),K=1,KL),L=1,2),
     3                    ((K2W(K,L),K=1,KL),L=1,2),
     4                    ((K11W(K,L),K=1,KL),L=1,2),
     5                    ((K21W(K,L),K=1,KL),L=1,2)
         WRITE (IU17,999) (AF11(M),M=1,ML+4),   (FKLAP(M),M=1,ML+4),
     1                    (FKLAP1(M),M=1,ML+4), (FKLAM(M),M=1,ML+4),
     2                    (FKLAM1(M),M=1,ML+4),
     3                    ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    5. WRITE COMMON COUPLE.
C        --------------------
C
 5000 CONTINUE
      IF (IFORM.NE.2) THEN
         WRITE (IU07) BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV
      ENDIF
      IF (IFORM.NE.1) THEN
         WRITE (IU17,999) BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    6. WRITE COMMON TABLE.
C        -------------------
C
 6000 CONTINUE
      IF (IFORM.NE.2) THEN
         WRITE (IU07) TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP
      ENDIF
      IF (IFORM.NE.1) THEN
         WRITE (IU17,999) TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    7. WRITE COMMON COUT.
C        ------------------
C
 7000 CONTINUE
      IF (IFORM.NE.2) THEN
         WRITE (IU07)  NGOUT, IGAR, IJAR
      ENDIF
      IF (IFORM.NE.1) THEN
         WRITE (IU17,998)  NGOUT, IGAR, IJAR
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    8. WRITE COMMON SHALLOW.
C        ---------------------
C
 8000 CONTINUE
      IF (IFORM.NE.2) THEN
         WRITE (IU07) DEPTH, DEPTHA, DEPTHD,
     1                ((TCGOND(L,M),L=1,NDEPTH),M=1,ML),
     2                ((TFAK(L,M),L=1,NDEPTH),M=1,ML),
     3                ((TSIHKD(L,M),L=1,NDEPTH),M=1,ML)
      ENDIF
      IF (IFORM.NE.1) THEN
         WRITE (IU17,999) DEPTH, DEPTHA, DEPTHD,
     1                    ((TCGOND(L,M),L=1,NDEPTH),M=1,ML),
     2                    ((TFAK(L,M),L=1,NDEPTH),M=1,ML),
     3                    ((TSIHKD(L,M),L=1,NDEPTH),M=1,ML)
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    9. WRITE COMMON CURRENT.
C        ---------------------
C
 9000 CONTINUE
      IF (IREFRA.EQ.2) THEN
         IF (IFORM.NE.2) WRITE (IU07) U, V
         IF (IFORM.NE.1) WRITE (IU17,999) U, V
      ENDIF
C
C ----------------------------------------------------------------------
C
C*   10. WRITE NAMELIST PARWAM.
C        ----------------------
C
10000 CONTINUE
 
      CALL OUTNAM
     . (NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER, MOUTP, MOUTT,
     .  NMAXC, NMAXF, NBINP, NIBL1, NIBLD, NBLD, NIBLC, NBLC ,
     .  ITAUMAX, JUMAX, IUSTAR, IALPHA, NDEPTH, IREFRA, iper)

      RETURN
      END
