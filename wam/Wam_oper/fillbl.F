      SUBROUTINE FILLBL (FL, IG, IU13)

C ----------------------------------------------------------------------
C
C**** *FILLBL* - ADDS FIRST AND LAST LATITUDE TO BLOCK
C****            AND SAVES SECOND TO LAST LAT FOR NEXT BLOCK.
C
C     HEINZ GUNTHER         ECMWF         SEPTEMBER 1989
C
C*    PURPOSE.
C     --------
C
C       THIS ROUTINE ADDS TO A BLOCK OF SPECTRA THE SEA POINTS IN
C       THE FIRST AND LAST LATITUDE. THE INFORMATION IS READ FROM
C       THE INTERMEDIATE STORAGE (SEE SUB SPLITBL) FOR THE LAST
C       LATITUDE AND TAKEN FROM COMMON OVER FOR THE FIRST LATITUDE.
C       THE COMMON IS FILLED WITH THE SECOND TO LAST LATITUDE TO
C       PASS THE INFORMATION TO THE NEXT BLOCK.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *FILLBL (FL, IG, IU13)*
C          *FL*      REAL     BLOCK OF SPECTRA.
C          *IG*      INTEGER  BLOCK NUMBER.
C          *IU13*    INTEGER  INPUT UNIT OF LAST LATITUDE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     METHOD.
C     -------
C
C          UNFORMATED READ FROM SPECIFIED UNITS OR COPY.
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
#include "comgrid.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
#include "comover.h"
C
#include "txtover.h"
C
C ----------------------------------------------------------------------
C
      DIMENSION FL(NINF-1:NSUP,NANG,NFRE)
C
C ----------------------------------------------------------------------
C
C*    1. ADD LATITUDE BEFORE LAST OF PREVIOUS BLOCK INTO FIRST
C        LATITUDE OF BLOCK.
C        ------------------------------------------------------
C
      IF(IG.NE.1) THEN
         IJA = 1
         IJE = IJS(IG)-1
         DO 1201 M=1,NFRE
            DO 1202 K=1,NANG
               DO 1203 IJ=IJA,IJE
                  FL(IJ,K,M) = OV(IJ,K,M)
 1203          CONTINUE
 1202       CONTINUE
 1201    CONTINUE
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. ADD SECOND LATITUDE OF NEXT BLOCK INTO LAST LATITUDE OF BLOCK.
C        --------------------------------------------------------------
C
      IF(IG.NE.IGL) THEN
         IJA = IJL(IG)+1
         IJE = IJLT(IG)
         READ(IU13) (((FL(IJ,K,M),IJ=IJA,IJE),K=1,NANG),M=1,NFRE)
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3. LATITUDE BEFORE LAST OF BLOCK IS SAVED FOR FIRST
C        LATITUDE OF NEXT BLOCK.
C        -------------------------------------------------
C
      IF (IG.NE.IGL) THEN
         IJA = IJLS(IG)
         IJH = 1-IJA
         IJE = IJL(IG)
         DO 3201 M=1,NFRE
            DO 3202 K=1,NANG
               DO 3203 IJ=IJA,IJE
                  OV(IJ+IJH,K,M) = FL(IJ,K,M)
 3203          CONTINUE
 3202       CONTINUE
 3201    CONTINUE
      ENDIF


      RETURN
      END
