CICE
C-----------------------------------------------------------------------
C
      SUBROUTINE SETICE (FL3,IJS,IJL,IG)
C
C-----------------------------------------------------------------------
C
C**** *SETICE* ROUTINE TO SET SPECTRA ON ICE TO ZERO.
C
C     R.PORTZ      MPI         OKT.1992
C     J. BIDLOT    ECMWF       JUNE 1996  MESSAGE PASSING
C     P. JANSSEN   ECMWF       OCTOBER 1997 SIMPLIFY BY INTRODUCING
C                                           ICEMASK
C
C     PURPOSE.
C     -------
C
C          *READICE* SET ICE SPECTRA (FL3) TO ZERO
C
C**   INTERFACE.
C     ----------
C
C         *CALL* *SETICE*
C
C     METHOD.
C     -------
C
C          NONE.
C
C
C     EXTERNALS.
C     ----------
C
C          .NONE.
C
C     REFERENCE.
C     ----------
C
C          NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "comice.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "txtmpp.h"
C
      DIMENSION FL3(NINF-1:NSUP,NANG,NFRE)
C
C ----------------------------------------------------------------------
C
C*    1. SET SPECTRA TO ZERO AT ICE POINTS.
C     ------------------------------------
C
      DO M = 1, NFRE
         DO K = 1, NANG
            DO IJ = IJS,IJL
               ZICE = FLOAT(ICEMASK(IJ,IG))
               FL3(IJ,K,M) = ZICE*FL3(IJ,K,M)
            ENDDO
         ENDDO
      ENDDO
C      
      RETURN
C      
      END
CICE
