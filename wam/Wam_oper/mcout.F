      SUBROUTINE MCOUT

C ----------------------------------------------------------------------
C
C**** *MCOUT* - ROUTINE TO COMPUTE OUTPUT INDICES (COMMON COUT).
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       *MCOUT* COMPUTES THE INDICES OF SPECTRA OUTPUT POINTS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MCOUT*
C
C     METHOD.
C     -------
C
C       THE LATITUDE AND LOGITUDE ARE CONVERTED TO INDICES.
C
C     EXTERNALS.
C     ----------
C
C       *FINDB*     - FIND BLOCK AND SEA POINT NUMBERS.
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
#include "parameter_topo.h"
C
#include "comcinp.h"
C
#include "comcout.h"
C
#include "comgrid.h"
C
#include "commap.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
C*    1. NO OUTPUT POINTS SPECIFIED.
C        ---------------------------
C
 1000 CONTINUE
      IF (NGOUT.EQ.0) THEN
         WRITE(IU06,'(1H1,'' SPECIAL OUTPUT POINTS FOR SPECTRA:'')')
         WRITE(IU06,*) 'OUTPUT POINTS ARE NOT DEFINED IN USER INPUT'
         RETURN
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
C        -----------------------------------------
C
 2000 CONTINUE
      CALL FINDB (NIBLO, NGOUT, OUTLAT, OUTLONG, IGAR, IJAR)
C
C ----------------------------------------------------------------------
C
C*    3. PRINTER PROTOCOL.
C        -----------------
C
 3000 CONTINUE
      WRITE(IU06,'(1H1,'' SPECIAL OUTPUT POINTS FOR SPECTRA:'')')
      WRITE(IU06,'(''    NUMBER OF OUTPUT POINTS IS NGOUT = '',I4)')
     1      NGOUT
      WRITE(IU06,'(4X,''     |-----INPUT-----|-NEAREST POINT-|'',
     1              ''-POINT INDEX--|'')')
      WRITE(IU06,'(4X,''  NO.    LAT.   LONG.    LAT.   LONG.  BLOCK.'',
     1             ''  POINT.'')')
      DO 3001 IO=1,NGOUT
         IF (IJAR(IO).GT.0) THEN
            KX  = KXLT(IJAR(IO),IGAR(IO))
            IX  = IXLG(IJAR(IO),IGAR(IO))
            ALONG = AMOWEP + (IX-1)*ZDELLO(KX)
            ALAT  = AMOSOP + (KX-1)*XDELLA
         ELSE
            ALONG = 9999999
            ALAT  = 9999999
         ENDIF
         WRITE(IU06,'(4X,I5,4F8.2,2I8)')
     1     IO, OUTLAT(IO), OUTLONG(IO), ALAT, ALONG, IGAR(IO), IJAR(IO)
 3001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    4. REMOVE OUTPUT POINTS WHICH ARE NOT IN GRID.
C        -------------------------------------------
C
 4000 CONTINUE
      NG = 0
      DO 4001 IO=1,NGOUT
         IF (NG.GT.0 .AND. IO-NG.GT.0) THEN
            IJAR(IO-NG) = IJAR(IO)
            IGAR(IO-NG) = IGAR(IO)
         ENDIF
         IF (IJAR(IO).EQ.0) NG = NG+1
 4001 CONTINUE
      NGOUT = NGOUT-NG
      IF (NG.GT.0) THEN
         WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
         WRITE (IU06,*) ' +                                       +'
         WRITE (IU06,*) ' +     WARNING ERROR FROM SUB. MCOUT     +'
         WRITE (IU06,*) ' +     =============================     +'
         WRITE (IU06,*) ' +                                       +'
         WRITE (IU06,*) ' + NO SEAPOINT FOUND FOR NG = ',NG
         WRITE (IU06,*) ' + OUTPUT POINT REQUESTS (SEE ABOVE LIST)+'
         WRITE (IU06,*) ' + THESE POINTS WILL NOT BE TAKEN.       +'
         WRITE (IU06,*) ' + NUMBER OF OUTPUT POINTS IS NGOUT = ', NGOUT
         WRITE (IU06,*) ' +                                       +'
         WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
      ENDIF

      RETURN
      END
