C------------------------------------------------------------------
C
      SUBROUTINE READT(IU06,FILNM,ITEST,TEMPG)
C
C----------------------------------------------------------------------
C
C**** *READT*  READ TEMPERATURE GRID.
C
C     R. PORTZ      MPI HAMBURG OCTOBER 1992
C
C*    PURPOSE.
C     --------
C
C       INPUT OF GLOBAL TEMPARATURE.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *READT* *(IU06,FILNM,ITEST,TEMPG)*
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
C      *FILNM*     CHARACTER*24 DATA INPUT FILENAME.
C      *ITEST*     INTEGER   OUTPUT FLAG IF LE 2 NO EXTRA OUTPUT IS
C                            GENERATED.
C      *TEMPG*     REAL      GLOBAL TEMPRERATURES (KELVIN).
C
C     METHOD.
C     -------
C
C       FORMATED READ FROM UNIT 1.
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
C-------------------------------------------------------------------
C
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.
C
#include "param.h"
#include "commap.h"
#include "comtemp.h"
#include "parcons.h"
C
      REAL,DIMENSION(NXT,NYT) :: TEMPG
      REAL,ALLOCATABLE :: WORK(:,:)
      INTEGER,ALLOCATABLE :: KLONRGG(:)
C
      CHARACTER*24 FILNM
      CHARACTER*40 MSG, TITL*100
C
C-----------------------------------------------------------------------
C
C*    1. INPUT OF TEMPERATURE DATA.
C     -----------------------------
C
 1000 CONTINUE
C
      ALLOCATE(KLONRGG(NYT))
      CALL INMARS  ( FILNM, NXT, NYT, ZMISS, CTDATE, IFORP, IPARAM,
     .               NGXT, NGYT, TMOWEP, TMOSOP, TMOEAP, TMONOP,
     .               KRGG,KLONRGG,KZLEV,TEMPG, IERR, MSG)
      IF ( IERR .EQ. 0 ) THEN
        WRITE(IU06,*)'            '
        WRITE(IU06,*)' SUB. READT: ICE BOUNDARY DATE=',CTDATE
      ELSE
        WRITE(IU06,*)'   SUB. READT: ICE BOUNDARY IERR=',IERR
        WRITE(IU06,*) MSG
      ENDIF
      CALL ADJUST (TMOWEP,TMOEAP)
      DELLA  = (TMONOP - TMOSOP)/REAL(NGYT-1)
      DELLO  = (TMOEAP - TMOWEP)/REAL(NGXT-1)
C
      IF (ITEST.GT.3) THEN
        ALLOCATE(WORK(NXT,NYT))
        WRITE(IU06,*) 'TEMPERATURE GRID:'
        WRITE(IU06,*) '-----------------'
        WRITE(IU06,*)
        WRITE(IU06,*) '  RESOLUTION LAT-LON  ',DELLA, DELLO
        WRITE(IU06,*) '   SOUTHERN LAT: ',TMOSOP,
     &                '   NORTHERN LAT: ',TMONOP
        WRITE(IU06,*) '   WESTERN LONG: ',TMOWEP,
     &                '   EASTERN LONG: ',TMOEAP
        WRITE(IU06,*)
        WORK=TEMPG-273.15
        TITL= '   TEMPERATURES IN CELCIUS         '
        CALL OUTPP (CTDATE, IU06, NXT, NYT, NGXT, NGYT,
     &              TITL, 1.,
     &              WORK, TMOWEP, TMOSOP, TMOEAP, TMONOP)
        DEALLOCATE(WORK)
      ENDIF
C
C*    1.1 ERROR MESSAGES
C     ------------------
C
 1100 CONTINUE
C
      IF (NGXT .GT. NXT) THEN
         WRITE(IU06,*)'********************************'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* FATAL ERROR IN SUB READT     *'
         WRITE(IU06,*)'* ========================     *'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* PARAMETER NXT IS TO SMALL    *'
         WRITE(IU06,*)'* SET NXT = ',NGXT
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'********************************'
         CALL ABORT1
      ENDIF
      IF (NGYT .GT. NYT) THEN
         WRITE(IU06,*)'********************************'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* FATAL ERROR IN SUB READT     *'
         WRITE(IU06,*)'* ========================     *'
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'* PARAMETER NXT IS TO SMALL    *'
         WRITE(IU06,*)'* SET NYT = ',NGYT
         WRITE(IU06,*)'*                              *'
         WRITE(IU06,*)'********************************'
         CALL ABORT1
      ENDIF
C
      IF(ALLOCATED(KLONRGG))DEALLOCATE(KLONRGG)

      RETURN
      END
