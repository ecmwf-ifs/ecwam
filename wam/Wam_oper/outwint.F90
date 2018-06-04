      SUBROUTINE OUTWINT

! ----------------------------------------------------------------------

!**** *OUTWINT* -

!*    PURPOSE.
!     --------


!**   INTERFACE.
!     ----------

!     METHOD.
!     -------

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU, LIFS_IO_SERV_ENABLED
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWCOUT  , ONLY : JPPFLAG ,NIPRMOUT , NINFOBOUT,              &
     &                      INFOBOUT,BOUT,LWAM_USE_IO_SERV
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEF   ,CDTPRO   ,              &
     &            CFDBSF   ,MARSTYPE ,NWFDBREF ,LFDBOPEN
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "outint_io_serv.intfb.h"
#include "abort1.intfb.h"
#include "difdate.intfb.h"
#include "incdate.intfb.h"
#include "outint.intfb.h"

      INTEGER(KIND=JWIM) :: IFCST, INHOUR, ISHIFT
      INTEGER(KIND=JWIM) :: IY,IM,ID,IH,IMN,ISS

      REAL(KIND=JWRB) :: ZHOOK_HANDLE

      CHARACTER(LEN=14) :: CDATE, CDATE1, CDATE2

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTWINT',0,ZHOOK_HANDLE)

!     DETERMINE THE FORECAST STEP INFORMATION
      IF(MARSTYPE.EQ.'4v') THEN
!*      THIS IS A 4V CASE
        CDATE=CDATEA
        READ(CDATEA,'(I4,5I2)')IY,IM,ID,IH,IMN,ISS
        CALL DIFDATE(CDATEA,CDTPRO,IFCST)
        IFCST=IFCST/3600
      ELSEIF(CDTPRO.LE.CDATEF) THEN
!*      THIS IS AN ANALYSIS DATE
        CDATE=CDTPRO
        READ(CDTPRO,'(I4,5I2)') IY,IM,ID,IH,IMN,ISS
        IFCST = 0
      ELSE
!*      THIS IS A  FORECAST DATE
        CDATE=CDATEF
        READ(CDATEF,'(I4,5I2)')IY,IM,ID,IH,IMN,ISS
!       FIND THE STEP IN HOURS ! 
        CDATE1=CDATEF
        CDATE2=CDATE1(1:8)//'010000'
!       MORE THAN A YEAR IN SECONDS 366*24*3600=31622400
        CALL INCDATE(CDATE2,31622400) 
        CDATE2=CDATE2(1:8)//'010000'
        INHOUR=0
!       FIND THE NUMBER OF FULL YEARS AND CONVERT THEM IN HOURS
        DO WHILE (CDATE2.LT.CDTPRO)
          CALL DIFDATE(CDATE1,CDATE2,ISHIFT)
          INHOUR=INHOUR+ISHIFT/3600
          CDATE1=CDATE2
          CDATE2=CDATE1(1:8)//'010000'
          CALL INCDATE(CDATE2,31622400) 
          CDATE2=CDATE2(1:8)//'010000'
        ENDDO
        CALL DIFDATE(CDATE1,CDTPRO,IFCST)
        IF(MOD(IFCST,3600) == 0 ) THEN
          IFCST=IFCST/3600
        ELSE
          WRITE(IU06,*) ' -----------------------------------------'
          WRITE(IU06,*) ' ERROR in routine OUTINT :'
          WRITE(IU06,*) ' forecast step must be  multiple of hours!'
          WRITE(IU06,*) ' IFCST =', IFCST
          WRITE(IU06,*) ' CDATEF=', CDATEF 
          WRITE(IU06,*) ' CDATE1=', CDATE1 
          WRITE(IU06,*) ' CDTPRO=', CDTPRO 
          WRITE(IU06,*) ' INHOUR=', INHOUR 
          WRITE(IU06,*) ' -----------------------------------------'
          CALL ABORT1
        ENDIF
        IFCST=INHOUR+IFCST
      ENDIF

      ! Use IFS IO server?
      IF (LWCOU .AND. LIFS_IO_SERV_ENABLED .AND. LWAM_USE_IO_SERV) THEN
          CALL OUTINT_IO_SERV(NIPRMOUT, IJSLOC, IJLLOC, BOUT, INFOBOUT, MARSTYPE, CDATE, IFCST)
      ELSE
          CALL OUTINT(CDATE, IFCST)
      ENDIF

    IF (LHOOK) CALL DR_HOOK('OUTWINT',1,ZHOOK_HANDLE)

END SUBROUTINE OUTWINT
