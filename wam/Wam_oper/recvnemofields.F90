SUBROUTINE RECVNEMOFIELDS(BLK2LOC, WVENVI, NEMO2WAM, FF_NOW, LREST, LINIT) 

!****  *RECVNEMOFIELDS* - UPDATE FIELDS WAVE FIELDS WITH NEMO INFORMATION

!      KRISTIAN MOGENSEN ECMWF    MARCH 2013

!      MODIFICATION.
!      -------------
!                                            

!     PURPOSE.                                                          
!     --------                                                          

!          THIS SUBROUTINE PASSES NEMO INFORMATION THROUGH TO
!          WAM VIA THE NEMO SINGLE EXECUTABLE COUPLING INTERFACE

!*    INTERFACE.                                                        
!     ----------                                                        


!     METHOD.                                                           
!     -------                                                           

!          PARALLEL INTERPOLATION BASED ON PREDETERMINED WEIGHTS

!     EXTERNALS.                                                        
!     ----------                                                        

!          NEMOGCMCOUP_WAM_GET  -  UPDATE NEMO FIELDS IN WAM

!     REFERENCES.                                                       
!     -----------                                                       

!          NONE                                                         

! -------------------------------------------------------------------   

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

! GRID POINTS CHUNKS
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, NTOTIJ, KIJL4CHNK
! MODULES NEEDED FOR LAKE MASK HANDLING
      USE YOWWIND  , ONLY : FIELDG   ,LLNEWCURR 
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_DATA_MODULE, ONLY : MPL_COMM
      USE MPL_MODULE
! COUPLING INFORMATION
      USE YOWCOUP  , ONLY : LWCOU, LWNEMOCOUCIC, LWNEMOCOUCIT, LWNEMOCOUCUR
! ICE AND CURRENT INFORMATION 
      USE YOWCURR  , ONLY : CURRENT_MAX
! OUTPUT FORTRAN UNIT
      USE YOWTEST  , ONLY : IU06
! NEMO FIELDS ON WAVE GRID
      USE YOWNEMOFLDS,ONLY: LNEMOCITHICK, LNEMOICEREST  
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

! -------------------------------------------------------------------   

      IMPLICIT NONE
#include "abort1.intfb.h"

      TYPE(WVGRIDLOC), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: WVENVI
      TYPE(OCEAN2WAVE), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: NEMO2WAM
      TYPE(FORCING_FIELDS), DIMENSION(NPROMA_WAM, NCHNK), INTENT(INOUT) :: FF_NOW ! FORCING FIELDS
      LOGICAL, INTENT(IN) :: LREST ! RESTART SO UPDATE FROM RESTART VALUES
      LOGICAL, INTENT(IN) :: LINIT ! UPDATE CICOVER, CITHICK, UCUR, VCUR AT INITIAL TIME
                                   ! IF NEEDED.

      INTEGER(KIND=JWIM), PARAMETER :: NFIELD = 5
      INTEGER(KIND=JWIM) :: IX, IY, IJ
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IFLD
      REAL(KIND=JWRB), DIMENSION(NTOTIJ, NFIELD) :: ZNEMOTOWAM
      REAL(KIND=JWRB) :: ZHOOK_HANDLE

! -------------------------------------------------------------------   

IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',0,ZHOOK_HANDLE)

ASSOCIATE(IFROMIJ => BLK2LOC%IFROMIJ, &
 &        JFROMIJ => BLK2LOC%JFROMIJ, &
 &        UCUR => WVENVI%UCUR, &
 &        VCUR => WVENVI%VCUR, &
 &        CICOVER => FF_NOW%CICOVER, &
 &        CITHICK => FF_NOW%CITHICK, &
 &        NEMOSST => NEMO2WAM%NEMOSST, &
 &        NEMOCICOVER => NEMO2WAM%NEMOCICOVER, &
 &        NEMOCITHICK => NEMO2WAM%NEMOCITHICK, &
 &        NEMOUCUR => NEMO2WAM%NEMOUCUR, &
 &        NEMOVCUR => NEMO2WAM%NEMOVCUR)


      ! IF WE ARE IN RESTART JUST COPY THE RESTART FILE INFO.
      IF (LREST) THEN
        NEMOCICOVER(IJS:IJL)=CICOVER(IJS:IJL)
        NEMOCITHICK(IJS:IJL)=CITHICK(IJS:IJL)
        NEMOUCUR(IJS:IJL)=UCUR(IJS:IJL)
        NEMOVCUR(IJS:IJL)=VCUR(IJS:IJL)
        LNEMOICEREST=.TRUE.
        LNEMOCITHICK=.TRUE.
      ELSE
#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_GET( IRANK-1, NPROC, MPL_COMM,     &
     &                            NTOTIJ, NFIELD, ZNEMOTOWAM,   &
     &                            LNEMOCITHICK )

        LLNEWCURR=.TRUE. 

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL, IFLD)
        DO ICHNK = 1, NCHNK
            KIJS = SUM(KIJL4CHNK(1:MAX(1,ICHNK-1))) 
            KIJL = KIJS + KIJL4CHNK(ICHNK) -1
!debile check
            IFLD = 1
            NEMOSST(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMOCICOVER(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMOCITHICK(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMOUCUR(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMOVCUR(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
        ENDDO
!$OMP   END PARALLEL DO

#endif
        LNEMOICEREST=.FALSE.

      ENDIF

!     UPDATE CICOVER, CITHICK UCUR AND VCUR AT INITIAL TIME ONLY !!!!
      IF (LINIT) THEN

        WRITE(IU06,*)' RECVNEMOFIELDS: INITIALISE OCEAN FIELDS'

        IF (LWCOU) THEN
          IF (LWNEMOCOUCIC) THEN
            DO IJ = IJS,IJL
              IX = IFROMIJ(IJ)
              IY = JFROMIJ(IJ)
!            if lake cover = 0, we assume open ocean point, then get sea ice directly from NEMO
              IF (FIELDG(IX,IY)%LKFR <= 0.0_JWRB ) THEN
                CICOVER(IJ)=NEMOCICOVER(IJ)
              ELSE
!            get ice information from atmopsheric model
                CICOVER(IJ)=FIELDG(IX,IY)%CICOVER 
              ENDIF
            ENDDO
          ENDIF

          IF (LWNEMOCOUCIT) THEN
            DO IJ = IJS,IJL
              IX = IFROMIJ(IJ)
              IY = JFROMIJ(IJ)
!            if lake cover = 0, we assume open ocean point, then get sea ice thickness directly from NEMO
              IF (FIELDG(IX,IY)%LKFR <= 0.0_JWRB ) THEN
                CITHICK(IJ)=NEMOCICOVER(IJ)*NEMOCITHICK(IJ)
              ELSE
                CICOVER(IJ)=0.5_JWRB*NEMOCICOVER(IJ)
              ENDIF
            ENDDO
          ENDIF

          IF (LWNEMOCOUCUR) THEN
              DO IJ = IJS,IJL
                IX = IFROMIJ(IJ)
                IY = JFROMIJ(IJ)
!              if lake cover = 0, we assume open ocean point, then get currents directly from NEMO
                IF (FIELDG(IX,IY)%LKFR <= 0.0_JWRB ) THEN
                  UCUR(IJ) = SIGN(MIN(ABS(NEMOUCUR(IJ)),CURRENT_MAX),NEMOUCUR(IJ))
                  VCUR(IJ) = SIGN(MIN(ABS(NEMOVCUR(IJ)),CURRENT_MAX),NEMOVCUR(IJ))
                ELSE
                  UCUR(IJ)=0.0_JWRB
                  VCUR(IJ)=0.0_JWRB
                ENDIF
              ENDDO
          ENDIF

        ELSE

          IF (LWNEMOCOUCIC) THEN
            CICOVER(IJS:IJL)=NEMOCICOVER(IJS:IJL)
          ENDIF
          IF (LWNEMOCOUCIT) THEN
            CITHICK(IJS:IJL)=NEMOCICOVER(IJS:IJL)*NEMOCITHICK(IJS:IJL)
          ENDIF
          IF (LWNEMOCOUCUR) THEN
             UCUR(IJS:IJL)=NEMOUCUR(IJS:IJL)
             VCUR(IJS:IJL)=NEMOVCUR(IJS:IJL)
          ENDIF
        ENDIF

      ENDIF

      IF (LWNEMOCOUCIT.AND.(.NOT.LNEMOCITHICK)) THEN
        WRITE(IU06,*) ' --------------------------------'
        WRITE(IU06,*) ' LWNEMOCOUCIT ONLY MAKES SENSES  '
        WRITE(IU06,*) ' IF LIM IS ACTIVATED IN NEMO     '
        WRITE(IU06,*) ' --------------------------------'
        CALL ABORT1
      ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',1,ZHOOK_HANDLE)

END SUBROUTINE RECVNEMOFIELDS
