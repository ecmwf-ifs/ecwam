! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE RECVNEMOFIELDS(BLK2LOC, WVENVI, NEMO2WAM,  &
 &                        NXS, NXE, NYS, NYE, FIELDG, &
 &                        FF_NOW, LREST, LINIT) 

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

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

! GRID POINTS CHUNKS
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, NTOTIJ, KIJL4CHNK, IJFROMCHNK
      USE YOWCOUT  , ONLY : NREAL
! MODULES NEEDED FOR LAKE MASK HANDLING
      USE YOWWIND  , ONLY : LLNEWCURR 
! MPP INFORMATION
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE MPL_MODULE, ONLY : MPL_COMM
! COUPLING INFORMATION
      USE YOWCOUP  , ONLY : LWCOU, LWNEMOCOUCIC, LWNEMOCOUCIT, LWNEMOCOUCUR, LWNEMOCOUDEBUG
! ICE AND CURRENT INFORMATION 
      USE YOWCURR  , ONLY : CURRENT_MAX
! OUTPUT FORTRAN UNIT
      USE YOWTEST  , ONLY : IU06
! NEMO FIELDS ON WAVE GRID
      USE YOWNEMOFLDS,ONLY: LNEMOCITHICK, LNEMOICEREST, &
     &                      LLAWPEND, IJINFLAW, IJSUPLAW, CFILELAW, CDLAWHDR, RFIELDLAW
! DR. HOOK
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! -------------------------------------------------------------------   

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "outmdldcp.intfb.h"
#include "iwam_get_unit.intfb.h"

      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(OCEAN2WAVE), INTENT(INOUT) :: NEMO2WAM
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG

      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW ! FORCING FIELDS
      LOGICAL, INTENT(IN) :: LREST ! RESTART SO UPDATE FROM RESTART VALUES
      LOGICAL, INTENT(IN) :: LINIT ! UPDATE CICOVER, CITHICK, UCUR, VCUR AT INITIAL TIME
                                   ! IF NEEDED.


      INTEGER(KIND=JWIM), PARAMETER :: NFIELD = 5
      INTEGER(KIND=JWIM) :: IX, JY, IJ
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IC, IFLD
      INTEGER(KIND=JWIM) :: IJSB, IJLB, IUNIT
      REAL(KIND=JWRO), DIMENSION(NTOTIJ, NFIELD) :: ZNEMOTOWAM
       REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      LOGICAL :: LLFLDUPDT
! -------------------------------------------------------------------   

IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',0,ZHOOK_HANDLE)


      ! IN RESTART, KEEP THE NEMO FIELDS GETSTRESS RESTORED FROM THE LAW FILE:
      ! THEY ARE THE ONES THE RUN THAT WROTE IT LAST RECEIVED, AND THIS FIRST
      ! UPDATE REPEATS THAT RUN'S LAST FORCING UPDATE. DO NOT ASK NEMO: ITS OWN
      ! RESTART IS ONE EXCHANGE AHEAD, SO ITS FIELDS ARE THOSE OF THE NEXT
      ! UPDATE, AND IT STILL FLAGS THEM AS NEW FOR THAT ONE. (COPYING THE ICE
      ! FROM FF_NOW, AS BEFORE, WAS ONE EXCHANGE BEHIND INSTEAD.)
      IF (LREST) THEN
        ! NEMOCITHICK IS NEMO'S RAW THICKNESS AGAIN, NOT SCALED BY ICE COVER.
        LNEMOICEREST=.FALSE.
        LNEMOCITHICK=.TRUE.
      ELSE
#ifdef WITH_NEMO
        CALL NEMOGCMCOUP_WAM_GET( IRANK-1, NPROC, MPL_COMM,     &
     &                            NTOTIJ, NFIELD, ZNEMOTOWAM,   &
     &                            LNEMOCITHICK, LLFLDUPDT, LWNEMOCOUDEBUG )

        LLNEWCURR=.TRUE.

        IF (LLFLDUPDT) THEN
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL, IC, IFLD)
          DO ICHNK = 1, NCHNK

            KIJS = 1
            DO IC = 1, ICHNK-1
               KIJS = KIJS + KIJL4CHNK(IC)
            ENDDO
            KIJL = KIJS + KIJL4CHNK(ICHNK) - 1

            IFLD = 1
            NEMO2WAM%NEMOSST(1:KIJL4CHNK(ICHNK),ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMO2WAM%NEMOCICOVER(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMO2WAM%NEMOCITHICK(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMO2WAM%NEMOUCUR(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)
            IFLD = IFLD + 1
            NEMO2WAM%NEMOVCUR(1:KIJL4CHNK(ICHNK), ICHNK) = ZNEMOTOWAM(KIJS:KIJL, IFLD)

            IF ( KIJL4CHNK(ICHNK) < NPROMA_WAM ) THEN
!             values for fictious points
              NEMO2WAM%NEMOSST(KIJL4CHNK(ICHNK)+1:NPROMA_WAM,ICHNK)      = NEMO2WAM%NEMOSST(1, ICHNK)
              NEMO2WAM%NEMOCICOVER(KIJL4CHNK(ICHNK)+1:NPROMA_WAM, ICHNK) = NEMO2WAM%NEMOCICOVER(1, ICHNK)
              NEMO2WAM%NEMOCITHICK(KIJL4CHNK(ICHNK)+1:NPROMA_WAM, ICHNK) = NEMO2WAM%NEMOCITHICK(1, ICHNK)
              NEMO2WAM%NEMOUCUR(KIJL4CHNK(ICHNK)+1:NPROMA_WAM,ICHNK)     = NEMO2WAM%NEMOUCUR(1, ICHNK)
              NEMO2WAM%NEMOVCUR(KIJL4CHNK(ICHNK)+1:NPROMA_WAM,ICHNK)     = NEMO2WAM%NEMOVCUR(1, ICHNK)
            ENDIF
          ENDDO
!$OMP     END PARALLEL DO

          IF (LWNEMOCOUDEBUG) CALL OUTMDLDCP(NEMO2WAM=NEMO2WAM,NTYPE=1)

       ENDIF


#endif
        LNEMOICEREST=.FALSE.

      ENDIF

!     A LAW RESTART FILE WRITTEN SINCE THE LAST RECEIVE IS STILL WAITING FOR
!     THE NEMO FIELDS OF THIS ONE (SEE SAVSTRESS): REWRITE IT WITH THEM.
      IF (LLAWPEND .AND. .NOT.LREST) THEN
        DO ICHNK = 1, NCHNK
          KIJL = KIJL4CHNK(ICHNK)
          IJSB = IJFROMCHNK(1, ICHNK)
          IJLB = IJFROMCHNK(KIJL, ICHNK)
          RFIELDLAW(IJSB:IJLB,17) = NEMO2WAM%NEMOSST(1:KIJL,ICHNK)
          RFIELDLAW(IJSB:IJLB,18) = NEMO2WAM%NEMOCICOVER(1:KIJL,ICHNK)
          RFIELDLAW(IJSB:IJLB,19) = NEMO2WAM%NEMOCITHICK(1:KIJL,ICHNK)
          RFIELDLAW(IJSB:IJLB,20) = NEMO2WAM%NEMOUCUR(1:KIJL,ICHNK)
          RFIELDLAW(IJSB:IJLB,21) = NEMO2WAM%NEMOVCUR(1:KIJL,ICHNK)
        ENDDO
        IUNIT = IWAM_GET_UNIT(IU06, CFILELAW(1:LEN_TRIM(CFILELAW)), 'w', 'u', 0, 'READWRITE')
        WRITE(IUNIT) CDLAWHDR
        DO IFLD = 1, NREAL
          WRITE(IUNIT) (RFIELDLAW(IJ,IFLD), IJ=IJINFLAW, IJSUPLAW)
        ENDDO
        CLOSE(IUNIT)
        DEALLOCATE(RFIELDLAW)
        LLAWPEND = .FALSE.
        WRITE(IU06,*) ' RECVNEMOFIELDS: LAW RESTART REWRITTEN WITH NEMO FIELDS: ', TRIM(CFILELAW)
      ENDIF

!     UPDATE CICOVER, CITHICK UCUR AND VCUR AT INITIAL TIME ONLY !!!!
      IF (LINIT) THEN

        WRITE(IU06,*)' RECVNEMOFIELDS: INITIALISE OCEAN FIELDS'

        IF (LWCOU) THEN
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, IJ, IX, JY)
          DO ICHNK = 1, NCHNK

           IF (LWNEMOCOUCIC) THEN
              DO IJ = 1, NPROMA_WAM
                IX = BLK2LOC%IFROMIJ(IJ,ICHNK)
                JY = BLK2LOC%JFROMIJ(IJ,ICHNK)
!              if lake cover = 0, we assume open ocean point, then get sea ice directly from NEMO
                IF (FIELDG%LKFR(IX,JY) <= 0.0_JWRB ) THEN
                  FF_NOW%CICOVER(IJ,ICHNK)=NEMO2WAM%NEMOCICOVER(IJ,ICHNK)
                ELSE
!              get ice information from atmopsheric model
                  FF_NOW%CICOVER(IJ,ICHNK)=FIELDG%CICOVER(IX,JY)
                ENDIF
              ENDDO
            ENDIF

            IF (LWNEMOCOUCIT) THEN
              DO IJ = 1, NPROMA_WAM
                IX = BLK2LOC%IFROMIJ(IJ,ICHNK)
                JY = BLK2LOC%JFROMIJ(IJ,ICHNK)
!              if lake cover = 0, we assume open ocean point, then get sea ice thickness directly from NEMO
                IF (FIELDG%LKFR(IX,JY) <= 0.0_JWRB ) THEN
                  FF_NOW%CITHICK(IJ,ICHNK)=NEMO2WAM%NEMOCICOVER(IJ,ICHNK)*NEMO2WAM%NEMOCITHICK(IJ,ICHNK)
                ELSE
                  FF_NOW%CICOVER(IJ,ICHNK)=0.5_JWRB*NEMO2WAM%NEMOCICOVER(IJ,ICHNK)
                ENDIF
              ENDDO
            ENDIF

            IF (LWNEMOCOUCUR) THEN
              DO IJ = 1, NPROMA_WAM
                IX = BLK2LOC%IFROMIJ(IJ,ICHNK)
                JY = BLK2LOC%JFROMIJ(IJ,ICHNK)
!              if lake cover = 0, we assume open ocean point, then get currents directly from NEMO
                IF (FIELDG%LKFR(IX,JY) <= 0.0_JWRB ) THEN
                  WVENVI%UCUR(IJ,ICHNK) = SIGN(MIN(ABS(NEMO2WAM%NEMOUCUR(IJ,ICHNK)),REAL(CURRENT_MAX,JWRO)), &
 &                                             NEMO2WAM%NEMOUCUR(IJ,ICHNK))
                  WVENVI%VCUR(IJ,ICHNK) = SIGN(MIN(ABS(NEMO2WAM%NEMOVCUR(IJ,ICHNK)),REAL(CURRENT_MAX,JWRO)), &
 &                                             NEMO2WAM%NEMOVCUR(IJ,ICHNK))
                ELSE
                  WVENVI%UCUR(IJ,ICHNK)=0.0_JWRB
                  WVENVI%VCUR(IJ,ICHNK)=0.0_JWRB
                ENDIF
              ENDDO
            ENDIF

          ENDDO
!$OMP   END PARALLEL DO

        ELSE

!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK)
          DO ICHNK = 1, NCHNK
            IF (LWNEMOCOUCIC) FF_NOW%CICOVER(:,ICHNK)=NEMO2WAM%NEMOCICOVER(:,ICHNK)
            IF (LWNEMOCOUCIT) FF_NOW%CITHICK(:,ICHNK)=NEMO2WAM%NEMOCICOVER(:,ICHNK)*NEMO2WAM%NEMOCITHICK(:,ICHNK)
            IF (LWNEMOCOUCUR) THEN
             WVENVI%UCUR(:,ICHNK)=NEMO2WAM%NEMOUCUR(:,ICHNK)
             WVENVI%VCUR(:,ICHNK)=NEMO2WAM%NEMOVCUR(:,ICHNK)
            ENDIF
          ENDDO
!$OMP     END PARALLEL DO

        ENDIF

      ENDIF

      IF (LWNEMOCOUCIT.AND.(.NOT.LNEMOCITHICK)) THEN
        WRITE(IU06,*) ' --------------------------------'
        WRITE(IU06,*) ' LWNEMOCOUCIT ONLY MAKES SENSES  '
        WRITE(IU06,*) ' IF LIM IS ACTIVATED IN NEMO     '
        WRITE(IU06,*) ' --------------------------------'
        CALL ABORT1
      ENDIF

IF (LHOOK) CALL DR_HOOK('RECVNEMOFIELDS',1,ZHOOK_HANDLE)

END SUBROUTINE RECVNEMOFIELDS
