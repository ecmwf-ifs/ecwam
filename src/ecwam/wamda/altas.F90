SUBROUTINE ALTAS(IJSG, IJLG, MIJ, BLK2GLO,           &
 &               WVENVI, FF_NOW, INTFLDS, WAM2NEMO,  &
 &               WVPRPT, FL1, XLLWS)

! ----------------------------------------------------------------------

!****  *ALTAS* - WAVE DATA ASSIMILATION USING ALTIMETER DATA ALONE. 

!      J. BIDLOT       ECMWF     MAY     1999  SPLIT FROM WAMASSI 
!      S. ABDALLA      ECMWF     NOV     2011  ADD ODB;  CONVERT TO F90

!     PURPOSE.                                                          
!     --------                                                          

!         TO ANALYSE THE WAVE FIELD SUBSTITUTING THE FIRST GUESS        
!         SPECTRA WITH ANALYSED SPECTRA AND FIRST GUESS VALUES          
!         WITH ANALYSED VALUES IN THE GLOBAL GRID                       

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL*ALTAS (IJSG, IJLG, MIJ, BLK2GLO,
!                  WVENVI, FF_NOW, INTFLDS, WAM2NEMO, 
!                  WVPRPT, FL1, XLLWS)
!     *IJSG,IJLG : BLOCK INDEXES (ALTAS has not yet been modified to use the NPROMA CHUNKS)
!     *MIJ*    - LAST FREQUENCY INDEX OF THE PROGNOSTIC RANGE.
!     *BLK2GLO*- BLOCK TO GRID TRANSFORMATION
!     *WVENVI* - WAVE ENVIRONMENT FIELDS
!     *FF_NOW* - FORCING FIELDS AT CURRENT TIME.
!     *INTFLDS*- INTEGRATED/DERIVED PARAMETERS
!     *WAM2NEMO- WAVE FIELDS PASSED TO NEMO
!     *WVPRPT* - WAVE PROPERTIES
!     *FL1*    - INPUT SPECTRUM.
!     *XLLWS*  - WINDSEA MASK FROM INPUT SOURCE TERM

!     METHOD.                                                           
!     -------                                                           

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS, &
     &                         INTGT_PARAM_FIELDS, WAVE2OCEAN

      USE YOWABORT , ONLY : WAM_ABORT
      USE YOWALTAS , ONLY : INTLMAX  ,KMINLMAX ,KMAXLMAX  ,NOBSPE,     &
     &                      IJALT    ,DIFFALTFG,LODBRALT  ,ALTDATA,    &
     &                      ALTEXDATA,ALTUNDATA,NIJALT    ,CDATEOBS,   &
     &                      LECWAMDA
      USE YOWGRID  , ONLY : COSPH    ,NPROMA_WAM, NCHNK, IJFROMCHNK, KIJL4CHNK,   &
     &                      ICHNKFROMIJ, IPRMFROMIJ
      USE YOWMAP   , ONLY : XDELLA   ,AMOSOP   ,NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE   ,LLUNSTR
      USE YOWPCONS , ONLY : DEG      ,R
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06
#ifdef WAM_HAVE_UNWAM
      USE YOWPD    , ONLY : RANK
#endif

      USE MPL_MODULE
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "grfield.intfb.h"
#include "mpmdl4alt.intfb.h"
#include "mppewithindist.intfb.h"
#include "oifield.intfb.h"
#include "semean.intfb.h"
#include "updatewfld.intfb.h"

#ifdef WITH_ODB
#include "wam2odb.intfb.h"
#endif

      INTEGER(KIND=JWIM), INTENT(IN)                                           :: IJSG, IJLG
      INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK), INTENT(OUT)            :: MIJ
      TYPE(WVGRIDGLO), INTENT(IN)                                              :: BLK2GLO
      TYPE(ENVIRONMENT), INTENT(IN)                                            :: WVENVI
      TYPE(FORCING_FIELDS), INTENT(INOUT)                                      :: FF_NOW
      TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT)                                  :: INTFLDS
      TYPE(WAVE2OCEAN), INTENT(INOUT)                                          :: WAM2NEMO
      TYPE(FREQUENCY), INTENT(IN)                                              :: WVPRPT
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(OUT)   :: XLLWS


      INTEGER(KIND=JWIM) :: IJ, KJ, IR, IC, IP
      INTEGER(KIND=JWIM) :: MINIJS, MAXIJL
      INTEGER(KIND=JWIM) :: IJSB, IJLB, KIJS, KIJL, ICHNK, IPRM

      REAL(KIND=JWRB) :: XLMIN, XLMAX, SIGH, SIGMOD, REDUC
      REAL(KIND=JWRB) :: DISTMAX, EXTENDMAX, DMAX
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK) :: WORK, HSOIBCH, HSANCH, U10FGCH 
      REAL(KIND=JWRB), DIMENSION(IJSG:IJLG) :: DIST, HSMOD, CICVR
      REAL(KIND=JWRB), DIMENSION(IJSG:IJLG) :: HSOIB, HSAN, U10FG, U10AN
      REAL(KIND=JWRB), ALLOCATABLE,DIMENSION(:) :: HSMOD_H, CICVR_H

      LOGICAL :: LLEPSMIN

!     -------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ALTAS',0,ZHOOK_HANDLE)


      WRITE(IU06,*)' '
      WRITE(IU06,*)'  -----ALT DATA ASSIMILATION BEGINS--------'

      LLEPSMIN=.TRUE.

!*    1. MEASUREMENTS AND MODEL ARE MERGED BY OPTIMUM INTERPOLATION.
!        -----------------------------------------------------------

!     specify the correlation distance (DIST) and the maximum spreading
!     factor  EXTENDMAX and get first guess wave height (HSMOD_H)

!!! if DIST were to be reset during the run, then MPPEWITHINDIST has to be called
!!! again !!!!

!!! for now, we will disactivate the geographically varying DIST
!!! by setting XLMAX to 0 and XLMIN to maximum required value.

      IF ( LECWAMDA ) THEN
        IF (XDELLA < 0.36_JWRB) THEN
          XLMIN = 150000._JWRB
          XLMAX = 0.0_JWRB
        ELSE
          XLMIN = 300000._JWRB
          XLMAX = 0.0_JWRB
        ENDIF
      ELSE
!       METEOFRANCE CHOICE
        IF(XDELLA < 0.05_JWRB) THEN
          XLMIN = 50000._JWRB
          XLMAX = 0._JWRB
        ELSE IF(XDELLA < 0.3_JWRB) THEN
          XLMIN = 170000._JWRB
          XLMAX = 0._JWRB
        ELSE
          XLMIN = 300000._JWRB
          XLMAX = 0._JWRB
        ENDIF
      ENDIF

      SIGH=0.75_JWRB


!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, IPRM, KIJS, IJSB, KIJL, IJLB, IJ, KJ)
      DO ICHNK = 1, NCHNK

        CALL SEMEAN (FL1(:,:,:,ICHNK), 1, NPROMA_WAM, WORK(:,ICHNK), LLEPSMIN)

        DO IPRM = 1, KIJL4CHNK(ICHNK)
          IF (WORK(IPRM, ICHNK) > 0.0_JWRB) THEN
            WORK(IPRM, ICHNK) = 4.0_JWRB * SQRT(WORK(IPRM, ICHNK))
          ELSE
            WORK(IPRM, ICHNK) = 0.0_JWRB
          ENDIF
        ENDDO

        KIJS = 1
        IJSB = IJFROMCHNK(KIJS, ICHNK)
        KIJL = KIJL4CHNK(ICHNK)
        IJLB = IJFROMCHNK(KIJL, ICHNK)

        HSMOD(IJSB:IJLB) = WORK(KIJS:KIJL, ICHNK)
        CICVR(IJSB:IJLB) = FF_NOW%CICOVER(KIJS:KIJL, ICHNK)

        DO IJ = IJSB, IJLB
          KJ = BLK2GLO%KXLT(IJ)
          DIST(IJ) = (XLMIN + XLMAX*COSPH(KJ))/R
        ENDDO

      ENDDO
!$OMP END PARALLEL DO


      EXTENDMAX = SQRT(2.0_JWRB) * 2.0_JWRB

      DISTMAX = MAXVAL(DIST(:))

!!    DMAX IS THE ABSOLUTE MAXIMUM OF DISTMAX
      DMAX = EXTENDMAX*(XLMIN+XLMAX)/R

!     DETERMINE THE TABLE WHICH TELLS WHICH PE POTENTIALLY SHARES 
!     ALTIMETER DATA WITH OTHER PE'S.
!     IT IS BASED ON THE VALUES OF DIST BEFORE REDUCTION FOR LOW WAVE HEIGHT
      IF (.NOT. ALLOCATED(INTLMAX) ) THEN
        ALLOCATE(INTLMAX(NPROC))
        ALLOCATE(KMINLMAX(NPROC))
        ALLOCATE(KMAXLMAX(NPROC))

        CALL MPPEWITHINDIST(BLK2GLO, DMAX, INTLMAX, KMINLMAX, KMAXLMAX)
      ENDIF

!     SPECIFY THE MODEL ERROR:
!     -----------------------

!     note that the observation error is now specified in rfl4wam
      SIGMOD = 0.5_JWRB

!     REDUCE DIST FOR LOW WAVE HEIGHT:
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB, IJ, REDUC)
      DO ICHNK = 1, NCHNK
        KIJS = 1
        IJSB = IJFROMCHNK(KIJS,ICHNK)
        KIJL = KIJL4CHNK(ICHNK)
        IJLB = IJFROMCHNK(KIJL,ICHNK)

        DO IJ = IJSB, IJLB
          REDUC = 1.0_JWRB - EXP(-MIN((HSMOD(IJ)/SIGH)**2, 50.0_JWRB)) 
          DIST(IJ) = MAX(0.1_JWRB, REDUC) * DIST(IJ)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO


!     DISTRIBUTE MODEL VALUES TO WHERE THEY ARE NEEDED
!     ------------------------------------------------

      IF (NPROC > 1) THEN
        IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
          MINIJS=NIBLO
          MAXIJL=1
          DO IR = 1,NPROC
            IF (INTLMAX(IR) == 1) THEN
              DO IP=1,RANK(IR)%NP
                MINIJS=MIN(RANK(IR)%IPLG(IP), MINIJS)
                MAXIJL=MAX(RANK(IR)%IPLG(IP), MAXIJL)
              ENDDO
            ENDIF
          ENDDO
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
        ELSE
          MINIJS=IJSG
          MAXIJL=IJLG
          DO IR = 1, NPROC
            IF (INTLMAX(IR) == 1) THEN
              MINIJS=MIN(NSTART(IR), MINIJS)
              MAXIJL=MAX(NEND(IR), MAXIJL)
            ENDIF
          ENDDO
        ENDIF
      ELSE
        MINIJS=IJSG
        MAXIJL=IJLG
      ENDIF

!!!   HSMOD_H and CICVR_H have global indexing !!!
!!! it is not a very smart way to do the localisation !!!!
!!!! will need to be revised !!!
      ALLOCATE(HSMOD_H(MINIJS:MAXIJL))
      ALLOCATE(CICVR_H(MINIJS:MAXIJL))

      CALL MPMDL4ALT(IJSG, IJLG, HSMOD, CICVR,                      &
     &               MINIJS , MAXIJL , HSMOD_H, CICVR_H)


!*    2.2 THINNED MEASUREMENTS ARE READ IN
!         --------------------------------

      ALLOCATE (NOBSPE(NPROC))

      CALL GRFIELD(BLK2GLO, FL1, INTFLDS, MINIJS, MAXIJL, HSMOD_H,  &
     &             CICVR_H, SIGMOD, DMAX)

!*    2.3 OI.
!         ---
      CALL OIFIELD (IJSG, IJLG, MINIJS, MAXIJL,          &
     &              BLK2GLO,                             &
     &              SIGMOD, DISTMAX, EXTENDMAX,          &
     &              DIST,                                &
     &              HSMOD_H, CICVR_H,                    &
     &              HSOIB)

      DEALLOCATE(HSMOD_H)
      DEALLOCATE(CICVR_H)

!*    3. ANALYSING THE SPECTRA. 
!        ----------------------

!*    3.1 LOOP OVER THE BLOCKS.
!         ---------------------

!*    3.1.1 TRANSFERS RESULTS OF OPTIMAL INT TO BLOCKS.
!           -------------------------------------------

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, KIJS, IJSB, KIJL, IJLB)
      DO ICHNK = 1, NCHNK

        KIJS = 1
        IJSB = IJFROMCHNK(KIJS,ICHNK)
        KIJL = KIJL4CHNK(ICHNK)
        IJLB = IJFROMCHNK(KIJL,ICHNK)

        HSOIBCH(KIJS:KIJL, ICHNK) = HSOIB(IJSB:IJLB)
        IF ( KIJL < NPROMA_WAM ) THEN
        ! fictitious points
          HSOIBCH(KIJL+1:NPROMA_WAM, ICHNK) = HSOIBCH(1, ICHNK)
        ENDIF
 
        CALL UPDATEWFLD (1, NPROMA_WAM, MIJ(:,ICHNK),                         &
     &                   HSOIBCH(:,ICHNK), HSANCH(:,ICHNK), U10FGCH(:,ICHNK), &
     &                   WVENVI%DEPTH(:,ICHNK), FF_NOW%WSWAVE(:,ICHNK), FF_NOW%WDWAVE(:,ICHNK), &
     &                   FF_NOW%UFRIC(:,ICHNK), FF_NOW%CICOVER(:,ICHNK), FF_NOW%CITHICK(:,ICHNK), &
     &                   FF_NOW%Z0M(:,ICHNK), FF_NOW%Z0B(:,ICHNK), FF_NOW%CHRNCK(:,ICHNK), &
     &                   FF_NOW%TAUW(:,ICHNK), FF_NOW%TAUWDIR(:,ICHNK), FF_NOW%AIRD(:,ICHNK), &
     &                   FF_NOW%WSTAR(:,ICHNK), INTFLDS%USTOKES(:,ICHNK), INTFLDS%VSTOKES(:,ICHNK), &
     &                   INTFLDS%STRNMS(:,ICHNK), WAM2NEMO%NEMOUSTOKES(:,ICHNK), WAM2NEMO%NEMOVSTOKES(:,ICHNK), &
     &                   WAM2NEMO%NEMOSTRN(:,ICHNK), WVPRPT%WAVNUM(:,:,ICHNK), WVPRPT%STOKFAC(:,:,ICHNK), &
     &                   WVPRPT%CINV(:,:,ICHNK), WVPRPT%XK2CG(:,:,ICHNK), FL1(:,:,:,ICHNK), XLLWS(:,:,:,ICHNK) )

        U10AN(IJSB:IJLB) = FF_NOW%WSWAVE(KIJS:KIJL, ICHNK)
        HSAN(IJSB:IJLB) = HSANCH(KIJS:KIJL, ICHNK)
        U10FG(IJSB:IJLB) = U10FGCH(KIJS:KIJL, ICHNK)

      ENDDO
!$OMP END PARALLEL DO


!     OUTPUT TO ODB:
!     -------------
      IF (LODBRALT) THEN
#ifdef WITH_ODB
        CALL WAM2ODB(IJSG, IJLG, HSOIB, HSAN, U10FG, U10AN)
#endif
      ENDIF

!     CLEANING UP
!     -----------
      DEALLOCATE (NOBSPE)
      IF (ALLOCATED(DIFFALTFG)) DEALLOCATE(DIFFALTFG)
      IF (ALLOCATED(IJALT)) DEALLOCATE(IJALT)
      IF (ALLOCATED(ALTDATA)) DEALLOCATE(ALTDATA)
      IF (ALLOCATED(ALTEXDATA)) DEALLOCATE(ALTEXDATA)
      IF (ALLOCATED(ALTUNDATA)) DEALLOCATE(ALTUNDATA)
      IF (ALLOCATED(CDATEOBS)) DEALLOCATE(CDATEOBS)

      WRITE(IU06,*) '  ALTAS ENDS NORMALLY'
      WRITE(IU06,*) ' '
      CALL FLUSH(IU06)

IF (LHOOK) CALL DR_HOOK('ALTAS',1,ZHOOK_HANDLE)

END SUBROUTINE ALTAS
