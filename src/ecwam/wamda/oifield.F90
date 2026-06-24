SUBROUTINE OIFIELD (IJS, IJL, MINIJS, MAXIJL,                     &
 &                  BLK2GLO,                                      &
 &                  SIGMOD, DISTMAX, EXTENDMAX,                   &
 &                  DIST,                                         &
 &                  XMO, CICVR,                                   &
 &                  XOI)

!**** *OIFIELD* - OPTIMUM INTERPOLATION.                                

!     P.LIONELLO     ECMWF       APRIL 1990                             
!     J. BIDLOT      ECMWF       SEPTEMBER 96  MESSAGE PASSING
!     J. BIDLOT      ECMWF       OCTOBER 2009 : INTRODUCE VARIABLE DIST

!     PURPOSE.                                                          
!     --------                                                          

!       TO PRODUCE A MAP OF THE FIELD X, MERGING MEASUREMENT            
!       AND MODEL , BY OPTIMUM INTERPOLATION. THE ARRAY XOI
!       AT THE END OF THE SUBROUTINE CONTAINS THE VALUES                
!       TO BE USED TO ANALYSE THE SPECTRA, HAVING NEGATIVE RETURN       
!       CODES WHERE O.I. PRODUCED NO RESULTS.                           
!       THE FIRST GUESS FIELD XMO IS NOT MODIFIED                       
!       IN THIS SUBROUTINE.                                             

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *OIFIELD (IJS, IJL, MINIJS,MAXIJL,XMO,CICVR,
!    &                   BLK2GLO,
!    &                   XOI,SIGMOD,DIST,DISTMAX,EXTENDMAX)

!         *IJS*     INTEGER  MINIMUM INDEX OF XOI AND DIST (local indexing !)
!         *IJL*     INTEGER  MAXIMUM INDEX OF XOI AND DIST (local indexing !)
!         *MINIJS*  INTEGER  MINIMUM INDEX OF XMO AND CICVR (global indexing !)
!         *MAXIJL*  INTEGER  MAXIMUM INDEX OF XMO AND CICVR (global indexing !)
!         *XMO*     REAL     FIELD FROM MODEL (FIRST GUESS)(global indexing !)
!         *CICVR*   REAL     MODEL SEA ICE (global indexing !)
!         *BLK2GLO*-         BLOCK TO GRID TRANSFORMATION
!         *XOI*     REAL     FIELD FROM O.I. ON OUTPUT (local indexing).
!         *SIGMOD*  REAL     MODEL ERROR ESTIMATE.                             
!         *DIST*    REAL     CORRELATION DISTANCE (in radian)
!         *DISTMAX* REAL     MAXIMUM OF DIST
!         *EXTENDMAX* REAL   EXTENDMAX*DIST = MAXIMUM SPREADING DISTANCE

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWABORT , ONLY : WAM_ABORT
      USE YOWALTAS , ONLY : NALTAVLB ,IJALT    ,LALTPASSIV,             &
     &                      NOBSPE   ,ALTDATA  ,ALTUNDATA,              &
     &                      XLONOBS  ,SIGRATIO2,DIFFALTFG
      USE YOWICE   , ONLY : CITHRSH_SAT 
      USE YOWGRID  , ONLY : SINPH    ,COSPH    , NPROMA_WAM, NCHNK
      USE YOWMAP   , ONLY : IPER     ,XDELLA   ,ZDELLO   ,NLONRGG,     &
     &                      AMOWEP   ,NGX      ,NGY
      USE YOWMPP   , ONLY : IRANK
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWPCONS , ONLY : RAD      ,DEG      ,R
      USE YOWTEST  , ONLY : IU06

#ifdef WAM_HAVE_UNWAM
      USE YOWPD,     ONLY : NODES=>NODES_GLOBAL, RANK
#endif
      USE YOWSPHERE, ONLY : SPHERICAL_COORDINATE_DISTANCE

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR

! ----------------------------------------------------------------------
      IMPLICIT NONE

#include "wamoi.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IJS, IJL, MINIJS, MAXIJL 
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      REAL(KIND=JWRB), INTENT(IN) :: SIGMOD, DISTMAX, EXTENDMAX
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(IN) :: DIST
      REAL(KIND=JWRB), DIMENSION(MINIJS:MAXIJL), INTENT(IN) :: XMO, CICVR
      REAL(KIND=JWRB), DIMENSION(IJS:IJL), INTENT(OUT) :: XOI


      INTEGER(KIND=JWIM) :: IJ, IJG, K, KK, I, II, IC, IOBS, IALT, KALT, IIALT, IP
      INTEGER(KIND=JWIM) :: IIMIN, IIMAX, KKMIN, KKMAX, NBLK
      INTEGER(KIND=JWIM) :: LMAX, NDIM2, NOBSMAX, KMIN, KMAX
      INTEGER(KIND=JWIM) :: JKGLO, KIJS, KIJL, NPROMA
      INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: NOBS
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:) :: KALTMIN, KALTMAX
      INTEGER(KIND=JWIM), ALLOCATABLE, DIMENSION(:,:) :: IOBS4IJ

      REAL(KIND=JWRB) :: DIMAX, DOBS, DOBS2, DELLON, COSLON, XMODLON 
      REAL(KIND=JWRB) :: XII, XIIALT, ALTLON
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJS:IJL) :: HS_LOC
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: W

      REAL(KIND=JWRU)  :: XLONJD, XLATJD, DISTD, DISTMAXD
      REAL(KIND=JWRU),ALLOCATABLE,DIMENSION(:) :: XLONID, XLATID

      LOGICAL :: LLINVIEW
      LOGICAL :: LLINLATBAND
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LLSIMASK

! ----------------------------------------------------------------------
 IF (LHOOK) CALL DR_HOOK('OIFIELD',0,ZHOOK_HANDLE)

!*    1. INITIALISE SEARCH DISTANCE AND CORRELATION LENGTH.
!     -----------------------------------------------------

      DIMAX = EXTENDMAX*DISTMAX

      LMAX = INT(DEG*DIMAX/XDELLA)+1

!!! ??? may not be optimal for unstructured grid !!!
      NDIM2 = (5+2*(NALTAVLB-1))*(2*LMAX+1)

!      LOOP OVER GRID POINTS.                                          
!      ----------------------                                          

!     FIND THE INDEX OF ALL OBSERVATIONS INFLUENCING A GIVEN GRID POINT

      IF (NOBSPE(IRANK) > 0) THEN
        ALLOCATE(KALTMIN(NOBSPE(IRANK)))
        ALLOCATE(KALTMAX(NOBSPE(IRANK)))
        ALLOCATE(XLONOBS(NOBSPE(IRANK)))
        ALLOCATE(SIGRATIO2(NOBSPE(IRANK)))
        ALLOCATE(DIFFALTFG(NOBSPE(IRANK)))
      ENDIF

1000  CONTINUE
      ALLOCATE(IOBS4IJ(IJS:IJL,NDIM2))
      ALLOCATE(W(IJS:IJL,NDIM2))

      IOBS4IJ(IJS:IJL,1:NDIM2)=0
      NOBS(IJS:IJL)=0
      W(IJS:IJL,1:NDIM2)=1._JWRB

      NOBSMAX=0

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        ALLOCATE(XLONID(IJS:IJL))
        ALLOCATE(XLATID(IJS:IJL))
        DO IJ=IJS,IJL
!!! should be replaced with a permanent array that contains lat lon of local grid points
          IJG=RANK(IRANK)%IPLG(IJ)
          HS_LOC(IJ) = XMO(IJG)
          XLONID(IJ) = NODES(IJG)%X
          XLATID(IJ) = NODES(IJG)%Y
        ENDDO
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ELSE
        DO IJ=IJS,IJL
          HS_LOC(IJ) = XMO(IJ)
        ENDDO
        KMIN=NGY
        KMAX=1
        DO IOBS=1,NOBSPE(IRANK)
          IALT = BLK2GLO%IXLG(IJALT(IOBS,1))
          KALT = BLK2GLO%KXLT(IJALT(IOBS,1))
          KALTMIN(IOBS) = KALT-LMAX
          KMIN=MIN(KMIN,KALTMIN(IOBS))
          KALTMAX(IOBS) = KALT+LMAX
          KMAX=MAX(KMAX,KALTMAX(IOBS))
          XLONOBS(IOBS) = REAL(IALT-1)*ZDELLO(KALT)
        ENDDO
        KMIN=MAX(1,KMIN)
        KMAX=MIN(NGY,KMAX)
      ENDIF

      DO IOBS=1,NOBSPE(IRANK)
        SIGRATIO2(IOBS) = (ALTDATA(IOBS,2)/SIGMOD)**2
        IF (IJALT(IOBS,3) /= -1) THEN
          DIFFALTFG(IOBS) =  ALTDATA(IOBS,1)-XMO(IJALT(IOBS,1))
        ELSE
          DIFFALTFG(IOBS) = 0.0_JWRB 
        ENDIF
      ENDDO

!     SET LAND SEA ICE MASK (LLSIMASK=.true. if over land or sea ice)
!     BASED ON FIRST GUESS HS AND SEA ICE COVER WITH THRESHOLD CITHRSH_SAT
      IF (NOBSPE(IRANK) > 0) THEN
        IF (LLUNSTR) THEN
!!! we will need to code it !!!!
        ELSE
          IF (.NOT.ALLOCATED(LLSIMASK))THEN
             ALLOCATE(LLSIMASK(-LMAX:NGX+LMAX+1,KMIN:KMAX))
          ENDIF
          DO K = KMIN,KMAX
            DO I = 1,NLONRGG(K)
              LLSIMASK(I,K) = .TRUE. 
            ENDDO
          ENDDO
          DO IJ = MINIJS, MAXIJL
            K = BLK2GLO%KXLT(IJ)
            IF ( K >= KMIN .AND. K <= KMAX ) THEN
              I = BLK2GLO%IXLG(IJ)
              IF (XMO(IJ) > 0.0_JWRB .AND. CICVR(IJ) <= CITHRSH_SAT) THEN
                LLSIMASK(I,K) = .FALSE. 
              ENDIF
            ENDIF
          ENDDO
!         periodicity
          IF (IPER == 1) THEN
            DO K = KMIN,KMAX
              DO I = -LMAX,0 
                LLSIMASK(I,K) = LLSIMASK(NLONRGG(K)+I-1,K)
              ENDDO
              DO I = NLONRGG(K)+1,NLONRGG(K)+LMAX+1
               LLSIMASK(I,K) = LLSIMASK(I-NLONRGG(K),K)
              ENDDO
            ENDDO
          ENDIF

        ENDIF

      ENDIF

!     FIND HOW MANY VALID OBSERVATIONS CAN INFLUENCE EACH GRID POINT
      DO IOBS=1,NOBSPE(IRANK)
        IF (ALTDATA(IOBS,1) > 0._JWRB .AND. IJALT(IOBS,3) == 1) THEN
          IF (LLUNSTR) THEN
            XLATJD = ALTUNDATA(IOBS,1)
            XLONJD = ALTUNDATA(IOBS,2)
          ELSE
            IALT = BLK2GLO%IXLG(IJALT(IOBS,1))
            KALT = BLK2GLO%KXLT(IJALT(IOBS,1))
!           ALTIMETER LONGITUDE WITH RESPECT TO AMOWEP
            ALTLON=(IALT-1)*ZDELLO(KALT)
          ENDIF

          DO IJ = IJS, IJL
!           TAKE ONLY POSITIVE ACTIVE HS
!           OVER AREAS WITH SEA ICE COVER <= CITHRSH_SAT
            IF (HS_LOC(IJ) >= 0.01_JWRB .AND. .NOT.LALTPASSIV(IJALT(IOBS,2)))THEN

              IF (LLUNSTR) THEN
                CALL SPHERICAL_COORDINATE_DISTANCE(XLONID(IJ),XLONJD,XLATID(IJ),XLATJD,DISTD)
!!!! will need to add the test for latitude band
                LLINLATBAND=.TRUE.
                DOBS = REAL(DISTD,JWRB)
              ELSE
                I = BLK2GLO%IXLG(IJ)
                K = BLK2GLO%KXLT(IJ)
                LLINLATBAND=(K.GE.KALTMIN(IOBS).AND.K.LE.KALTMAX(IOBS))
                IF (IJ == IJALT(IOBS,1)) THEN
                  DOBS = 0._JWRB
                ELSE
                  DELLON = XLONOBS(IOBS) - REAL(I-1)*ZDELLO(K)
                  COSLON = COS(DELLON*RAD)
                  DOBS2  = COSLON*COSPH(KALT)*COSPH(K) + SINPH(KALT)*SINPH(K)
                  DOBS = ACOS(MAX(MIN(DOBS2,1.),-1.))
                ENDIF
              ENDIF

              IF (LLINLATBAND) THEN
!               FIND WHETHER THE OBSERVATION IS WITHIN THE MAXIMUM
!               SPREADING DISTANCE FROM THE GRID POINT.

                IF (DOBS <= EXTENDMAX*DIST(IJ)) THEN
!                 FIND IF LAND OR ICE IS NOT BLOCKING THE
!                 LINE OF VIEW BETWEEN OBSERVATION AND GRID POINT.
!!! put all this in a subroutine !!!1
                  IF (LLUNSTR) THEN

!!!                no land mask yet !!!
                    LLINVIEW=.TRUE.
                  ELSE
!                   MODEL LONGITUDE WITH RESPECT TO AMOWEP
                    XMODLON=(I-1)*ZDELLO(K)

                    LLINVIEW=.TRUE.

                    IF (K == KALT) THEN

                      ! periodicity ?
                      IF (IPER == 1 .AND. ABS(IALT-I) > NLONRGG(K)/2) THEN
                        IF (IALT > I) THEN
                          IIALT=IALT-NLONRGG(K)
                        ELSE
                          IIALT=IALT+NLONRGG(K)
                        ENDIF
                      ELSE
                        IIALT=IALT
                      ENDIF
                      IIMIN=MAX(MIN(I,IIALT),-LMAX)
                      IIMAX=MIN(MAX(I,IIALT),NGX+LMAX+1)
                      IC=IIMIN
                      DO WHILE (IC <= IIMAX .AND. LLINVIEW)
                        IF (LLSIMASK(IC,K)) THEN
                           LLINVIEW=.FALSE.
                           EXIT
                        ENDIF
                        IC=IC+1
                      ENDDO

                    ELSE
                      KKMIN=MIN(K+1,KALT)
                      KKMAX=MAX(K-1,KALT)
                      NBLK=0
                      DO KK=KKMIN,KKMAX
                        XII = XMODLON/ZDELLO(KK)+1.
                        II = MAX(1,INT(XII))
                        XIIALT = ALTLON/ZDELLO(KK)+1.
                        IIALT = MAX(1,INT(XIIALT))
                        ! periodicity ?
                        IF (IPER == 1 .AND. ABS(IIALT-II) > NLONRGG(KK)/2) THEN
                          IF (IIALT > II) THEN
                            IIALT=IIALT-NLONRGG(KK)
                          ELSE
                            IIALT=IIALT+NLONRGG(KK)
                          ENDIF
                        ENDIF
                        IIMIN=MAX(MIN(II,IIALT),-LMAX)
                        IIMAX=MIN(MAX(II,IIALT),NGX+LMAX+1)
                        IC=IIMIN
                        DO WHILE (IC <= IIMAX)
                          IF (LLSIMASK(IC,KK)) NBLK=NBLK+1
                          IC=IC+1
                        ENDDO
                      ENDDO
                      IF (NBLK >= (KKMAX-KKMIN+1)) LLINVIEW=.FALSE.

                    ENDIF
                  ENDIF

                  IF (LLINVIEW) THEN
                    NOBS(IJ)=NOBS(IJ)+1
                    IF (NOBS(IJ) <= NDIM2) THEN
                      IOBS4IJ(IJ,NOBS(IJ))=IOBS
                      W(IJ,NOBS(IJ)) = EXP(-DOBS/DIST(IJ))
                    ELSE
                      NOBSMAX=MAX(NOBSMAX,NOBS(IJ))
                    ENDIF
                  ENDIF

                ENDIF

              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO

!     IN CASE THE MAXIMUM NUMBER OF OBSERVATIONS INFLUENCING A POINT IS
!     TOO SMALL, RESET IT AND TRY AGAIN.
      IF (NOBSMAX > 0) THEN
        WRITE(NULERR,*) ''
        WRITE(NULERR,*) 'WARNING IN WAVE MODEL :'
        WRITE(NULERR,*) '**********************'
        WRITE(NULERR,*) 'IN OIFIELD : THE MAXIMUM NUMBER OF OBSERVATIONS'
        WRITE(NULERR,*) '             INFLUENCING A POINT WAS RESET !!!'
        WRITE(NULERR,*) 'IF THAT HAPPENS OFTEN, YOU SHOULD CHANGE NDIM2 !!'
        WRITE(NULERR,*) 'IRANK = ',IRANK
        WRITE(NULERR,*) 'NDIM2 = ',NDIM2
        WRITE(NULERR,*) 'LMAX = ',LMAX
        WRITE(NULERR,*) 'NOBSMAX = ',NOBSMAX
        WRITE(NULERR,*) ''
        WRITE(IU06,*) ''
        WRITE(IU06,*) 'IN OIFIELD : THE MAXIMUM NUMBER OF OBSERVATIONS'
        WRITE(IU06,*) '             INFLUENCING A POINT WAS RESET !!!'
        WRITE(IU06,*) 'IF THAT HAPPENS OFTEN, YOU SHOULD CHANGE NDIM2 !'
        WRITE(IU06,*) 'NDIM2 = ',NDIM2
        WRITE(IU06,*) 'LMAX = ',LMAX
        WRITE(IU06,*) 'NOBSMAX = ',NOBSMAX
        WRITE(IU06,*) ''
        CALL FLUSH(IU06)
        NDIM2=MIN(2*NDIM2,NOBSPE(IRANK))
        IF (ALLOCATED(IOBS4IJ)) DEALLOCATE(IOBS4IJ)
        IF (ALLOCATED(W)) DEALLOCATE(W)
        GOTO 1000
      ENDIF

      IF (ALLOCATED(KALTMIN)) DEALLOCATE(KALTMIN)
      IF (ALLOCATED(KALTMAX)) DEALLOCATE(KALTMAX)

!     LOOP ON EACH GRID POINT

      NPROMA=NPROMA_WAM
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO, KIJS, KIJL)
      DO JKGLO=IJS,IJL,NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1, IJL)
        CALL WAMOI(NOBS(KIJS),IOBS4IJ(KIJS:KIJL,:),W(KIJS:KIJL,:),      &
     &             KIJS,KIJL,NDIM2,                                     &
     &             BLK2GLO,                                             &
     &             DIST(KIJS), HS_LOC(KIJS), XOI(KIJS))
      ENDDO
!$OMP END PARALLEL DO

      IF (ALLOCATED(XLONOBS)) DEALLOCATE(XLONOBS)
      IF (ALLOCATED(SIGRATIO2)) DEALLOCATE(SIGRATIO2)
      IF (ALLOCATED(IOBS4IJ)) DEALLOCATE(IOBS4IJ)
      IF (ALLOCATED(W)) DEALLOCATE(W)
      IF (ALLOCATED(LLSIMASK)) DEALLOCATE(LLSIMASK)
      IF (ALLOCATED(XLONID)) DEALLOCATE(XLONID)
      IF (ALLOCATED(XLATID)) DEALLOCATE(XLATID)

IF (LHOOK) CALL DR_HOOK('OIFIELD',1,ZHOOK_HANDLE)

END SUBROUTINE OIFIELD
