      SUBROUTINE OUTSPP (FL3, IJS, IJL, IJ_OFFSET, IG, IU25, IU26)

! ----------------------------------------------------------------------

!**** *OUTSPP* -  OUTPUT OF SPECTRA , WAVE HEIGHT , WAVE DIRECTION,
!****             USTAR AND WIND DIRECTION AT GIVEN OUTPUT POINTS.

!      S.D.HASSELMANN.
!      P.JANSSEN       KNMI         FEBRUARY 1986

!      P. LIONELLO                  FEBRUARY 1987
!                  OUTPUT OF SWELL 2-D DISTRIBUTION , SWELL WAVE HEIGHT
!                  MEAN SWELL FREQUENCY , MEAN SWELL DIRECTION ,WIND-SE
!                  WAVE HEIGTH AND WIND SEA WAVE DIRECTION.

!      L. ZAMBRESKY   GKSS/ECMWF    JULY 89
!                 VARIANCE ENERGY SUMMED OVER FREQUENCY AND DIRECTION.

!      H. GUNTHER     ECMWF         OCTOBER 1990
!                 EXTENDED HEADER AND OUTPUT OF 1D SPECTRA IN PRINT.
!      J. BIDLOT      ECMWF         APRIL 1996 MESSAGE PASSING
!      J. BIDLOT      ECMWF         JANURARY 97 MODIFY DEFINITION OF TH0

!**   INTERFACE.
!     ----------

!        *CALL  OUTSPP (FL3, IG, IU25, IU26)
!           *FL3*    -   BLOCK OF SPECTRA.
!           *IJS*    - INDEX OF FIRST GRIDPOINT
!           *IJL*    - INDEX OF LAST GRIDPOINT
!           *IJ_OFFSET* 
!           *IG*     -   BLOCK NUMBER.
!           *IU25*   -   OUTPUT UNIT FOR SPECTRA.
!           *IU26*   -   OUTPUT UNIT FOR SEA AND SWELL SPECTRA.

!     EXTERNALS.
!     ----------

!       *PRSPP*     - PRINT A SPECTRUM.
!       *PRSPPS*    - PRINT A SWELL SPECTRUM.

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : FFLAG    ,PFLAG    ,NGOUT    ,IGAR     ,
     &            IJAR     ,LOUTSPP
      USE YOWFRED  , ONLY : FR       ,DELTH    ,TH       ,FRATIO
      USE YOWICE   , ONLY : CITHRSH
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,SMEAN
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : DEG      ,EPSMIN
      USE YOWSPEC, ONLY   : NSTART   ,NEND     , 
     &            U10NEW   ,THWNEW   ,USNEW
      USE YOWSTAT  , ONLY : CDTPRO   ,CDTSPT   ,CDTSPS
      USE YOWSWEL  , ONLY : ESWELL   ,FSWELL   ,THSWELL  ,THWISEA
      USE YOWTEST  , ONLY : IU06
      USE YOWUNPOOL, ONLY : LLUNSTR

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IG
      INTEGER, INTENT(IN) :: IJS, IJL, IJ_OFFSET, IU25, IU26
      REAL, DIMENSION(IJS:IJL,NANG,NFRE), INTENT(IN) :: FL3


      INTEGER :: IJ, M, K, NGOU
      INTEGER :: IX, KX
      INTEGER :: IRECV, NSPFLD,  NSCFLD

      REAL :: XLAT, XLON
      REAL :: ZHS, ZTHQ, ZFM, ZUST, ZTHW, ZU10
      REAL :: XANG, XFRE, TH0
      REAL :: EWSEA, FWSEA

      REAL, DIMENSION(NANG,NFRE) :: SPEC
      REAL,ALLOCATABLE :: FLPTS(:,:,:,:)
      REAL,ALLOCATABLE,DIMENSION(:) :: U10_G,THW_G,US_G
      REAL,ALLOCATABLE,DIMENSION(:) :: EMEAN_G, FMEAN_G, THQ_G,
     &                                 SMEAN_G, ESWELL_G, FSWELL_G,
     &                                 THSWELL_G, THWISEA_G

      CHARACTER(LEN=20) :: NAME
      CHARACTER(LEN=72) :: CPATH

! ----------------------------------------------------------------------

!     0. GATHER THE FIELDS ON PROCESS 1 
!        ------------------------------

      LOUTSPP=.TRUE.

      IF(NGOUT.GE.1) THEN 

        IF (LLUNSTR) THEN
         write(iu06,*) ' this bit is NOT ready for unstructured grid !'
         write(*,*) ' this bit is NOT ready for unstructured grid !'
!!! see how it was done in outgrid. we might not need ij_offset
         call abort1
        ENDIF



        ALLOCATE(U10_G(NIBLO),THW_G(NIBLO),US_G(NIBLO))
        DO IJ=IJS,IJL
          U10_G(IJ)=U10NEW(IJ)
          THW_G(IJ)=THWNEW(IJ)
          US_G(IJ)=USNEW(IJ)
        ENDDO

        ALLOCATE(EMEAN_G(NIBLO))
        CALL LOC2GLO(EMEAN(IJS),EMEAN_G,IJS,IJL,IJ_OFFSET,IG,CITHRSH)

        ALLOCATE(FMEAN_G(NIBLO))
        CALL LOC2GLO(FMEAN(IJS),FMEAN_G,IJS,IJL,IJ_OFFSET,IG,CITHRSH)

        ALLOCATE(THQ_G(NIBLO))
        CALL LOC2GLO(THQ(IJS),THQ_G,IJS,IJL,IJ_OFFSET,IG,CITHRSH)

        ALLOCATE(SMEAN_G(NIBLO))
        IF(.NOT.ALLOCATED(SMEAN)) THEN
          SMEAN_G(:)=0.
        ELSE
          CALL LOC2GLO(SMEAN(IJS),SMEAN_G,IJS,IJL,IJ_OFFSET,IG,CITHRSH)
        ENDIF

        ALLOCATE(ESWELL_G(NIBLO))
        IF(.NOT.ALLOCATED(ESWELL)) THEN
          ESWELL_G(:)=0.
        ELSE
          CALL LOC2GLO(ESWELL,ESWELL_G,IJS,IJL,IJ_OFFSET,IG,CITHRSH)
        ENDIF

        ALLOCATE(FSWELL_G(NIBLO))
        IF(.NOT.ALLOCATED(FSWELL)) THEN
          FSWELL_G(:)=0.
        ELSE
          CALL LOC2GLO(FSWELL,FSWELL_G,IJS,IJL,IJ_OFFSET,IG,CITHRSH)
        ENDIF

        ALLOCATE(THSWELL_G(NIBLO))
        IF(.NOT.ALLOCATED(THSWELL)) THEN
          THSWELL_G(:)=0.
        ELSE
          CALL LOC2GLO(THSWELL,THSWELL_G,IJS,IJL,IJ_OFFSET,IG,CITHRSH)
        ENDIF

        ALLOCATE(THWISEA_G(NIBLO))
        IF(.NOT.ALLOCATED(THWISEA)) THEN
          THWISEA_G(:)=0.
        ELSE
          CALL LOC2GLO(THWISEA,THWISEA_G,IJS,IJL,IJ_OFFSET,IG,CITHRSH)
        ENDIF

        IRECV=1
        NSPFLD=1 
        NSCFLD=11 
        ALLOCATE(FLPTS(NSPFLD,NGOUT,NANG,NFRE))
!!!!! not ready for unstructured grid !!!
        CALL MPGATHERSPP(IRECV,25,NSTART,NEND,NSPFLD,NSCFLD,
     &                   IJS, IJL, FL3, FLPTS,
     &                   EMEAN_G,FMEAN_G,THQ_G,SMEAN_G,ESWELL_G,
     &                   FSWELL_G,THSWELL_G,THWISEA_G,
     &                   U10_G,THW_G,US_G)

      ENDIF

!*    1. LOOP OVER OUTPUT POINTS.
!        ------------------------

      IF(IRANK.EQ.1) THEN

        DO NGOU=1,NGOUT
          IF (IG.EQ.IGAR(NGOU)) THEN
            IJ = IJAR(NGOU)
            IX = IXLG(IJ,IG)
            KX = KXLT(IJ,IG)
            XLON = AMOWEP + REAL(IX-1)*ZDELLO(KX)
            XLAT = AMOSOP + REAL(KX-1)*XDELLA

!*    1.1 TOTAL SPECTRUM.
!         ---------------

            IF (CDTSPT.EQ.CDTPRO) THEN

!*    1.1.1 FILE OUTPUT TO IU25.
!           --------------------

              IF (FFLAG(17)) THEN
!
!               USE SIMPLIFIED OUTPUT FOR SWAMP/SINGLE GRIDPOINT CASE
!
                IF(CLDOMAIN.EQ.'s') THEN
                  DO M=1,NFRE
                    DO K=1,NANG
                      SPEC(K,M) = FLPTS(1,NGOU,K,M)
                    ENDDO
                  ENDDO
                  ZHS = 4.*SQRT(MAX(EMEAN_G(IJ),0.0)) 
                  ZTHQ = THQ_G(IJ)
                  ZFM = FMEAN_G(IJ)
                  ZUST = US_G(IJ)
                  ZTHW = THW_G(IJ)
                  ZU10 = U10_G(IJ)
                  XFRE = REAL(NFRE)
                  TH0 = TH(1)*DEG

                  CALL OUT_ONEGRDPT_SP(SPEC,FR,TH,XLON,XLAT,CDTPRO,
     &              ZHS,ZTHQ,ZFM,ZUST,ZTHW,ZU10,IU06,IU25)

                ELSE
                  XANG = REAL(NANG)
                  XFRE = REAL(NFRE)
                  TH0 = TH(1)*DEG
                  WRITE(IU25) XLON, XLAT, CDTPRO, XANG, XFRE,
     &             TH0, FR(1), FRATIO 
                  WRITE(IU25) 4.*SQRT(MAX(EMEAN_G(IJ),0.0)),
     &             DEG*THQ_G(IJ),
     &             FMEAN_G(IJ), US_G(IJ), DEG*THW_G(IJ),
     &             U10_G(IJ)
                  WRITE(IU25) ((FLPTS(1,NGOU,K,M),K=1,NANG),M=1,NFRE)
                ENDIF
              ENDIF

!*    1.1.2 PRINT SPECTRUM.
!           ---------------

              IF (PFLAG(17)) THEN
                DO M=1,NFRE
                  DO K=1,NANG
                    SPEC(K,M) = FLPTS(1,NGOU,K,M)
                  ENDDO
                ENDDO
                WRITE (NAME,'(''BL= '',I3,'' POINT= '',I4)') IG,IJ
                CALL PRSPP (IU06, CDTPRO, XLON, XLAT, NAME,
     &           4.*SQRT(MAX(EMEAN_G(IJ),0.0)), DEG*THQ_G(IJ),
     &           FMEAN_G(IJ), US_G(IJ), DEG*THW_G(IJ),
     &           U10_G(IJ),
     &           FR, TH, SPEC, NANG, NFRE, NANG, NFRE)
              ENDIF
            ENDIF
!PL****
!*    1.2 OUTPUT SWELL INFORMATION.
!         -------------------------


            IF (CDTSPS.EQ.CDTPRO) THEN
              EWSEA=EMEAN_G(IJ)-ESWELL_G(IJ)
              IF(FMEAN_G(IJ).GT.EPSMIN .AND.
     &           FSWELL_G(IJ).GT.EPSMIN .AND.
     &           EWSEA.GT.EPSMIN ) THEN
                FWSEA=EWSEA/(EMEAN_G(IJ)/FMEAN_G(IJ)
     &               -ESWELL_G(IJ)/FSWELL_G(IJ))
              ELSE
                FWSEA=0.
              ENDIF
              EWSEA = 4.*SQRT(MAX(EWSEA,0.0))

!*    1.2.1 FILE OUTPUT TO IU26.
!           --------------------

              IF (FFLAG(18)) THEN
                XANG = REAL(NANG)
                XFRE = REAL(NFRE)
                TH0 = TH(1)*DEG
                WRITE(IU26) XLON, XLAT, CDTPRO, XANG, XFRE,
     &           TH0, FR(1), FRATIO 
                WRITE(IU26) EWSEA, 4.*SQRT(MAX(ESWELL_G(IJ),0.0)),
     &           DEG*THWISEA_G(IJ), DEG*THSWELL_G(IJ),
     &           FWSEA, FSWELL_G(IJ)

!       writing of swell specrtra was disabled
                SPEC=0.0
              ENDIF

!*    1.2.2 PRINT OUTPUT TO IU26.
!           ---------------------

              IF (PFLAG(18)) THEN
                WRITE (NAME,'(''BL= '',I3,'' POINT= '',I4)') IG,IJ
                CALL PRSPPS (IU06, CDTPRO, XLON, XLAT, NAME,
     &           EWSEA, DEG*THWISEA_G(IJ), FWSEA,
     &           4.*SQRT(MAX(ESWELL_G(IJ),0.0)), DEG*THSWELL_G(IJ),
     &           FSWELL_G(IJ), FR, TH, SPEC, NANG, NFRE, NANG, NFRE)
              ENDIF
            ENDIF
          ENDIF
        ENDDO

      ENDIF
      IF(NGOUT.GE.1) DEALLOCATE(FLPTS)
      IF(NGOUT.GE.1) DEALLOCATE(U10_G,THW_G,US_G)
      IF(NGOUT.GE.1) DEALLOCATE(EMEAN_G,FMEAN_G,THQ_G)
      IF(NGOUT.GE.1) DEALLOCATE(SMEAN_G,ESWELL_G,FSWELL_G)
      IF(NGOUT.GE.1) DEALLOCATE(THSWELL_G,THWISEA_G)

      END SUBROUTINE OUTSPP
