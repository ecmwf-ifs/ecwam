      SUBROUTINE OUTSPP (FL3, FL1, IJS, IJL,IG, IU25, IU26,NSTART,NEND,
     &                   U10NEW,THWNEW,USNEW)

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

!        *CALL  OUTSPP (FL3, FL1, IG, IU25, IU26,,NSTART,NEND,
!                       U10NEW,THWNEW,USNEW)*
!           *FL3*    -   BLOCK OF SPECTRA.
!           *FL1*    -   SWELL SPECTRUM.
!           *IJS*    - INDEX OF FIRST GRIDPOINT
!           *IJL*    - INDEX OF LAST GRIDPOINT
!           *IG*     -   BLOCK NUMBER.
!           *IU25*   -   OUTPUT UNIT FOR SPECTRA.
!           *IU26*   -   OUTPUT UNIT FOR SEA AND SWELL SPECTRA.
!           *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!           *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!           *U10NEW*    NEW WIND SPEED IN M/S.
!           *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                       NOTATION (POINTING ANGLE OF WIND VECTOR,
!                       CLOCKWISE FROM NORTH).
!           *USNEW*     NEW FRICTION VELOCITY IN M/S.

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

      USE YOWCOUT  , ONLY : FFLAG    ,PFLAG    ,NGOUT    ,IGAR     ,IJAR
      USE YOWFRED  , ONLY : FR       ,TH       ,FRATIO
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,SMEAN
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : DEG      ,EPSMIN
      USE YOWSTAT  , ONLY : CDTPRO   ,CDTSPT   ,CDTSPS
      USE YOWSWEL  , ONLY : ESWELL   ,FSWELL   ,THSWELL  ,THWISEA
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL3,FL1
      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW
      CHARACTER NAME*20
      REAL,ALLOCATABLE :: SPEC(:,:)

! ----------------------------------------------------------------------

!     WORK ARRAYS
      REAL,ALLOCATABLE :: FLPTS(:,:,:,:)
      REAL,ALLOCATABLE,DIMENSION(:) :: U10,THW,US
      REAL,ALLOCATABLE,DIMENSION(:) :: A1, A2, A3, A4, A5, A6, A7, A8

!     0. GATHER THE FIELDS ON PROCESS 1 
!        ------------------------------

      ALLOCATE(SPEC(NANG,NFRE))
      IF(NGOUT.GE.1) THEN 
        ALLOCATE(U10(NIBLO),THW(NIBLO),US(NIBLO))
        DO IJ=NSTART(IRANK),NEND(IRANK)
          U10(IJ)=U10NEW(IJ)
          THW(IJ)=THWNEW(IJ)
          US(IJ)=USNEW(IJ)
        ENDDO
      ENDIF 
      IF((LMESSPASS).AND.(NGOUT.GE.1)) THEN 
        ALLOCATE(A1(NIBLO),A2(NIBLO),A3(NIBLO))
        ALLOCATE(A4(NIBLO),A5(NIBLO),A6(NIBLO),A7(NIBLO),A8(NIBLO))
        IRECV=1
        NSPFLD=2 
        NSCFLD=11 

        CALL LOC2GLO(EMEAN(IJS),A1,IJS,IJL,IG)
        CALL LOC2GLO(FMEAN(IJS),A2,IJS,IJL,IG)
        CALL LOC2GLO(THQ(IJS),A3,IJS,IJL,IG)

        IF(.NOT.ALLOCATED(SMEAN)) THEN
          ALLOCATE(SMEAN(IJS:IJL))
          SMEAN=0.
        ENDIF
        CALL LOC2GLO(SMEAN(IJS),A4,IJS,IJL,IG)

        IF(.NOT.ALLOCATED(ESWELL)) THEN
          ALLOCATE(ESWELL(IJS:IJL))
          ESWELL=0.
        ENDIF
        CALL LOC2GLO(ESWELL,A5,IJS,IJL,IG)

        IF(.NOT.ALLOCATED(FSWELL)) THEN
          ALLOCATE(FSWELL(IJS:IJL))
          FSWELL=0.
        ENDIF
        CALL LOC2GLO(FSWELL,A6,IJS,IJL,IG)

        IF(.NOT.ALLOCATED(THSWELL)) THEN
          ALLOCATE(THSWELL(IJS:IJL))
          THSWELL=0.
        ENDIF
        CALL LOC2GLO(THSWELL,A7,IJS,IJL,IG)

        IF(.NOT.ALLOCATED(THWISEA)) THEN
          ALLOCATE(THWISEA(IJS:IJL))
          THWISEA=0.
        ENDIF
        CALL LOC2GLO(THWISEA,A8,IJS,IJL,IG)

        ALLOCATE(FLPTS(NSPFLD,NGOUT,NANG,NFRE))
        CALL MPGATHERSPP(IRECV,25,NSTART,NEND,NSPFLD,NSCFLD,
     &                   FL3,FL1,FLPTS,
     &                   A1,A2,A3,A4,A5,A6,A7,
     &                   A8,U10,THW,US)

      ENDIF

!*    1. LOOP OVER OUTPUT POINTS.
!        ------------------------

      IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
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
!               use a formatted output for swamp cases
                IF(CLDOMAIN.EQ.'s') THEN
                  WRITE(IU25,'(1X,f7.3,1x,f8.3,1x,a12,3(1x,F8.4))') 
     &              XLAT,XLON,CDTPRO,U10NEW(IJ),THWNEW(IJ),USNEW(IJ)
                  dum=0.
                  WRITE(IU25,110) dum,(TH(K)*DEG,K=1,NANG)
                  DO M=1,NFRE
                    WRITE(IU25,110) FR(M),
     &                 (FL3(IJ,K,M),K=1,NANG)
                  ENDDO
110               FORMAT(25(1X,F8.4))
                ELSE
                  XANG = REAL(NANG)
                  XFRE = REAL(NFRE)
                  TH0 = TH(1)*DEG
                  WRITE(IU25) XLON, XLAT, CDTPRO, XANG, XFRE,
     &             TH0, FR(1), FRATIO 
                  WRITE(IU25) 4.*SQRT(EMEAN(IJ)), DEG*THQ(IJ),
     &             FMEAN(IJ), US(IJ), DEG*THW(IJ),
     &             U10(IJ)
                  IF(LMESSPASS) THEN
                    WRITE(IU25) ((FLPTS(1,NGOU,K,M),K=1,NANG),M=1,NFRE)
                  ELSE
                    WRITE(IU25) ((FL3(IJ,K,M),K=1,NANG),M=1,NFRE)
                  ENDIF
                ENDIF
              ENDIF

!*    1.1.2 PRINT SPECTRUM.
!           ---------------

              IF (PFLAG(17)) THEN
                DO M=1,NFRE
                  DO K=1,NANG
                    IF(LMESSPASS) THEN
                      SPEC(K,M) = FLPTS(1,NGOU,K,M)
                    ELSE
                      SPEC(K,M) = FL3(IJ,K,M)
                    ENDIF
                  ENDDO
                ENDDO
                WRITE (NAME,'(''BL= '',I3,'' POINT= '',I4)') IG,IJ
                CALL PRSPP (IU06, CDTPRO, XLON, XLAT, NAME,
     &           4.*SQRT(EMEAN(IJ)), DEG*THQ(IJ),
     &           FMEAN(IJ), US(IJ), DEG*THW(IJ),
     &           U10(IJ),
     &           FR, TH, SPEC, NANG, NFRE, NANG, NFRE)
              ENDIF
            ENDIF
!PL****
!*    1.2 OUTPUT SWELL INFORMATION.
!         -------------------------

            IF (CDTSPS.EQ.CDTPRO) THEN
              EWSEA=EMEAN(IJ)-ESWELL(IJ)
              IF(FMEAN(IJ).GT.EPSMIN .AND. FSWELL(IJ).GT.EPSMIN .AND.
     &           EWSEA.GT.EPSMIN ) THEN
                FWSEA=EWSEA/(EMEAN(IJ)/FMEAN(IJ)-ESWELL(IJ)/FSWELL(IJ))
              ELSE
                FWSEA=0.
              ENDIF
              EWSEA = 4.*SQRT(EWSEA)
              IF(LMESSPASS) THEN
                DO M=1,NFRE
                  DO K=1,NANG
                    SPEC(K,M) = FLPTS(2,NGOU,K,M)*FLPTS(1,NGOU,K,M)
                  ENDDO
                ENDDO
              ELSE
                DO M=1,NFRE
                  DO K=1,NANG
                    SPEC(K,M) = FL1(IJ,K,M)*FL3(IJ,K,M)
                  ENDDO
                ENDDO
              ENDIF

!*    1.2.1 FILE OUTPUT TO IU26.
!           --------------------

              IF (FFLAG(18)) THEN
                XANG = REAL(NANG)
                XFRE = REAL(NFRE)
                TH0 = TH(1)*DEG
                WRITE(IU26) XLON, XLAT, CDTPRO, XANG, XFRE,
     &           TH0, FR(1), FRATIO 
                WRITE(IU26) EWSEA, 4.*SQRT(ESWELL(IJ)),
     &           DEG*THWISEA(IJ), DEG*THSWELL(IJ),
     &           FWSEA, FSWELL(IJ)
                WRITE(IU26) ((SPEC(K,M),K=1,NANG),M=1,NFRE)
              ENDIF

!*    1.2.2 PRINT OUTPUT TO IU26.
!           ---------------------

              IF (PFLAG(18)) THEN
                WRITE (NAME,'(''BL= '',I3,'' POINT= '',I4)') IG,IJ
                CALL PRSPPS (IU06, CDTPRO, XLON, XLAT, NAME,
     &           EWSEA, DEG*THWISEA(IJ), FWSEA,
     &           4.*SQRT(ESWELL(IJ)), DEG*THSWELL(IJ),
     &           FSWELL(IJ), FR, TH, SPEC, NANG, NFRE, NANG, NFRE)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      IF((LMESSPASS).AND.(NGOUT.GE.1)) DEALLOCATE(FLPTS)
      IF(NGOUT.GE.1) DEALLOCATE(U10,THW,US)
      IF((LMESSPASS).AND.(NGOUT.GE.1)) DEALLOCATE(A1,A2,A3) 
      IF((LMESSPASS).AND.(NGOUT.GE.1)) DEALLOCATE(A4,A5,A6,A7,A8) 
      DEALLOCATE(SPEC)

      RETURN
      END SUBROUTINE OUTSPP
