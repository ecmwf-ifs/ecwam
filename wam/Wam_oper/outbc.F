      SUBROUTINE OUTBC (FL3, IJS, IJL, IG, IU19, THWNEW, USNEW)

! ----------------------------------------------------------------------

!**** *OUTBC* - OUTPUT OF THE COARSE GRID BOUNDARY VALUES.

!     R. PORTZ     MPI          JANUARY 1991
!     J. BIDLOT    ECMWF        FEBRARY 1996  MESSAGE PASSING

!*    PURPOSE.
!     --------

!        WRITE THE BOUNDARY VALUE OUTPUT FILE.

!**   INTERFACE.
!     ----------

!    *CALL* *OUTBC (FL3, IJS, IJL, IG, IU19, THWNEW, USNEW)*
!      *FL3*     - BLOCK OF SPECTRA.
!      *IJS*     - INDEX OF FIRST GRIDPOINT.
!      *IJL*     - INDEX OF LAST GRIDPOINT.
!      *IG*      - BLOCK NUMBER.
!      *IU19*    - OUTPUT UNIT OF BOUNDARY VALUES.
!      *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!      *USNEW*     NEW FRICTION VELOCITY IN M/S.


!     METHOD.
!     -------

!       SEQUENCIAL UNFORMATED WRITE TO UNIT.

!     EXTERNALS.
!     ----------

!       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
!       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
!       *STHQ*      - COMPUTATION OF MEAN WAVE DIRECTION AT EACH
!                     GRID POINT.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWCBOU  , ONLY : NBOUNC   ,IJARC    ,IGARC
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ      ,AKMEAN
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB

! ----------------------------------------------------------------------

      REAL :: FL3(NINF-1:NSUP,NANG,NFRE)

      REAL,DIMENSION(NINF:NSUP) :: THWNEW,USNEW
      LOGICAL :: LLAK

! ----------------------------------------------------------------------

!*    1. COMPUTE MEAN PARAMETERS.
!        ------------------------

      LLAK=.FALSE.
      CALL FEMEAN (FL3, IJS, IJL, EMEAN(IJS), FMEAN(IJS),
     &             AKMEAN(IJS), LLAK)
      CALL STHQ   (FL3, IJS, IJL, THQ(IJS))
      IF (ITEST.GE.3) THEN
        IF (ITESTB.GE.IG) THEN
          WRITE(IU06,*) '      SUB. OUTBC: INTEGRATED',
     &     ' PARAMETERS COMPUTED FOR BOUNDARY POINTS OUTPUT'
        ENDIF
      ENDIF

!*    2. LOOP OVER BOUNDARY POINTS.
!        --------------------------

      DO NGOU=1,NBOUNC
        IF (IG .LT. IGARC(NGOU)) EXIT 
        IF (IG.EQ.IGARC(NGOU)) THEN
          IJ = IJARC(NGOU)
          IX = IXLG(IJ,IG)
          KX = KXLT(IJ,IG)
          XLON = AMOWEP + REAL(IX-1)*ZDELLO(KX)
          XLAT = AMOSOP + REAL(KX-1)*XDELLA

!*    2.1 WRITE INFORMATION TO FILE IU19.
!         -------------------------------

          XIJ   = REAL(IJ)
          XIG   = REAL(IG)
          WRITE(IU19) XLON, XLAT, CDTPRO, EMEAN(IJ),
     &     THQ(IJ), FMEAN(IJ), USNEW(IJ), THWNEW(IJ)
          WRITE(IU19) ((FL3(IJ,K,M),K=1,NANG),M=1,NFRE)
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE OUTBC
