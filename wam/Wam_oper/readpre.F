      SUBROUTINE READPRE (IU07)

! ----------------------------------------------------------------------

!**** *READPRE*  READ GRID OUTPUT FROM PREPROC.

!     H. GUNTHER      GKSS/ECMWF     MAY 1990
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT            ECMWF       11/2003
!                          IF YOUR ARE RUNNING AT ECMWF:
!                          BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                          STRUCTURE OF THE OUTPUT FILE YOU WILL HAVE TO
!                          MAKE SURE THAT IT IS CREATED FOR YOUR RUN,
!                          OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                          THAT IS ALREADY ON DISK.
!                          YOU ALSO HAVE TO CHANGE OUTCOM.
! ALSO CHECK
!            READPREB in Alt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*    PURPOSE.
!     --------

!       INPUT OF PREPROC GRID OUTPUT.

!**   INTERFACE.
!     ----------

!       *CALL* *READPRE (IU07)*
!          *IU07 *  - INPUT UNIT OF PREPROC GRID FILE.

!     METHOD.
!     -------

!       UNFORMATED READ FROM UNIT.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWALTAS , ONLY : EGRCRV   ,AGRCRV   ,BGRCRV   ,AFCRV    ,
     &            BFCRV    ,ESH      ,ASH      ,BSH      ,ASWKM    ,
     &            BSWKM
      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,
     &            XKAPPA   ,XNLEV    ,TAUWSHELTER, ITSHELT,
     &            TAILFACTOR, TAILFACTOR_PM
      USE YOWCOUT  , ONLY : NGOUT    ,IGAR     ,IJAR
      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,
     &            DELTH    ,DELTR    ,TH       ,COSTH    ,SINTH
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,
     &            NLONRGG  ,IGL      ,IJS      ,IJL2     ,IJLS     ,
     &            IJL      ,IJLT
      USE YOWINDN  , ONLY : KFRH     ,IKP      ,IKP1     ,IKM      ,
     &            IKM1     ,K1W      ,K2W      ,K11W     ,K21W     ,
     &            AF11     ,FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,
     &            ACL1     ,ACL2     ,CL11     ,CL21     ,DAL1     ,
     &            DAL2     ,FRH      ,MFRSTLW  ,MLSTHG
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,NX       ,NY       ,
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,IQGAUSS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      , 
     &            NBLO     ,NIBLO    ,NOVER    ,NIBL1    ,CLDOMAIN ,
     &            IMDLGRDID 
      USE YOWPHYS  , ONLY : ALPHAPMAX
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,
     &            TCGOND   ,TFAK     ,TSIHKD   ,TFAC_ST  ,TOOSHALLOW
      USE YOWSTAT  , ONLY : IPHYS
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,IUSTAR   ,IALPHA   , 
     &            JPLEVT   ,
     &            FAC0     ,FAC1     ,FAC2     ,
     &            FAC3     ,FAK      ,FRHF     ,DFIMHF   ,NFREHF   ,
     &            MR       ,XMR      ,MA       ,XMA      ,NFREH    , 
     &            NANGH    ,NMAX     ,OMEGA    ,DFDTH    ,THH      ,
     &            DELTHH   ,IM_P     ,IM_M     ,TA       ,TB       ,
     &            TC_QL    ,TT_4M    ,TT_4P    ,TFAKH
      USE YOWTEST  , ONLY : IU06
      USE YOWUNPOOL, ONLY : LLUNSTR  ,LPREPROC
      USE UNWAM
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IU07
      INTEGER :: IREAD
      INTEGER :: IDUM, KIBLD, KBLD, KIBLC, KBLC
      INTEGER :: KMDLGRDID, KMDLGRBID_G, KMDLGRBID_M

      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('READPRE',0,ZHOOK_HANDLE)
#endif

      IREAD=1
      KTAG=1

      CALL GSTATS(1771,0)
      IF(IRANK.EQ.IREAD) THEN
!       READ MODEL IDENTIFIERS
        READ(IU07,ERR=1000) KMDLGRDID, KMDLGRBID_G, KMDLGRBID_M
        IF(KMDLGRDID.NE.IMDLGRDID) THEN
          WRITE(IU06,*) '*****************************************'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '*  FATAL ERROR(S) IN SUB. READPRE       *'
          WRITE(IU06,*) '*  ==============================       *'
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* THE PROGRAM HAS DETECTED DIFFERENT    *'
          WRITE(IU06,*) '* MODEL GRID IDENTIFIER.                *' 
          WRITE(IU06,*) '* MAKE SURE YOU HAVE RUN PREPROC !!!!   *'
          WRITE(IU06,*)    KMDLGRDID, IMDLGRDID
          WRITE(IU06,*) '*                                       *'
          WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.     *'
          WRITE(IU06,*) '* ---------------   --------------      *'
          WRITE(IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF

!*    0. READ YOWPARAM (BLOCK SIZES). 
!        ----------------------------

        READ(IU07) NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     &             KFRH, MFRSTLW, MLSTHG,
     &             NIBL1, IDUM, KIBLD, KBLD, KIBLC, KBLC, CLDOMAIN 

!*    1. READ MODULE YOWFRED (FREQUENCY DIRECTION GRID).
!        ----------------------------------------------

        IF(.NOT.ALLOCATED(FR)) ALLOCATE(FR(NFRE))
        IF(.NOT.ALLOCATED(DFIM)) ALLOCATE(DFIM(NFRE))
        IF(.NOT.ALLOCATED(GOM)) ALLOCATE(GOM(NFRE))
        IF(.NOT.ALLOCATED(C)) ALLOCATE(C(NFRE))
        IF(.NOT.ALLOCATED(TH)) ALLOCATE(TH(NANG))
        IF(.NOT.ALLOCATED(COSTH)) ALLOCATE(COSTH(NANG))
        IF(.NOT.ALLOCATED(SINTH)) ALLOCATE(SINTH(NANG))

        READ(IU07) FR, DFIM, GOM, C,
     &             DELTH, DELTR, TH, COSTH, SINTH

!*    2. READ MODULE YOWGRID (GENERAL GRID ORGANISATION).
!        ------------------------------------------------

        IF(.NOT.ALLOCATED(DELLAM)) ALLOCATE(DELLAM(NGY))
        IF(.NOT.ALLOCATED(NLONRGG)) ALLOCATE(NLONRGG(NGY))
        IF(.NOT.ALLOCATED(SINPH)) ALLOCATE(SINPH(NGY))
        IF(.NOT.ALLOCATED(COSPH)) ALLOCATE(COSPH(NGY))
        IF(.NOT.ALLOCATED(IJS)) ALLOCATE(IJS(NBLO))
        IF(.NOT.ALLOCATED(IJL2)) ALLOCATE(IJL2(NBLO))
        IF(.NOT.ALLOCATED(IJLS)) ALLOCATE(IJLS(NBLO))
        IF(.NOT.ALLOCATED(IJL)) ALLOCATE(IJL(NBLO))
        IF(.NOT.ALLOCATED(IJLT)) ALLOCATE(IJLT(NBLO))

        READ(IU07) DELPHI, DELLAM, NLONRGG, SINPH, COSPH,
     &             IGL, IJS, IJL2, IJLS, IJL, IJLT

!*    3. READ MODULE YOWMAP (LONG. AND LAT. INDICES OF GRID POINTS).
!        --------------------------------------------------------

        IF(ALLOCATED(IXLG)) DEALLOCATE(IXLG)
        ALLOCATE(IXLG(NIBLO,NBLO))
        IF(ALLOCATED(KXLT)) DEALLOCATE(KXLT)
        ALLOCATE(KXLT(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(ZDELLO)) ALLOCATE(ZDELLO(NGY))

        READ(IU07) IXLG, KXLT, NX, NY, IPER,
     &             AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,
     &             ZDELLO, IRGG

!     DETERMINE IF WE ARE USING A QUASI GAUSSIAN GRID OR 
!     LAT-LONG GRID (REGULAR OR IRREGULAR).

        IF(IPER.EQ.1 .AND. AMONOP.EQ.ABS(AMOSOP) .AND.
     &     MOD(NY,2).EQ.0 .AND. IRGG.EQ.1 ) THEN
          IQGAUSS=1
        ELSE
          IQGAUSS=0
        ENDIF

!*    4. READ MODULE YOWINDNL (NON-LINEAR INTERACTION).
!       -----------------------------------------------

        IF(.NOT.ALLOCATED(IKP)) ALLOCATE(IKP(MFRSTLW:MLSTHG))
        IF(.NOT.ALLOCATED(IKP1)) ALLOCATE(IKP1(MFRSTLW:MLSTHG))
        IF(.NOT.ALLOCATED(IKM)) ALLOCATE(IKM(MFRSTLW:MLSTHG))
        IF(.NOT.ALLOCATED(IKM1)) ALLOCATE(IKM1(MFRSTLW:MLSTHG))
        IF(.NOT.ALLOCATED(K1W)) ALLOCATE(K1W(NANG,2))
        IF(.NOT.ALLOCATED(K2W)) ALLOCATE(K2W(NANG,2))
        IF(.NOT.ALLOCATED(K11W)) ALLOCATE(K11W(NANG,2))
        IF(.NOT.ALLOCATED(K21W)) ALLOCATE(K21W(NANG,2))
        IF(.NOT.ALLOCATED(AF11)) ALLOCATE(AF11(MFRSTLW:MLSTHG))
        IF(.NOT.ALLOCATED(FKLAP)) ALLOCATE(FKLAP(MFRSTLW:MLSTHG))
        IF(.NOT.ALLOCATED(FKLAP1)) ALLOCATE(FKLAP1(MFRSTLW:MLSTHG))
        IF(.NOT.ALLOCATED(FKLAM)) ALLOCATE(FKLAM(MFRSTLW:MLSTHG))
        IF(.NOT.ALLOCATED(FKLAM1)) ALLOCATE(FKLAM1(MFRSTLW:MLSTHG))
        IF(.NOT.ALLOCATED(FRH)) ALLOCATE(FRH(KFRH))

        READ(IU07) IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W,
     &             AF11, FKLAP, FKLAP1, FKLAM, FKLAM1,
     &             ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH

!*    5. READ MODULE YOUCOUP.
!        -------------------

        READ(IU07) IPHYS,BETAMAX,ZALP,ALPHA,
     &     ALPHAPMAX,
     &     TAUWSHELTER,ITSHELT,
     &     TAILFACTOR, TAILFACTOR_PM, XKAPPA,XNLEV

!*    5.1 READ MODULE YOWALTAS
!         --------------------

        READ(IU07) EGRCRV,AGRCRV,BGRCRV,AFCRV,BFCRV,
     &             ESH,ASH,BSH,ASWKM,BSWKM


!*    7. READ MODULE YOWCOUT (INDICES OF OUTPUT POINTS).
!        --------------------------------------------

        READ(IU07)  NGOUT
        IF(NGOUT.GT.0) THEN
          IF(.NOT.ALLOCATED(IGAR)) ALLOCATE(IGAR(NGOUT))
          IF(.NOT.ALLOCATED(IJAR)) ALLOCATE(IJAR(NGOUT))
          READ(IU07)  IGAR, IJAR
        ENDIF

!*    8. READ MODULE YOWSHAL (DEPTH AND SHALLOW WATER TABLES).
!        ----------------------------------------------------

        READ(IU07)  NDEPTH, DEPTHA, DEPTHD
!!!   note that the size of DEPTH will be readjusted in mpdecomp !!!
        IF(ALLOCATED(DEPTH)) DEALLOCATE(DEPTH)
        ALLOCATE(DEPTH(NIBLO,NBLO))
        IF(.NOT.ALLOCATED(TCGOND)) ALLOCATE(TCGOND(NDEPTH,NFRE))
        IF(.NOT.ALLOCATED(TFAK)) ALLOCATE(TFAK(NDEPTH,NFRE))
        IF(.NOT.ALLOCATED(TSIHKD)) ALLOCATE(TSIHKD(NDEPTH,NFRE))
        IF(.NOT.ALLOCATED(TFAC_ST)) ALLOCATE(TFAC_ST(NDEPTH,NFRE))

        READ(IU07)  DEPTH, TCGOND, TFAK, TSIHKD, TFAC_ST

!*    9. READ MODULE YOWTABL (2ND AND 3RD part).
!        -------------------

        IF(.NOT.ALLOCATED(FAC0)) ALLOCATE(FAC0(NANG,NANG,NFREHF,NFREHF))
        IF(.NOT.ALLOCATED(FAC1)) ALLOCATE(FAC1(NANG,NANG,NFREHF,NFREHF))
        IF(.NOT.ALLOCATED(FAC2)) ALLOCATE(FAC2(NANG,NANG,NFREHF,NFREHF))
        IF(.NOT.ALLOCATED(FAC3)) ALLOCATE(FAC3(NANG,NANG,NFREHF,NFREHF))
        IF(.NOT.ALLOCATED(FAK)) ALLOCATE(FAK(NFREHF))
        IF(.NOT.ALLOCATED(FRHF)) ALLOCATE(FRHF(NFREHF))
        IF(.NOT.ALLOCATED(DFIMHF)) ALLOCATE(DFIMHF(NFREHF))

        READ(IU07) FAC0,FAC1,FAC2,FAC3,FAK,FRHF,DFIMHF

        READ(IU07) MR, XMR, MA, XMA, NFREH, NANGH, NMAX

        IF(.NOT.ALLOCATED(OMEGA)) ALLOCATE(OMEGA(NFREH))
        IF(.NOT.ALLOCATED(THH))   ALLOCATE(THH(NANGH))
        IF(.NOT.ALLOCATED(DFDTH)) ALLOCATE(DFDTH(NFREH))
        IF(.NOT.ALLOCATED(TA)) ALLOCATE(TA(NDEPTH,NANGH,NFREH,NFREH))
        IF(.NOT.ALLOCATED(TB)) ALLOCATE(TB(NDEPTH,NANGH,NFREH,NFREH))
        IF(.NOT.ALLOCATED(TC_QL)) 
     &         ALLOCATE(TC_QL(NDEPTH,NANGH,NFREH,NFREH))
        IF(.NOT.ALLOCATED(TT_4M))
     &         ALLOCATE(TT_4M(NDEPTH,NANGH,NFREH,NFREH))
        IF(.NOT.ALLOCATED(TT_4P))
     &         ALLOCATE(TT_4P(NDEPTH,NANGH,NFREH,NFREH))
        IF(.NOT.ALLOCATED(IM_P))
     &         ALLOCATE(IM_P(NFREH,NFREH))
        IF(.NOT.ALLOCATED(IM_M)) 
     &         ALLOCATE(IM_M(NFREH,NFREH))
        IF(.NOT.ALLOCATED(TFAKH)) ALLOCATE(TFAKH(NFREH,NDEPTH))

        READ(IU07) OMEGA, DFDTH, THH, DELTHH, IM_P, IM_M, 
     &             TA, TB, TC_QL, TT_4M, TT_4P, TFAKH


!*      10. READ MODULE YOWCURR.
!           --------------------

!       THE OCEAN CURRENTS ARE INPUT IN INITMDL

!       THE UNSTRUCTURED BITS (if pre-computed by PREPROC)
        IF (LLUNSTR .AND. LPREPROC) THEN
          CALL SET_UNWAM_HANDLES
          CALL UNWAM_IN(IU07)
        END IF

      ENDIF

      CLOSE (UNIT=IU07)

      CALL GSTATS(1771,1)

!     SEND INFORMATION FROM READPRE TO ALL PE's
      CALL GSTATS(694,0)
      CALL MPBCASTGRID(IU06,IREAD,KTAG)
      CALL GSTATS(694,1)

      TOOSHALLOW=0.1*DEPTHA

!     RECALCULATE THE UNSTRUCTURED BITS (if not read in)

      IF (LLUNSTR .AND. .NOT. LPREPROC) THEN
        CALL INIT_UNWAM
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('READPRE',1,ZHOOK_HANDLE)
#endif
      RETURN

1000  CONTINUE
      WRITE(IU06,*) '*************************************'
      WRITE(IU06,*) '*                                   *'
      WRITE(IU06,*) '*  READ ERROR IN SUB. READPRE       *'
      WRITE(IU06,*) '*  ==========================       *'
      WRITE(IU06,*) '*                                   *'
      WRITE(IU06,*) '*  READ ERROR TO UNIT fort.',IU07
      WRITE(IU06,*) '*  IS THE FILE PRESENT ????         *' 
      WRITE(IU06,*) '*                                   *'
      WRITE(IU06,*) '*************************************'
      CALL ABORT1

      END SUBROUTINE READPRE
