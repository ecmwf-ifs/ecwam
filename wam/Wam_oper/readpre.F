      SUBROUTINE READPRE (IU07, IREFRA)

! ----------------------------------------------------------------------

!**** *READPRE*  READ GRID OUTPUT FROM PREPROC.

!     H. GUNTHER      GKSS/ECMWF     MAY 1990

!*    PURPOSE.
!     --------

!       INPUT OF PREPROC GRID OUTPUT.

!**   INTERFACE.
!     ----------

!       *CALL* *READPRE (IU07, IREFRA)*
!          *IU07 *  - INPUT UNIT OF PREPROC GRID FILE.
!          *IREFRA *- REFRACTION OPTION.

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

      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA   ,
     &            XNLEV
      USE YOWCOUT  , ONLY : NGOUT    ,IGAR     ,IJAR
      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,
     &            DELTH    ,DELTR    ,TH       ,COSTH    ,SINTH
      USE YOWGRIBHD, ONLY : IMDLGRBID_G,IMDLGRBID_M
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,
     &            NLONRGG  ,IGL      ,IJS      ,IJL2     ,IJLS     ,
     &            IJL      ,IJLT
      USE YOWINDN  , ONLY : KFRH     ,IKP      ,IKP1     ,IKM      ,
     &            IKM1     ,K1W      ,K2W      ,K11W     ,K21W     ,
     &            AF11     ,FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,
     &            ACL1     ,ACL2     ,CL11     ,CL21     ,DAL1     ,
     &            DAL2     ,FRH
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,NX       ,NY       ,
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      , 
     &            NBLO     ,NIBLO    ,NOVER    ,NIBL1    ,NIBLD    ,
     &            NBLD     ,NIBLC    ,NBLC     ,CLDOMAIN ,IMDLGRDID 
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,
     &            TCGOND   ,TFAK     ,TSIHKD
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,IUSTAR   ,IALPHA   , 
     &            JPLEVT   ,TAUT     ,DELTAUW  ,DELU     ,TAUHFT   ,
     &            DELUST   ,DELALP   ,FAC0     ,FAC1     ,FAC2     ,
     &            FAC3     ,FAK      ,FRHF     ,DFIMHF   ,NFREHF
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------
      INTEGER :: KMDLGRDID, KMDLGRBID_G, KMDLGRBID_M

!     READ MODEL IDENTIFIERS
      READ(IU07) KMDLGRDID, KMDLGRBID_G, KMDLGRBID_M
      IF(KMDLGRDID.NE.IMDLGRDID .OR.
     &   KMDLGRBID_G.NE.IMDLGRBID_G .OR.
     &   KMDLGRBID_M.NE.IMDLGRBID_M) THEN
         WRITE(IU06,*) '*****************************************'
         WRITE(IU06,*) '*                                       *'
         WRITE(IU06,*) '*  FATAL ERROR(S) IN SUB. READPRE       *'
         WRITE(IU06,*) '*  ==============================       *'
         WRITE(IU06,*) '*                                       *'
         WRITE(IU06,*) '* THE PROGRAM HAS DETECTED DIFFERENT    *'
         WRITE(IU06,*) '* MODEL IDENTIFIERS. MAKE SURE YOU HAVE *' 
         WRITE(IU06,*) '* RUN PREPROC !!!!                      *'
         WRITE(IU06,*) KMDLGRDID, IMDLGRDID
         WRITE(IU06,*) KMDLGRBID_G, IMDLGRBID_G
         WRITE(IU06,*) KMDLGRBID_M, IMDLGRBID_M
         WRITE(IU06,*) '*                                       *'
         WRITE(IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.     *'
         WRITE(IU06,*) '* ---------------   --------------      *'
         WRITE(IU06,*) '*****************************************'
         CALL ABORT1
      ENDIF

!*    0. READ YOWPARAM (BLOCK SIZES). 
!        ----------------------------

      READ(IU07) NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     &           NIBL1, IREFRA, NIBLD, NBLD, NIBLC, NBLC, CLDOMAIN 

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
     &           DELTH, DELTR, TH, COSTH, SINTH

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
     &           IGL, IJS, IJL2, IJLS, IJL, IJLT

!*    3. READ MODULE YOWMAP (LONG. AND LAT. INDICES OF GRID POINTS).
!        --------------------------------------------------------

      IF(.NOT.ALLOCATED(IXLG)) ALLOCATE(IXLG(NIBLO,NBLO))
      IF(.NOT.ALLOCATED(KXLT)) ALLOCATE(KXLT(NIBLO,NBLO))
      IF(.NOT.ALLOCATED(ZDELLO)) ALLOCATE(ZDELLO(NGY))

      READ(IU07) IXLG, KXLT, NX, NY, IPER,
     &           AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,
     &           ZDELLO, IRGG

!*    4. READ MODULE YOWINDNL (NON-LINEAR INTERACTION).
!       -----------------------------------------------

      IF(.NOT.ALLOCATED(IKP)) ALLOCATE(IKP(NFRE+4))
      IF(.NOT.ALLOCATED(IKP1)) ALLOCATE(IKP1(NFRE+4))
      IF(.NOT.ALLOCATED(IKM)) ALLOCATE(IKM(NFRE+4))
      IF(.NOT.ALLOCATED(IKM1)) ALLOCATE(IKM1(NFRE+4))
      IF(.NOT.ALLOCATED(K1W)) ALLOCATE(K1W(NANG,2))
      IF(.NOT.ALLOCATED(K2W)) ALLOCATE(K2W(NANG,2))
      IF(.NOT.ALLOCATED(K11W)) ALLOCATE(K11W(NANG,2))
      IF(.NOT.ALLOCATED(K21W)) ALLOCATE(K21W(NANG,2))
      IF(.NOT.ALLOCATED(AF11)) ALLOCATE(AF11(NFRE+4))
      IF(.NOT.ALLOCATED(FKLAP)) ALLOCATE(FKLAP(NFRE+4))
      IF(.NOT.ALLOCATED(FKLAP1)) ALLOCATE(FKLAP1(NFRE+4))
      IF(.NOT.ALLOCATED(FKLAM)) ALLOCATE(FKLAM(NFRE+4))
      IF(.NOT.ALLOCATED(FKLAM1)) ALLOCATE(FKLAM1(NFRE+4))
      IF(.NOT.ALLOCATED(FRH)) ALLOCATE(FRH(KFRH))

      READ(IU07) IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W,
     &           AF11, FKLAP, FKLAP1, FKLAM, FKLAM1,
     &           ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH

!*    5. READ MODULE YOUCOUP.
!        -------------------

      READ(IU07) BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV

!*    6. READ MODULE YOWTABL (INDICES OF OUTPUT POINTS).
!        ---------------------------------------------

      IF(.NOT.ALLOCATED(TAUT)) ALLOCATE(TAUT(0:ITAUMAX,0:JUMAX,JPLEVT))
      IF(.NOT.ALLOCATED(TAUHFT)) ALLOCATE(TAUHFT(0:IUSTAR,0:IALPHA))

      READ(IU07) TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP

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

      IF(.NOT.ALLOCATED(DEPTH)) ALLOCATE(DEPTH(NIBLO,NBLO))
      IF(.NOT.ALLOCATED(TCGOND)) ALLOCATE(TCGOND(NDEPTH,NFRE))
      IF(.NOT.ALLOCATED(TFAK)) ALLOCATE(TFAK(NDEPTH,NFRE))
      IF(.NOT.ALLOCATED(TSIHKD)) ALLOCATE(TSIHKD(NDEPTH,NFRE))

      READ(IU07)  DEPTH, DEPTHA, DEPTHD, TCGOND, TFAK, TSIHKD

!*    9. READ MODULE YOWTABL (2nd part).
!        -------------------

      IF(.NOT.ALLOCATED(FAC0)) ALLOCATE(FAC0(NANG,NANG,NFREHF,NFREHF))
      IF(.NOT.ALLOCATED(FAC1)) ALLOCATE(FAC1(NANG,NANG,NFREHF,NFREHF))
      IF(.NOT.ALLOCATED(FAC2)) ALLOCATE(FAC2(NANG,NANG,NFREHF,NFREHF))
      IF(.NOT.ALLOCATED(FAC3)) ALLOCATE(FAC3(NANG,NANG,NFREHF,NFREHF))
      IF(.NOT.ALLOCATED(FAK)) ALLOCATE(FAK(NFREHF))
      IF(.NOT.ALLOCATED(FRHF)) ALLOCATE(FRHF(NFREHF))
      IF(.NOT.ALLOCATED(DFIMHF)) ALLOCATE(DFIMHF(NFREHF))

      READ(IU07) FAC0,FAC1,FAC2,FAC3,FAK,FRHF,DFIMHF

!*    10. READ MODULE YOWCURR.
!         --------------------

      IF (IREFRA.EQ.2) THEN 
         IF(.NOT.ALLOCATED(U)) ALLOCATE(U(0:NIBLC,NBLC))
         IF(.NOT.ALLOCATED(V)) ALLOCATE(V(0:NIBLC,NBLC))
         READ(IU07)  U, V
      ENDIF

      RETURN
      END SUBROUTINE READPRE
