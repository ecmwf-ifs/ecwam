      SUBROUTINE OUTCOM (IU07, IU17, IFORM, IREFRA, ML, KL)

! ----------------------------------------------------------------------

!**** *OUTCOM* - ROUTINE TO WRITE MODULES TO DISK

!     H.GUNTHER            ECMWF       04/04/1990

!     J. BIDLOT            ECMWF       10/1998.
!!!!!!                     COMMON BLOCKS HAVE BEEN CONVERTED TO MODULES

!*    PURPOSE.
!     -------

!       TO WRITE OUT THE COMPUTED MODULES
!       (MODULE UBUF IS WRITTEN IN OUTUBUF)

!**   INTERFACE.
!     ----------

!       *CALL* *OUTCOM (IU07, IU17, IFORM, IREFRA, ML, KL)*
!          *IU07*   - LOGICAL UNIT FOR  UNFORMATED WRITE.
!          *IU17*   - LOGICAL UNIT FOR    FORMATED WRITE.
!          *IFORM*   - FORMAT OPTION  = 1  UNFORMATED WRITE.
!                                     = 2  FORMATED WRITE.
!                                     OTHERWISE BOTH.
!          *IREFRA*  - REFRACTION OPTION.
!          *ML*      - NUMBER OF FREQUENCIES.
!          *KL*      - NUMBER OF DIRECTIONS.

!     METHOD.
!     -------

!       MODULES YOWPARAM, YOWCOUPL, YOWCURR, YOWFRED, YOWINDNL, YOWGRID,
!       YOWMAP, YOWCOUT, YOWTABL, AND YOWSHAL ARE WRITTEN TO UNIT.
!       ALL FREQUENCY AND DIRECTION DEPENDENT ARRAYS
!       ARE WRITTEN FROM 1 TO THE USED NUMBER OF FREQUENCIES (ML),
!       AND THE USED NUMBER OF DIRECTIONS (KL). OTHER ARRAYS ARE
!       WRITTEN ACCORDING TO THEIR DIMENSIONS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,NOVER    ,NIBL1    ,NIBLD    ,
     &            NBLD     ,NIBLC    ,NBLC     ,CLDOMAIN
      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA   ,
     &            XNLEV
      USE YOWCPBO  , ONLY : NBOUNC
      USE YOWCURR  , ONLY : U        ,V
      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,
     &            DELTH    ,DELTR    ,TH       ,COSTH    ,SINTH
      USE YOWFPBO  , ONLY : NBOUNF
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,
     &            NLONRGG  ,IGL      ,IJS      ,IJL2     ,IJLS     ,
     &            IJL      ,IJLT
      USE YOWINDN  , ONLY : IKP      ,IKP1     ,IKM      ,IKM1     ,
     &            K1W      ,K2W      ,K11W     ,K21W     ,AF11     ,
     &            FKLAP    ,FKLAP1   ,FKLAM    ,FKLAM1   ,ACL1     ,
     &            ACL2     ,CL11     ,CL21     ,DAL1     ,DAL2     ,
     &            FRH
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,NX       ,NY       ,
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWCOUT  , ONLY : NGOUT    ,IGAR     ,IJAR
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,
     &            TCGOND   ,TFAK     ,TSIHKD
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,IUSTAR   ,IALPHA   ,
     &            TAUT     ,DELTAUW  ,DELU     ,TAUHFT   ,DELUST   ,
     &            DELALP

! ----------------------------------------------------------------------

  997 FORMAT(12I8,1X,A1)
  998 FORMAT(10I8)
  999 FORMAT(5E16.7)

! ----------------------------------------------------------------------
!*    0. WRITE YOWPARAM (BLOCK SIZES).
!        ----------------------------

      IF (IFORM.NE.2) THEN
        WRITE(IU07) NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     &              NIBL1, IREFRA, NIBLD, NBLD, NIBLC, NBLC, CLDOMAIN
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE(IU17,997) NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     &              NIBL1, IREFRA, NIBLD, NBLD, NIBLC, NBLC, CLDOMAIN
      ENDIF


!*    1. WRITE MODULE YOWFRED 
!        --------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) (FR(M),M=1,ML), (DFIM(M),M=1,ML),
     &   (GOM(M),M=1,ML), (C(M),M=1,ML),
     &   DELTH, DELTR, (TH(K),K=1,KL),
     &   (COSTH(K),K=1,KL), (SINTH(K),K=1,KL)
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,999) (FR(M),M=1,ML), (DFIM(M),M=1,ML),
     &   (GOM(M),M=1,ML), (C(M),M=1,ML),
     &   DELTH, DELTR, (TH(K),K=1,KL),
     &   (COSTH(K),K=1,KL), (SINTH(K),K=1,KL)
      ENDIF

! ----------------------------------------------------------------------

!*    2. WRITE MODULE YOWGRID.
!        ---------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) DELPHI, (DELLAM(L),L=1,NY), (NLONRGG(L),L=1,NY),
     &   (SINPH(L),L=1,NY), (COSPH(L),L=1,NY),
     &   IGL, IJS, IJL2, IJLS, IJL, IJLT
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,999) DELPHI,(DELLAM(L),L=1,NY),(NLONRGG(L),L=1,NY),
     &   (SINPH(L),L=1,NY), (COSPH(L),L=1,NY)
        WRITE (IU17,998) IGL, IJS, IJL2, IJLS, IJL, IJLT
      ENDIF

! ----------------------------------------------------------------------

!*    3. WRITE MODULE YOWMAP.
!        --------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) IXLG, KXLT, NX, NY, IPER,
     &   AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO,
     &   ZDELLO, IRGG
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,998) IXLG, KXLT, NX, NY, IPER
        WRITE (IU17,999) AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &   XDELLA, XDELLO,ZDELLO
      ENDIF

! ----------------------------------------------------------------------

!*    4. WRITE MODULE YOWINDNL.
!        ---------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) (IKP(M),M=1,ML+4), (IKP1(M),M=1,ML+4),
     &   (IKM(M),M=1,ML+4), (IKM1(M),M=1,ML+4),
     &   ((K1W(K,L),K=1,KL),L=1,2),
     &   ((K2W(K,L),K=1,KL),L=1,2),
     &   ((K11W(K,L),K=1,KL),L=1,2),
     &   ((K21W(K,L),K=1,KL),L=1,2),
     &   (AF11(M),M=1,ML+4), (FKLAP(M),M=1,ML+4),
     &   (FKLAP1(M),M=1,ML+4), (FKLAM(M),M=1,ML+4),
     &   (FKLAM1(M),M=1,ML+4),
     &   ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,998) (IKP(M),M=1,ML+4), (IKP1(M),M=1,ML+4),
     &   (IKM(M),M=1,ML+4), (IKM1(M),M=1,ML+4),
     &   ((K1W(K,L),K=1,KL),L=1,2),
     &   ((K2W(K,L),K=1,KL),L=1,2),
     &   ((K11W(K,L),K=1,KL),L=1,2),
     &   ((K21W(K,L),K=1,KL),L=1,2)
        WRITE (IU17,999) (AF11(M),M=1,ML+4),   (FKLAP(M),M=1,ML+4),
     &   (FKLAP1(M),M=1,ML+4), (FKLAM(M),M=1,ML+4),
     &   (FKLAM1(M),M=1,ML+4),
     &   ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
      ENDIF

! ----------------------------------------------------------------------

!*    5. WRITE MODULE YOWCOUPL.
!        ----------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,999) BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV
      ENDIF

! ----------------------------------------------------------------------

!*    6. WRITE MODULE YOWTABL.
!        ---------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,999) TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP
      ENDIF

! ----------------------------------------------------------------------

!*    7. WRITE MODULE YOWCOUT.
!        ---------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07)  NGOUT
        IF(NGOUT.GT.0) WRITE (IU07)  IGAR, IJAR
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,998)  NGOUT
        IF(NGOUT.GT.0) WRITE (IU17,998)  IGAR, IJAR
      ENDIF

! ----------------------------------------------------------------------

!*    8. WRITE MODULE YOWSHAL.
!        --------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) DEPTH, DEPTHA, DEPTHD,
     &   ((TCGOND(L,M),L=1,NDEPTH),M=1,ML),
     &   ((TFAK(L,M),L=1,NDEPTH),M=1,ML),
     &   ((TSIHKD(L,M),L=1,NDEPTH),M=1,ML)
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,999) DEPTH, DEPTHA, DEPTHD,
     &   ((TCGOND(L,M),L=1,NDEPTH),M=1,ML),
     &   ((TFAK(L,M),L=1,NDEPTH),M=1,ML),
     &   ((TSIHKD(L,M),L=1,NDEPTH),M=1,ML)
      ENDIF

! ----------------------------------------------------------------------

!*    9. WRITE MODULE YOWCURR.
!        ---------------------

      IF (IREFRA.EQ.2) THEN
        IF (IFORM.NE.2) WRITE (IU07) U, V
        IF (IFORM.NE.1) WRITE (IU17,999) U, V
      ENDIF

! ----------------------------------------------------------------------

!*   10. WRITE NAMELIST PARWAM.
!        ----------------------
!     NBINP NOUTT have no meaning here since the size of the arrays
!     will be determined at input.
      NBINP=-1
      NOUTT=-1
      
 
      CALL OUTNAM
     & (NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER, NGOUT, NOUTT,
     &  NBOUNC, NBOUNF, NBINP, NIBL1, NIBLD, NBLD, NIBLC, NBLC ,
     &  ITAUMAX, JUMAX, IUSTAR, IALPHA, NDEPTH, IREFRA, iper)

      RETURN
      END SUBROUTINE OUTCOM
