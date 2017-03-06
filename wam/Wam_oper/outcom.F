      SUBROUTINE OUTCOM (IU07, IU17, IFORM)

! ----------------------------------------------------------------------

!**** *OUTCOM* - ROUTINE TO WRITE MODULES TO DISK

!     H.GUNTHER            ECMWF       04/04/1990

!     J. BIDLOT            ECMWF       10/1998.
!!!!!!                     COMMON BLOCKS HAVE BEEN CONVERTED TO MODULES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     J. BIDLOT            ECMWF       11/2003
!                          IF YOUR ARE RUNNING AT ECMWF:
!                          BE AWARE THAT IF YOU CHANGE ANYTHING TO THE
!                          STRUCTURE OF THE OUTPUT FILE YOU WILL HAVE TO
!                          MAKE SURE THAT IT IS CREATED FOR YOUR RUN, 
!                          OTHERWISE IT MIGHT PICK UP THE DEFAULT ONE
!                          THAT IS ALREADY ON DISK.
!                          YOU ALSO HAVE TO CHANGE READPRE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*    PURPOSE.
!     -------

!       TO WRITE OUT THE COMPUTED MODULES
!       (MODULE UBUF IS WRITTEN IN MUBUF)

!**   INTERFACE.
!     ----------

!       *CALL* *OUTCOM (IU07, IU17, IFORM)*
!          *IU07*   - LOGICAL UNIT FOR  UNFORMATED WRITE.
!          *IU17*   - LOGICAL UNIT FOR    FORMATED WRITE.
!          *IFORM*   - FORMAT OPTION  = 1  UNFORMATED WRITE.
!                                     = 2  FORMATED WRITE.
!                                     OTHERWISE BOTH.

!     METHOD.
!     -------

!       MODULES YOWPARAM, YOWCOUPL, YOWCURR, YOWFRED, YOWINDNL, YOWGRID,
!       YOWMAP, YOWCOUT, YOWTABL, AND YOWSHAL ARE WRITTEN TO UNIT.
!       ALL FREQUENCY AND DIRECTION DEPENDENT ARRAYS
!       ARE WRITTEN FROM 1 TO THE USED NUMBER OF FREQUENCIES (NFRE),
!       AND THE USED NUMBER OF DIRECTIONS (NANG). OTHER ARRAYS ARE
!       WRITTEN ACCORDING TO THEIR DIMENSIONS.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWGRIBHD, ONLY : IMDLGRBID_G,IMDLGRBID_M 
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,
     &            NBLO     ,NIBLO    ,NOVER    ,NIBL1    ,NIBLD    ,
     &            NBLD     ,NIBLC    ,NBLC     ,CLDOMAIN ,IMDLGRDID
      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA   ,
     &            XNLEV    ,TAUWSHELTER, ITSHELT
      USE YOWCPBO  , ONLY : NBOUNC
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
     &            FRH      ,KFRH     ,MFRSTLW  ,MLSTHG
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,NX       ,NY       ,
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWCOUT  , ONLY : NGOUT    ,IGAR     ,IJAR
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTH    ,DEPTHA   ,DEPTHD   ,
     &            TCGOND   ,TFAK     ,TSIHKD   ,TFAC_ST
      USE YOWSTAT  , ONLY : IPHYS
      USE YOWTABL  , ONLY : ITAUMAX  ,JUMAX    ,IUSTAR   ,IALPHA   ,
     &            TAUT     ,DELTAUW  ,DELU     ,TAUHFT   ,DELUST   ,
     &            TAUHFT2  ,DELTAIL  ,ILEVTAIL ,
     &            DELALP   ,FAC0     ,FAC1     ,FAC2     ,FAC3     ,
     &            FAK      ,FRHF     ,DFIMHF   ,
     &            MR       ,XMR      ,MA       ,XMA      ,NFREH    , 
     &            NANGH    ,NMAX     ,OMEGA    ,DFDTH    ,THH      ,
     &            DELTHH   ,IM_P     ,IM_M     ,TA       ,TB       ,
     &            TC_QL    ,TT_4M    ,TT_4P    ,TFAKH
      USE YOWUNPOOL ,ONLY : LLUNSTR  ,LPREPROC

      USE UNWAM
! ----------------------------------------------------------------------

      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IU07, IU17, IFORM
      INTEGER :: IDUM, K, M, L
      INTEGER :: NBINP, NOUTT
      INTEGER :: NKIND !Precision used when writing

! ----------------------------------------------------------------------

  995 FORMAT(5I8)
  996 FORMAT(I8,1X,2E16.7)
  997 FORMAT(12I8,1X,A1)
  998 FORMAT(10I8)
  999 FORMAT(5E16.7)

! ----------------------------------------------------------------------
!        WRITE IDENTIFIERS (MAKE SURE TO UPDATE THE VALUES IF YOU CHANGE
!        ANYTHING TO THE MODEL). 
      
      NKIND = KIND(DELPHI)

      IF (IFORM.NE.2) THEN
        WRITE(IU07) NKIND,IMDLGRDID, IMDLGRBID_G, IMDLGRBID_M
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE(IU17,998) NKIND,IMDLGRDID, IMDLGRBID_G, IMDLGRBID_M
      ENDIF
!*    0. WRITE YOWPARAM (BLOCK SIZES).
!        ----------------------------

!     IDUM REPLACES IREFRA WHICH IS NO LONGER USED IN PREPROC.
      IDUM=0

      IF (IFORM.NE.2) THEN
        WRITE(IU07) NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     &              KFRH, MFRSTLW, MLSTHG,
     &              NIBL1, IDUM, NIBLD, NBLD, NIBLC, NBLC, CLDOMAIN
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE(IU17,997) NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     &              KFRH, MFRSTLW, MLSTHG,
     &              NIBL1, IDUM, NIBLD, NBLD, NIBLC, NBLC, CLDOMAIN
      ENDIF


!*    1. WRITE MODULE YOWFRED 
!        --------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) (FR(M),M=1,NFRE), (DFIM(M),M=1,NFRE),
     &   (GOM(M),M=1,NFRE), (C(M),M=1,NFRE),
     &   DELTH, DELTR, (TH(K),K=1,NANG),
     &   (COSTH(K),K=1,NANG), (SINTH(K),K=1,NANG)
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,999) (FR(M),M=1,NFRE), (DFIM(M),M=1,NFRE),
     &   (GOM(M),M=1,NFRE), (C(M),M=1,NFRE),
     &   DELTH, DELTR, (TH(K),K=1,NANG),
     &   (COSTH(K),K=1,NANG), (SINTH(K),K=1,NANG)
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
        WRITE(IU07)(IKP(M),M=MFRSTLW,MLSTHG),
     &   (IKP1(M),M=MFRSTLW,MLSTHG),
     &   (IKM(M),M=MFRSTLW,MLSTHG), (IKM1(M),M=MFRSTLW,MLSTHG),
     &   ((K1W(K,L),K=1,NANG),L=1,2),
     &   ((K2W(K,L),K=1,NANG),L=1,2),
     &   ((K11W(K,L),K=1,NANG),L=1,2),
     &   ((K21W(K,L),K=1,NANG),L=1,2),
     &   (AF11(M),M=MFRSTLW,MLSTHG), (FKLAP(M),M=MFRSTLW,MLSTHG),
     &   (FKLAP1(M),M=MFRSTLW,MLSTHG), (FKLAM(M),M=MFRSTLW,MLSTHG),
     &   (FKLAM1(M),M=MFRSTLW,MLSTHG),
     &   ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE(IU17,998)(IKP(M),M=MFRSTLW,MLSTHG),
     &   (IKP1(M),M=MFRSTLW,MLSTHG),
     &   (IKM(M),M=MFRSTLW,MLSTHG), (IKM1(M),M=MFRSTLW,MLSTHG),
     &   ((K1W(K,L),K=1,NANG),L=1,2),
     &   ((K2W(K,L),K=1,NANG),L=1,2),
     &   ((K11W(K,L),K=1,NANG),L=1,2),
     &   ((K21W(K,L),K=1,NANG),L=1,2)
        WRITE(IU17,999)(AF11(M),M=MFRSTLW,MLSTHG),
     &    (FKLAP(M),M=MFRSTLW,MLSTHG),
     &   (FKLAP1(M),M=MFRSTLW,MLSTHG), (FKLAM(M),M=MFRSTLW,MLSTHG),
     &   (FKLAM1(M),M=MFRSTLW,MLSTHG),
     &   ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH
      ENDIF

! ----------------------------------------------------------------------

!*    5. WRITE MODULE YOWCOUPL.
!        ----------------------

      IF (ABS(TAUWSHELTER).LE.0.0) THEN
        ITSHELT=0
      ELSE
        ITSHELT=1
      ENDIF

      IF (IFORM.NE.2) THEN
        WRITE (IU07) IPHYS,BETAMAX,ZALP,ALPHA,TAUWSHELTER,ITSHELT,
     &    XKAPPA,XNLEV
      ENDIF
      IF (IFORM.NE.1) THEN
       WRITE (IU17,999)IPHYS,BETAMAX,ZALP,ALPHA,TAUWSHELTER,ITSHELT,
     &    XKAPPA,XNLEV
      ENDIF

! ----------------------------------------------------------------------

!*    6. WRITE MODULE YOWTABL (1st part).
!        ---------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) DELTAUW, DELU, DELUST, DELALP, DELTAIL 
        WRITE (IU07) TAUT
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,999) DELTAUW, DELU, DELUST, DELALP, DELTAIL 
        WRITE (IU17,999) TAUT
      ENDIF
      IF (ABS(TAUWSHELTER).LE.0.0) THEN
        IF (IFORM.NE.2) WRITE (IU07) TAUHFT
        IF (IFORM.NE.1) WRITE (IU17,999) TAUHFT
      ELSE
        IF (IFORM.NE.2) WRITE (IU07) TAUHFT2
        IF (IFORM.NE.1) WRITE (IU17,999) TAUHFT2
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
        WRITE (IU07) NDEPTH, DEPTHA, DEPTHD
        WRITE (IU07) DEPTH,
     &   ((TCGOND(L,M),L=1,NDEPTH),M=1,NFRE),
     &   ((TFAK(L,M),L=1,NDEPTH),M=1,NFRE),
     &   ((TSIHKD(L,M),L=1,NDEPTH),M=1,NFRE),
     &   ((TFAC_ST(L,M),L=1,NDEPTH),M=1,NFRE)
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,996) NDEPTH, DEPTHA, DEPTHD
        WRITE (IU17,999) DEPTH,
     &   ((TCGOND(L,M),L=1,NDEPTH),M=1,NFRE),
     &   ((TFAK(L,M),L=1,NDEPTH),M=1,NFRE),
     &   ((TSIHKD(L,M),L=1,NDEPTH),M=1,NFRE),
     &   ((TFAC_ST(L,M),L=1,NDEPTH),M=1,NFRE)
      ENDIF

! ----------------------------------------------------------------------

!*    9. WRITE MODULE YOWTABL (2ND AND 3RD PART).
!        ---------------------

      IF (IFORM.NE.2) THEN
        WRITE (IU07) FAC0,FAC1,FAC2,FAC3,FAK,FRHF,DFIMHF 
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,999) FAC0,FAC1,FAC2,FAC3,FAK,FRHF,DFIMHF 
      ENDIF

      IF (IFORM.NE.2) THEN
        WRITE (IU07) MR, XMR, MA, XMA, NFREH, NANGH, NMAX
 
        WRITE (IU07) OMEGA, DFDTH, THH, DELTHH, IM_P, IM_M, 
     &               TA, TB, TC_QL, TT_4M, TT_4P, TFAKH
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE (IU17,995) MR, MA, NFREH, NANGH, NMAX 
        WRITE (IU17,998) IM_P, IM_M 
        WRITE (IU17,999) XMR, XMA, OMEGA, DFDTH, THH, DELTHH,
     &                   TA, TB, TC_QL, TT_4M, TT_4P, TFAKH 
      ENDIF
! ----------------------------------------------------------------------

!*    10. WRITE MODULE YOWCURR.
!         ---------------------

!      THE CURRENTS ARE NO LONGER PART OF THE CONSTANT FILES.

! ----------------------------------------------------------------------

!*   11. WRITE NAMELIST PARWAM.
!        ----------------------
!     NBINP NOUTT have no meaning here since the size of the arrays
!     will be determined at input.
      NBINP=-1
      NOUTT=-1
      
 
      CALL OUTNAM
     & (NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER, NGOUT, NOUTT,
     &  KFRH, MFRSTLW, MLSTHG,
     &  NBOUNC, NBOUNF, NBINP, NIBL1, NIBLD, NBLD, NIBLC, NBLC ,
     &  ITAUMAX, JUMAX, IUSTAR, IALPHA, NDEPTH, IDUM, IPER)

      IF (LLUNSTR .AND. LPREPROC) THEN
        CALL UNWAM_OUT(IU07)
      END IF

      RETURN
      END SUBROUTINE OUTCOM
