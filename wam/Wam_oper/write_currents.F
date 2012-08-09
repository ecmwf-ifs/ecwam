! ----------------------------------------------------------------------

      PROGRAM WRITE_CURRENTS

! ----------------------------------------------------------------------

!**** *WRITE_CURRENTS* - 

!**   INTERFACE.
!     ----------

!          ---- FORMAL PARAMETERS ----

!     *WRITE_CURRENTS*

!          ---- INPUT/OUTPUT ---

!           *IU06*   - PRINTER OUTPUT.
!           *currents*   - INPUT CURRENTS
!           *wam_grid_tables*   - INPUT  UNIT OF PRECOMPUTED COMMON BLOCKS.
!                      (OUTPUT OF PREPROC). fort.7 must be there !
!           *wam_subgrid*   - INPUT  UNIT OF COMMON UBUF.
!                      (OUTPUT OF PREPROC). fort.8 must be there !
!           *IUOUT*  - THE OUTPUT IS DIRECTED TO fort.30 

!          MESSAGE PASSING NOT CATERED FOR YET.
!          IF MESSAGE PASSING ENABLED THE UNITS WILL BE REASSIGNED.

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE
!     ---------


! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWCURR  , ONLY : U        ,V
      USE YOWGRIBHD, ONLY : NTENCODE ,HOPERI
      USE YOWMAP   , ONLY : IXLG     ,KXLT
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &          MPMAXLENGTH,KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO    ,
     &            CLDOMAIN ,LL1D
      USE YOWPCONS , ONLY : DEG      ,ZMISS
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWSTAT  , ONLY : CFDBSF   ,NENSFNB  ,NTOTENS  , 
     &            MARSTYPE ,YCLASS   ,YEXPVER  ,NLOCGRB  ,
     &            IPROPAGS ,LSUBGRID
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH

! -------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IREAD, LIU, IJ, I, J, NPR, MAXLEN, IG
      INTEGER :: IU30, IPARAM, KFDB
      INTEGER :: IY, IM, ID, IH, IMN, ISS, IFCST
      INTEGER :: NLON,NLAT

      INTEGER I4(2)

      REAL Z4(2)
      REAL, ALLOCATABLE, DIMENSION(:,:) :: CURMAG, CURDIR
      REAL :: RSOUTW,RNORTW
      REAL :: PRPLRADI

      CHARACTER*14 :: CDATEIN
      CHARACTER*28 FILNM

      LOGICAL :: LLCURRENT, LLIRANK, LFDB

! ----------------------------------------------------------------------

!*    1. INITIAL VALUES
!        --------------

      PRPLRADI=1.0
      CALL INIWCST(PRPLRADI)

      IU06 = 6
      IU30=30

      NPROC = 1
      IRANK = 1
      LMESSPASS = .FALSE.

      IPROPAGS = 0
      LSUBGRID = .FALSE.

      LL1D = .FALSE.
      LLIRANK = .FALSE.

      LFDB = .FALSE.

      MARSTYPE  ='an'
      YCLASS    = 'rd'
      YEXPVER   = 'xxxx' 
      CLDOMAIN  = 'g'
      NENSFNB   = 0
      NTOTENS   = 0
      CFDBSF    = ""
      NLOCGRB   = 0
      CPATH  = '.' 
      ICPLEN=LEN_TRIM(CPATH)

      ITEST = 3


! ---------------------------------------------------------------------

!     2. GLOBAL SETTINGS.
!        ----------------

      Z4=1.
      NPRECR = KIND(Z4)
      I4=1
      NPRECI = KIND(I4)

      IREAD=1

      LWCOU = .FALSE.
      NLON=1
      NLAT=1

      CALL WVWAMINIT (LWCOU,IU06,NLON,NLAT,RSOUTW,RNORTW)

      NPR=NPROC

      CALL MPDECOMP(NPR,MAXLEN,LLIRANK)

      MPMAXLENGTH=MAXLEN

!     SET GRIB HEADERS FOR INPUTS/OUTPUTS
      CALL SETWGRIBHD

! ----------------------------------------------------------------------

!     3. GETTING CURRENTS ONTO WAM GRID.
!        -------------------------------

      ALLOCATE(U(NINF:NSUP,NBLO))
      ALLOCATE(V(NINF:NSUP,NBLO))

      LLCURRENT=.FALSE.

      FILNM = 'currents'
      LIU   = LEN_TRIM(FILNM)
      FILNM=FILNM(1:LIU)
      INQUIRE(FILE=FILNM,EXIST=LLCURRENT)

      IREAD=1
      KTAG=1

      IF (LLCURRENT) THEN
        CALL CURRENT2WAM (FILNM,IREAD,CDATEIN)
      ELSE
        U=0.
        V=0.
        WRITE(IU06,*)' '
        WRITE(IU06,*)'    ****************************'
        WRITE(IU06,*)'     FILE ',FILNM,' NOT FOUND '
        WRITE(IU06,*)'     CURRENTS WILL BE SET TO 0. '
        WRITE(IU06,*)'    ****************************'
        WRITE(IU06,*)' '
        CALL FLUSH(IU06)
      ENDIF

!     4. OUTPUT CURRENTS IN GRIB AS WAVE HEIGHT AND DIRECTION 
!        ----------------------------------------------------

      ALLOCATE(CURMAG(NGX,NGY))
      ALLOCATE(CURDIR(NGX,NGY))

      CURMAG=ZMISS
      CURDIR=ZMISS

      IG=1
      DO IJ = 1,NIBLO
        I = IXLG(IJ,IG)
        J = NGY- KXLT(IJ,IG) +1
        CURMAG(I,J) = SQRT(U(IJ,IG)**2+V(IJ,IG)**2)
        IF(CURMAG(I,J).GT.0.) THEN
          CURDIR(I,J) = ATAN2(U(IJ,IG),V(IJ,IG))
        ELSE
          CURDIR(I,J) = 0. 
        ENDIF
      ENDDO

!     convert directions to degrees and meteorological convention
      DO J=1,NGY
        DO I=1,NGX
          IF(CURDIR(I,J).NE.ZMISS)
     &      CURDIR(I,J) = MOD(DEG * CURDIR(I,J) +180.,360.)
        ENDDO
      ENDDO

      READ(CDATEIN,'(I4,5I2)') IY, IM, ID, IH, IMN, ISS
      IFCST = 0

      IPARAM=229
      CALL GRIBPAC(IU06,IU30,NTENCODE,CURMAG,NGX,NGY,IPARAM,
     &  IY,IM,ID,IH,IMN,
     &  IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,HOPERI,
     &  CFDBSF,NENSFNB,NTOTENS,KFDB)

      IPARAM=230
      CALL GRIBPAC(IU06,IU30,NTENCODE,CURDIR,NGX,NGY,IPARAM,
     &  IY,IM,ID,IH,IMN,
     &  IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,HOPERI,
     &  CFDBSF,NENSFNB,NTOTENS,KFDB)


      END
