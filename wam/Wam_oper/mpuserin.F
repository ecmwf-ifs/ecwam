! ----------------------------------------------------------------------

      SUBROUTINE MPUSERIN (IU05, LMESSPASS, LFDB, LWCOU)

! ----------------------------------------------------------------------

!**** *MPUSERIN* - ROUTINE TO READ LMESSPASS AND LFDB FROM NAMELIST
!****              INPUT.

!     J. BIDLOT      ECMWF          JULY    1996
!     B. HANSEN      ECMWF          APRIL   1997
!     B. HANSEN      ECMWF          JANUARY 1998 NAMELIST INPUT.


!*    PURPOSE.
!     --------
!     READ LMESSPASS AND LFDB FROM USER INPUT SKIPPING ALL THE OTHER
!     INPUT PARAMETER WHICH WILL BE READ IN WITH USERIN.

!**   INTERFACE.
!     ----------

!      *CALL* *MPUSERIN (IU05, LMESSPASS, LFDB) 
!        *IU05*      -  INPUT UNIT USER INPUT FILE.
!        *LMESSPASS* - LOGICAL USED TO SET PROGRAM TO COMPUTER
!                      ARCHITECHTURE
!        *LFDB*      - LOGICAL USED TO TOGGLE ACCESS TO FDB

!     METHOD.
!     -------

!     SEE USERIN

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : FFLAG    ,PFLAG    ,GFLAG    ,NFLAG
      USE YOWTEST  , ONLY : ITEST    ,ITESTB

! ----------------------------------------------------------------------

      CHARACTER *2 YCLASS
      CHARACTER *3 USERID, RUNID
      CHARACTER *4 YEXPVER
      CHARACTER*12 CDATEF, CDATER, CDATES, CBPLTDT, CEPLTDT
      CHARACTER*60 PATH
      CHARACTER*70 CLHEADER
      CHARACTER*80 CPATH
      CHARACTER*256 CFDBSF, CFDB2DSP
      LOGICAL LPBIOIN, LGRIBIN, LGRIBOUT, LPBIOOUT, LMESSPASS, LFDB,
     &        LWCOU, LFDBIOOUT, LNOCDIN, LALTCOR, L4VTYPE, LFRSTFLD,
     &        LALTAS, LSARAS,LALTLRGR


      CHARACTER*1 CLMTSU(4), CLOTSU(7)

      NAMELIST /NALINE/ CLHEADER,
     &   CBPLTDT, CEPLTDT, CDATEF,
     &   IDELPRO, IDELT, IDELWO, IDELWI, CLMTSU, IDELALT,
     &   IDELINT, IDELINS, IDELAGI, IDELSPT, IDELSPS, IDELAGS, IDELRES,
     &   CLOTSU, CDATER, CDATES,
     &   PFLAG,  FFLAG,  GFLAG, NFLAG,
     &   LFDB, LPBIOIN, LGRIBIN, LPBIOOUT, LGRIBOUT, LFDBIOOUT,
     &   ICASE, ISHALLO, ITEST, ITESTB, IREST, IASSI, IFORCA,
     &   IBOUNC, IBOUNF,
     &   USERID, RUNID,  PATH, YCLASS, YEXPVER, CPATH,
     &   NENSFNB, NTOTENS,
     &   LMESSPASS, LWCOU, NTASKS, NSIZE, CFDBSF, CFDB2DSP, LNOCDIN,
     &   LALTCOR, L4VTYPE, LFRSTFLD, LALTAS, LSARAS, XKAPPA2, NCALL,
     $   SWAMPWIND,LALTLRGR,HSCOEFCOR,HSCONSCOR

! ----------------------------------------------------------------------

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      LFDB =.FALSE.
      LMESSPASS =.TRUE.
      LWCOU=.FALSE.

! ----------------------------------------------------------------------

!*    1. READ NAMELIST NALINE AND PRINT SETTINGS.
!        ----------------------------------------

      REWIND(IU05)
      READ (IU05, NALINE)
      REWIND(IU05)

      WRITE(6,*) '==============================================='
      WRITE(6,*) '*** MPUSERIN has read the following settings'
      WRITE(6,*) '*** LMESSPASS = ',LMESSPASS
      WRITE(6,*) '*** LFDB = ',LFDB
      WRITE(6,*) '*** LWCOU= ',LWCOU
      WRITE(6,*) '==============================================='

      RETURN
      END SUBROUTINE MPUSERIN
