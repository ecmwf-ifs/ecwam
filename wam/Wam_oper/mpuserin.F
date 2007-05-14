! ----------------------------------------------------------------------

      SUBROUTINE MPUSERIN (IU05, LMESSPASS, LFDB, LWCOU, LL1D)

! ----------------------------------------------------------------------

!**** *MPUSERIN* - ROUTINE TO READ LMESSPASS AND LFDB FROM NAMELIST
!****              INPUT.

!     J. BIDLOT      ECMWF          JULY    1996
!     B. HANSEN      ECMWF          APRIL   1997
!     B. HANSEN      ECMWF          JANUARY 1998 NAMELIST INPUT.
!     S. ABDALLA     ECMWF          OCTOBER 2001 LGUST & LADEN ADDED

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

      USE YOWALTAS , ONLY : NUMALT   ,IBUFRSAT  ,ALTSDTHRSH,ALTBGTHRSH,
     &            HSALTCUT, LALTGRDOUT,
     &            XKAPPA2  ,HSCOEFCOR,HSCONSCOR ,LALTCOR   ,LALTLRGR
      USE YOWCOUT  , ONLY : FFLAG    ,PFLAG    ,GFLAG    ,NFLAG
      USE YOWTEST  , ONLY : ITEST    ,ITESTB
      USE YOWCOUP  , ONLY : LWCOUNORMS

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
     &        LWCOU, LFDBIOOUT, LNOCDIN, L4VTYPE, LFRSTFLD,
     &        LALTAS, LSARAS, LSARINV, LTURN90,
     &        LGUST, LADEN  ,LL1D, LGRHDIFS ,LNEWLVTP


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
     &   NENSFNB, NTOTENS, NSYSNB, NMETNB,
     &   LMESSPASS, LWCOU, NTASKS, NSIZE, CFDBSF, CFDB2DSP, LNOCDIN,
     &   LALTCOR, L4VTYPE, LFRSTFLD, LALTAS, LSARAS, LSARINV, XKAPPA2,
     &   IBUFRSAT,SWAMPWIND, SWAMPWIND2, DTNEWWIND, LTURN90,
     &   LALTLRGR, HSCOEFCOR, HSCONSCOR,ALTSDTHRSH,ALTBGTHRSH,HSALTCUT,
     &   ISTREAM, NLOCGRB, IREFDATE,
     &   NCONSENSUS, NDWD, NMFR, NNCEP, NUKM,
     &   LGUST, LADEN, LALTGRDOUT,
     &   NWAM_BLKS, LL1D, LWCOUNORMS, LGRHDIFS ,LNEWLVTP

! ----------------------------------------------------------------------

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      LFDB =.FALSE.
      LMESSPASS =.TRUE.
      LWCOU=.FALSE.
      LL1D=.TRUE.
      LWCOUNORMS = .FALSE.

! ----------------------------------------------------------------------

!*    1. READ NAMELIST NALINE AND PRINT SETTINGS.
!        ----------------------------------------

      REWIND(IU05)
      READ (IU05, NALINE)
      CLOSE(IU05)


      WRITE(*,*) '==============================================='
      WRITE(*,*) '*** MPUSERIN has read the following settings'
      WRITE(*,*) '*** LMESSPASS = ',LMESSPASS
      WRITE(*,*) '*** LFDB = ',LFDB
      WRITE(*,*) '*** LWCOU= ',LWCOU
      WRITE(*,*) '*** LWCOUNORMS= ',LWCOUNORMS
      WRITE(*,*) '==============================================='

      RETURN
      END SUBROUTINE MPUSERIN
