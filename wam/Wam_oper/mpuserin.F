! ----------------------------------------------------------------------

      SUBROUTINE MPUSERIN (IU05, LMESSPASS, LFDB, LWCOU, LL1D, LSUBGRID,
     &                     IPROPAGS)

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
!        *IPROPAGS*  - CONTROLS WHICH PROPAGATION SCHEME IS USED.

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
     &            HSALTCUT, LALTGRDOUT, LALTPAS,
     &            XKAPPA2  ,HSCOEFCOR,HSCONSCOR ,LALTCOR   ,LALTLRGR
      USE YOWCOUT  , ONLY : FFLAG    ,PFLAG    ,GFLAG    ,NFLAG
      USE YOWTEST  , ONLY : ITEST    ,ITESTB
      USE YOWCOUP  , ONLY : LWCOUNORMS

! ----------------------------------------------------------------------

      CHARACTER *2 YCLASS
      CHARACTER *3 USERID, RUNID
      CHARACTER *4 YEXPVER
      CHARACTER*14 CDATEF,CDATER,CDATES,CBPLTDT,CEPLTDT,CDATECURA
      CHARACTER*60 PATH
      CHARACTER*70 CLHEADER
      CHARACTER*80 CPATH
      CHARACTER*256 CFDBSF, CFDB2DSP
      LOGICAL LPBIOIN, LGRIBIN, LGRIBOUT, LPBIOOUT, LMESSPASS, LFDB,
     &        LWCOU, LFDBIOOUT, LNOCDIN, L4VTYPE, LFRSTFLD,
     &        LALTAS, LSARAS, LSARINV, LTURN90,
     &        LGUST, LADEN  ,LL1D, LSUBGRID, LGRHDIFS ,LNEWLVTP, LICERUN
      LOGICAL LBIWBK 
      LOGICAL LLSOURCE
      LOGICAL LNSESTART

      CHARACTER*1 CLMTSU(4), CLOTSU(7)

      NAMELIST /NALINE/ CLHEADER,
     &   CBPLTDT, CEPLTDT, CDATEF,
     &   IDELPRO, IDELT, IDELWO, IDELWI, CLMTSU, IDELALT,
     &   IDELINT, IDELINS, IDELAGI, IDELSPT, IDELSPS, IDELAGS, IDELRES,
     &   IDELCUR, CDATECURA,
     &   CLOTSU, CDATER, CDATES,
     &   PFLAG,  FFLAG,  GFLAG, NFLAG,
     &   LFDB, LPBIOIN, LGRIBIN, LPBIOOUT, LGRIBOUT, LFDBIOOUT,
     &   ICASE, ISHALLO, ITEST, ITESTB, IREST, IASSI, IFORCA,
     &   IPROPAGS,
     &   IREFRA  ,
     &   ISNONLIN,
     &   LBIWBK  ,
     &   IBOUNC, IBOUNF,
     &   USERID, RUNID,  PATH, YCLASS, YEXPVER, CPATH,
     &   NENSFNB, NTOTENS, NSYSNB, NMETNB,
     &   LMESSPASS, LWCOU, NTASKS, NSIZE, CFDBSF, CFDB2DSP, LNOCDIN,
     &   LALTCOR, L4VTYPE, LFRSTFLD, LALTAS, LSARAS, LSARINV, XKAPPA2,
     &   IBUFRSAT,SWAMPWIND, SWAMPWIND2, DTNEWWIND, LTURN90,
     &   LALTLRGR, HSCOEFCOR, HSCONSCOR,ALTSDTHRSH,ALTBGTHRSH,HSALTCUT,
     &   ISTREAM, NLOCGRB, IREFDATE,
     &   NCONSENSUS, NDWD, NMFR, NNCEP, NUKM,
     &   LGUST, LADEN, LALTGRDOUT, LSUBGRID, LALTPAS,
     &   LLSOURCE,
     &   LNSESTART,
     &   NWAM_BLKS, LL1D, LWCOUNORMS, LGRHDIFS ,LNEWLVTP, LICERUN

! ----------------------------------------------------------------------

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      LFDB =.FALSE.
      LMESSPASS =.TRUE.
      LWCOU=.FALSE.
      LL1D=.TRUE.
      LSUBGRID=.TRUE.
      LWCOUNORMS = .FALSE.

      IPROPAGS=0

! ----------------------------------------------------------------------

!*    1. READ NAMELIST NALINE AND PRINT SETTINGS.
!        ----------------------------------------

      REWIND(IU05)
      READ (IU05, NALINE)
      CLOSE(IU05)


      WRITE(6,*) '==============================================='
      WRITE(6,*) '*** MPUSERIN has read the following settings'
      WRITE(6,*) '*** LMESSPASS = ',LMESSPASS
      WRITE(6,*) '*** LFDB = ',LFDB
      WRITE(6,*) '*** LWCOU= ',LWCOU
      WRITE(6,*) '*** LSUBGRID= ',LSUBGRID
      WRITE(6,*) '*** IPROPAGS= ',IPROPAGS
      WRITE(6,*) '*** LWCOUNORMS= ',LWCOUNORMS
      WRITE(6,*) '==============================================='

      RETURN
      END SUBROUTINE MPUSERIN
