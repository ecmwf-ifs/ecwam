C ----------------------------------------------------------------------
C
      SUBROUTINE MPUSERIN (IU05, LMESSPASS, LFDB, LWCOU)
C
C ----------------------------------------------------------------------
C
C**** *MPUSERIN* - ROUTINE TO READ LMESSPASS AND LFDB FROM NAMELIST
C****              INPUT.
C
C     J. BIDLOT      ECMWF          JULY    1996
C     B. HANSEN      ECMWF          APRIL   1997
C     B. HANSEN      ECMWF          JANUARY 1998 NAMELIST INPUT.
C
C
C*    PURPOSE.
C     --------
C     READ LMESSPASS AND LFDB FROM USER INPUT SKIPPING ALL THE OTHER
C     INPUT PARAMETER WHICH WILL BE READ IN WITH USERIN.
C
C**   INTERFACE.
C     ----------
C
C      *CALL* *MPUSERIN (IU05, LMESSPASS, LFDB) 
C        *IU05*      -  INPUT UNIT USER INPUT FILE.
C        *LMESSPASS* - LOGICAL USED TO SET PROGRAM TO COMPUTER
C                      ARCHITECHTURE
C        *LFDB*      - LOGICAL USED TO TOGGLE ACCESS TO FDB
C
C     METHOD.
C     -------
C
C     SEE USERIN
C
C     EXTERNALS.
C     ----------
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
      CHARACTER *2 YCLASS
      CHARACTER *3 USERID, RUNID
      CHARACTER *4 YEXPVER
      CHARACTER*12 CDATEA, CDATEE, CDATEF, CDTPRO, CDTRES, CDATER,
     1             CDATES, CDTINTT, CDTINTS, CDTSPT, CDTSPS,
     2             CBPLTDT, CEPLTDT
      CHARACTER*60 PATH
      CHARACTER*70 CLHEADER
      CHARACTER*80 CPATH
      CHARACTER*256 CFDBSF, CFDB2DSP
      PARAMETER (JPPFLAG=21)
      LOGICAL FFLAG(JPPFLAG), PFLAG(JPPFLAG), GFLAG(JPPFLAG),
     1        LPBIOIN, LGRIBIN, LPBIOOUT, LMESSPASS, LFDB, LWCOU,
     2        LFDBIOIN, LFDBIOOUT, LNOCDIN

C
      CHARACTER*1 CLMTSU(4), CLOTSU(7)
      NAMELIST /NALINE/ CLHEADER,
     1   CBPLTDT, CEPLTDT, CDATEF,
     2   IDELPRO, IDELT, IDELWO, IDELWI, CLMTSU, IDELALT,
     3   IDELINT, IDELINS, IDELAGI, IDELSPT, IDELSPS, IDELAGS, IDELRES,
     4   CLOTSU, CDATER, CDATES,
     5   PFLAG,  FFLAG,  GFLAG,
     6   LFDB,  LPBIOIN, LGRIBIN, LPBIOOUT, LFDBIOIN, LFDBIOOUT,
     7   ICASE, ISHALLO, IREFRA, ITEST, ITESTB, IREST, IASSI, IFORCA,
     8   IBOUNC, IBOUNF,
     9   USERID, RUNID,  PATH,   YCLASS, YEXPVER, CPATH,
     x   NENSFNB, NTOTENS,
     1   LMESSPASS, LWCOU, NTASKS, NSIZE, CFDBSF, CFDB2DSP, LNOCDIN
C
C ----------------------------------------------------------------------
C
C*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
C        ---------------------------------------------
C
      LFDB =.FALSE.
      LMESSPASS =.TRUE.
      LWCOU=.FALSE.

C ----------------------------------------------------------------------
C
C*    1. READ NAMELIST NALINE AND PRINT SETTINGS.
C        ----------------------------------------
C
      REWIND(IU05)
      READ (IU05, NALINE)
      REWIND(IU05)
C
      WRITE(6,*) '==============================================='
      WRITE(6,*) '*** MPUSERIN has read the following settings'
      WRITE(6,*) '*** LMESSPASS = ',LMESSPASS
      WRITE(6,*) '*** LFDB = ',LFDB
      WRITE(6,*) '*** LWCOU= ',LWCOU
      WRITE(6,*) '==============================================='

      RETURN
      END
