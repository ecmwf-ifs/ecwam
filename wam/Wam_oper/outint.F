      SUBROUTINE OUTINT (IU06, IU20, IU21, IU30, ITEST, KFDB)

! ----------------------------------------------------------------------

!**** *OUTINT* - OUTPUT OF INTEGRATED FIELDS IS WRITTEN.

!     H. GUNTHER       ECMWF  JUNE     1990
!     A. GUILLAUME     ECMWF  JUNE     1992 HANDLING OF MISSING VALUES,
!                                           AND FDB OUTPUT.
!     A. GUILLAUME     ECMWF  JULY     1992 DIRECTIONS.

!     P.A.E.M. JANSSEN ECMWF  AUGUST   1995

!     J. BIDLOT        ECMWF  APRIL    1996 MESSAGE PASSING.
!     J. BIDLOT        ECMWF  OCTOBER  1997 ADD OUTPUT OF U10G.
!     B. HANSEN        ECMWF  FEBRUARY 1998 INTERFACE TO FDB.
!     P. TOWERS        ECMWF  APRIL    2003 ISETFIELDCOUNTFDB

!*    PURPOSE.
!     --------

!       OUTPUT OF WAVE AND WIND FIELDS IN GRIDDED FORMAT
!       TO DISK ,PRINTER AND FDB.


!**   INTERFACE.
!     ----------
!       *CALL* *OUTINT (IU06, IU20, IU21, IU30, ITEST, KFDB)
!          *IU06*   - PRINTER OUTPUT UNIT.
!          *IU20*   - OUTPUT UNIT OF PARAMETERS OF TOTAL SPECTRUM.
!          *IU21    - OUTPUT UNIT OF SEA AND SWELL PARAMETER.
!          *IU30    - OUTPUT UNIT OF GRIBBED PARAMETERS OF TOTAL
!                     SPECTRUM.
!          *ITEST*  - TESTOUTPUT GIVEN IF ITEST LT 2.
!          *KFDB*   - REFERENCE TO THE FDB DATABASE.

!     EXTERNALS.
!     ----------

!       *OUTPP*     - PRINTER OUTPUT OF AN ARRAY.
!       *GRIBPAC*   -
!       *MPGATHERGRDFLD* - GATHER GRID FIELD ONTO ONE PROCESS

!     METHOD.
!     -------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU    ,ALPHA    ,XKAPPA   ,XNLEV
      USE YOWCOUT  , ONLY : FFLAG    ,FFLAG20  ,FFLAG21  ,PFLAG    ,
     &            PFLAG20  ,PFLAG21  ,GFLAG    ,GFLAG20  ,GFLAG21  ,
     &            LFDB     ,IPFGTBL
      USE YOWGRIBHD, ONLY : NTENCODE ,HOPERI
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,WHGTAG   ,CWHTAG   ,RANGAG   ,
     &            MWP1G    ,MWP2G    ,WSPRDG
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WMNFSG   ,WHGTWG   ,
     &            WDIRWG   ,WMNFWG   ,MWP1SG   ,MWP2SG   ,WSPRDSG  ,
     &            MWP1WG   ,MWP2WG   ,WSPRDWG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY      ,CLDOMAIN
      USE YOWPCONS , ONLY : G        ,DEG      ,ZMISS    ,EGRCRV   ,
     &            AGRCRV   ,BGRCRV   ,AFCRV    ,BFCRV    ,EPSUS    ,
     &            EPSU10
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEF   ,CDTPRO   ,CDTINTT  ,
     &            CDTINTS  ,CFDBSF   ,IASSI    ,NENSFNB  ,NTOTENS  ,
     &            MARSTYPE ,YCLASS   ,YEXPVER

! ----------------------------------------------------------------------

      CHARACTER TITL*100

! ----------------------------------------------------------------------
!
      INTEGER IGLOBAL,ILOCAL      ! FDB field counters
      INTEGER ISETFIELDCOUNTFDB
      INTEGER IERR
      EXTERNAL ISETFIELDCOUNTFDB

      IGLOBAL=0
      ILOCAL=0

      IF(CLDOMAIN.EQ.'s') THEN
!     the swamp case has been reduced to a one grid point model
!     the necessary parameters for the study of the non dimensional
!     growth curves are output
!!!!  MAKE SURE that whgttg(i,j) is the wave wave height.


        call difdate(cdatea,cdtpro,itime)
        if(allocated(whgttg).and.
     &     allocated(cdg).and.
     &     allocated(u10g).and.
     &     allocated(mwp1g).and.
     &     allocated(tauwg).and.
     &     allocated(wpkftg).and.
     &     allocated(wmnftg) ) then
          do j=1,ngy
          do i=1,ngx
            if(whgttg(i,j).ne.zmiss) then
!             take the definiton from the model for the relation
!             between ustar and u10g (see buildstress)
              ustar2=max(cdg(i,j)*max(u10g(i,j)**2,epsu10),epsus)
              ustar=sqrt(ustar2)
              tstar=itime*g/ustar
              estar=(whgttg(i,j)**2/16)*g**2/(ustar**4)
              anfit=EGRCRV*(tstar/(AGRCRV+tstar))**BGRCRV
              fmstar=ustar/(mwp1g(i,j)*g)
              fmstar1=ustar*wmnftg(i,j)/g
              fmfit=AFCRV*fmstar**BFCRV
              cdsqrtinv = min(1./sqrt(cdg(i,j)),100.0)
              z0=xnlev(1)*exp(-xkappa*cdsqrtinv)
              charnock=max(g*z0/ustar2,alpha)
              X=min(tauwg(i,j),0.99999999)
              charnock2=alpha/sqrt(1.-X)
              write(iu06,*) 'grow ',tstar,estar,ustar,whgttg(i,j)**2/16,
     &                      cdg(i,j),u10g(i,j),anfit,fmstar,fmfit,
     &                      wpkftg(i,j),charnock,charnock2
              write(iu06,*) ' '
            endif
          enddo
          enddo
        endif

      ENDIF

      IF (CDTINTT.EQ.CDTPRO) THEN

!*      0.2 CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
!           -----------------------------------------------------------

        IF(FFLAG(2).OR.PFLAG(2).OR.GFLAG(2)) THEN
          IF ((LMESSPASS.AND.(IRANK.EQ.IPFGTBL(2)))
     &     .OR.(.NOT.LMESSPASS)) THEN
            DO J=1,NGY
              DO I=1,NGX
                IF(WDIRTG(I,J).NE.ZMISS)
     &            WDIRTG(I,J) = MOD(DEG * WDIRTG(I,J) +180.,360.)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

        IF(FFLAG(5).OR.PFLAG(5).OR.GFLAG(5)) THEN
          IF ((LMESSPASS.AND.(IRANK.EQ.IPFGTBL(5)))
     &     .OR.(.NOT.LMESSPASS)) THEN
            DO J=1,NGY
              DO I=1,NGX
                IF(UDIRG(I,J).NE.ZMISS)
     &            UDIRG(I,J) = MOD(DEG * UDIRG(I,J) +180.,360.)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

      ENDIF

      IF (CDTINTS.EQ.CDTPRO) THEN

!*      0.4 CONVERT DIRECTIONS TO DEGREES AND METEOROLOGICAL CONVENTION
!           -----------------------------------------------------------
        IF(FFLAG(14).OR.PFLAG(14).OR.GFLAG(14)) THEN
          IF ((LMESSPASS.AND.(IRANK.EQ.IPFGTBL(14)))
     &     .OR.(.NOT.LMESSPASS)) THEN
            DO J=1,NGY
              DO I=1,NGX
              IF(WDIRSG(I,J).NE.ZMISS)
     &          WDIRSG(I,J) = MOD(DEG * WDIRSG(I,J) + 180.,360.)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

        IF(FFLAG(13).OR.PFLAG(13).OR.GFLAG(13)) THEN
          IF ((LMESSPASS.AND.(IRANK.EQ.IPFGTBL(13)))
     &     .OR.(.NOT.LMESSPASS)) THEN
            DO J=1,NGY
              DO I=1,NGX
                IF(WDIRWG(I,J).NE.ZMISS)
     &            WDIRWG(I,J) = MOD(DEG * WDIRWG(I,J) + 180.,360.)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

      ENDIF

!     1. COLLECT INTEGRATED PARAMETERS FOR OUTPUT ON A SPECIFIC PROCESS
!        --------------------------------------------------------------

      IF(LMESSPASS.AND.(FFLAG20.OR.PFLAG20.OR.GFLAG20)) THEN
!       1.1 COLLECT TOTAL SEA PARAMETER ARRAYS TO PROCESS NPROC
!           ---------------------------------------------------
        IRECV=NPROC

        IF (FFLAG(1)) CALL MPGATHERGRDFLD(IPFGTBL(1),IRECV,1,WHGTTG)
        IF (FFLAG(2)) CALL MPGATHERGRDFLD(IPFGTBL(2),IRECV,2,WDIRTG)
        IF (FFLAG(3)) CALL MPGATHERGRDFLD(IPFGTBL(3),IRECV,3,WMNFTG)
        IF (FFLAG(4)) CALL MPGATHERGRDFLD(IPFGTBL(4),IRECV,4,USTARG)
        IF (FFLAG(5)) CALL MPGATHERGRDFLD(IPFGTBL(5),IRECV,5,UDIRG)
        IF (FFLAG(6)) CALL MPGATHERGRDFLD(IPFGTBL(6),IRECV,6,WPKFTG)
        IF (FFLAG(7)) CALL MPGATHERGRDFLD(IPFGTBL(7),IRECV,7,CDG)
        IF (FFLAG(8)) CALL MPGATHERGRDFLD(IPFGTBL(8),IRECV,8,TAUWG)
        IF (FFLAG(9)) CALL MPGATHERGRDFLD(IPFGTBL(9),IRECV,9,SMEANG)
        IF (FFLAG(10))CALL MPGATHERGRDFLD(IPFGTBL(10),IRECV,10,U10G)
        IF (FFLAG(25))CALL MPGATHERGRDFLD(IPFGTBL(25),IRECV,25,MWP1G)
        IF (FFLAG(26))CALL MPGATHERGRDFLD(IPFGTBL(26),IRECV,26,MWP2G)
        IF (FFLAG(27))CALL MPGATHERGRDFLD(IPFGTBL(27),IRECV,27,WSPRDG)

      ENDIF

      IF(LMESSPASS) THEN
!       1.2 COLLECT INTEGRATED PARAMETER ARRAYS TO PROCESS 1
!           ONLY FOR PRINT OUTPUT OR GRIB IF OUTPUT TO FILE
!           ------------------------------------------------
        IRECV=1

        IF (PFLAG(1).OR.(.NOT.LFDB.AND.GFLAG(1)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(1),IRECV,101,WHGTTG)
        IF (PFLAG(2).OR.(.NOT.LFDB.AND.GFLAG(2)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(2),IRECV,102,WDIRTG)
        IF (PFLAG(3).OR.(.NOT.LFDB.AND.GFLAG(3)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(3),IRECV,103,WMNFTG)
        IF (PFLAG(4).OR.(.NOT.LFDB.AND.GFLAG(4)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(4),IRECV,104,USTARG)
        IF (PFLAG(5).OR.(.NOT.LFDB.AND.GFLAG(5)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(5),IRECV,105,UDIRG)
        IF (PFLAG(6).OR.(.NOT.LFDB.AND.GFLAG(6)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(6),IRECV,106,WPKFTG)
        IF (PFLAG(7).OR.(.NOT.LFDB.AND.GFLAG(7)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(7),IRECV,107,CDG)
        IF (PFLAG(8).OR.(.NOT.LFDB.AND.GFLAG(8)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(8),IRECV,108,TAUWG)
        IF (PFLAG(9).OR.(.NOT.LFDB.AND.GFLAG(9)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(9),IRECV,109,SMEANG)
        IF (PFLAG(10).OR.(.NOT.LFDB.AND.GFLAG(10)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(10),IRECV,110,U10G)
        IF (PFLAG(11).OR.(.NOT.LFDB.AND.GFLAG(11)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(11),IRECV,111,WHGTWG)
        IF (PFLAG(12).OR.(.NOT.LFDB.AND.GFLAG(12)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(12),IRECV,112,WHGTSG)
        IF (PFLAG(13).OR.(.NOT.LFDB.AND.GFLAG(13)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(13),IRECV,113,WDIRWG)
        IF (PFLAG(14).OR.(.NOT.LFDB.AND.GFLAG(14)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(14),IRECV,114,WDIRSG)
        IF (PFLAG(15).OR.(.NOT.LFDB.AND.GFLAG(15)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(15),IRECV,115,WMNFWG)
        IF (PFLAG(16).OR.(.NOT.LFDB.AND.GFLAG(16)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(16),IRECV,116,WMNFSG)
        IF (PFLAG(25).OR.(.NOT.LFDB.AND.GFLAG(25)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(25),IRECV,125,MWP1G)
        IF (PFLAG(26).OR.(.NOT.LFDB.AND.GFLAG(26)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(26),IRECV,126,MWP2G)
        IF (PFLAG(27).OR.(.NOT.LFDB.AND.GFLAG(27)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(27),IRECV,127,WSPRDG)
        IF (PFLAG(28).OR.(.NOT.LFDB.AND.GFLAG(28)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(28),IRECV,128,MWP1WG)
        IF (PFLAG(29).OR.(.NOT.LFDB.AND.GFLAG(29)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(29),IRECV,129,MWP1SG)
        IF (PFLAG(30).OR.(.NOT.LFDB.AND.GFLAG(30)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(30),IRECV,130,MWP2WG)
        IF (PFLAG(31).OR.(.NOT.LFDB.AND.GFLAG(31)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(31),IRECV,131,MWP2SG)
        IF (PFLAG(32).OR.(.NOT.LFDB.AND.GFLAG(32)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(32),IRECV,132,WSPRDWG)
        IF (PFLAG(33).OR.(.NOT.LFDB.AND.GFLAG(33)))
     &      CALL MPGATHERGRDFLD(IPFGTBL(33),IRECV,133,WSPRDSG)

      ENDIF

      IF(LMESSPASS.AND.(FFLAG21.OR.PFLAG21.OR.GFLAG21)) THEN
!       1.3 COLLECT SWELL AND WIND SEA PARAMETER ARRAYS TO PROCESS IRECV
!           ------------------------------------------------------------
        IRECV=MIN(2,NPROC)

        IF (FFLAG(11)) CALL MPGATHERGRDFLD(IPFGTBL(11),IRECV,11,WHGTWG)
        IF (FFLAG(12)) CALL MPGATHERGRDFLD(IPFGTBL(12),IRECV,12,WHGTSG)
        IF (FFLAG(13)) CALL MPGATHERGRDFLD(IPFGTBL(13),IRECV,13,WDIRWG)
        IF (FFLAG(14)) CALL MPGATHERGRDFLD(IPFGTBL(14),IRECV,14,WDIRSG)
        IF (FFLAG(15)) CALL MPGATHERGRDFLD(IPFGTBL(15),IRECV,15,WMNFWG)
        IF (FFLAG(16)) CALL MPGATHERGRDFLD(IPFGTBL(16),IRECV,16,WMNFSG)
        IF (FFLAG(28)) CALL MPGATHERGRDFLD(IPFGTBL(28),IRECV,28,MWP1WG)
        IF (FFLAG(29)) CALL MPGATHERGRDFLD(IPFGTBL(29),IRECV,29,MWP1SG)
        IF (FFLAG(30)) CALL MPGATHERGRDFLD(IPFGTBL(30),IRECV,30,MWP2WG)
        IF (FFLAG(31)) CALL MPGATHERGRDFLD(IPFGTBL(31),IRECV,31,MWP2SG)
        IF (FFLAG(32)) CALL MPGATHERGRDFLD(IPFGTBL(32),IRECV,32,WSPRDWG)
        IF (FFLAG(33)) CALL MPGATHERGRDFLD(IPFGTBL(33),IRECV,33,WSPRDSG)

      ENDIF


!*    2. OUTPUT OF INTEGRATED PARAMETERS OF TOTAL SEA.
!        ---------------------------------------------

      IF (CDTINTT.EQ.CDTPRO) THEN

!*      2.1 WRITE PARAMETERS OF TOTAL SEA AND WIND TO FILE.
!           -----------------------------------------------

        IF(FFLAG20) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.NPROC)).OR.(.NOT.LMESSPASS)) THEN

            DNGX = REAL(NGX)
            DNGY = REAL(NGY)

            WRITE (IU20) CDTPRO, DNGX, DNGY, AMOWEP,AMOSOP,AMOEAP,
     &                   AMONOP, ZDELLO
            WRITE (IU20) (FFLAG(J),J=1,9)
            IF (FFLAG(1)) WRITE (IU20) WHGTTG
            IF (FFLAG(2)) WRITE (IU20) WDIRTG
            IF (FFLAG(3)) WRITE (IU20) WMNFTG
            IF (FFLAG(4)) WRITE (IU20) USTARG
            IF (FFLAG(5)) WRITE (IU20) UDIRG
            IF (FFLAG(6)) WRITE (IU20) WPKFTG
            IF (FFLAG(7)) WRITE (IU20) CDG
            IF (FFLAG(8)) WRITE (IU20) TAUWG
            IF (FFLAG(9)) WRITE (IU20) SMEANG
            IF (FFLAG(10))WRITE (IU20) U10G
            IF (FFLAG(25))WRITE (IU20) MWP1G
            IF (FFLAG(26))WRITE (IU20) MWP2G
            IF (FFLAG(27))WRITE (IU20) WSPRDG
          ENDIF
        ENDIF

!*      2.2 PRINT PARAMETERS OF TOTAL SEA AND WIND.
!           ---------------------------------------

        IF(CLDOMAIN.NE.'s') THEN
        IF (PFLAG(1)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' SIGNIFICANT WAVE HEIGHT ( 0.1 METRES )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     &                 WHGTTG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(2)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WAVE DIRECTION ( DEGREE FROM NORTH TO )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     &                 WDIRTG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(3)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WAVE FREQUENCY ( 0.01 HERTZ )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 WMNFTG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(4)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' FRICTION VELOCITY ( 0.01 METRES/SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 USTARG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(5)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' FRICTION DIRECTION ( DEGREE FROM NORTH TO )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     &                 UDIRG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(6)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' PEAK WAVE FREQUENCY ( 0.01 HERTZ )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 WPKFTG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(7)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' DRAG COEFFICIENT (IN  0.1 *PROMILLE )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10000.,
     &                 CDG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(8)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' NORMALISED WAVE STRESS (IN % )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 TAUWG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(9)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SQUARE SLOPE (IN % )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 SMEANG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(10)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' WIND SPEED ( METRES/SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     &                 U10G, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(22).AND.MARSTYPE.EQ.'an'.AND.IASSI.EQ.1) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' ALTIMETER WAVE HEIGHT ( 0.1 METRES )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     &                 WHGTAG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(23).AND.MARSTYPE.EQ.'an'.AND.IASSI.EQ.1) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' CORRECTED ALTIMETER WAVE HEIGHT (0.1 METRES)'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     &                 CWHTAG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(24).AND.MARSTYPE.EQ.'an'.AND.IASSI.EQ.1) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' ALTIMETER RANGE RELATIVE CORRECTION (1000)'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 0.001,
     &                 RANGAG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(25)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN PERIOD (1) ( SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     &                 MWP1G, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(26)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN PERIOD (2) ( SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     &                 MWP2G, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(27)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN DIRECTIONAL SPREAD (x 0.01)'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 WSPRDG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    3. OUTPUT OF INTEGRATED PARAMETERS OF SEA AND SWELL.
!        -------------------------------------------------

      IF (CDTINTS.EQ.CDTPRO) THEN

!*      3.1 WRITE PARAMETERS OF SWELL AND WIND DIRECTION TO FILE.
!           -----------------------------------------------------

        IF(FFLAG21) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.MIN(2,NPROC)))
     &     .OR.(.NOT.LMESSPASS)) THEN

            DNGX = REAL(NGX)
            DNGY = REAL(NGY)

            WRITE (IU21) CDTPRO, DNGX, DNGY, AMOWEP,AMOSOP,AMOEAP,
     &                   AMONOP, ZDELLO
            WRITE (IU21) (FFLAG(J),J=11,16)
            IF (FFLAG(11)) WRITE (IU21) WHGTWG
            IF (FFLAG(13)) WRITE (IU21) WDIRWG
            IF (FFLAG(15)) WRITE (IU21) WMNFWG
            IF (FFLAG(12)) WRITE (IU21) WHGTSG
            IF (FFLAG(14)) WRITE (IU21) WDIRSG
            IF (FFLAG(16)) WRITE (IU21) WMNFSG
            IF (FFLAG(28)) WRITE (IU21) MWP1WG
            IF (FFLAG(29)) WRITE (IU21) MWP1SG
            IF (FFLAG(30)) WRITE (IU21) MWP2WG
            IF (FFLAG(31)) WRITE (IU21) MWP2SG
            IF (FFLAG(32)) WRITE (IU21) WSPRDWG
            IF (FFLAG(33)) WRITE (IU21) WSPRDSG
          ENDIF
        ENDIF


!*      3.2 PRINT PARAMETERS OF SWELL AND WIND WAVE DIRECTION.
!           --------------------------------------------------
        IF(CLDOMAIN.NE.'s') THEN
        IF (PFLAG(11)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' WIND SEA WAVE HEIGHT ( 0.1 METRES )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     &                 WHGTWG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(13)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WIND SEA DIRECTION ( DEGREE FROM NORTH FROM )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     &                 WDIRWG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(15)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WIND SEA FREQUENCY ( 0.01 HERTZ )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 WMNFWG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(12)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' SWELL WAVE HEIGHT ( 0.1 METRES )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     &                 WHGTSG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(14)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SWELL DIRECTION ( DEGREE FROM NORTH FROM )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     &                 WDIRSG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(16)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SWELL FREQUENCY ( 0.01 HERTZ )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 WMNFSG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(28)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WIND SEA PERIOD (1) ( 0.1 SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     &                 MWP1WG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(29)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SWELL PERIOD (1) ( 0.1 SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     &                 MWP1SG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(30)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WIND SEA PERIOD (2) ( 0.1 SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     &                 MWP2WG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(31)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SWELL PERIOD (2) ( 0.1 SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     &                 MWP2SG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(32)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WIND SEA DIRECTIONAL SPREAD (x 0.01)'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 WSPRDWG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        IF (PFLAG(33)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SWELL DIRECTIONAL SPREAD (x 0.01 )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     &                 WSPRDSG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
        ENDIF
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    4.  PACK INTEGRATED PARAMETERS INTO GRIB FM94.
!         ------------------------------------------

      IF(MARSTYPE.EQ.'4v') THEN
!*      THIS IS A 4V CASE
        READ(CDATEA,'(I4,4I2)')IY,IM,ID,IH,IMN
        CALL DIFDATE(CDATEA,CDTPRO,IFCST)
        IFCST=IFCST/3600
      ELSEIF(CDTPRO.LE.CDATEF) THEN
!*      THIS IS AN ANALYSIS DATE
        READ(CDTPRO,'(I4,4I2)') IY, IM, ID, IH, IMN
        IFCST = 0
      ELSE
!*      THIS IS A  FORECAST DATE
        READ(CDATEF,'(I4,4I2)')IY,IM,ID,IH,IMN
        CALL DIFDATE(CDATEF,CDTPRO,IFCST)
        IFCST=IFCST/3600
      ENDIF

      IF (GFLAG20.AND.CDTINTT.EQ.CDTPRO) THEN

!*      WAVE HEIGHTS.
!       -------------

        IF (GFLAG(1)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(1)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=229
            CALL GRIBPAC(IU06,IU30,NTENCODE,WHGTTG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,HOPERI,
     &       CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      WAVE DIRECTION.
!       ---------------

        IF (GFLAG(2)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(2)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=230
            CALL GRIBPAC(IU06,IU30,NTENCODE,WDIRTG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,HOPERI,
     &       CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      WAVE MEAN PERIOD.
!       -----------------

        IF (GFLAG(3)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(3)))
     &      .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &      .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=232
            DO JY=1,NGY
              DO JX=1,NGX
                IF(WMNFTG(JX,JY).NE.0.)THEN
                  IF(WMNFTG(JX,JY).NE.ZMISS)THEN
                    WMNFTG(JX,JY)=1./WMNFTG(JX,JY)
                  ELSE
                    WMNFTG(JX,JY)=ZMISS
                  ENDIF
                ELSE
                  WMNFTG(JX,JY)=0.
                ENDIF
              ENDDO
            ENDDO
            CALL GRIBPAC(IU06,IU30,NTENCODE,WMNFTG,NGX,NGY,IPARAM,IY,IM,
     &                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      WIND DIRECTION.
!       ---------------

        IF (GFLAG(5)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(5)))
     &      .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &      .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=249
            CALL GRIBPAC(IU06,IU30,NTENCODE,UDIRG,NGX,NGY,IPARAM,IY,IM,
     &                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      WAVE PEAK PERIOD.
!       -----------------

        IF (GFLAG(6)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(6)))
     &      .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &      .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=231
            DO JY=1,NGY
              DO JX=1,NGX
                IF(WPKFTG(JX,JY).NE.0.)THEN
                  IF(WPKFTG(JX,JY).NE.ZMISS)THEN
                    WPKFTG(JX,JY)=1./WPKFTG(JX,JY)
                  ELSE
                    WPKFTG(JX,JY)=ZMISS
                  ENDIF
                ELSE
                  WPKFTG(JX,JY)=0.
                ENDIF
              ENDDO
            ENDDO
            CALL GRIBPAC(IU06,IU30,NTENCODE,WPKFTG,NGX,NGY,IPARAM,IY,IM,
     &                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      DRAG COEFFICIENT
!       ----------------

        IF (GFLAG(7)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(7)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=233
            IZLEV=NINT(XNLEV(1))
            CALL GRIBPAC(IU06,IU30,NTENCODE,CDG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,IZLEV,LFDB,ITEST,MARSTYPE,YCLASS,
     &       YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      MEAN SQUARE SLOPE.
!       ------------------

        IF (GFLAG(9)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(9)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=244
            CALL GRIBPAC(IU06,IU30,NTENCODE,SMEANG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,HOPERI,
     &       CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      WIND SPEED.
!       ----------

        IF (GFLAG(10)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(10)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=245
            IZLEV=NINT(XNLEV(1))
            CALL GRIBPAC(IU06,IU30,NTENCODE,U10G,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,IZLEV,LFDB,ITEST,MARSTYPE,YCLASS,
     &       YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      ALTIMETER WAVE HEIGHT.
!       ---------------------

        IF (GFLAG(22).AND.MARSTYPE.EQ.'an'.AND.IASSI.EQ.1) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(22)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=246
            CALL GRIBPAC(IU06,IU30,NTENCODE,WHGTAG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &       YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      CORRECTED ALTIMETER WAVE HEIGHT
!       -------------------------------

        IF (GFLAG(23).AND.MARSTYPE.EQ.'an'.AND.IASSI.EQ.1) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(23)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=247
            CALL GRIBPAC(IU06,IU30,NTENCODE,CWHTAG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &       YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      ALTIMETER RANGE RELATIVE CORRECTION.
!       -----------------------------------

        IF (GFLAG(24).AND.MARSTYPE.EQ.'an'.AND.IASSI.EQ.1) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(24)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=248
            CALL GRIBPAC(IU06,IU30,NTENCODE,RANGAG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &       YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      MEAN PERIOD 1.
!       -------------

        IF (GFLAG(25)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(25)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=220
            CALL GRIBPAC(IU06,IU30,NTENCODE,MWP1G,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &       YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      MEAN PERIOD 2.
!       -------------

        IF (GFLAG(26)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(26)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=221
            CALL GRIBPAC(IU06,IU30,NTENCODE,MWP2G,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &       YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      MEAN DIRECTIONAL SPREAD.
!       -----------------------

        IF (GFLAG(27)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(27)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=222
            CALL GRIBPAC(IU06,IU30,NTENCODE,WSPRDG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &       YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

      ENDIF


! ----------------------------------------------------------------------

!*    5.  PACK WIND SEA AND SWELL PARAMETERS INTO GRIB FM94.
!         -------------------------------------------------

      IF (GFLAG21.AND.CDTINTS.EQ.CDTPRO) THEN

!*      WIND SEA WAVE HEIGHTS.
!       ---------------------

        IF (GFLAG(11)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(11)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=234
            CALL GRIBPAC(IU06,IU30,NTENCODE,WHGTWG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,HOPERI,
     &       CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      WIND SEA WAVE DIRECTION.
!       ------------------------


        IF (GFLAG(13)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(13)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=235
            CALL GRIBPAC(IU06,IU30,NTENCODE,WDIRWG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,HOPERI,
     &       CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      WIND SEA MEAN PERIOD.
!       --------------------

        IF (GFLAG(15)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(15)))
     &      .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &      .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=236
            DO JY=1,NGY
              DO JX=1,NGX
                IF(WMNFWG(JX,JY).NE.0.)THEN
                  IF(WMNFWG(JX,JY).NE.ZMISS)THEN
                    WMNFWG(JX,JY)=1./WMNFWG(JX,JY)
                  ELSE
                    WMNFWG(JX,JY)=ZMISS
                  ENDIF
                ELSE
                  WMNFWG(JX,JY)=0.
                ENDIF
              ENDDO
            ENDDO
            CALL GRIBPAC(IU06,IU30,NTENCODE,WMNFWG,NGX,NGY,IPARAM,IY,IM,
     &                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      SWELL WAVE HEIGHTS.
!       -------------------

        IF (GFLAG(12)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(12)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=237
            CALL GRIBPAC(IU06,IU30,NTENCODE,WHGTSG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,HOPERI,
     &       CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      SWELL WAVE DIRECTION.
!       --------------------


        IF (GFLAG(14)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(14)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=238
            CALL GRIBPAC(IU06,IU30,NTENCODE,WDIRSG,NGX,NGY,IPARAM,IY,IM,
     &       ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,HOPERI,
     &       CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      SWELL SEA MEAN PERIOD.
!       ---------------------

        IF (GFLAG(16)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(16)))
     &      .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &      .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=239
            DO JY=1,NGY
              DO JX=1,NGX
                IF(WMNFSG(JX,JY).NE.0.)THEN
                  IF(WMNFSG(JX,JY).NE.ZMISS)THEN
                    WMNFSG(JX,JY)=1./WMNFSG(JX,JY)
                  ELSE
                    WMNFSG(JX,JY)=ZMISS
                  ENDIF
                ELSE
                  WMNFSG(JX,JY)=0.
                ENDIF
              ENDDO
            ENDDO
            CALL GRIBPAC(IU06,IU30,NTENCODE,WMNFSG,NGX,NGY,IPARAM,IY,IM,
     &                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      MEAN WIND SEA PERIOD 1.
!       -----------------------

        IF (GFLAG(28)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(28)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=223
            CALL GRIBPAC(IU06,IU30,NTENCODE,MWP1WG,NGX,NGY,IPARAM,IY,IM,
     &                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      MEAN WIND SEA PERIOD 2.
!       -----------------------

        IF (GFLAG(30)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(30)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=224
            CALL GRIBPAC(IU06,IU30,NTENCODE,MWP2WG,NGX,NGY,IPARAM,IY,IM,
     &                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      MEAN WIND SEA DIRECTIONAL SPREAD.
!       ---------------------------------

        IF (GFLAG(32)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(32)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=225
            CALL GRIBPAC(IU06,IU30,NTENCODE,WSPRDWG,NGX,NGY,IPARAM,
     &                   IY,IM,ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

      ENDIF

!*      MEAN SWELL PERIOD 1.
!       --------------------

        IF (GFLAG(29)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(29)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=226
            CALL GRIBPAC(IU06,IU30,NTENCODE,MWP1SG,NGX,NGY,IPARAM,IY,IM,
     &                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      MEAN SWELL PERIOD 2.
!       --------------------

        IF (GFLAG(31)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(31)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=227
            CALL GRIBPAC(IU06,IU30,NTENCODE,MWP2SG,NGX,NGY,IPARAM,IY,IM,
     &                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

!*      MEAN SWELL DIRECTIONAL SPREAD.
!       ------------------------------

        IF (GFLAG(33)) THEN
          iglobal=iglobal+1
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(33)))
     &     .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     &     .OR.(.NOT.LMESSPASS)) THEN
            ilocal=ilocal+1
            IPARAM=228
            CALL GRIBPAC(IU06,IU30,NTENCODE,WSPRDSG,NGX,NGY,IPARAM,
     &                   IY,IM,ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     &                   YEXPVER,HOPERI,CFDBSF,NENSFNB,NTOTENS,KFDB)
          ENDIF
        ENDIF

      IF(LFDB.and.ILOCAL.gt.0) THEN
         IERR = ISETFIELDCOUNTFDB(KFDB,IGLOBAL,ILOCAL)
         IF(IERR.NE.0)THEN
            WRITE(IU06,*) ' ------------------------'
            WRITE(IU06,*) ' ERROR setting fdb field count '
            WRITE(IU06,*) ' in routine OUTINT '
            WRITE(IU06,*) ' FDB ERROR CODE IS ',IERR
            WRITE(IU06,*) ' ------------------------'
            CALL ABORT1
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE OUTINT
