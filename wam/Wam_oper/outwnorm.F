      SUBROUTINE OUTWNORM

! ----------------------------------------------------------------------

!**** *WAMNORM* - PRINTS AVERAGE, MINIMUM AND MAXIMUM VALUES OF
!                 INTEGRATED PARAMETER FIELDS 

!     J. BIDLOT       ECMWF   JULY  1998 

!*    PURPOSE.
!     --------


!**   INTERFACE.
!     ----------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : FFLAG    ,PFLAG    ,GFLAG    ,NFLAG    ,
     &            NFLAGALL ,IPFGTBL
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G
      USE YOWINTS  , ONLY : WHGTSG   ,WDIRSG   ,WHGTWG   ,WDIRWG   ,
     &            WMNFWG
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,KTAG
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWPCONS , ONLY : EPSMIN   ,ZMISS
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IF(NFLAGALL) THEN

        WRITE(IU06,*) '   WAMNORM ON ',CDTPRO
        WRITE(IU06,*) '          AVERAGE            MINIMUM  ',
     &   '           MAXIMUM,   NUMBER OF SEA POINTS' 
        WRITE(IU06,*) ' '

        IF(NFLAG(1).AND.(GFLAG(1).OR.FFLAG(1).OR.PFLAG(1))) THEN
          CALL WAMNORM(WHGTTG,NGX,NGY,ZMISS,EPSMIN,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(1),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE ENERGY'
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(2).AND.(GFLAG(2).OR.FFLAG(2).OR.PFLAG(2))) THEN
          CALL WAMNORM(WDIRTG,NGX,NGY,ZMISS,EPSMIN,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(2),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE MEAN DIRECTION'
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(3).AND.(GFLAG(3).OR.FFLAG(3).OR.PFLAG(3))) THEN
          CALL WAMNORM(WMNFTG,NGX,NGY,ZMISS,EPSMIN,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(3),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE MEAN PERIOD'
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(4).AND.(GFLAG(4).OR.FFLAG(4).OR.PFLAG(4))) THEN
          CALL WAMNORM(USTARG,NGX,NGY,ZMISS,EPSMIN,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(4),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR FRICTION VELOCITY'
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(5).AND.(GFLAG(5).OR.FFLAG(5).OR.PFLAG(5))) THEN
          CALL WAMNORM(UDIRG,NGX,NGY,ZMISS,EPSMIN,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(5),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR MEAN WIND DIRECTION'
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(6).AND.(GFLAG(6).OR.FFLAG(6).OR.PFLAG(6))) THEN
          CALL WAMNORM(WPKFTG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(6),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE PEAK FREQUENCY: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(7).AND.(GFLAG(7).OR.FFLAG(7).OR.PFLAG(7))) THEN
          CALL WAMNORM(CDG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(7),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR DRAG COEFFICIENT: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(8).AND.(GFLAG(8).OR.FFLAG(8).OR.PFLAG(8))) THEN
          CALL WAMNORM(TAUWG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(8),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE INDUCED STRESS: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(9).AND.(GFLAG(9).OR.FFLAG(9).OR.PFLAG(9))) THEN
          CALL WAMNORM(SMEANG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(9),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR MEAN SQUARE SLOPE: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(10).AND.(GFLAG(10).OR.FFLAG(10).OR.PFLAG(10))) THEN
          CALL WAMNORM(U10G,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(10),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WAVE WIND SPEED: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(11).AND.(GFLAG(11).OR.FFLAG(11).OR.PFLAG(11))) THEN
          CALL WAMNORM(WHGTWG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(11),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR WIND SEA WAVE HEIGHTS: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(12).AND.(GFLAG(12).OR.FFLAG(12).OR.PFLAG(12))) THEN
          CALL WAMNORM(WHGTSG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(12),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM FOR SWELL WAVE HEIGHTS: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(13).AND.(GFLAG(13).OR.FFLAG(13).OR.PFLAG(13))) THEN
          CALL WAMNORM(WDIRWG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(13),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM WIND SEA WAVE DIRECTION: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(14).AND.(GFLAG(14).OR.FFLAG(14).OR.PFLAG(14))) THEN
          CALL WAMNORM(WDIRSG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(14),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM SWELL WAVE DIRECTION: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(15).AND.(GFLAG(15).OR.FFLAG(15).OR.PFLAG(15))) THEN
          CALL WAMNORM(WMNFWG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(15),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM WIND SEA MEAN PERIOD: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

        IF(NFLAG(16).AND.(GFLAG(16).OR.FFLAG(16).OR.PFLAG(16))) THEN
          CALL WAMNORM(WMNFWG,NGX,NGY,ZMISS,0.,FAVG,FMIN,FMAX,
     &     NAVG,LMESSPASS,IRANK,IPFGTBL(16),1,KTAG)
          KTAG=KTAG+1
          IF(IRANK.EQ.1) THEN
            WRITE(IU06,*) '   WAMNORM SWELL MEAN PERIOD: '
            WRITE(IU06,*) '  ',FAVG,FMIN,FMAX,NAVG
            WRITE(IU06,*) ' '
          ENDIF
        ENDIF

      ENDIF

      RETURN
      END SUBROUTINE OUTWNORM
