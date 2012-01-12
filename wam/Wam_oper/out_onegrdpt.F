
      SUBROUTINE OUT_ONEGRDPT(IU06)
!
!--------------------------------------------------------------------
!
!*****OUT_ONEGRDPT** OUTPUT OF ENERGY, PEAK FREQUENCY, DIMENSIONLESS
!                    ENERGY, ROUGHNESS, FETCH  FOR ONE GRID POINT.
!
!     P.JANSSEN JUNE 2005
!
!     PURPOSE
!     -------
!             TO OUTPUT INFORMATIONS AT A FIXED GRID POINT.
!
!     INTERFACE
!     ---------
!             *CALL* *OUT_ONEGRDPT*
!
!
!     METHOD
!     ------
!             NONE
!
!     EXTERNALS
!     ---------
!             NONE
!
!     REFERENCES
!     ----------
!             NONE
!
!--------------------------------------------------------------------
!
      USE YOWCOUP  , ONLY : LWCOU    ,ALPHA    ,XKAPPA   ,XNLEV
      USE YOWCOUT  , ONLY : FFLAG    ,FFLAG20  ,FFLAG21  ,PFLAG    ,
     &            PFLAG20  ,PFLAG21  ,GFLAG    ,GFLAG20  ,GFLAG21  ,
     &            LFDB     ,IPFGTBL
      USE YOWGRID  , ONLY : NLONRGG  ,DELPHI
      USE YOWINTP  , ONLY : WHGTTG   ,WDIRTG   ,WPKFTG   ,WMNFTG   ,
     &            USTARG   ,UDIRG    ,TAUWG    ,CDG      ,SMEANG   ,
     &            U10G     ,WHGTAG   ,CWHTAG   ,RANGAG   ,MWP1G    ,
     &            MWP2G    ,WSPRDG   ,C4G      ,BFG      ,QPG      ,
     &            DEPTHG   ,PHIEPSG  ,PHIAWG   ,TAUOCG
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
      USE YOWSHAL  , ONLY : DEPTH    ,TFAK     ,INDEP
    
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEF   ,CDTPRO   ,CDTINTT  ,
     &            CDTINTS  ,CFDBSF   ,ISHALLO  ,IASSI    ,NENSFNB  ,
     &            NTOTENS  ,MARSTYPE ,YCLASS   ,YEXPVER
      USE FDBSUBS_MOD

! ----------------------------------------------------------------------
      
      CHARACTER TITL*100, CPATH*72
      LOGICAL LPARAM, FRSTIME
      DATA FRSTIME/.TRUE./   
!
! -------------------------------------------------------------------
!
!       THE SWAMP CASE HAS BEEN REDUCED TO A ONE GRID POINT MODEL
!       NECESSARY PARAMETERS ARE OUTPUT HERE
!
      IU_INTP = 999

      IF (FRSTIME) THEN
        CPATH = 'out_intp'
        OPEN (IU_INTP,FILE=CPATH,STATUS='UNKNOWN', FORM='FORMATTED')
        FRSTIME = .FALSE.
      ENDIF
         
      CALL DIFDATE(CDATEA,CDTPRO,ITIME)

      LPARAM = ( ALLOCATED(WHGTTG).AND.ALLOCATED(CDG).AND.
     &           ALLOCATED(U10G).AND.ALLOCATED(MWP1G).AND.
     &           ALLOCATED(TAUWG).AND.ALLOCATED(WPKFTG).AND.
     &           ALLOCATED(WMNFTG) )

      IF (LPARAM) THEN
        I = MAX(1,NGX/2)
        DO J = 1,NGY
          IF (WHGTTG(I,J).NE.ZMISS) THEN
!
!            TAKE RELATION BETWEEN USTAR AND U10G FROM BUILDSTRESS
!
            CD     = CDG(I,J)
            U10    = U10G(I,J)
            HS     = WHGTTG(I,J) 
            USTAR2 = MAX(CD*MAX(U10**2,EPSU10),EPSUS)
            USTAR  = SQRT(USTAR2)
            SQTAUW   = SQRT(TAUWG(I,J)*USTAR2)

            TSTAR = G*ITIME/USTAR
            DSTAR = G*DEPTH(I,J)/(USTAR**2)
            E     = HS**2/16.
            ESTAR = G**2*E/(USTAR**4)

            FMSTAR=USTAR/(MWP1G(I,J)*G)

            TSTAR_0 = 4.26*10.**5
            XP      = 1.5
            BETA_K  = 100.
            ALPHA_K = 2.25/10000.
              
            E_LIM = BETA_K**2/16.
            E_STAR_OBS = E_LIM/(1.+TSTAR_0/TSTAR)**XP 

            FP     = WPKFTG(I,J)
            XNUSTAR = USTAR*FP/G
            XNU_OBS = (ALPHA_K/E_STAR_OBS)**(1./3.)

            CDSQRTINV = MIN(1./SQRT(CD),100.0)
            Z0        = XNLEV(1)*EXP(-XKAPPA*CDSQRTINV)
            BETA      = MAX(G*Z0/USTAR2,ALPHA)

            DFETCH = (NGY-J+1)*DELPHI
            FETCHSTAR = G*DFETCH/USTAR2

            T10 = G*ITIME/U10
            E10 = G**2*E/U10**4
            FP10 = U10*FP/G
            FETCH10 = G*DFETCH/U10**2

            T_0 = 1.93*10.**4
            XP      = 1.5
            BETA_K  = 0.22
            E_LIM = BETA_K**2/16.
            E_OBS = E_LIM/(1.+T_0/T10)**XP

            WRITE(IU_INTP,60) NGY-J+1,DEPTH(1,1),ITIME/3600.,
     &                       LOG10(TSTAR),LOG10(FETCHSTAR),HS,
     &                       FP,LOG10(ESTAR),LOG10(E_STAR_OBS),
     &                       LOG10(XNUSTAR),LOG10(XNU_OBS),U10,
     &                       USTAR,CD,BETA,T10,FETCH10,E10,
     &                       TSTAR, FMSTAR, DSTAR, ESTAR

            WRITE(IU06,61) NGY-J+1,DEPTH(1,1),ITIME/3600.,
     &                     TSTAR,FETCHSTAR,HS,FP,ESTAR,
     &                     E_STAR_OBS,XNUSTAR,XNU_OBS,U10,
     &                     USTAR,CD,BETA,T10,FETCH10,E10,
     &                     TSTAR, FMSTAR, DSTAR, ESTAR

          ENDIF
        ENDDO
        CALL FLUSH(IU_INTP)
      ENDIF
   60 FORMAT(I4,F7.2,F7.2,2E10.3,2F8.3,4E12.3,F6.1,F7.3,
     V       2F12.5,2E10.3,7(1x,F14.5))
   61 FORMAT(I4,F7.2,F7.2,2E10.3,2F8.3,2F10.3,2F8.5,F6.1,F7.3,
     V       2F12.5,2E10.3,7(1xF14.5))
!!
!     YOUNG-VERHAGEN LIMITS
! 
!     XLOGE_YV = 1.3*XD+LOG10(1.06/1000.)
!     XLOGF_YV = -0.375*XD+LOG10(0.2)
!
!     USTAR-SCALING
!
!     XLOGE_YV = 1.3*XD+LOG10(111.*1.06/1000.)
!     XLOGF_YV = -0.375*XD+LOG10(0.2*0.43)
!
!     CERC
!
!      XLOGE_YV = 1.5*XD+LOG10(0.04)
!      XLOGF_YV = -0.375*XD+LOG10(0.2*0.43)
!
! -------------------------------------------------------------------
!
      RETURN
      END
