      SUBROUTINE OUT_ONEGRDPT_SP(SPEC,FR,TH,XLON,XLAT,CDTPRO,HS,
     &              THQ,FMEAN,USTAR,THW,U10,IU06,IU25)
!
!--------------------------------------------------------------------
!
!*****OUT_ONEGRDPT_SP** OUTPUT OF A NUMBER OF INTEGRATED PARAMETERS AND
!                       OF THE ONE AND TWO-DIMENSIONAL SPECTRUM AT
!                       ONE GRID POINT.
!
!     P.JANSSEN JUNE 2005
!
!     PURPOSE
!     -------
!             TO OUTPUT INFORMATIONS AT A FIXED GRID POINT.
!
!     INTERFACE
!     ---------
!             *CALL* *OUT_ONEGRDPT_SP*
!
!           *SPEC*   -   TWO-DIMENSIONAL SPECTRA.
!           *FR*     -   NFRE FREQUENCIES.
!           *TH*     -   NANG DIRECTIONS.
!           *XLON*   -   LONGITUDE.
!           *XLAT*   -   LATITUDE.
!           *CDTPRO* -   DATE TIME GROUP
!           *HS*     -   SIGNIFICANT WAVE HEIGHT
!           *TH_WA*  -   MEAN WAVE DIRECTION (DEGREES)
!           *FMEAN*  -   MEAN PERIOD (F**-1 MOMENT)
!           *USTAR*  -   FRICTION VELOCITY
!           *THW*    -   WIND DIRECTION (DEGREES)
!           *U10*    -   WINDSPEED AT 10 M HEIGHT
!           *IU06*   -   FORMATTED WRITE UNIT 
!           *IU25*   -   UNFORMATTED WRITE UNIT
!
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
      USE YOWFRED  , ONLY : DELTH    ,FRATIO
      USE YOWPARAM , ONLY : NANG     ,NFRE     
      USE YOWPCONS , ONLY : G        ,ZPI      ,DEG      ,EPSMIN

! ----------------------------------------------------------------------

      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IU06,IU25

      REAL, INTENT(IN) :: XLON, XLAT, HS
      REAL, INTENT(IN) :: THQ, FMEAN, USTAR, THW, U10
      REAL, DIMENSION(NANG,NFRE), INTENT(IN) :: SPEC
      REAL, DIMENSION(NFRE), INTENT(IN) :: FR
      REAL, DIMENSION(NANG), INTENT(IN) :: TH

      CHARACTER(LEN=12), INTENT(IN) :: CDTPRO


      INTEGER ::  K, M
      INTEGER :: IU_SPEC
      INTEGER :: MMAX
      REAL :: DUM, SUMMAX 
      REAL :: FPK, X1, X2, X3, X4, X5, X6, AF, ASPF, AF2, ASPF2
      REAL, DIMENSION(NFRE) :: SPF

      CHARACTER(LEN=72) :: CPATH
      CHARACTER(LEN=100) :: TITL

      LOGICAL :: LPARAM
      LOGICAL, SAVE :: FRSTIME

      DATA FRSTIME/.TRUE./   
!
! -------------------------------------------------------------------
!
!       THE SWAMP CASE HAS BEEN REDUCED TO A ONE GRID POINT MODEL
!       NECESSARY PARAMETERS ARE OUTPUT HERE
!

      IU_SPEC = 998

      IF (FRSTIME) THEN
        CPATH = 'out_spec'
        OPEN (IU_SPEC,FILE=CPATH,STATUS='UNKNOWN', FORM='FORMATTED')
        FRSTIME = .FALSE.
      ENDIF
         
      WRITE(IU25,60) XLAT,XLON,CDTPRO,U10,THW,USTAR

      DUM=0.
      WRITE(IU25,61) DUM,(TH(K)*DEG,K=1,NANG)

      DO M=1,NFRE
        WRITE(IU25,61) FR(M),(MIN(SPEC(K,M),999.),K=1,NANG)
      ENDDO

      DO M=1,NFRE
        SPF(M) = 0.
        DO K=1,NANG 
          SPF(M) = SPF(M) + SPEC(K,M)*DELTH
        ENDDO
      ENDDO 

      SUMMAX = 0.
      DO M=1,NFRE
        IF(SPF(M).GT.SUMMAX)THEN
          SUMMAX = SPF(M)
          MMAX = M
        END IF
      ENDDO


      IF(SUMMAX.EQ.0.)THEN
         FPK=0.
      ELSEIF(MMAX.EQ.1) THEN
         FPK=FR(MMAX)
      ELSEIF(MMAX.EQ.NFRE) THEN
         FPK=FR(MMAX)
      ELSE
         X1 = FR(MMAX)**2-FR(MMAX+1)**2
         X2 = FR(MMAX)-FR(MMAX+1)
         X3 = SPF(MMAX)-SPF(MMAX+1)
         X4 = FR(MMAX-1)**2-FR(MMAX+1)**2
         X5 = FR(MMAX-1)-FR(MMAX+1)
         X6 = SPF(MMAX-1)-SPF(MMAX+1)
         FPK =0.5*X1*(X3*X4-X1*X6)/(X1*X3*X5-X1*X2*X6)
      END IF 

      DO M=1,NFRE
         AF = FR(M)/FPK
         ASPF = (ZPI*FPK)**3/(G*USTAR)*SPF(M)
         AF2 = USTAR*FR(M)/G
         ASPF2 = ZPI**3*G**3/USTAR**5*SPF(M)
!!         WRITE(IU_SPEC,62) XLON,XLAT,FR(M),SPF(M),AF,ASPF,AF2,ASPF2
         WRITE(IU_SPEC,63) XLON,XLAT,FR(M),SPF(M)
      ENDDO

 60   FORMAT(1X,F7.3,1X,F8.3,1X,A12,3(1X,F8.4))
 61   FORMAT(25(1X,F8.4))
 62   FORMAT(4(1X,F8.4),1X,F10.6,E10.3)
 63   FORMAT(2(1X,F8.4),2(1x,F14.10))
!
! -------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE OUT_ONEGRDPT_SP
