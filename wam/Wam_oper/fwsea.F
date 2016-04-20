      SUBROUTINE FWSEA (F, EMN, EWFG, FMWFG, ETOI, USMO, THMO, IJS, IJL)

! ----------------------------------------------------------------------

!**** *FWSEA* - ANALYSE THE MODEL SPECTRA TO PROVIDE THE WIND-SEA PART  
!               OF THE SPECTRUM, ITS ENERGY, ITS MEAN FREQUENCY.        

!      P.LIONELLO     ECMWF    FEBRUARY 1989                       
!        (MODIFICATION OF A CODE BY P.JANSSEN                           
!         AND P.LIONELLO - SUMMER '87 - )                               
!      J.BIDLOT       ECMWF    JANUARY 1997 : CORRECTION FOR INITIAL
!                                             WAVE DIRECTION NOT BEING 0

!     J.BIDLOT       ECMWF    MARCH 2001 : SIMPLIFY THE SEACH BY USING
!                                          THE SAME CRITERIA AS IN
!                                          SEPWISW

!      METHOD.                                                          
!      -------                                                          

!      A SPECTRAL COMPONENT IS WINDSEA IF
!      1.2 * 28 * U*/C COS(THETA-THMO) > 1

!      WIND SEA ENERGY AND MEAN FREQUENCY ARE COMPUTED BY THE EXTERNAL  
!      WSMFEN                                                           

!     F        :     MODEL SPECTRUM (FIRST GUESS)                       
!     EWFG     :     FIRST GUESS WIND-SEA ENERGY                        
!     FMWFG    :     WINDSEA MEAN FREQUENCY (FIRST GUESS)               
!     ETOI     :     TOTAL ENERGY FROM O.I. (ANALYSIS)                  
!     USMO     :     USTAR (FIRST GUESS)                                
!     THMO     :     WIND DIRECTION                                     

! ----------------------------------------------------------------------

      USE YOWFRED  , ONLY : TH       ,C      ,FRIC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJS,IJL,IJ,K,M,NWP,IG

      REAL :: THRESHOLD,CHECKTA
      REAL :: ZHOOK_HANDLE
      REAL :: CM(NFRE)
      REAL, DIMENSION(IJS:IJL) :: EMN
      REAL, DIMENSION(IJS:IJL) :: EWFG, FMWFG, ETOI, THMO, USMO
      REAL, DIMENSION(NANG,NFRE) :: FSEA
      REAL :: F(IJS:IJL,NANG,NFRE)

      LOGICAL :: LLFOUND

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FWSEA',0,ZHOOK_HANDLE)
#endif

      IG=1

!*   1. PARAMETERS AND INITIALIZATION.                                  
!       ------------------------------                                  

      NWP=0                                                             
      DO IJ=IJS,IJL                                                
        EWFG(IJ)  = -999.                                              
        FMWFG(IJ) = -9999.                                             
      ENDDO


!*    2. LOOP ON THE POINTS INSIDE THE BLOCK.                           
!        ------------------------------------                           

      DO M=1,NFRE
        CM(M) = FRIC/C(M)
      ENDDO

      DO IJ=IJS,IJL
        LLFOUND=.FALSE.
!       SKIP THE LAND POINTS.
        IF (EMN(IJ).GT.0.01 .AND. ETOI(IJ).GT.0.01) THEN
          DO K=1,NANG
            THRESHOLD = 1.2*USMO(IJ)*COS(TH(K)-THMO(IJ))
            DO M=1,NFRE
              CHECKTA = CM(M)*THRESHOLD
              IF (CHECKTA.LT.1.) THEN
                FSEA(K,M) = 0.
              ELSE
                FSEA(K,M) = F(IJ,K,M)
                LLFOUND = .TRUE.
              ENDIF
            ENDDO
          ENDDO

!         COMPUTATION OF THE WINDSEA ENERGY WINDSEA MEAN FREQUENCY.
          IF(LLFOUND) THEN
            CALL WSMFEN (FSEA,EWFG(IJ),FMWFG(IJ),USMO(IJ),DEPTH(IJ,IG))
            NWP = NWP + 1
          ENDIF
        ENDIF
      ENDDO

      IF (ITEST.GT.1)  WRITE(IU06,*)
     &   '          SUB. FWSEA: WIND SEA FOUND AT ',NWP, ' POINTS'

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('FWSEA',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE FWSEA
