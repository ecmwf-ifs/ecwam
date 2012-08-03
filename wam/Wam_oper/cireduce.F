      SUBROUTINE CIREDUCE (CICOVER,CITHICK,CIWA)

! ----------------------------------------------------------------------

!**** *CIREDUCE* - COMPUTE SEA ICE REDUCTION FACTOR FOR SORUCE TERMS 
!                  AND THE SEA ICE WAVE ATTENUATION FACTORS

!!!! currently also setting parametric sea ice thickness !!!!

!*    PURPOSE.
!     --------

!       CIREDUCE COMPUTES SEA ICE SOURCE TERM REDUCTION FACTOR.

!**   INTERFACE.
!     ----------

!       *CALL* *CIREDUCE (CICOVER,CITHICK,CIWA)

!          *CICOVER*-  SEA ICE COVER.
!          *CITHICK*-  SEA ICE THICKNESS. 
!          *CIWA*-     SEA ICE WAVE ATTENUATION FACTOR. 

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------


! ----------------------------------------------------------------------

      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,CITHRSH   ,CIBLOCK   ,
     &            NICH     ,HICMIN   , DHIC
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL       ,SINPH
      USE YOWMAP   , ONLY : KXLT
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NBLO     ,NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWSTAT  , ONLY : NWAM_BLKS 

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
!     CONSTANTS FROM PARAMETRISATION OF SEA ICE THICKNESS FOLLOWING
!     KRINNER ET AL. 2010, Clim Dyn DOI 10.1007/s00382-009-0587-7 FOR THE ARTIC
!     EXCEPT A CONSTANT IS ADDED (C0) FOR THE ANTARCTICA
!     AND THE MINIMUM SEA ICE COVER IS SET TO 0. BECAUSE WE WILL NOT TREAT AREAS
!     THAT ARE EXTEMSIVELY COVERED WITh SEA ICE (SEE CIBLOCK VALUE) 
      REAL, PARAMETER :: C0=0.2
      REAL, PARAMETER :: C1=0.2
      REAL, PARAMETER :: C2=2.8
      REAL, PARAMETER :: C3=1.0
      REAL, PARAMETER :: C4=2.0
      REAL, PARAMETER :: CICOVERMIN=0.0

      INTEGER :: IG, IJ, JH, M 
      INTEGER :: JKGLO,KIJS,KIJL,NPROMA

      REAL :: HICMAX 
      REAL :: DICMEAN
      REAL,DIMENSION(NINF:NSUP,NFRE,NBLO) :: CIWA
      REAL,DIMENSION(NINF:NSUP,NBLO) :: CICOVER,CITHICK
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('CIREDUCE',0,ZHOOK_HANDLE)

      HICMAX=HICMIN+(NICH-1)*DHIC

      DICMEAN=36.0


      DO IG=1,IGL
        NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1

        IF( .NOT. LICERUN .OR. LMASKICE ) THEN
!         NO REDUCTION, EITHER THERE IS NO SEA ICE INFORMATION OR
!         ALL SEA ICE COVER POINTS WILL BE MASKED
          CALL GSTATS(1493,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC)
!$OMP+    PRIVATE(JKGLO,KIJS,KIJL,IJ,M)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            DO IJ=KIJS,KIJL
              CITHICK(IJ,IG)=0.0
            ENDDO
            DO M=1,NFRE
              DO IJ=KIJS,KIJL
                CIWA(IJ,M,IG)=1.0
              ENDDO
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO

          CALL GSTATS(1493,1)
        ELSE


          CALL GSTATS(1493,0)

!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            DO IJ=KIJS,KIJL
!!!!!! we do not have ice thickness so try a parametrisation based on sea ice cover???
!!!!   using Kinner et al. in the Artic (NH)
!!!!   and a constant C0+Krinner in the Antarctica (SH)
              IF(CICOVER(IJ,IG).GT.0.0) THEN
                JH=KXLT(IJ,IG)
                IF(SINPH(JH).GT.0.) THEN
                  CITHICK(IJ,IG)=MAX((C1+C2*CICOVERMIN**2)*
     &                   (C3+C4*(CICOVER(IJ,IG)-CICOVERMIN)),0.)
                ELSE
                  CITHICK(IJ,IG)=C0+MAX((C1+C2*CICOVERMIN**2)*
     &                   (C3+C4*(CICOVER(IJ,IG)-CICOVERMIN)),0.)
                ENDIF
              ELSE
                CITHICK(IJ,IG)=0.0
              ENDIF
            ENDDO

!           DETERMINE THE WAVE ATTENUATION FACTOR
            CALL CIWAF(KIJS,KIJL,CICOVER(KIJS:KIJL,IG),
     &                 CITHICK(KIJS:KIJL,IG), CIWA(KIJS:KIJL,:,IG))

          ENDDO
!$OMP     END PARALLEL DO

          CALL GSTATS(1493,1)
        ENDIF

      ENDDO

      IF (LHOOK) CALL DR_HOOK('CIREDUCE',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE CIREDUCE
