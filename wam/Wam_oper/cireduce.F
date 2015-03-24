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

      USE YOWICE   , ONLY : LICERUN  ,LMASKICE 
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NBLO     ,NFRE
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWSTAT  , ONLY : NWAM_BLKS 

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: IG, IJ, JH, M 
      INTEGER :: JKGLO,KIJS,KIJL,NPROMA

      REAL,DIMENSION(NINF:NSUP,NFRE,NBLO) :: CIWA
      REAL,DIMENSION(NINF:NSUP,NBLO) :: CICOVER,CITHICK
      REAL :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CIREDUCE',0,ZHOOK_HANDLE)
#endif

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
!         DETERMINE THE WAVE ATTENUATION FACTOR
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO=IJS(IG),IJL(IG),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
            CALL CIWAF(KIJS,KIJL,CICOVER(KIJS:KIJL,IG),
     &                 CITHICK(KIJS:KIJL,IG), CIWA(KIJS:KIJL,:,IG))
          ENDDO
!$OMP     END PARALLEL DO

          CALL GSTATS(1493,1)
        ENDIF

      ENDDO

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CIREDUCE',1,ZHOOK_HANDLE)
#endif
      RETURN
      END SUBROUTINE CIREDUCE
