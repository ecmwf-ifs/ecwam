SUBROUTINE MPPEWITHINDIST(BLK2GLO,DISTMAX,INTLMAX,KMINLMAX,KMAXLMAX)

! ----------------------------------------------------------------------

!****  *MPPEWITHINDIST* -  DETERMINES A TABLE WHICH TELLS WHETHER
!                          A PE IS WITHIN A DISTANCE DISTMAX FROM 
!                          THE OTHERS.

!      J. BIDLOT       ECMWF     MAY 2002 

!     PURPOSE.                                                          
!     --------                                                          

!     DETERMINES A TABLE WHICH TELLS WHETHER A PE IS WITHIN A DISTANCE
!     DISTMAX FROM ANOTHER PE.

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *MPPEWITHINDIST(BLK2GLO,DISTMAX,INTLMAX,KMINLMAX,KMAXLMAX)

!     *BLK2GLO*  BLOCK TO GRID TRANSFORMATION
!     *DISTMAX*  REAL     DISTANCE IN RADIAN !!!!
!     *INTLMAX*  INTEGER  TABLE INDICATING WHETHER A PE IS WITHIN 
!                         DISTMAX OF ANOTHER PE.
!     *ONLY MEANINGFUL FOR STRUCTURED GRIDS:
!     *KMINLMAX* INTEGER  SOUTHERN LATITUDE INDEX OF THE MOST SOUTHERN
!                         LATITUDE DISTMAX AWAY FROM EACH PE
!     *KMAXLMAX* INTEGER  NORTHERN LATITUDE INDEX OF THE MOST NORTHERN
!                         LATITUDE DISTMAX AWAY FROM EACH PE

!     METHOD.                                                           
!     -------                                                           

!     SEARCH THE DISTANCE BETWEEN GRID POINTS FROM ONE PE AND THE OTHER
!     ONES UNTIL THE DISTANCE IS LESS THE DISTMAX OR THERE IS NO POINT
!     LEFT.

!     EXTERNALS.                                                        
!     ----------                                                        

!     REFERENCES.                                                       
!     -----------                                                       

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWABORT , ONLY : WAM_ABORT
      USE YOWGRID  , ONLY : SINPH    ,COSPH   ,IJSLOC   ,IJLLOC
      USE YOWMAP   , ONLY : XDELLA   ,ZDELLO  ,NGY
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : LL1D     ,LLUNSTR
      USE YOWPCONS , ONLY : DEG      ,RAD
      USE YOWSPEC,   ONLY : NSTART   ,NEND
#ifdef WAM_HAVE_UNWAM
      USE YOWPD,     ONLY : NODES=>NODES_GLOBAL, RANK
#endif
      USE YOWSPHERE, ONLY : SPHERICAL_COORDINATE_DISTANCE

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!-----------------------------------------------------------------------
      IMPLICIT NONE

      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      INTEGER(KIND=JWIM), DIMENSION(NPROC),INTENT(OUT) :: INTLMAX, KMINLMAX, KMAXLMAX
      REAL(KIND=JWRB), INTENT(IN) :: DISTMAX


      INTEGER(KIND=JWIM) :: IJ, IJJ, IR, JR, KI, KJ, IC, IP, IPP
      INTEGER(KIND=JWIM) :: LMAX

      REAL(KIND=JWRB) :: COSDISTMAX
      REAL(KIND=JWRB) :: COSLON, DIS2 
      REAL(KIND=JWRB)  :: XLONJ
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(IJSLOC:IJLLOC) :: XLONI
      REAL(KIND=JWRB), ALLOCATABLE,DIMENSION(:,:) :: CPH2, SPH2

      REAL(KIND=JWRU)  :: XLONJD, XLATJD, DISTD, DISTMAXD
      REAL(KIND=JWRU), DIMENSION(IJSLOC:IJLLOC) :: XLONID, XLATID

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MPPEWITHINDIST',0,ZHOOK_HANDLE)

      LMAX = INT(DEG*DISTMAX/XDELLA)+1
!     add a small number to DISTMAX as a safety precaution
!     i.e each PE will get slightly more data that they should
!     strictly have.
      COSDISTMAX=COS(DISTMAX+0.1_JWRB)
      DISTMAXD=REAL(DISTMAX+0.1_JWRB,JWRU)

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
!       GENERAL CASE FOR UNSTRUCTURED GRID

        DO JR = 1, NPROC
          INTLMAX(JR)=0
        ENDDO

        IR = IRANK

        DO IP=1,RANK(IR)%NP
          XLONID(IP) = NODES(RANK(IR)%IPLG(IP))%X
          XLATID(IP) = NODES(RANK(IR)%IPLG(IP))%Y
        ENDDO

        DO JR = 1, NPROC
          OUTERU : DO IPP=1,RANK(JR)%NP
                    XLONJD = NODES(RANK(JR)%IPLG(IPP))%X
                    XLATJD = NODES(RANK(JR)%IPLG(IPP))%Y
!!! limit on lat separation ???
                  DO IP=1,RANK(IR)%NP
                      CALL SPHERICAL_COORDINATE_DISTANCE(XLONID(IP),XLONJD,XLATID(IP),XLATJD,DISTD)
                      IF (DISTD <= DISTMAXD) THEN
                        INTLMAX(JR)=1
                        EXIT OUTERU 
                      ENDIF
                  ENDDO
!!!!                ENDIF 
          ENDDO OUTERU 

        ENDDO
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif

      ELSEIF (LL1D) THEN
!       FOR THE SPECIAL CASE OF 1-D DECOMPOSITION ON A STRUCTURED GRID
        IR = IRANK
        IJ=NSTART(IR)
        KMINLMAX(IR) = MAX(1,BLK2GLO%KXLT(IJ)-LMAX)

        IJ=NEND(IR)
        KMAXLMAX(IR) = MIN(NGY,BLK2GLO%KXLT(IJ)+LMAX)

        DO JR = 1, NPROC
          INTLMAX(JR)=0
          IJ=NSTART(JR)
          IF (BLK2GLO%KXLT(IJ) >= KMINLMAX(IR) .AND. BLK2GLO%KXLT(IJ) <= KMAXLMAX(IR) ) THEN 
            INTLMAX(JR)=1
          ENDIF

          IJ=NEND(JR)
          IF (BLK2GLO%KXLT(IJ) >= KMINLMAX(IR) .AND. BLK2GLO%KXLT(IJ) <= KMAXLMAX(IR) ) THEN
            INTLMAX(JR)=1
          ENDIF
        ENDDO

      ELSE
!       FOR THE MORE GENERAL CASE OF 2-D DECOMPOSITION ON A STRUCTURED GRID
        IR = IRANK 
        KMINLMAX(IR) = NGY
        KMAXLMAX(IR) = 1
        DO IJ=NSTART(IR),NEND(IR)
          KMINLMAX(IR) = MIN(KMINLMAX(IR),BLK2GLO%KXLT(IJ)-LMAX)
          KMAXLMAX(IR) = MAX(KMAXLMAX(IR),BLK2GLO%KXLT(IJ)+LMAX)
        ENDDO
        KMINLMAX(IR) = MAX(1,KMINLMAX(IR))
        KMAXLMAX(IR) = MIN(NGY,KMAXLMAX(IR))

        ALLOCATE(CPH2(NGY,NGY))
        ALLOCATE(SPH2(NGY,NGY))
        DO KJ=1,NGY
          DO KI=1,NGY
            CPH2(KI,KJ) = COSPH(KJ)*COSPH(KI)
            SPH2(KI,KJ) = SINPH(KJ)*SINPH(KI)
          ENDDO
        ENDDO

        DO JR = 1, NPROC
          INTLMAX(JR)=0
        ENDDO

        DO IJ=NSTART(IR), NEND(IR)
          XLONI(IJ) = REAL(BLK2GLO%IXLG(IJ)-1)*ZDELLO(BLK2GLO%KXLT(IJ))*RAD
        ENDDO

        DO JR = 1, NPROC
          OUTER : DO IJJ=NSTART(JR),NEND(JR)
              XLONJ = REAL(BLK2GLO%IXLG(IJJ)-1)*ZDELLO(BLK2GLO%KXLT(IJJ))*RAD
              KJ = BLK2GLO%KXLT(IJJ)
              IF (KJ >= KMINLMAX(IR) .AND. KJ <= KMAXLMAX(IR) ) THEN
                DO IJ=NEND(IR), NSTART(IR), -1
                    KI = BLK2GLO%KXLT(IJ)
                    COSLON = COS(XLONJ-XLONI(IJ))
                    DIS2 = COSLON*CPH2(KJ,KI)+SPH2(KJ,KI)

                    IF (DIS2 >= COSDISTMAX) THEN
                      INTLMAX(JR)=1
                      EXIT OUTER 
                    ENDIF
                ENDDO
              ENDIF 
          ENDDO OUTER 
        ENDDO

        DEALLOCATE(CPH2)
        DEALLOCATE(SPH2)

      ENDIF

IF (LHOOK) CALL DR_HOOK('MPPEWITHINDIST',1,ZHOOK_HANDLE)

END SUBROUTINE MPPEWITHINDIST
