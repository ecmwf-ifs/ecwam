      SUBROUTINE MBLOCK (BATHY, KA, KE, IPP)

! ----------------------------------------------------------------------

!**** *MBLOCK* - ROUTINE TO ARRANGE WAMODEL GRID FOR ONE BLOCK.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       *MBLOCK* ARRANGES WAMODEL GRID FOR A BLOCK AND
!                COMPUTES VARIOUS MODEL CONSTANTS

!**   INTERFACE.
!     ----------

!       *CALL* *MBLOCK (BATHY, KA, KE, IPP)*
!          *BATHY*   - BATHYMETRY DATA.
!          *KA*      - NUMBER OF FIRST LAT IN BLOCK.
!          *KE*      - NUMBER OF LAST LAT IN BLOCK.
!          *IPP*     - NUMBER OF SEA POINTS PER LAT.

!     METHOD.
!     -------

!       THE LAND POINTS ARE REMOVED. ALL MODEL CONSTANTS WHICH ARE
!       GRID DEPENDENT ARE YOWPUTED AND STORED IN THE YOWMON BLOCKS.

!     EXTERNALS.
!     ----------

!       *ABORT1*     - TERMINATES PROCESSING.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NGX      ,NGY      ,NBLO     ,NIBLO
      USE YOWGRID  , ONLY : NLONRGG  ,IGL      ,IJS      ,IJL2     ,
     &            IJLS     ,IJL      ,IJLT
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,NY       ,AMOSOP   ,
     &            XDELLA
      USE YOWSHAL  , ONLY : DEPTH
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      INTEGER, ALLOCATABLE :: IDUM(:), JDUM(:,:)
      REAL, ALLOCATABLE :: ZDUM(:,:)

      INTEGER :: IPP(NGY)
      REAL :: BATHY(NGX, NGY)

! ----------------------------------------------------------------------

!*    1. UPDATE BLOCK NUMBER AND INITIALIZES ARRAYS.
!        -------------------------------------------

      IGL = IGL + 1

      IF(ALLOCATED(IJS)) THEN
        ALLOCATE(IDUM(IGL-1))
        DO IC=1,IGL-1
          IDUM(IC)=IJS(IC)
        ENDDO
        DEALLOCATE(IJS)
        ALLOCATE(IJS(IGL))
        DO IC=1,IGL-1
          IJS(IC)=IDUM(IC)
        ENDDO
        DEALLOCATE(IDUM)
      ELSE
        ALLOCATE(IJS(IGL))
      ENDIF

      IF(ALLOCATED(IJL2)) THEN
        ALLOCATE(IDUM(IGL-1))
        DO IC=1,IGL-1
          IDUM(IC)=IJL2(IC)
        ENDDO
        DEALLOCATE(IJL2)
        ALLOCATE(IJL2(IGL))
        DO IC=1,IGL-1
          IJL2(IC)=IDUM(IC)
        ENDDO
        DEALLOCATE(IDUM)
      ELSE
        ALLOCATE(IJL2(IGL))
      ENDIF

      IF(ALLOCATED(IJLS)) THEN
        ALLOCATE(IDUM(IGL-1))
        DO IC=1,IGL-1
          IDUM(IC)=IJLS(IC)
        ENDDO
        DEALLOCATE(IJLS)
        ALLOCATE(IJLS(IGL))
        DO IC=1,IGL-1
          IJLS(IC)=IDUM(IC)
        ENDDO
        DEALLOCATE(IDUM)
      ELSE
        ALLOCATE(IJLS(IGL))
      ENDIF

      IF(ALLOCATED(IJL)) THEN
        ALLOCATE(IDUM(IGL-1))
        DO IC=1,IGL-1
          IDUM(IC)=IJL(IC)
        ENDDO
        DEALLOCATE(IJL)
        ALLOCATE(IJL(IGL))
        DO IC=1,IGL-1
          IJL(IC)=IDUM(IC)
        ENDDO
        DEALLOCATE(IDUM)
      ELSE
        ALLOCATE(IJL(IGL))
      ENDIF

      IF(ALLOCATED(IJLT)) THEN
        ALLOCATE(IDUM(IGL-1))
        DO IC=1,IGL-1
          IDUM(IC)=IJLT(IC)
        ENDDO
        DEALLOCATE(IJLT)
        ALLOCATE(IJLT(IGL))
        DO IC=1,IGL-1
          IJLT(IC)=IDUM(IC)
        ENDDO
        DEALLOCATE(IDUM)
      ELSE
        ALLOCATE(IJLT(IGL))
      ENDIF

      IF(ALLOCATED(DEPTH)) THEN
        ALLOCATE(ZDUM(NIBLO,IGL-1))
        DO IC=1,IGL-1
          DO IJ=1,NIBLO
            ZDUM(IJ,IC)=DEPTH(IJ,IC)
          ENDDO
        ENDDO
        DEALLOCATE(DEPTH)
        ALLOCATE(DEPTH(NIBLO,IGL))
        DO IC=1,IGL-1
          DO IJ=1,NIBLO
            DEPTH(IJ,IC)=ZDUM(IJ,IC)
          ENDDO
        ENDDO
        DEALLOCATE(ZDUM)
      ELSE
        ALLOCATE(DEPTH(NIBLO,IGL))
      ENDIF

      IF(ALLOCATED(IXLG)) THEN
        ALLOCATE(JDUM(NIBLO,IGL-1))
        DO IC=1,IGL-1
          DO IJ=1,NIBLO
            JDUM(IJ,IC)=IXLG(IJ,IC)
          ENDDO
        ENDDO
        DEALLOCATE(IXLG)
        ALLOCATE(IXLG(NIBLO,IGL))
        DO IC=1,IGL-1
          DO IJ=1,NIBLO
            IXLG(IJ,IC)=JDUM(IJ,IC)
          ENDDO
        ENDDO
        DEALLOCATE(JDUM)
      ELSE
        ALLOCATE(IXLG(NIBLO,IGL))
      ENDIF

      IF(ALLOCATED(KXLT)) THEN
        ALLOCATE(JDUM(NIBLO,IGL-1))
        DO IC=1,IGL-1
          DO IJ=1,NIBLO
            JDUM(IJ,IC)=KXLT(IJ,IC)
          ENDDO
        ENDDO
        DEALLOCATE(KXLT)
        ALLOCATE(KXLT(NIBLO,IGL))
        DO IC=1,IGL-1
          DO IJ=1,NIBLO
            KXLT(IJ,IC)=JDUM(IJ,IC)
          ENDDO
        ENDDO
        DEALLOCATE(JDUM)
      ELSE
        ALLOCATE(KXLT(NIBLO,IGL))
      ENDIF

      IF(NBLO.LE.0) NBLO = IGL

      IF (IGL.GT.NBLO) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
        WRITE (IU06,*) '*        ==========================          *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '* MORE BLOCKS THAN DIMENSION ALLOWS.         *'
        WRITE (IU06,*) '* BLOCK NUMBER IS                 IGL = ', IGL
        WRITE (IU06,*) '* DIMENSION IS                   NBLO = ', NBLO
        WRITE (IU06,*) '* NUMBER OF FIRST LATITUDE IS      KA = ', KA
        WRITE (IU06,*) '* NUMBER OF LAST  LATITUDE IS      KE = ', KE
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '**********************************************'
        CALL ABORT1
      ENDIF

      DO IJ=1,NIBLO
        DEPTH(IJ,IGL) = 0.
        IXLG(IJ,IGL) = 0
        KXLT(IJ,IGL) = 0
      ENDDO

! ----------------------------------------------------------------------

!*    2. THE FIRST AND LAST BLOCK MUST CONTAIN MORE THAN 2
!*       ALL OTHER BLOCKS MORE  THAN 3 LATITUDES.
!        -------------------------------------------------

      IF( KA.EQ.1 .AND. KE.EQ.1 .AND. NY.EQ.1 ) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '* ALLOWS FOR THE 1 GRID POINT MODEL          *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '**********************************************'
      ELSE IF ((KE.EQ.1) .OR. (KA.EQ.NY) .OR.
     &    ((KA.NE.1) .AND. (KE.EQ.NY) .AND. (KE-KA.LT.2))) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
        WRITE (IU06,*) '*        ==========================          *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '* BLOCK LENGTH IS TOO SHORT.                 *'
        WRITE (IU06,*) '* LESS THAN 2 LATITUDES IN FIRST OR LAST, OR *'
        WRITE (IU06,*) '* LESS THAN 3 LATITUDES IN OTHER BLOCKS.     *'
        WRITE (IU06,*) '* BLOCK NUMBER IS               IGL = ', IGL
        WRITE (IU06,*) '* BLOCK LENGTH IS             NIBLO = ', NIBLO
        WRITE (IU06,*) '* NUMBER OF FIRST LATITUDE IS    KA = ', KA
        WRITE (IU06,*) '* NUMBER OF LAST  LATITUDE IS    KE = ', KE
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '**********************************************'
        CALL ABORT1
      ENDIF

! ----------------------------------------------------------------------

!*    3. COMPUTE INDICES OF FIRST, SECOND, BEFORE LAST, AND LAST LAT.
!        -----------------------------------------------------------

      IF (KA.EQ.1) THEN
        IJS (IGL) = 1
        IJL2(IGL) = IPP(1)
      ELSE
        IJS (IGL) = IPP(KA)+1
        IJL2(IGL) = IPP(KA)+IPP(KA+1)
      ENDIF
      IJLT(IGL) = 0
      DO K=KA,KE
        IJLT(IGL) = IJLT(IGL)+IPP(K)
      ENDDO
      IF (KE.EQ.NY) THEN
        IJL (IGL) = IJLT(IGL)
      ELSE
        IJL (IGL) = IJLT(IGL)-IPP(KE)
      ENDIF
      IJLS(IGL) = IJL(IGL)-IPP(KE-1)+1

! ----------------------------------------------------------------------

!*    4. REMOVE LAND POINTS AND STORE COS AND SIN OF LAT.
!        ------------------------------------------------

      IP = 0
      DO K=KA,KE
        DO I=1,NLONRGG(K)
          IF (BATHY(I,K).GT.-990.) THEN
            IP = IP+1
            DEPTH(IP,IGL) = BATHY(I,K)
            IXLG(IP,IGL) = I
            KXLT(IP,IGL) = K
          ENDIF
        ENDDO
      ENDDO
      IF (IP.NE.IJLT(IGL)) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*        FATAL ERROR IN SUB. MBLOCK          *'
        WRITE (IU06,*) '*        ==========================          *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '* TOTAL NUMBER OF SEAPOINTS DO NOT MATCH.    *'
        WRITE (IU06,*) '* BLOCK NUMBER                    IGL = ', IGL
        WRITE (IU06,*) '* NO. OF SEAPOINTS COUNTED         IP = ', IP
        WRITE (IU06,*) '* NO. OF SEAPOINTS EXPECTED IJLT(IGL) = ',
     &   IJLT(IGL)
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*  PROGRAM WILL BE ABORTED                   *'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '**********************************************'
        CALL ABORT1
      ENDIF

! ----------------------------------------------------------------------

!*    5. PRINTER PROTOCOL OF BLOCK.
!        --------------------------

      IF (IGL.EQ.1) THEN
        WRITE (IU06,'(1H0,'' BLOCKING INFORMATION:'')')
        WRITE (IU06,'(1H ,''            LATITUDES   '',
     &   ''   SECOND LAT. INDEX '',
     &   '' SECOND TO LAST LAT  '',
     &   ''   TOTAL'')')
        WRITE (IU06,'(1H ,''  NO     SOUTH     NORTH'',
     &   ''     START       END'',
     &   ''     START       END'',
     &   ''    POINTS'')')
      ENDIF
      WRITE (IU06,'(1X,I4,2F10.2,5I10)')
     &        IGL, AMOSOP+(KA-1)*XDELLA, AMOSOP+(KE-1)*XDELLA,
     &        IJS(IGL), IJL2(IGL), IJLS(IGL), IJL(IGL), IJLT(IGL)

      RETURN
      END SUBROUTINE MBLOCK
