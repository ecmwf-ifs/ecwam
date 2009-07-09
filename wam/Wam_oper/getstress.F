      SUBROUTINE GETSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                     ROAIRO, ZIDLOLD, ICEMASK,
     &                     NSTART, NEND, IREAD)
! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997 

!     S. ABDALLA   ECMWF      OCTOBER 2001  AIR DENSITY & Zi/L

!*    PURPOSE.
!     --------
!     READS RESTART WIND AND STRESS FIELDS FROM DISK OR USE GRIB WIND
!     AND DRAG COEFFICIENT TO RECONSTRUCT ALL WIND AND STRESS FIELDS.

!**   INTERFACE.
!     ----------
!     *CALL**GETSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,NSTART,NEND,IREAD)
!     *U10OLD*    WIND SPEED.
!     *THWOLD*    WIND DIRECTION (RADIANS).
!     *USOLD*     FRICTION VELOCITY.
!     *TAUW*      WAVE STRESS.
!     *Z0OLD*     ROUGHNESS LENGTH IN M.
!     *ROAIRO*    AIR DENSITY IN KG/M3.
!     *ZIDLOLD*   Zi/L (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!     *ICEMASK*   SEA ICE MASK
!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK

!     METHOD.
!     -------
!     USE READSTRESS TO READ IN THE FILE DEPENDING ON THE PBIO
!     SOFTWARE OR NOT. 
!     THE READING IS ONLY DONE ON PE 1, THEREFORE
!     THE RELEVANT INFORMATION IS SENT TO THE OTHER PE'S USING 
!     MPDISTRIBSCFLD

!     EXTERNALS.
!     ----------
!     BUILDSTRESS
!     GETENV
!     GRSTNAME
!     MPDISTRIBSCFLD
!     MPL_BARRIER
!     READSTRESS

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : NREAL    ,NINTGR
      USE YOWCURR  , ONLY : U        ,V
      USE YOWGRID  , ONLY : IGL
      USE YOWMESPAS, ONLY : LMESSPASS,LPBIOIN  ,LGRIBIN  ,LWAVEWIND
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,KTAG
      USE YOWPARAM , ONLY : NBLO     ,NIBLO
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEF   ,CDTPRO   ,IREFRA   ,
     &            NWAM_BLKS,LNSESTART
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,LRESTARTED
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE

      CHARACTER*14 ZERO
      CHARACTER*120 FILENAME

! ----------------------------------------------------------------------
      INTEGER,DIMENSION(NPROC) :: NSTART, NEND
      REAL,DIMENSION(NINF:NSUP,NBLO) :: U10OLD,THWOLD,USOLD,Z0OLD,TAUW,
     &                                  ROAIRO, ZIDLOLD 
      INTEGER,DIMENSION(NINF:NSUP,NBLO) :: ICEMASK 
! ----------------------------------------------------------------------
!     define local allocatable arrays
      INTEGER,ALLOCATABLE,DIMENSION(:) :: IBUF
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: RFIELD
      INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: IFIELD 

      INTEGER :: IFLD
! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA
! End Mod for OPENMP
      REAL ZHOOK_HANDLE


      IF (LHOOK) CALL DR_HOOK('GETSTRESS',0,ZHOOK_HANDLE)

      ZERO = ' '

      IF (ITEST > 3) THEN
        WRITE(IU06,*) ' SUB: GETSTRESS. '
        WRITE(IU06,*) ' ABOUT TO READ WIND AND STRESS FILE '
        WRITE(IU06,*) ' CDTPRO    =', CDTPRO
        WRITE(IU06,*) ' LMESSPASS =', LMESSPASS
        WRITE(IU06,*) ' LPBIOIN   =', LPBIOIN
        WRITE(IU06,*) ' LGRIBIN   =', LGRIBIN
        WRITE(IU06,*) ' LRESTARTED=', LRESTARTED
        CALL FLUSH (IU06)
      ENDIF

      ALLOCATE(RFIELD(NIBLO,NBLO,NREAL))
      ALLOCATE(IFIELD(NIBLO,NBLO,NINTGR))

      IF(LMESSPASS) THEN


!       READ RESTART FILE FROM PE IREAD

        CALL GSTATS(1882,0)

        IF(LGRIBIN.AND..NOT.LRESTARTED) THEN
!         CREATES WIND AND STRESS FIELDS FROM GRIB WINDS AND DRAG COEFFICIENT.
          CALL BUILDSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                     ROAIRO,ZIDLOLD,ICEMASK,
     &                     IREAD)
          IF (ITEST.GE.1)
     &     WRITE(IU06,*)
     &     'SUB. GETSTRESS: BUILDSTRESS CALLED'
        ELSE


          IF (IRANK.EQ.IREAD) THEN
            CALL GRSTNAME(CDATEA,CDATEF,'LAW',ICPLEN,CPATH,FILENAME)
            CALL READSTRESS(NREAL, RFIELD, NINTGR, IFIELD,
     &                      FILENAME)
          ENDIF

!         DISTRIBUTE THE DIFFERENT CONTRIBUTIONS TO THE OTHER PE'S

          DO IFLD=1,NREAL
            CALL MPDISTRIBSCFLD(IREAD,KTAG,NSTART,NEND,
     &                          RFIELD(:,:,IFLD))
            KTAG=KTAG+1
          ENDDO

          DO IFLD=1,NINTGR
            CALL MPDISTRIBINTFLD(IREAD,KTAG,NSTART,NEND,
     &                           IFIELD(:,:,IFLD))
            KTAG=KTAG+1
          ENDDO

          IF (ITEST.GE.2)
     &     WRITE(IU06,*)
     &     'SUB. GETSTRESS: RESTART WIND AND STRESS COLLECTED'


!         KEEP CORRESPONDING CONTRIBUTION 
          DO IG=1,NBLO
            NPROMA=(NEND(IRANK)-NSTART(IRANK)+1)/NWAM_BLKS+1
!$OMP       PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
            DO JKGLO=NSTART(IRANK),NEND(IRANK),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,NEND(IRANK))
              DO IJ=KIJS,KIJL
                U10OLD(IJ,IG)=RFIELD(IJ,IG,1)
                THWOLD(IJ,IG)=RFIELD(IJ,IG,2)
                USOLD(IJ,IG)=RFIELD(IJ,IG,3)
                TAUW(IJ,IG)=RFIELD(IJ,IG,4)
                Z0OLD(IJ,IG)=RFIELD(IJ,IG,5)
                ROAIRO(IJ,IG)=RFIELD(IJ,IG,6)
                ZIDLOLD(IJ,IG)=RFIELD(IJ,IG,7)
!               for U and V see below
                ICEMASK(IJ,IG)=IFIELD(IJ,IG,1)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
          ENDDO

          IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
            DO IG=1,NBLO
              NPROMA=(NEND(IRANK)-NSTART(IRANK)+1)/NWAM_BLKS+1
!$OMP         PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
              DO JKGLO=NSTART(IRANK),NEND(IRANK),NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,NEND(IRANK))
                DO IJ=KIJS,KIJL
                  U(IJ,IG)=RFIELD(IJ,IG,8)
                  V(IJ,IG)=RFIELD(IJ,IG,9)
                ENDDO
              ENDDO
!$OMP         END PARALLEL DO
              U(NINF-1,IG)=0.
              V(NINF-1,IG)=0.

!!!           U AND V MUST ALSO BE DEFINED OVER THE HALO
              CALL MPEXCHNG(U(:,IG),1,1)
              CALL MPEXCHNG(V(:,IG),1,1)
            ENDDO

!           RE-COMPUTE REFRACTION TERMS
            CALL PROPDOT

!           SET LOGICAL TO RECOMPUTE THE WEIGHTS IN CTUW.
            LUPDTWGHT=.TRUE.

          ENDIF

        ENDIF
        CALL GSTATS(1882,1)

!       BROADCAST THE 4 DATES FROM RESTART FILE TO THE OTHER PE'S
!       AS WELL AS LWAVEWIND
        IF((IREAD.EQ.0).OR.(NPROC.EQ.1)) THEN
      
        ELSE
          IBUFLENGTH=57+1
          ALLOCATE(IBUF(IBUFLENGTH))
          IF(IRANK.EQ.IREAD) THEN
            IF(LWAVEWIND) THEN
              IBUF(1)=1
            ELSE
              IBUF(1)=0
            ENDIF
            KCOUNT=1
            DO IC=1,14
              KCOUNT=KCOUNT+1
              READ(CDTPRO(IC:IC),'(I1)') IBUF(KCOUNT)
            ENDDO
            DO IC=1,14
              KCOUNT=KCOUNT+1
              READ(CDATEWO(IC:IC),'(I1)') IBUF(KCOUNT)
            ENDDO
            DO IC=1,14
              KCOUNT=KCOUNT+1
              READ(CDAWIFL(IC:IC),'(I1)') IBUF(KCOUNT)
            ENDDO
            DO IC=1,14
              KCOUNT=KCOUNT+1
              READ(CDATEFL(IC:IC),'(I1)') IBUF(KCOUNT)
              IF(KCOUNT.GT.IBUFLENGTH) THEN
                WRITE (IU06,*) ' '
                WRITE (IU06,*) ' **************************************'
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' *   FATAL ERROR IN SUB. GETSTRESS:   *'
                WRITE (IU06,*) ' *   ==============================   *'
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' * IBUFLENGTH TOO SMALL !!!           *'
                WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS   *'
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' **************************************'
                CALL ABORT1
              ENDIF
            ENDDO
          ENDIF

          CALL GSTATS(621,0)
          CALL MPL_BROADCAST(IBUF(1:IBUFLENGTH),KROOT=IREAD,
     &     KTAG=KTAG,CDSTRING='GETSTRESS:')
          CALL GSTATS(621,1)

          KTAG=KTAG+1

          IF(IRANK.NE.IREAD) THEN
            IF(IBUF(1).EQ.1) THEN
              LWAVEWIND=.TRUE. 
            ELSE
              LWAVEWIND=.FALSE. 
            ENDIF
            KCOUNT=1
            DO IC=1,14
              KCOUNT=KCOUNT+1
              WRITE(CDTPRO(IC:IC),'(I1)') IBUF(KCOUNT)
            ENDDO
            DO IC=1,14
              KCOUNT=KCOUNT+1
              WRITE(CDATEWO(IC:IC),'(I1)') IBUF(KCOUNT)
            ENDDO
            DO IC=1,14
              KCOUNT=KCOUNT+1
              WRITE(CDAWIFL(IC:IC),'(I1)') IBUF(KCOUNT)
            ENDDO
            DO IC=1,14
              KCOUNT=KCOUNT+1
              WRITE(CDATEFL(IC:IC),'(I1)') IBUF(KCOUNT)
              IF(KCOUNT.GT.IBUFLENGTH) THEN
                WRITE (IU06,*) ' '
                WRITE (IU06,*) ' **************************************'
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' *   FATAL ERROR IN SUB. GETSTRESS:   *'
                WRITE (IU06,*) ' *   ==============================   *'
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' * IBUFLENGTH TOO SMALL !!!           *'
                WRITE (IU06,*) ' * PROGRAM ABORTS    PROGRAM ABORTS   *'
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' **************************************'
                CALL ABORT1
              ENDIF
            ENDDO
          ENDIF

          DEALLOCATE(IBUF)

        ENDIF

        IF(CDTPRO.EQ.'00000000000000') CDTPRO = ZERO
        IF(CDATEWO.EQ.'00000000000000')CDATEWO = ZERO
        IF(CDAWIFL.EQ.'00000000000000')CDAWIFL = ZERO
        IF(CDATEFL.EQ.'00000000000000') CDATEFL= ZERO

      ELSE
        IF(LGRIBIN.AND..NOT.LRESTARTED) THEN
          CALL BUILDSTRESS(U10OLD,THWOLD,USOLD,TAUW,Z0OLD,
     &                     ROAIRO,ZIDLOLD,ICEMASK,
     &                     IREAD)
        ELSE
          CALL GRSTNAME(CDATEA,CDATEF,'LAW',ICPLEN,CPATH,FILENAME)
          CALL READSTRESS(NREAL, RFIELD, NINTGR, IFIELD,
     &                    FILENAME)

!         KEEP CORRESPONDING CONTRIBUTION 
          DO IG=1,NBLO
            NPROMA=(NEND(IRANK)-NSTART(IRANK)+1)/NWAM_BLKS+1
!$OMP       PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
            DO JKGLO=NSTART(IRANK),NEND(IRANK),NPROMA
              KIJS=JKGLO
              KIJL=MIN(KIJS+NPROMA-1,NEND(IRANK))
              DO IJ=KIJS,KIJL
                U10OLD(IJ,IG)=RFIELD(IJ,IG,1)
                THWOLD(IJ,IG)=RFIELD(IJ,IG,2)
                USOLD(IJ,IG)=RFIELD(IJ,IG,3)
                TAUW(IJ,IG)=RFIELD(IJ,IG,4)
                Z0OLD(IJ,IG)=RFIELD(IJ,IG,5)
                ROAIRO(IJ,IG)=RFIELD(IJ,IG,6)
                ZIDLOLD(IJ,IG)=RFIELD(IJ,IG,7)
!               for U and V see below
                ICEMASK(IJ,IG)=IFIELD(IJ,IG,1)
              ENDDO
            ENDDO
!$OMP       END PARALLEL DO
          ENDDO

          IF(ALLOCATED(U) .AND. ALLOCATED(V) ) THEN
            DO IG=1,NBLO
              NPROMA=(NEND(IRANK)-NSTART(IRANK)+1)/NWAM_BLKS+1
!$OMP         PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
              DO JKGLO=NSTART(IRANK),NEND(IRANK),NPROMA
                KIJS=JKGLO
                KIJL=MIN(KIJS+NPROMA-1,NEND(IRANK))
                DO IJ=KIJS,KIJL
                  U(IJ,IG)=RFIELD(IJ,IG,8)
                  V(IJ,IG)=RFIELD(IJ,IG,9)
                ENDDO
              ENDDO
!$OMP         END PARALLEL DO
              U(NINF-1,IG)=0.
              V(NINF-1,IG)=0.
            ENDDO

!           RE-COMPUTE REFRACTION TERMS
            CALL PROPDOT

!           SET LOGICAL TO RECOMPUTE THE WEIGHTS IN CTUW.
            LUPDTWGHT=.TRUE.

          ENDIF
        ENDIF
      ENDIF

      IF(LNSESTART .AND. .NOT.LRESTARTED) THEN
!       RESET WAVE INDUCED STRESS TO ZERO 
!       BECAUSE THE SPECTRA WILL BE RESET TO NOISE
        DO IG=1,NBLO
          NPROMA=(NEND(IRANK)-NSTART(IRANK)+1)/NWAM_BLKS+1
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
          DO JKGLO=NSTART(IRANK),NEND(IRANK),NPROMA
            KIJS=JKGLO
            KIJL=MIN(KIJS+NPROMA-1,NEND(IRANK))
            DO IJ=KIJS,KIJL
              TAUW(IJ,IG)=0.
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
        ENDDO
      ENDIF

      DEALLOCATE(RFIELD, IFIELD)

      WRITE(IU06,*) ''
      WRITE(IU06,*) ' WIND AND STRESS FILES READ IN.......',
     &              ' CDTPRO  = ', CDTPRO
      CALL FLUSH (IU06)

      IF (LHOOK) CALL DR_HOOK('GETSTRESS',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE GETSTRESS
