      SUBROUTINE GETCURR(LWCUR, IREAD)

!****  *GETCURR* - READS SURFACE CURRENTS FROM FILE (IF UNCOUPLED)
!                  OR EXTRACT THEM FROM FORCING_FIELDS DATA STRUCTURE (IF COUPLED)

!     PURPOSE.
!     --------

!       GET SURFACE CURRENTS AND COMPUTE THE ASSOCIATED
!       REFRACTION TERMS.

!*    INTERFACE.
!     ----------

!     CALL *GETCURR*(LWCUR, IREAD)

!     *LWCUR*   -  LOGICAL CONTROLLING THE PRESENCE OF MEANINGFUL
!                  SURFACE CURRENTS IN ARRAY FORCING_FIELDS DATA STRUCTURE   
!     *IREAD*      INTEGER   PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                            (IF NEEDED).


!     METHOD.
!     -------

!     IF UNCOUPLED
!     READS CURRENTS FROM FILE "currents" WHEN IT IS NEEDED
!     IF COUPLED
!     GET CURRENTS AS PROCESSED BY *IFSTOWAM*
!     THEN 
!     CALL PROPDOT FOR THE CALCULATION OF THE REFRACTION TERMS.
!     ALSO SET LLCHKCFLA=.TRUE. TO CHECK ON THE CFL CRITERIA.

!     EXTERNALS.
!     ----------

!       *INCDATE*
!       *ABORT1*
!       *CURRENT2WAM*
!       *PROPDOT*

!     REFERENCES.
!     ----------- 

!       NONE

! ------------------------------------------------------------------- 

      USE YOWCOUP  , ONLY : LWCOU, LWNEMOCOUCUR
      USE YOWNEMOFLDS,ONLY: NEMOUCUR, NEMOVCUR
      USE YOWCURR  , ONLY : U        ,V        ,CDTCUR   ,IDELCUR  ,
     &             LLCHKCFL,LLCHKCFLA
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : NINF     ,NSUP
      USE YOWPARAM , ONLY : NIBLO    ,NBLO
      USE YOWSTAT  , ONLY : CDTPRO   ,IREFRA   ,NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE YOWWIND  , ONLY : LLNEWCURR 

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! --------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IG
      INTEGER :: LIU, IREAD
      INTEGER :: JKGLO,KIJS,KIJL,NPROMA,IJ

      REAL :: ZHOOK_HANDLE

      CHARACTER(LEN=14) :: CDATEIN, CDTNEWCUR
      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: LWCUR
      LOGICAL :: LLCURRENT

! --------------------------------------------------------------------- 
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('GETCURR',0,ZHOOK_HANDLE)
#endif

!     1.0  GET NEW CURRENTS IF IT IS REQUIRED.
!          ----------------------------------

      CALL GSTATS(1984,0)

      IF (LLNEWCURR) THEN
        IF ( (LWCOU .AND. LWCUR ) .OR.
     &        IREFRA.EQ.2 .OR. IREFRA.EQ.3) THEN

          IF(.NOT.ALLOCATED(U)) ALLOCATE(U(NINF-1:NSUP,NBLO))
          IF(.NOT.ALLOCATED(V)) ALLOCATE(V(NINF-1:NSUP,NBLO))

          CDTNEWCUR=CDTCUR

          IF(.NOT.LWCOU) CALL INCDATE(CDTNEWCUR,IDELCUR/2)

          IF(CDTPRO.GE.CDTNEWCUR) THEN

            LLCURRENT=.FALSE.

            CALL INCDATE(CDTCUR,IDELCUR)

            IF(LWCOU) THEN
!             CURRENTS FROM COUPLING INTERFACE 
!             --------------------------------
              IF(LWCUR.AND..NOT.LWNEMOCOUCUR) THEN

                DO IG=1,IGL
!!!               from NINF to NSUP as the halo has to be included
!!!               for the calculation of the gradients !!!
! Mod for OPENMP
                  CALL GSTATS(1444,0)
                  NPROMA=(NSUP-NINF+1)/NWAM_BLKS+1
!$OMP           PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
                  DO JKGLO=NINF,NSUP,NPROMA
                    KIJS=JKGLO
                    KIJL=MIN(KIJS+NPROMA-1,NSUP)
                    CALL WAMCUR (U, V, NINF, NSUP, KIJS, KIJL, IG)
                  ENDDO
!$OMP           END PARALLEL DO
                  CALL GSTATS(1444,1)
                  U(NINF-1,IG)=0.
                  V(NINF-1,IG)=0.
                ENDDO ! END LOOP ON IG

                IF (ITEST.GE.1) THEN
                   WRITE (IU06,*) ' '
                   WRITE (IU06,*) '  SUB. GETCURR:',
     &             ' INPUT CURRENTS FIELDS CONVERTED TO BLOCKS'
                   CALL FLUSH(IU06)
                ENDIF
!             CURRENTS FROM NEMO
!             ------------------
              ELSEIF(LWCUR.AND.LWNEMOCOUCUR) THEN
                WRITE(IU06,*)'NEMO CURRENTS'!
                IG=1
                DO IJ = IJS(IG),IJL(IG)
                  U(IJ,IG) = NEMOUCUR(IJ)
                  V(IJ,IG) = NEMOVCUR(IJ)
                ENDDO
                U(NINF-1,IG)=0.
                V(NINF-1,IG)=0.
                CALL MPEXCHNG(U,1,1)
                CALL MPEXCHNG(V,1,1)
              ELSE
                U=0.
                V=0.
                WRITE(IU06,*)' '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)'     IREFRA = ',IREFRA
                WRITE(IU06,*)'     LWCUR IS FALSE !!!!' 
                WRITE(IU06,*)'     CURRENTS ARE SET TO 0. '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)' '
                CALL FLUSH(IU06)
              ENDIF
            ELSE
!             CURRENTS FROM INPUT FILE
!             ------------------------
              FILNM = 'currents'
              LIU   = LEN_TRIM(FILNM)
              FILNM=FILNM(1:LIU)
              INQUIRE(FILE=FILNM,EXIST=LLCURRENT)

              IF (LLCURRENT) THEN
                WRITE(IU06,*) ' '
                WRITE(IU06,*) ' SUB. GETCURR: GETTING OCEAN CURRENTS',
     &                        ' FOR DATE ',CDTCUR
                CALL FLUSH(IU06)

                CALL CURRENT2WAM (FILNM,IREAD,CDATEIN)

                IF(CDATEIN.NE.CDTCUR) THEN
                WRITE (IU06,*) ' **************************************'
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' * PROBLEM IN GETCURR :               *'
                WRITE (IU06,*) ' * THE REQUESTED DATE FOR THE CURRENTS*'
                WRITE (IU06,*) ' * DOES NOT CORRESPOND TO THE DECODED *'
                WRITE (IU06,*) ' * DATE !!!!                          *'
                WRITE (IU06,*) ' * CDTCUR =',CDTCUR 
                WRITE (IU06,*) ' * CDATEIN=',CDATEIN
                WRITE (IU06,*) ' *                                    *'
                WRITE (IU06,*) ' **************************************'
                CALL ABORT1
                ENDIF
              ELSE
                U=0.
                V=0.
                WRITE(IU06,*)' '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)'     FILE ',FILNM,' NOT FOUND '
                WRITE(IU06,*)'     CURRENTS ARE SET TO 0. '
                WRITE(IU06,*)'    ****************************'
                WRITE(IU06,*)' '
                CALL FLUSH(IU06)
              ENDIF

            ENDIF


!           COMPUTE REFRACTION TERMS
            CALL PROPDOT
            IF (ITEST.GE.2)
     &       WRITE(IU06,*) ' SUB. GETCURR: REFRACTION TERMS',
     &       ' INITIALIZED'
             CALL FLUSH(IU06)

            LLCHKCFLA=.TRUE.

!           SET LOGICAL TO RECOMPUTE THE WEIGHTS IN CTUW.
            LUPDTWGHT=.TRUE.

          ELSE
            LLCHKCFLA=.FALSE.
          ENDIF
        ENDIF
      ELSE
        LLCHKCFLA=.FALSE.
      ENDIF

      CALL GSTATS(1984,1)

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('GETCURR',1,ZHOOK_HANDLE)
#endif

      RETURN                                                            
      END SUBROUTINE GETCURR 
