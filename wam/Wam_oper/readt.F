!------------------------------------------------------------------

      SUBROUTINE READT(IU06,FILNM,ITEST,TEMPG)

!----------------------------------------------------------------------

!**** *READT*  READ TEMPERATURE GRID.

!     R. PORTZ      MPI HAMBURG OCTOBER 1992

!*    PURPOSE.
!     --------

!       INPUT OF GLOBAL TEMPARATURE.

!**   INTERFACE.
!     ----------

!       *CALL* *READT* *(IU06,FILNM,ITEST,TEMPG)*

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
!      *FILNM*     CHARACTER*24 DATA INPUT FILENAME.
!      *ITEST*     INTEGER   OUTPUT FLAG IF LE 2 NO EXTRA OUTPUT IS
!                            GENERATED.
!      *TEMPG*     REAL      GLOBAL TEMPRERATURES (KELVIN).

!     METHOD.
!     -------

!       FORMATED READ FROM UNIT 1.

!     EXTERNALS.
!     ----------

!       ABORT1
!       ADJUST
!       INMARS 
!       OUTPP

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3
      USE YOWPCONS , ONLY : ZMISS
      USE YOWTEMP  , ONLY : NXT      ,NYT      ,CTDATE   ,NGXT     ,
     &            NGYT     ,TMOWEP   ,TMOSOP   ,TMOEAP   ,TMONOP   ,
     &            DELLA    ,DELLO

!-----------------------------------------------------------------------

      REAL,DIMENSION(NXT,NYT) :: TEMPG
      REAL,ALLOCATABLE :: WORK(:,:)
      INTEGER,ALLOCATABLE :: KLONRGG(:)

      CHARACTER*24 FILNM
      CHARACTER*40 MSG, TITL*100

!-----------------------------------------------------------------------

!*    1. INPUT OF TEMPERATURE DATA.
!     -----------------------------

      ALLOCATE(KLONRGG(NYT))
      CALL INMARS  ( FILNM, NXT, NYT, ZMISS, NKSEK0, NKSEK1, NKSEK2,
     &               NKSEK3, NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &               CTDATE, IFORP, IPARAM, NGXT, NGYT,
     &               TMOWEP, TMOSOP, TMOEAP, TMONOP, KRGG, KLONRGG,
     &               KZLEV,TEMPG, IERR, MSG, IU06)
      IF ( IERR .EQ. 0 ) THEN
        WRITE(IU06,*)' SUB. READT: ICE BOUNDARY DATE=',CTDATE
        CALL FLUSH(IU06)
      ELSE
        WRITE(IU06,*)'   SUB. READT: ICE BOUNDARY IERR=',IERR
        WRITE(IU06,*) MSG
        CALL ABORT1
      ENDIF
      CALL ADJUST (TMOWEP,TMOEAP)
      DELLA  = (TMONOP - TMOSOP)/REAL(NGYT-1)
      DELLO  = (TMOEAP - TMOWEP)/REAL(NGXT-1)

      IF (ITEST.GT.3) THEN
        ALLOCATE(WORK(NXT,NYT))
        WRITE(IU06,*) 'TEMPERATURE GRID:'
        WRITE(IU06,*) '-----------------'
        WRITE(IU06,*)
        WRITE(IU06,*) '  RESOLUTION LAT-LON  ',DELLA, DELLO
        WRITE(IU06,*) '   SOUTHERN LAT: ',TMOSOP,
     &                '   NORTHERN LAT: ',TMONOP
        WRITE(IU06,*) '   WESTERN LONG: ',TMOWEP,
     &                '   EASTERN LONG: ',TMOEAP
        WRITE(IU06,*)
        WORK=TEMPG-273.15
        TITL= '   TEMPERATURES IN CELCIUS         '
        CALL OUTPP (CTDATE, IU06, NXT, NYT, NGXT, NGYT,
     &              TITL, 1.,
     &              WORK, TMOWEP, TMOSOP, TMOEAP, TMONOP)
        DEALLOCATE(WORK)
        CALL FLUSH(IU06)
      ENDIF

!*    1.1 ERROR MESSAGES
!     ------------------


      IF (NGXT .GT. NXT) THEN
        WRITE(IU06,*)'********************************'
        WRITE(IU06,*)'*                              *'
        WRITE(IU06,*)'* FATAL ERROR IN SUB READT     *'
        WRITE(IU06,*)'* ========================     *'
        WRITE(IU06,*)'*                              *'
        WRITE(IU06,*)'* PARAMETER NXT IS TO SMALL    *'
        WRITE(IU06,*)'* SET NXT = ',NGXT
        WRITE(IU06,*)'*                              *'
        WRITE(IU06,*)'********************************'
        CALL ABORT1
      ENDIF
      IF (NGYT .GT. NYT) THEN
        WRITE(IU06,*)'********************************'
        WRITE(IU06,*)'*                              *'
        WRITE(IU06,*)'* FATAL ERROR IN SUB READT     *'
        WRITE(IU06,*)'* ========================     *'
        WRITE(IU06,*)'*                              *'
        WRITE(IU06,*)'* PARAMETER NXT IS TO SMALL    *'
        WRITE(IU06,*)'* SET NYT = ',NGYT
        WRITE(IU06,*)'*                              *'
        WRITE(IU06,*)'********************************'
        CALL ABORT1
      ENDIF

      IF(ALLOCATED(KLONRGG))DEALLOCATE(KLONRGG)

      RETURN
      END SUBROUTINE READT
