      SUBROUTINE READICE (IU06,FILNM,ITEST,NSTART,NEND,IREAD,ITAG)

!--------------------------------------------------------------------

!**** *READICE* *DEFINES SEA POINTS AS ICE IF THE CORRESPONDING
!                SST IS LT 271.5 KELVIN (IF SST IS USED AS INPUT)
!                OR
!                SEA ICE FRACTION GT 0.3 (IF IT IS USED AS INPUT)

!     B. HANSEN    MPIFM/ECMWF   NOVEMBER 92
!     J. BIDLOT    ECMWF         JUNE 1996  MESSAGE PASSING
!     P. JANSSEN   ECMWF         OCTOBER 1997 SIMPLIFY BY INTRODUCING
!                                             ICEMASK
!     J. BIDLOT    ECMWF         JANUARY 1998 : CORRECT THRESHOLD FOR
!                                               ICE TO BE LT 271.50 
!     J. BIDLOT    ECMWF         FEBRUARY 2000 : USE OF SEA ICE FRACTION
 

!     J. BIDLOT    ECMWF
!!!! this is now obsolete as the information for the ice mask is now
!!!! read with the wind (see *READWIND*)


!     PURPOSE.
!     --------
!            DYNAMICALLY ADJUST NORTH AND SOUTH ICE BOUNDARIES
!            OF THE WAMODEL.

!**   INTERFACE
!     ---------
!     *CALL* *READICE (IU06,FILNM,ITEST,NSTART,NEND,IREAD,ITAG)*

!     *IU06*      OUTPUT UNIT FOR STANDARD OUTPUT.
!     *FILNM*     DATA INPUT FILENAME.
!     *ITEST*     OUTPUT FLAG IF LE 2 NO EXTRA OUTPUT IS GENERATED
!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *IREAD*     RANK OF THE PROCESS WHICH INPUTS THE SST FIELD


!     METHOD.
!     -------
!            READS GLOBAL TEMPRERATURES (UNIT 99).
!            DEFINES BLOCK(IG) AND INDEX(IJ) FOR EACH POINT
!            WITH TEMPERATURE < 271.50 K  -->  ICE.
!            EXTRACT WATER POINTS.
!            OR
!            DOES THE SAME BUT WITH THE SEA ICE FRACTION WHERE
!            ANY POINT WITH SEA ICE FRACTION > 0.3 IS CONSIDERED ICE.
!          


!    EXTERNAL.
!    ---------
!             *MICEP*     DEFINE BLOCK AND INDEX IN A BLOCK FOR
!                         EACH ICE POINT (T < 271.5 KELVIN).
!             *READT*     READ GLOBAL TEMPERATURES.
!             *MPBCASTGRDFLD* BROADCAST SST FIELD TO ALL OTHER PE's

!---------------------------------------------------------------------

      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL
      USE YOWICE   , ONLY : ICEMASK
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWSTAT  , ONLY : NWAM_BLKS

!-----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA
! End Mod for OPENMP


!-----------------------------------------------------------------------

      CHARACTER*24 FILNM
      REAL,ALLOCATABLE :: FIELD(:,:)

!*    1. READ SST FROM UNIT 90.
!        ----------------------

      ALLOCATE(FIELD(NGX,NGY))

      IF (LMESSPASS) THEN

!        READ SST FIELD ON PE IREAD

        IF (IRANK.EQ.IREAD) CALL READT(IU06,FILNM,ITEST,FIELD,IPARAM)

!       BROADCAST FIELD TO OTHER PE'S

        CALL MPBCASTGRDFLD(IREAD,ITAG,NGX,NGY,FIELD)
        ITAG=ITAG+1

!       BROADCAST IPARAM TO OTHER PE'S

        CALL MPBCASTINTFLD(IREAD,ITAG,1,1,IPARAM)
        ITAG=ITAG+1

      ELSE
        CALL READT(IU06,FILNM,ITEST,FIELD,IPARAM)
      ENDIF

!-------------------------------------------------------------------

!*    2. DEFINE ICE POINTS AND COMPUTE ICEMASK
!        -------------------------------------

      DO IG=1,IGL
        MIJS=IJS(IG)
        MIJL=IJL(IG)
! Mod for OPENMP
        NPROMA=(MIJL-MIJS+1)/NWAM_BLKS+1
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=MIJS,MIJL,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,MIJL)
          CALL MICEP (FIELD, IPARAM, ICEMASK, KIJS, KIJL, IG)
        ENDDO
!$OMP   END PARALLEL DO
      ENDDO

      DEALLOCATE(FIELD)

      RETURN
      END SUBROUTINE READICE
