      SUBROUTINE READICE (IU06,FILNM,ITEST,NSTART,NEND,IREAD,ITAG)

!--------------------------------------------------------------------

!**** *READICE* *DEFINES SEA POINTS AS ICE IF THE CORRESPONDING
!                SST IS LT 271.5 KELVIN.

!     B. HANSEN    MPIFM/ECMWF     NOVEMBER 92
!     J. BIDLOT    ECMWF           JUNE 1996  MESSAGE PASSING
!     P. JANSSEN   ECMWF           OCTOBER 1997 SIMPLIFY BY INTRODUCING
!                                               ICEMASK
!     J. BIDLOT    ECMWF           JANURAY 1998 : CORRECT THRESHOLD FOR
!                                                 ICE TO BE LT 271.50 

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
!            READ GLOBAL TEMPRERATURES (UNIT 99).
!            DEFINE BLOCK(IG) AND INDEX(IJ) FOR EACH POINT
!            WITH TEMPERATURE < 271.50 K  -->  ICE.
!            EXTRACT WATER POINTS.
!            SORT ICE POINTS ACCORDING TO BLOCK AND INDEX IN A BLOCK.
!            KEEP INFORMATION IN COMMON COMICE FOR LATER USE.


!    EXTERNAL.
!    ---------
!             *MICEP*     DEFINE BLOCK AND INDEX IN A BLOCK FOR
!                         EACH ICE POINT (T < 271.5 KELVIN).
!             *READT*     READ GLOBAL TEMPERATURES.
!             *MPBCASTGRDFLD* BROADCAST SST FIELD TO ALL OTHER PE's

!---------------------------------------------------------------------

      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWTEMP  , ONLY : NXT      ,NYT      ,TMOWEP   ,TMOSOP   ,
     &            TMOEAP   ,TMONOP   ,DELLA    ,DELLO

!-----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND

!-----------------------------------------------------------------------

      CHARACTER*24 FILNM
      REAL,ALLOCATABLE :: TEMPG(:,:), WORK(:)

!*    1. READ SST FROM UNIT 90.
!        ----------------------


      ALLOCATE(TEMPG(NXT,NYT),WORK(6))
      TEMPG=273.15

      IF (LMESSPASS) THEN

!        READ SST FIELD ON PE IREAD

        IF (IRANK.EQ.IREAD) CALL READT(IU06,FILNM,ITEST,TEMPG)

!        BROADCAST  SST FIELD ARRAY INFORMATION TO OTHER PE'S

        IF (IRANK.EQ.IREAD) THEN
          WORK(1)=TMOWEP
          WORK(2)=TMOSOP
          WORK(3)=TMOEAP
          WORK(4)=TMONOP
          WORK(5)=DELLA
          WORK(6)=DELLO
        ENDIF

        CALL MPBCASTGRDFLD(IREAD,ITAG,6,1,WORK)
        ITAG=ITAG+1

        IF (IRANK.NE.IREAD) THEN
          TMOWEP=WORK(1)
          TMOSOP=WORK(2)
          TMOEAP=WORK(3)
          TMONOP=WORK(4)
          DELLA=WORK(5)
          DELLO=WORK(6)
        ENDIF

!        BROADCAST  SST FIELD TO OTHER PE'S

        CALL MPBCASTGRDFLD(IREAD,ITAG,NXT,NYT,TEMPG)
        ITAG=ITAG+1

      ELSE
        CALL READT(IU06,FILNM,ITEST,TEMPG)
      ENDIF

!-------------------------------------------------------------------

!*    2. DEFINE ICE POINTS AND COMPUTE ICEFLAG AND WRITE COMMON BLOCK
!        COMICE.
!        ----------------------------------------------

      CALL MICEP (IU06, NSTART, NEND, TEMPG)

      DEALLOCATE(TEMPG,WORK)

      RETURN
      END SUBROUTINE READICE
