      SUBROUTINE DOTDC (IJS, IJL, IG, IGL, ISHALLO, IREFRA, SL, FCONST)

! ----------------------------------------------------------------------

!**** *DOTDC* - SIGMA DOT FOR CURRENT REFRACTION.

!     H. GUNTHER   GKSS/ECMWF   17/02/91
!     J. BIDLOT    ECMWF FEBRUARY 1997 : MESSAGE PASSING

!*    PURPOSE.
!     --------

!       READ (FROM IU16) FOR PROPAGATION ..
!           - DEPTH AND CURRENT PART OF THETA DOT.
!           - SIGMA DOT TERM.
!       SCATTER SIGMA/ SINH (2*K*D) TABLE.

!**   INTERFACE.
!     ----------

!       *CALL* *DOTDC (IJS, IJL, IG, IGL, ISHALLO, IREFRA, SL, FCONST)*
!          *IJS*     - INDEX OF FIRST GRID POINT.
!          *IJL*     - INDEX OF LAST  GRID POINT.
!          *IG*      - BLOCK NUMBER.
!          *IGL*     - NUMBER OF BLOCKS.
!          *ISHALLO* - SHALLOW WATER OPTION.
!          *IREFRA*  - REFRACTION OPTION.
!          *SL*      - SIGMA DOT
!          *FCONST * - SIGMA/SINH(2KD)

!     METHOD.
!     -------

!       SEQUENTIAL READ FROM UNIT.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO    ,NIBLC
      USE YOWREFD  , ONLY : THDD     ,THDC
      USE YOWSHAL  , ONLY : TSIHKD   ,INDEP
      USE YOWUNIT  , ONLY : IU16

! ----------------------------------------------------------------------

      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: SL
      REAL,DIMENSION(NINF:NSUP,NFRE):: FCONST

! ----------------------------------------------------------------------

      CHARACTER*30 FILENAME

!*    1. READ COMMON REFDOT.
!        -------------------

      IF (IG.EQ.1) REWIND IU16
      IF (IGL.NE.1) THEN
        IF (ISHALLO.NE.1) THEN
          IF (IREFRA.EQ.1) THEN
            READ (IU16) THDD
          ELSE
            READ (IU16) THDD, THDC, SL
          ENDIF
        ELSE
          READ (IU16) THDC, ((SL(IJ,K,NFRE),IJ=1,NIBLC),K=1,NANG)
        ENDIF
      ELSE
        IF(LMESSPASS) THEN
          FILENAME='fort.16.%p'
          CALL EXPAND_STRING(IRANK,NPROC,0,0,FILENAME,1)
          OPEN(IU16,FILE=FILENAME,STATUS='UNKNOWN')

          READ(IU16) SL
          CLOSE(IU16)

        ELSE
          IF (ISHALLO.NE.1) THEN
            IF (IREFRA.EQ.2) THEN
              READ (IU16) SL
            ENDIF
          ELSE
            READ (IU16) ((SL(IJ,K,NFRE),IJ=1,NIBLC),K=1,NANG)
          ENDIF
        ENDIF
      ENDIF

!*    2. IF DEEP WATER RETURN.
!        ---------------------

      IF (ISHALLO.EQ.1) RETURN

!*    3. GATHER SIGMA /SINH(2KD) FROM TABLE.
!        -----------------------------------

      DO M=1,NFRE
        DO IJ=IJS,IJL
          FCONST(IJ,M) = TSIHKD(INDEP(IJ),M)
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DOTDC
