      SUBROUTINE MBOUNC (IU09, IU19, IFORM)

! ----------------------------------------------------------------------

!**** *MBOUNC* - MAKE COARSE GRID BOUNDARY.

!     R. PORTZ     MPI         15/01/1991

!*    PURPOSE.
!     -------

!       COMPUTE ALL INFORMATION FOR COARSE GRID BOUNDARY VALUE
!       OUTPUT (COMMON CBOUND).

!**   INTERFACE.
!     ----------

!       *CALL* *MBOUNC (IU09, IU19, IFORM)*
!          *IU09*   - LOGICAL UNIT FOR  UNFORMATED WRITE.
!          *IU19*   - LOGICAL UNIT FOR    FORMATED WRITE.
!          *IFORM*  - FORMAT OPTION  = 1  UNFORMATED WRITE.
!                                    = 2  FORMATED WRITE.
!                                      OTHERWISE BOTH.

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       *FINDB*     - FIND BLOCK AND SEA POINT NUMBERS.
!       *MBOXB*     - FIND LAT AND LONG OF BOUNDARY POINTS.
!       *PACKI*     - PACKS AN INTEGER ARRAY.
!       *PACKC*     - PACKS A  REAL    ARRAY.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWCPBO  , ONLY : NBOUNC   ,IJARC    ,IGARC    ,DLAMAC   ,
     &            DPHIAC   ,AMOSOC   ,AMONOC   ,AMOEAC   ,AMOWEC   ,
     &            BLATC    ,BLNGC
      USE YOWMAP   , ONLY : XDELLA   ,XDELLO
      USE YOWPRPROC,ONLY : NBMAX
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

  998 FORMAT(10I8)
  999 FORMAT(5E16.7)

! ----------------------------------------------------------------------

!*    1. INITIAL.
!        --------
      
      DPHIAC = XDELLA
      DLAMAC = XDELLO

      NBMAX = (NGX+NGY)*2-4

      ALLOCATE(IJARC(NBMAX))
      ALLOCATE(IGARC(NBMAX))
      ALLOCATE(BLATC(NBMAX))
      ALLOCATE(BLNGC(NBMAX))

      DO I = 1,NBMAX
        IJARC(I) = 0
        IGARC(I) = 0
        BLATC(I) = 0.
        BLNGC(I) = 0.
      ENDDO

! ----------------------------------------------------------------------

!*    2. COMPUTED THE SQUARE BOX.
!        ------------------------

      CALL MBOXB (NBOUNC, AMOWEC, AMOSOC, AMOEAC, AMONOC, BLATC, BLNGC)

! ----------------------------------------------------------------------

!*    3. SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
!        -----------------------------------------

      write(*,*) 'the call the findb has changed. It must be adapted'
      write(IU06,*) 'the call the findb has changed. It must be adapted'
      call flush(IU06)
      call abort1
     
      CALL FINDB (NBMAX, NBOUNC, BLATC, BLNGC, IGARC, IJARC)

! ----------------------------------------------------------------------

!*    4. PACKED ALL ARRAYS.
!        -------------------

      CALL PACKR (NBOUNC, NBOUNEW, NBMAX, IGARC, BLATC)
      CALL PACKR (NBOUNC, NBOUNEW, NBMAX, IGARC, BLNGC)
      CALL PACKI (NBOUNC, NBOUNEW, NBMAX, IGARC, IJARC)
      CALL PACKI (NBOUNC, NBOUNEW, NBMAX, IGARC, IGARC)
      NBOUNC = NBOUNEW

! ----------------------------------------------------------------------

!*    5. PRINTER PROTOCOL.
!        -----------------

      WRITE (IU06,'(''1BOUNDARY OUTPUT POINTS OF THIS GRID:'')')
      WRITE (IU06,*) '===================================='
      WRITE (IU06,*) '  NUMBER OF BOUNDARY POINTS IS NBOUNC = ', NBOUNC
      WRITE (IU06,'(4X,''     |-----INPUT-----|'',
     &               ''-POINT INDEX--|'')')
      WRITE (IU06,'(4X,''  NO.    LAT.   LONG.  BLOCK.'',
     &              ''  POINT.'')')
      DO IO=1,NBOUNC
        WRITE (IU06,'(4X,I5,2F8.2,2I8)')
     &   IO, BLATC(IO), BLNGC(IO), IGARC(IO), IJARC(IO)
      ENDDO

! ----------------------------------------------------------------------

!*    6. WRITE COMMON CBOUND.
!        --------------------

      IF (IFORM.NE.2) THEN
        WRITE(IU09) NBOUNC
        WRITE(IU09) (IGARC(I),I=1,NBOUNC), (IJARC(I),I=1,NBOUNC)
        WRITE(IU09) DLAMAC, DPHIAC, AMOSOC, AMONOC, AMOEAC, AMOWEC,
     &   (BLNGC(I),I=1,NBOUNC), (BLATC(I),I=1,NBOUNC)
      ENDIF
      IF (IFORM.NE.1) THEN
        WRITE(IU19,998) NBOUNC
        WRITE(IU19,998) (IGARC(I),I=1,NBOUNC), (IJARC(I),I=1,NBOUNC)
        WRITE(IU19,999) DLAMAC,DPHIAC, AMOSOC, AMONOC, AMOEAC, AMOWEC,
     &   (BLNGC(I),I=1,NBOUNC), (BLATC(I),I=1,NBOUNC)
      ENDIF

! ----------------------------------------------------------------------

      RETURN
      END SUBROUTINE MBOUNC
