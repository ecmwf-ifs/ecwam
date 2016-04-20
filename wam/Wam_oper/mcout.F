      SUBROUTINE MCOUT

! ----------------------------------------------------------------------

!**** *MCOUT* - ROUTINE TO COMPUTE OUTPUT INDICES (COMMON COUT).

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       *MCOUT* COMPUTES THE INDICES OF SPECTRA OUTPUT POINTS.

!**   INTERFACE.
!     ----------

!       *CALL* *MCOUT*

!     METHOD.
!     -------

!       THE LATITUDE AND LOGITUDE ARE CONVERTED TO INDICES.

!     EXTERNALS.
!     ----------

!       *FINDB*     - FIND BLOCK AND SEA POINT NUMBERS.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE YOWCINP  , ONLY : OUTLONG  ,OUTLAT
      USE YOWCOUT  , ONLY : NGOUT    ,IGAR     ,IJAR
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWPARAM , ONLY : NIBLO
      USE YOWSPEC,   ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IO, KX, IX, NG
      INTEGER :: NGOUTNEW
      INTEGER, ALLOCATABLE :: IDUM(:)
      INTEGER :: ListSTART(1), ListEND(1)
      REAL :: ALONG, ALAT

!*    1. NO OUTPUT POINTS SPECIFIED.
!        ---------------------------

      IF (NGOUT.EQ.0) THEN
        WRITE(IU06,'(1H1,'' SPECIAL OUTPUT POINTS FOR SPECTRA:'')')
        WRITE(IU06,*) 'OUTPUT POINTS ARE NOT DEFINED IN USER INPUT'
        RETURN
      ELSE
        ALLOCATE(IGAR(NGOUT))
        ALLOCATE(IJAR(NGOUT))
      ENDIF

! ----------------------------------------------------------------------

!*    2. SEARCH BLOCK NUMBER AND SEA POINT NUMBER.
!        -----------------------------------------

      ListSTART(1)=1
      ListEND(1)=1
      CALL FINDB (NGOUT, NGOUT, OUTLAT, OUTLONG, IGAR, IJAR,
     &            1, ListSTART,ListEND,1)

! ----------------------------------------------------------------------

!*    3. PRINTER PROTOCOL.
!        -----------------

      WRITE(IU06,'(1H1,'' SPECIAL OUTPUT POINTS FOR SPECTRA:'')')
      WRITE(IU06,'(''    NUMBER OF OUTPUT POINTS IS NGOUT = '',I4)')
     &      NGOUT
      WRITE(IU06,'(4X,''     |-----INPUT-----|-NEAREST POINT-|'',
     &              ''-POINT INDEX--|'')')
      WRITE(IU06,'(4X,''  NO.    LAT.   LONG.    LAT.   LONG.  BLOCK.'',
     &             ''  POINT.'')')
      DO IO=1,NGOUT
        IF (IJAR(IO).GT.0) THEN
          KX  = KXLT(IJAR(IO),IGAR(IO))
          IX  = IXLG(IJAR(IO),IGAR(IO))
          ALONG = AMOWEP + (IX-1)*ZDELLO(KX)
          ALAT  = AMOSOP + (KX-1)*XDELLA
        ELSE
          ALONG = 9999999
          ALAT  = 9999999
        ENDIF
        WRITE(IU06,'(4X,I5,4F8.2,2I8)')
     &   IO, OUTLAT(IO), OUTLONG(IO), ALAT, ALONG, IGAR(IO), IJAR(IO)
      ENDDO

! ----------------------------------------------------------------------

!*    4. REMOVE OUTPUT POINTS WHICH ARE NOT IN GRID.
!        -------------------------------------------

      NG = 0
      DO IO=1,NGOUT
        IF (NG.GT.0 .AND. IO-NG.GT.0) THEN
          IJAR(IO-NG) = IJAR(IO)
          IGAR(IO-NG) = IGAR(IO)
        ENDIF
        IF (IJAR(IO).EQ.0) NG = NG+1
      ENDDO
      NGOUTNEW = NGOUT-NG
      IF (NG.GT.0) THEN
        WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
        WRITE (IU06,*) ' +                                       +'
        WRITE (IU06,*) ' +     WARNING ERROR FROM SUB. MCOUT     +'
        WRITE (IU06,*) ' +     =============================     +'
        WRITE (IU06,*) ' +                                       +'
        WRITE (IU06,*) ' + NO SEAPOINT FOUND FOR NG = ',NG
        WRITE (IU06,*) ' + OUTPUT POINT REQUESTS (SEE ABOVE LIST)+'
        WRITE (IU06,*) ' + THESE POINTS WILL NOT BE TAKEN.       +'
        WRITE (IU06,*) ' + NUMBER OF OUTPUT POINTS IS NGOUT = ',NGOUTNEW
        WRITE (IU06,*) ' +                                       +'
        WRITE (IU06,*) ' +++++++++++++++++++++++++++++++++++++++++'
        ALLOCATE(IDUM(NGOUTNEW))
        DO IO=1,NGOUTNEW
           IDUM(IO) = IJAR(IO)
        ENDDO
        DEALLOCATE(IJAR)
        ALLOCATE(IJAR(NGOUTNEW))
        DO IO=1,NGOUTNEW
           IJAR(IO) = IDUM(IO)
        ENDDO
        DO IO=1,NGOUTNEW
           IDUM(IO) = IGAR(IO)
        ENDDO
        DEALLOCATE(IGAR)
        ALLOCATE(IGAR(NGOUTNEW))
        DO IO=1,NGOUTNEW
           IGAR(IO) = IDUM(IO)
        ENDDO
        DEALLOCATE(IDUM)
        NGOUT = NGOUTNEW
      ENDIF

      RETURN
      END SUBROUTINE MCOUT
