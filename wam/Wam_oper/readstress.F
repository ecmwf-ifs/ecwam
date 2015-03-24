      SUBROUTINE READSTRESS(NREAL, RFIELD, FILENAME)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     READS THE RESTART WIND AND STRESS FIELDS

!**   INTERFACE.
!     ----------
!     CALL *READSTRESS*(NREAL, RFIELD, FILENAME)  
!     *NREAL*    NUMBER OF FIELDS IN RFIELD
!     *RFIELD*     REAL RESTART FIELDS
!     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE

!     METHOD.
!     -------
!     READS ARRAYS FROM FILE. FORTRAN READS

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : NIBLO    ,LL1D
      USE YOWSTAT  , ONLY : CDTPRO
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWIND  , ONLY : CDAWIFL  ,CDATEWO  ,CDATEFL
      USE YOWUNPOOL, ONLY : LLUNSTR
      USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: LFILE, IUNIT
      INTEGER :: IFLD, IJ
      INTEGER :: NREAL
      INTEGER :: I_GET_UNIT

      REAL :: ZHOOK_HANDLE
      REAL,DIMENSION(NIBLO) :: RFIELD_G 
      REAL,DIMENSION(NIBLO,NREAL) :: RFIELD 

      CHARACTER(LEN=120) :: FILENAME

      LOGICAL :: LLEXIST

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('READSTRESS',0,ZHOOK_HANDLE)
#endif

      IF (ITEST.GE.2) THEN
          WRITE(IU06,*) ' SUB. READSTRESS: ABOUT TO INPUT RESTART FILES'
          CALL FLUSH (IU06)
      ENDIF

      LFILE=0
      IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
      LLEXIST=.FALSE.
      INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
      IF(.NOT. LLEXIST) THEN
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (IU06,*) '*  IN READSTRESS :                  *'
        WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILENAME
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        WRITE (*,*) '*************************************'
        WRITE (*,*) '*                                   *'
        WRITE (*,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (*,*) '*  IN READSTRESS :                  *'
        WRITE (*,*) '*  COULD NOT FIND FILE ',FILENAME
        WRITE (*,*) '*                                   *'
        WRITE (*,*) '*************************************'
        CALL ABORT1
      ENDIF
      IUNIT=I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0)

      READ(IUNIT) CDTPRO,CDATEWO,CDAWIFL,CDATEFL

      IF(LLUNSTR) THEN
        DO IFLD=1,NREAL
          READ(IUNIT) RFIELD_G 
          CALL UNBLKRORD(-1,1,1,1,1,RFIELD(1,IFLD),RFIELD_G(1))
        ENDDO
      ELSEIF(LL1D .OR. NPROC.EQ.1 ) THEN
        DO IFLD=1,NREAL
          READ(IUNIT) RFIELD(:,IFLD) 
        ENDDO
      ELSE
!     WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE
!     RE-LABELLED BUT THE BINARY INPUT FILES ARE IN THE OLD MAPPING
        DO IFLD=1,NREAL
          READ(IUNIT)(RFIELD(IJ2NEWIJ(IJ),IFLD),IJ=1,NIBLO)
        ENDDO
      ENDIF

      CLOSE(IUNIT)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) ' SUB. READSTRESS: INPUT OF RESTART FILES DONE'
        CALL FLUSH(IU06)
      ENDIF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('READSTRESS',1,ZHOOK_HANDLE)
#endif
      RETURN

      END SUBROUTINE READSTRESS
