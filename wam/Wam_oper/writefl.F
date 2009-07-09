      SUBROUTINE WRITEFL(FL, KINF, KSUP, MINF, MSUP,
     &                   FILENAME,IUNIT,LPBOPEN,LPBCLOSE)

! ---------------------------------------------------------------------
!     J. BIDLOT    ECMWF      MARCH 1997

!*    PURPOSE.
!     --------
!     WRITES ARRAY FL TO FILE.

!**   INTERFACE.
!     ----------
!     CALL *WRITEFL(FL, KINF, KSUP, MINF, MSUP,
!    &              FILENAME,IUNIT,LPBOPEN,LPBCLOSE)
!     *FL*        ARRAY TO BE WRITTEN TO FILE
!     *KINF*      SECOND LOWER DIMENSION BOUND OF FL
!     *KSUP*      SECOND UPPER DIMENSION BOUND OF FL
!     *MINF*      THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*      THIRD UPPER DIMENSION BOUND OF FL
!     *FILENAME*  FILENAME (INCLUDING PATH) OF TARGET FILE
!     *IUNIT*     PBIO UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LPBOPEN*   LOGICAL, TRUE IF PBOPEN HAS TO BE CALLED
!     *LPBCLOSE*  LOGICAL, TRUE IF PBCLOSE HAS TO BE CALLED

!     METHOD.
!     -------
!     WRITES ARRAY FL TO FILE. A FORTRAN WRITE, OR A PBIO
!     WRITE IS USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED
!     BY LPBIOOUT

!     IF PBIO IS USED THEN FL WILL BE SAVED IN FILENAME USING
!     PBWRITE ON CONDITION THAT FILENAME WAS NOT OPENED BEFORE AND
!     CONNECTED TO IUNIT.
!     IT WILL BE APPENDED TO THE BOTTOM OF THE CONTENT OF FILENAME IF
!     FILENAME WAS ALREADY OPENED AND CONNECTED TO UNIT IUNIT.

!     IF FORTRAN WRITE IS USED THEN IN THE CASE OF MESSAGE PASSING
!     FL IS WRITTEN AS UNFORMATTED BINARY TO FILENAME CONNECTED
!     TO UNIT IUNIT. IN THE NON MESSAGE PASSING CASE, FL IS WRITTEN AS
!     UNFORMATTED BINARY TO FORTRAN UNIT SPECIFIED BY IUNIT.


!     EXTERNALS.
!     ----------
!     *PBOPEN*
!     *PBWRITE*
!     *PBCLOSE*
!     *ABORT1*

!     REFERENCE.
!     ----------
!     NONE

! ----------------------------------------------------------------------

      USE YOWGRID  , ONLY : IGL
      USE YOWMESPAS, ONLY : LMESSPASS,LPBIOOUT
      USE YOWMPP   , ONLY : NPROC, NPRECR
      USE YOWPARAM , ONLY : NIBLO    ,CLDOMAIN ,LL1D
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      CHARACTER*  1 CLMODE
      DATA          CLMODE   / "w" /
      CHARACTER*120 FILENAME
      LOGICAL LPBOPEN,LPBCLOSE

! ----------------------------------------------------------------------

      REAL,DIMENSION(0:NIBLO,KINF:KSUP,MINF:MSUP) :: FL
      REAL,ALLOCATABLE :: RFL(:,:,:)

! ----------------------------------------------------------------------

      IF(LPBIOOUT) THEN
        IF(LPBOPEN) THEN
          CALL PBOPEN(IUNIT,FILENAME,CLMODE,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '******************************************'
            WRITE (IU06,*) '*                                        *'
            WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
            WRITE (IU06,*) '*   IN WRITEFL                           *'
            IF(KRET.EQ.-1)
     &       WRITE (IU06,*) 'COULD NOT OPEN FILE ',FILENAME
            IF(KRET.EQ.-2) WRITE (IU06,*) 'INVALID FILENAME ',FILENAME
            IF(KRET.EQ.-3) WRITE (IU06,*) 'INVALID OPEN MODE SPECIFIED'
            WRITE (IU06,*) '*                                        *'
            WRITE (IU06,*) '******************************************'
            CALL ABORT1
          ENDIF
        ENDIF
        KOUNT = (NIBLO+1)*(KSUP-KINF+1)*(MSUP-MINF+1)*NPRECR
        IF(LL1D .OR. NPROC.EQ.1) THEN
          CALL PBWRITE(IUNIT,FL,KOUNT,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '*****************************************'
            WRITE (IU06,*) '*                                       *'
            WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBWRITE FL *'
            WRITE (IU06,*) '*    IN WRITEFL                         *'
            WRITE (IU06,*) '*    FILE ',FILENAME
            WRITE (IU06,*) '*                                       *'
            WRITE (IU06,*) '*****************************************'
            CALL ABORT1
          ENDIF
        ELSE
!         WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!         BUT THE BINARY INPUT FILES SHOULD BE IN THE OLD MAPPING
          ALLOCATE(RFL(0:NIBLO,KINF:KSUP,MINF:MSUP))
          DO M=MINF,MSUP
            DO K=KINF,KSUP
              DO IJ=0,NIBLO
                RFL(IJ,K,M)=FL(IJ2NEWIJ(IJ),K,M)
              ENDDO
            ENDDO
          ENDDO
          CALL PBWRITE(IUNIT,RFL,KOUNT,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '*****************************************'
            WRITE (IU06,*) '*                                       *'
            WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBWRITE RFL*'
            WRITE (IU06,*) '*    IN WRITEFL                         *'
            WRITE (IU06,*) '*    FILE ',FILENAME
            WRITE (IU06,*) '*                                       *'
            WRITE (IU06,*) '*****************************************'
            CALL ABORT1
          ENDIF
          DEALLOCATE(RFL)
        ENDIF
        IF(LPBCLOSE) THEN
          CALL PBCLOSE(IUNIT,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '************************************'
            WRITE (IU06,*) '*                                  *'
            WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
            WRITE (IU06,*) '* IN WRITEFL                       *'
            WRITE (IU06,*) '* FILE ',FILENAME
            WRITE (IU06,*) '*                                  *'
            WRITE (IU06,*) '************************************'
            CALL ABORT1
          ENDIF
        ENDIF
      ELSE
        LFILE=0
        IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
        IF(LPBOPEN) THEN
          IUNIT=I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'w', 'u', 0)
        ELSE
          IUNIT=I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'a', 'u', 0)
        ENDIF

        IF(LL1D .OR. NPROC.EQ.1) THEN
          WRITE(IUNIT) (((FL(J1,J2,J3),
     &                    J1=0,NIBLO),
     &                    J2=KINF,KSUP),
     &                    J3=MINF,MSUP)
        ELSE
!       WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!       BUT THE BINARY INPUT FILES SHOULD BE IN THE OLD MAPPING
          WRITE(IUNIT) (((FL(IJ2NEWIJ(J1),J2,J3),
     &                    J1=0,NIBLO),
     &                    J2=KINF,KSUP),
     &                    J3=MINF,MSUP)

        ENDIF

        IF (ITEST.GE.2)
     &    WRITE(IU06,*) ' SUB. WRITEFL: OUTPUT OF FL DONE',
     &                   KINF,KSUP,MINF,MSUP

        CLOSE(IUNIT)
      ENDIF

      RETURN
      END SUBROUTINE WRITEFL
