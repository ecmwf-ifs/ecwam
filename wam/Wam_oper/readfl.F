      SUBROUTINE READFL(FL, KINF, KSUP, MINF, MSUP,
     &                  FILENAME, IUNIT, LPBOPEN, LPBCLOSE)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997

!*    PURPOSE.
!     --------
!     READS ARRAY FL

!**   INTERFACE.
!     ----------
!     CALL *READFL*(FL, KINF, KSUP, MINF, MSUP,
!                   FILENAME, IUNIT, LPBOPEN, LPBCLOSE)
!     *FL*       ARRAY TO BE WRITTEN TO FILE
!     *KINF*     SECOND LOWER DIMENSION BOUND OF FL
!     *KSUP*     SECOND UPPER DIMENSION BOUND OF FL
!     *MINF*     THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*     THIRD UPPER DIMENSION BOUND OF FL
!     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE
!     *IUNIT*    FILE UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LPBOPEN*  LOGICAL, TRUE IF PBOPEN HAS TO BE CALLED
!     *LPBCLOSE* LOGICAL, TRUE IF PBCLOSE HAS TO BE CALLED

!     METHOD.
!     -------
!     READS ARRAY FL FROM FILE. A FORTRAN READ OR PBIO READ ARE
!     USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED BY LPBIOIN

!     IF PBIO IS USED, THEN PBREAD IS CALLED TO READ (MSUP-MINF+1) 
!     FREQUENCY CONTRIBUTION TO FL.


!     EXTERNALS.
!     ----------
!     *PBOPEN*
!     *PBREAD*
!     *PBCLOSE*
!     *ABORT1*

!     REFERENCE.
!     ----------
!     NONE
! ----------------------------------------------------------------------

      USE YOWMESPAS, ONLY : LPBIOIN
      USE YOWMPP   , ONLY : NPROC    ,NPRECR
      USE YOWPARAM , ONLY : NIBLO     ,LL1D
      USE YOWSPEC  , ONLY : IJ2NEWIJ
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------

      CHARACTER*1 MODE
      CHARACTER*120 FILENAME
      LOGICAL LPBOPEN,LPBCLOSE

      REAL,DIMENSION(0:NIBLO,KINF:KSUP,MINF:MSUP) :: FL
      REAL,ALLOCATABLE :: RFL(:,:,:)
! ----------------------------------------------------------------------

      IF(LPBIOIN) THEN

        IF(LPBOPEN) THEN
          MODE = 'r'
          CALL PBOPEN(IUNIT,FILENAME,MODE,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '******************************************'
            WRITE (IU06,*) '*                                        *'
            WRITE (IU06,*) '*   ERROR FOLLOWING CALL TO PBOPEN       *'
            WRITE (IU06,*) '*   IN READFL                            *' 
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
          CALL PBREAD(IUNIT,FL,KOUNT,KRET)
        ELSE
!         WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!         BUT THE BINARY INPUT FILES ARE IN THE OLD MAPPING
          ALLOCATE(RFL(0:NIBLO,KINF:KSUP,MINF:MSUP))
          CALL PBREAD(IUNIT,RFL,KOUNT,KRET)
          DO M=MINF,MSUP
            DO K=KINF,KSUP
              DO IJ=0,NIBLO
                FL(IJ2NEWIJ(IJ),K,M)=RFL(IJ,K,M)
              ENDDO
            ENDDO
          ENDDO
          DEALLOCATE(RFL)
        ENDIF
        IF(KRET.LT.0) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBREAD     *'
          WRITE (IU06,*) '*    IN READFL                          *' 
          WRITE (IU06,*) '*    FILE ',FILENAME 
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    KRET = ', KRET 
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ELSEIF(KRET.NE.KOUNT) THEN
          WRITE (IU06,*) '*****************************************'
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*    ERROR FOLLOWING CALL TO PBREAD     *'
          WRITE (IU06,*) '*    IN READFL                          *' 
          WRITE (IU06,*) '*    FILE ',FILENAME 
          WRITE (IU06,*) '*    NUMBER OF BYTES READ : ',KRET
          WRITE (IU06,*) '*    IS NOT EQUAL TO THE NUMBER OF BYTES*'
          WRITE (IU06,*) '*    REQUESTED : ',KOUNT 
          WRITE (IU06,*) '*                                       *'
          WRITE (IU06,*) '*****************************************'
          CALL ABORT1
        ENDIF

        IF(LPBCLOSE) THEN
          CALL PBCLOSE(IUNIT,KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '************************************'
            WRITE (IU06,*) '*                                  *'
            WRITE (IU06,*) '* ERROR FOLLOWING CALL TO PBCLOSE   '
            WRITE (IU06,*) '* IN READFL                        *'
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
          IUNIT=I_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0)
        ENDIF

        IF(LL1D .OR. NPROC.EQ.1) THEN
          READ(IUNIT) (((FL(J1,J2,J3),
     &                    J1=0,NIBLO),
     &                    J2=KINF,KSUP),
     &                    J3=MINF,MSUP)
        ELSE
!         WHEN 2-D DECOMPOSITION IS USED THEN THE INDEXES IJ ARE RE-LABELLED
!         BUT THE BINARY INPUT FILES ARE IN THE OLD MAPPING
          READ(IUNIT) (((FL(IJ2NEWIJ(J1),J2,J3),
     &                    J1=0,NIBLO),
     &                    J2=KINF,KSUP),
     &                    J3=MINF,MSUP)

        ENDIF
        IF(LPBCLOSE) THEN
          CLOSE(IUNIT)
        ENDIF

      ENDIF

      IF (ITEST.GE.2)
     &            WRITE(IU06,*) ' SUB. READFL: INPUT OF FL DONE'

      RETURN
      END SUBROUTINE READFL
