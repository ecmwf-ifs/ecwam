      SUBROUTINE READFL(FL,MINF,MSUP,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)

! ----------------------------------------------------------------------
!     J. BIDLOT    ECMWF      SEPTEMBER 1997

!*    PURPOSE.
!     --------
!     READS ARRAY FL

!**   INTERFACE.
!     ----------
!     CALL *READFL*(FL,MINF,MSUP,FILENAME,IUNIT,LPBOPEN,LPBCLOSE)
!     *FL*       ARRAY TO BE WRITTEN TO FILE
!     *MINF*     THIRD LOWER DIMENSION BOUND OF FL
!     *MSUP*     THIRD UPPER DIMENSION BOUND OF FL
!     note the first and second dimension of FL are given 
!     by 0:NIBLO and NANG
!     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE
!     *IUNIT*    PBIO UNIT (ONLY ACTIVE IF PBIO OUTPUT)
!     *LPBOPEN*  LOGICAL, TRUE IF PBOPEN HAS TO BE CALLED
!     *LPBCLOSE* LOGICAL, TRUE IF PBCLOSE HAS TO BE CALLED
!     *FILENAME* FILENAME (INCLUDING PATH) OF INPUT FILE

!     METHOD.
!     -------
!     READS ARRAY FL FROM FILE. A FORTRAN READ OR PBIO READ ARE
!     USED DEPENDING ON THE CHOICE OF OUTPUT MODE AS SELECTED BY LPBIOIN

!     IF PBIO IS USED, THEN PBREAD IS CALLED TO READ (MSUP-MINF+1) 
!     FREQUENCY CONTRIBUTION TO FL.
!     IF FORTRAN READ IS USED THEN FILNAME IS READ AS UNFORMATTED BINARY
!     MINF SHOULD BE SET TO 1 AND MSUP TO NFRE


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
      USE YOWMPP   , ONLY : NPRECR
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWUNIT  , ONLY : IU12

! ----------------------------------------------------------------------

      CHARACTER*1 MODE
      CHARACTER*80 FILENAME
      LOGICAL LPBOPEN,LPBCLOSE

      REAL,DIMENSION(0:NIBLO,NANG,MINF:MSUP) :: FL
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

        KOUNT = (NIBLO+1)*NANG*(MSUP-MINF+1)*NPRECR
        CALL PBREAD(IUNIT,FL,KOUNT,KRET)
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

        IF(MINF.NE.1.OR.MSUP.NE.NFRE) THEN
          WRITE (IU06,*) '************************************'
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '* ERROR IN READFL                  *' 
          WRITE (IU06,*) '* MINF MUST BE = 1                 *'
          WRITE (IU06,*) '* MSUP MUST BE = NFRE              *'
          WRITE (IU06,*) '* MINF WAS ', MINF
          WRITE (IU06,*) '* MSUP WAS ', MSUP
          WRITE (IU06,*) '* NFRE WAS ', NFRE 
          WRITE (IU06,*) '*                                  *'
          WRITE (IU06,*) '************************************'
          CALL ABORT1
        ENDIF

        READ(IU12) FL
        CLOSE(IU12)

      ENDIF

      IF (ITEST.GE.2)
     &            WRITE(IU06,*) ' SUB. READFL: INPUT OF FL DONE'

      RETURN
      END SUBROUTINE READFL
