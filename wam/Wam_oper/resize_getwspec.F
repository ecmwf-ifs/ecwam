      SUBROUTINE RESIZE_GETWSPEC(IMOD, MAXJPSIW, DEFAULT_DIST, NSIW1D0)

! ----------------------------------------------------------------------

!**** *RESIZE_GETWSPEC* - CHANGES ARRAY DIMENSIONS USED IN GETWSPEC 

!     J. BIDLOT     ECMWF    MAY 1999.

!*    PURPOSE.
!     --------
!     TO CHANGE THE DIMENSIONS OF THE ALLOCATABLE ARRAYS
!     WITHOUT LOSING THE INFORMATION ALREADY CONTAINED IN THE ARRAY. 

!**   INTERFACE.
!     ----------
!     *CALL RESIZE_GETWSPEC(IMOD, MAXJPSIW, DEFAULT_DIST, NSIW1D0)*

!     METHOD.
!     -------
!     USE ALLOCATABLE STATEMENT AND DUMMY ARRAY TO EXCHANGE VALUES AND
!     TO CHANGE SIZE.

!     EXTERNALS.
!     ----------
!       NONE.

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------
!     MODULES:
!     --------

      USE YOWSARAS , ONLY : NSPECW   ,JPSIW    ,NTOSIW1D ,NSIW1D   ,
     &            DIST
      USE YOWPARAM , ONLY : NIBLO
      USE YOWTEST  , ONLY : IU06     ,ITEST

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, ALLOCATABLE :: I2DUM(:,:)
      INTEGER, ALLOCATABLE :: I1DUM(:)
      REAL, ALLOCATABLE :: R2DUM(:,:)

      INTEGER :: IMOD,MAXJPSIW
!                FUTURE FIRST AND SECOND INDEX OF RESIZABLE ARRAYS
      INTEGER :: I, J, NEWNSPECW, NEWJPSIW

      INTEGER :: NSIW1D0
!             DEFAULT VALUE FOR ARRAY NSIW1D 

      REAL :: DEFAULT_DIST
!             DEFAULT VALUE FOR ARRAY DIST.


      IF(IMOD.EQ.NSPECW.AND.MAXJPSIW.EQ.JPSIW) RETURN

!     NEW FIRST DIMENSION

      IF(IMOD.EQ.NSPECW) THEN
        GOTO 2000
      ELSEIF(IMOD.GT.NSPECW) THEN
        NEWNSPECW = IMOD + 500
        IF(NEWNSPECW.GT.NIBLO) THEN
          WRITE(IU06,*)' *******************************************'
          WRITE(IU06,*)' IN RESIZE_GETWSPEC : ' 
          WRITE(IU06,*)' NEWNSPECW IS TOO BIG. PROGRAM WILL ABORT '
          WRITE(IU06,*)' *******************************************'
          CALL ABORT1
        ENDIF
        WRITE(IU06,*)' *********************************************'
        WRITE(IU06,*)' * THE DIMENSION OF ALLOCATABLES ARRAYS IN   *'
        WRITE(IU06,*)' * GETWSPEC ARE RESIZED. IF IT HAPPENS OFTEN *'
        WRITE(IU06,*)' * CHANGE THE DEFAULT VALUE IN GETWSPEC !!!  *' 
        WRITE(IU06,*)' * OLD NSPECW = ',NSPECW
        WRITE(IU06,*)' * NEW NSPEC = ',NEWNSPECW
        WRITE(IU06,*)' *********************************************'
        CALL FLUSH(IU06)
      ELSE
        NEWNSPECW = IMOD
      ENDIF

!     DO THE SWITCH:

      ALLOCATE(R2DUM(NSPECW,JPSIW))

      R2DUM=DIST
      DEALLOCATE(DIST)
      ALLOCATE(DIST(NEWNSPECW,JPSIW))
      DO J = 1,JPSIW
        DO I = 1,MIN(NSPECW,NEWNSPECW)
          DIST(I,J)=R2DUM(I,J)
        ENDDO
        DO I = NSPECW+1,NEWNSPECW
          DIST(I,J) = DEFAULT_DIST
        ENDDO
      ENDDO

      DEALLOCATE(R2DUM)

      ALLOCATE(I2DUM(NSPECW,JPSIW))

      I2DUM=NSIW1D
      DEALLOCATE(NSIW1D)
      ALLOCATE(NSIW1D(NEWNSPECW,JPSIW))
      DO J = 1,JPSIW
        DO I = 1,MIN(NSPECW,NEWNSPECW)
          NSIW1D(I,J)=I2DUM(I,J)
        ENDDO
        DO I = NSPECW+1,NEWNSPECW
          NSIW1D(I,J)=NSIW1D0
        ENDDO
      ENDDO

      DEALLOCATE(I2DUM)


      ALLOCATE(I1DUM(NSPECW))

      I1DUM=NTOSIW1D
      DEALLOCATE(NTOSIW1D)
      ALLOCATE(NTOSIW1D(NEWNSPECW))
      DO I = 1,MIN(NSPECW,NEWNSPECW)
        NTOSIW1D(I)=I1DUM(I)
      ENDDO

      DEALLOCATE(I1DUM)

!     UPDATE NEW DIMENSION

      NSPECW = NEWNSPECW


2000  CONTINUE

!     NEW SECOND DIMENSION

      IF(MAXJPSIW.EQ.JPSIW) THEN
        GOTO 3000
      ELSEIF(MAXJPSIW.GT.JPSIW) THEN
        NEWJPSIW = MAXJPSIW
        IF(NEWJPSIW.GT.NIBLO) THEN
          WRITE(IU06,*)' *******************************************'
          WRITE(IU06,*)' IN RESIZE_GETWSPEC : ' 
          WRITE(IU06,*)' NEWJPSIW IS TOO BIG. PROGRAM WILL ABORT '
          WRITE(IU06,*)' *******************************************'
          CALL ABORT1
        ENDIF
        WRITE(IU06,*)' *********************************************'
        WRITE(IU06,*)' * THE DIMENSION OF ALLOCATABLES ARRAYS IN   *'
        WRITE(IU06,*)' * GETWSPEC ARE RESIZED. IF IT HAPPENS OFTEN *'
        WRITE(IU06,*)' * CHANGE THE DEFAULT VALUE IN GETWSPEC !!!  *' 
        WRITE(IU06,*)' * OLD JPSIW= ', JPSIW
        WRITE(IU06,*)' * NEW JPSIW= ',NEWJPSIW
        WRITE(IU06,*)' *********************************************'
        CALL FLUSH(IU06)
      ELSE
        NEWJPSIW = MAXJPSIW
      ENDIF

!     DO THE SWITCH:

      ALLOCATE(R2DUM(NSPECW,JPSIW))

      R2DUM=DIST
      DEALLOCATE(DIST)
      ALLOCATE(DIST(NSPECW,NEWJPSIW))
      DO J = 1,MIN(JPSIW,NEWJPSIW)
        DO I = 1,NSPECW
          DIST(I,J)=R2DUM(I,J)
        ENDDO
      ENDDO
      DO J = JPSIW+1,NEWJPSIW
        DO I = 1,NSPECW
          DIST(I,J) = DEFAULT_DIST
        ENDDO
      ENDDO

      DEALLOCATE(R2DUM)

      ALLOCATE(I2DUM(NSPECW,JPSIW))

      I2DUM=NSIW1D
      DEALLOCATE(NSIW1D)
      ALLOCATE(NSIW1D(NSPECW,NEWJPSIW))
      DO J = 1,MIN(JPSIW,NEWJPSIW)
        DO I = 1,NSPECW
          NSIW1D(I,J)=I2DUM(I,J)
        ENDDO
      ENDDO
      DO J = JPSIW+1,NEWJPSIW
        DO I = 1,NSPECW
          NSIW1D(I,J)=NSIW1D0
        ENDDO
      ENDDO

      DEALLOCATE(I2DUM)

      ! UPDATE NEW DIMENSION

      JPSIW = NEWJPSIW

3000  CONTINUE
      IF(ITEST.GT.0) THEN
        WRITE(IU06,*)' CALL TO RESIZE_GETWSPEC COMPLETE '
        WRITE(IU06,*)' NEW JPSIW ',JPSIW 
        WRITE(IU06,*)' NEW NSPECW ', NSPECW
        WRITE(IU06,*)'' 
        CALL FLUSH(IU06)
      ENDIF

      RETURN
      END SUBROUTINE RESIZE_GETWSPEC 
