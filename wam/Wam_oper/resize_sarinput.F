      SUBROUTINE RESIZE_SARINPUT(KSPEC)

! ----------------------------------------------------------------------

!**** *RESIZE_SARINPUT* - CHANGES ARRAY DIMENSIONS USED IN READSARSPEC.

!     J. BIDLOT     ECMWF    MAY 1999.

!*    PURPOSE.
!     --------
!     TO CHANGE THE DIMENSIONS OF THE ALLOCATABLE ARRAYS
!     WITHOUT LOSING THE INFORMATION ALREADY CONTAINED IN THE ARRAY. 

!**   INTERFACE.
!     ----------
!     *CALL RESIZE_SARINPUT*

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
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWSARAS , ONLY : NSPEC   ,SPEC     ,LONG    ,LAT
      USE YOWTEST  , ONLY : IU06    ,ITEST

! ----------------------------------------------------------------------
      IMPLICIT NONE

      REAL, ALLOCATABLE :: R4DUM(:,:,:,:)
      REAL, ALLOCATABLE :: R2DUM(:,:)

      INTEGER :: KSPEC
!                FUTURE DIMENSION OF RESIZABLE ARRAYS  
      INTEGER :: I,ISPEC, IANG, IFRE, NEWNSPEC

!     NEW DIMENSION

      IF(KSPEC.EQ.NSPEC) THEN
        RETURN
      ELSE IF(KSPEC.GT.NSPEC) THEN
        NEWNSPEC = NSPEC + 100

        WRITE(IU06,*)' ************************************************'
        WRITE(IU06,*)' * THE DIMENSION OF ALLOCATABLES ARRAYS IN      *'
        WRITE(IU06,*)' * READSARSPEC ARE RESIZED. IF IT HAPPENS OFTEN *'
        WRITE(IU06,*)' * CHANGE THE DEFAULT VALUE IN READSARSPEC !!!  *' 
        WRITE(IU06,*)' * OLD NSPEC = ',NSPEC 
        WRITE(IU06,*)' * NEW NSPEC = ',NEWNSPEC 
        WRITE(IU06,*)' ************************************************'
        CALL FLUSH(IU06)
      ELSE
        NEWNSPEC = KSPEC
      ENDIF

!     DO THE SWITCH:

      ALLOCATE(R4DUM(0:NSPEC,NANG,NFRE,2))
      R4DUM=SPEC
      DEALLOCATE(SPEC)
      ALLOCATE(SPEC(0:NEWNSPEC,NANG,NFRE,2))
      DO I = 1,2
        DO IANG=1,NANG
          DO IFRE=1,NFRE
            DO ISPEC=0,MIN(NSPEC,NEWNSPEC)
               SPEC(ISPEC,IANG,IFRE,I)=R4DUM(ISPEC,IANG,IFRE,I)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(R4DUM)

      ALLOCATE(R2DUM(NSPEC,2))

      R2DUM=LONG
      DEALLOCATE(LONG)
      ALLOCATE(LONG(NEWNSPEC,2))
      DO I = 1,2
        DO ISPEC=1,MIN(NSPEC,NEWNSPEC)
          LONG(ISPEC,I)=R2DUM(ISPEC,I)
        ENDDO
      ENDDO

      R2DUM=LAT
      DEALLOCATE(LAT)
      ALLOCATE(LAT(NEWNSPEC,2))
      DO I = 1,2
        DO ISPEC=1,MIN(NSPEC,NEWNSPEC)
          LAT(ISPEC,I)=R2DUM(ISPEC,I)
        ENDDO
      ENDDO

      DEALLOCATE(R2DUM)

!     RESET SIZE

      NSPEC = NEWNSPEC

      IF(ITEST.GT.0) THEN
        WRITE(IU06,*)' CALL TO RESIZE_SARINPUT COMPLETE '
        CALL FLUSH(IU06)
      ENDIF

      RETURN
      END SUBROUTINE RESIZE_SARINPUT 
