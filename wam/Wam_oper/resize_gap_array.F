      SUBROUTINE RESIZE_GAP_ARRAY

! ----------------------------------------------------------------------

!**** *RESIZE_GAP_ARRAY* - CHANGES ARRAY DIMENSIONS USED IN MAKEFRAMES.

!     J. BIDLOT     ECMWF    MAY 1999.

!*    PURPOSE.
!     --------
!     TO CHANGE THE DIMENSIONS OF THE ALLOCATABLE ARRAYS
!     WITHOUT LOSING THE INFORMATION ALREADY CONTAINED IN THE ARRAY. 

!**   INTERFACE.
!     ----------
!     *CALL RESIZE_GAP_ARRAY*

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
      USE YOWGAP   , ONLY : MGAP     ,LANG     ,HANG     ,LFRE     ,
     &            HFRE     ,EQUIGAP
      USE YOWPARAM , ONLY : NANG    ,NFRE
      USE YOWSARAS , ONLY : NSPEC
      USE YOWTEST  , ONLY : IU06    ,ITEST

! ----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, ALLOCATABLE :: I2DUM(:,:)
      LOGICAL, ALLOCATABLE :: L3DUM(:,:,:)

      INTEGER :: ISPEC, IGAP, JGAP, NEWMGAP, DEL


!     NEW DIMENSION

      DEL = 10
      NEWMGAP = MGAP + DEL 

      WRITE(IU06,*)' ************************************************'
      WRITE(IU06,*)' * THE DIMENSION OF ALLOCATABLES ARRAYS IN      *'
      WRITE(IU06,*)' * MAKEFRAMES ARE RESIZED. IF IT HAPPENS OFTEN  *'
      WRITE(IU06,*)' * CHANGE THE DEFAULT VALUE IN READSARSPEC !!!  *'
      WRITE(IU06,*)' * OLD MGAP= ',MGAP
      WRITE(IU06,*)' * NEW MGAP= ',NEWMGAP
      WRITE(IU06,*)' ************************************************'
      CALL FLUSH(IU06)

      IF(NEWMGAP.GT.NFRE*NANG+DEL) THEN
        WRITE(IU06,*)' ****************************************'
        WRITE(IU06,*)' PROBLEM WITH RESIZING GAP ARRAYS ' 
        WRITE(IU06,*)' THE REQUESTED SIZE IS LARGER THEN THE'
        WRITE(IU06,*)' NUMBER OF SPECTRAL COMPONENTS.' 
        WRITE(IU06,*)' THE PROGRAM WILL ABORT.'
        WRITE(IU06,*)' ****************************************'
        CALL ABORT1
      ENDIF

!     DO THE SWITCH:

      ALLOCATE(I2DUM(NSPEC,MGAP))

      I2DUM=LANG
      DEALLOCATE(LANG)
      ALLOCATE(LANG(NSPEC,NEWMGAP))
      DO IGAP = 1,MGAP
        DO ISPEC=1,NSPEC
          LANG(ISPEC,IGAP)=I2DUM(ISPEC,IGAP)
        ENDDO
      ENDDO

      I2DUM=HANG
      DEALLOCATE(HANG)
      ALLOCATE(HANG(NSPEC,NEWMGAP))
      DO IGAP = 1,MGAP
        DO ISPEC=1,NSPEC
          HANG(ISPEC,IGAP)=I2DUM(ISPEC,IGAP)
        ENDDO
      ENDDO

      I2DUM=LFRE
      DEALLOCATE(LFRE)
      ALLOCATE(LFRE(NSPEC,NEWMGAP))
      DO IGAP = 1,MGAP
        DO ISPEC=1,NSPEC
          LFRE(ISPEC,IGAP)=I2DUM(ISPEC,IGAP)
        ENDDO
      ENDDO

      I2DUM=HFRE
      DEALLOCATE(HFRE)
      ALLOCATE(HFRE(NSPEC,NEWMGAP))
      DO IGAP = 1,MGAP
        DO ISPEC=1,NSPEC
          HFRE(ISPEC,IGAP)=I2DUM(ISPEC,IGAP)
        ENDDO
      ENDDO

      DEALLOCATE(I2DUM)


      ALLOCATE(L3DUM(NSPEC,MGAP,MGAP))

      L3DUM=EQUIGAP
      DEALLOCATE(EQUIGAP)
      ALLOCATE(EQUIGAP(NSPEC,NEWMGAP,NEWMGAP))
      EQUIGAP=.FALSE.
      DO JGAP = 1,MGAP
        DO IGAP = 1,MGAP
          DO ISPEC=1,NSPEC
            EQUIGAP(ISPEC,IGAP,JGAP)=L3DUM(ISPEC,IGAP,JGAP)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(L3DUM)


!     RESET SIZE

      MGAP = NEWMGAP

      IF(ITEST.GT.0) THEN
        WRITE(IU06,*)' CALL TO RESIZE_GAP_ARRAY COMPLETE '
        CALL FLUSH(IU06)
      ENDIF

      RETURN
      END SUBROUTINE RESIZE_GAP_ARRAY
