      SUBROUTINE MPBCASTINTFLD(ISEND,ITAG,NXF,NYF,IFIELD) 

!****  *MPBCASTINTFLD* - BROADCAST 2-D INTEGER DATA FIELD FROM 1 PE 
!                        TO ALL OTHERS.

!     J. BIDLOT    ECMWF   FEBRUARY 1999. 

!     PURPOSE.
!     --------

!     BROADCAST DATA CONTAINED IN 2-D ARRAY IFIELD FROM PE ISEND 
!     TO ALL OTHER PE's.

!*    INTERFACE.
!     ----------

!     CALL *MPBCASTINTFLD*(ITAG,NSTART,NEND,IFIELD) 

!     *ISEND*     RANK OF THE PROCESS FROM WHICH THE FIELD IS
!                 BROADCASTED
!     *ITAG*      TAG ASSOCIATED WITH A PARTICULAR CALL TO SUBROUTINE
!     *NXF*       FIRST DIMENSION OF IFIELD
!     *NYF*       SECOND DIMENSION OF IFIELD
!     *IFIELD*    INPUT/OUTPUT ARRAY CONTAINING THE FIELD

!     METHOD.
!     -------
!     MPL BROADCAST OF ARRAY FIELD FROM ISEND
!     AND RECEIVE OF THE CONTRIBUTION ON THE OTHER PE's.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_BROADCAST

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------
 
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE MPL_MODULE

!----------------------------------------------------------------------

      INTEGER,DIMENSION(NXF,NYF) :: IFIELD

      INTEGER,ALLOCATABLE :: ICOMBUF(:)
      MPLENGTH=NXF*NYF

      IF((ISEND.EQ.0).OR.(NPROC.EQ.1)) THEN
        RETURN
      ELSE
        ALLOCATE(ICOMBUF(MPLENGTH))

        IF(IRANK.EQ.ISEND) THEN
          KCOUNT=0
          DO J=1,NYF
            DO I=1,NXF
              KCOUNT=KCOUNT+1
              ICOMBUF(KCOUNT)=IFIELD(I,J)
            ENDDO
          ENDDO
        ENDIF

        CALL MPL_BROADCAST(ICOMBUF(1:MPLENGTH),KROOT=ISEND,
     &   KTAG=ITAG,CDSTRING='MPBCASTINTFLD:')


        IF(IRANK.NE.ISEND) THEN
          KCOUNT=0
          DO J=1,NYF
            DO I=1,NXF
              KCOUNT=KCOUNT+1
              IFIELD(I,J)=ICOMBUF(KCOUNT)
            ENDDO
          ENDDO
        ENDIF

        DEALLOCATE(ICOMBUF)

      ENDIF

      RETURN
      END SUBROUTINE MPBCASTINTFLD
