      SUBROUTINE MPBCASTGRDFLD(ISEND,ITAG,NXF,NYF,FIELD) 

!****  *MPBCASTSCFLD* - BROADCAST 2-D DATA FIELD FROM 1 PE TO ALL OTHERS

!     J. BIDLOT    ECMWF   SEPTEMBER 1997

!     PURPOSE.
!     --------

!     BROADCAST DATA CONTAINED IN 2-D ARRAY FIELD FROM PE ISEND 
!     TO ALL OTHER PE's.

!*    INTERFACE.
!     ----------

!     CALL *MPGATHERSCFLD*(ITAG,NSTART,NEND,FIELD) 

!     *ISEND*     RANK OF THE PROCESS FROM WHICH THE FIELD IS
!                 BROADCASTED
!     *ITAG*      TAG ASSOCIATED WITH A PARTICULAR CALL TO SUBROUTINE
!     *NXF*       FIRST DIMENSION OF FIELD
!     *NYF*       SECOND DIMENSION OF FIELD
!     *FIELD*     INPUT/OUTPUT ARRAY CONTAINING THE FIELD

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

      REAL,DIMENSION(NXF,NYF) :: FIELD

      REAL,ALLOCATABLE :: ZCOMBUF(:)
      MPLENGTH=NXF*NYF

      IF((ISEND.EQ.0).OR.(NPROC.EQ.1)) THEN
        RETURN
      ELSE
!     1.1 SEND TO ALL PROCESSORS OTHER THAN ISEND
!         ----------------------------------------
        ALLOCATE(ZCOMBUF(MPLENGTH))

        IF(IRANK.EQ.ISEND) THEN
          KCOUNT=0
          DO J=1,NYF
            DO I=1,NXF
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=FIELD(I,J)
            ENDDO
          ENDDO
        ENDIF

        CALL MPL_BROADCAST(ZCOMBUF(1:MPLENGTH),KROOT=ISEND,
     &   KTAG=ITAG,CDSTRING='MPBCASTGRDFLD:')

        IF(IRANK.NE.ISEND) THEN
          KCOUNT=0
          DO J=1,NYF
            DO I=1,NXF
              KCOUNT=KCOUNT+1
              FIELD(I,J)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO
        ENDIF

        DEALLOCATE(ZCOMBUF)

      ENDIF

      RETURN
      END SUBROUTINE MPBCASTGRDFLD
