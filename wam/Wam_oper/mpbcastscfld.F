      SUBROUTINE MPBCASTSCFLD(NSTART,NEND,FIELD)

!****  *MPBCASTSCFLD* - BROADCAST SCALAR BLOCK DATA FIELD ONTO ALL
!****                   PROCESS

!     J. BIDLOT    ECMWF   SEPTEMBER 1996  MESSAGE PASSING

!     PURPOSE.
!     --------

!     BROADCAST SCALAR BLOCK DATA CONTAINED IN ARRAY FIELD DISTRIBUTED
!     ACROSS THE DIFFERENT PROCESSES ONTO ALL PE's

!*    INTERFACE.
!     ----------

!     CALL *MPBCASTSCFLD*(NSTART,NEND,FIELD)

!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *FIELD*     INPUT/OUTPUT ARRAY CONTAINING THE FIELD

!     METHOD.
!     -------
!     MPL ALLGATHER OF ARRAY FIELD TO ALL PE's

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_ALLGATHER

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE YOWPARAM , ONLY : NIBLO
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE MPL_MODULE
!----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER,DIMENSION(NPROC), INTENT(IN) :: NSTART,NEND
      REAL,DIMENSION(NIBLO), INTENT(INOUT) :: FIELD

      INTEGER  :: IP
      INTEGER,DIMENSION(NPROC) :: KRECVCOUNTS
      REAL,DIMENSION(NIBLO) :: ZTMP
!----------------------------------------------------------------------


      IF(NPROC.EQ.1) THEN
        RETURN
      ELSE

        DO IP=1,NPROC
          KRECVCOUNTS(IP)=NEND(IP) - NSTART(IP)+1
        ENDDO

        CALL GSTATS(616,0)
        CALL MPL_ALLGATHERV(FIELD(NSTART(IRANK):NEND(IRANK)),
     &                      ZTMP,KRECVCOUNTS,
     &                      CDSTRING='MPBCASTSCFLD:')
        FIELD(:) = ZTMP(:)
        CALL GSTATS(616,1)

      ENDIF

      RETURN
      END SUBROUTINE MPBCASTSCFLD
