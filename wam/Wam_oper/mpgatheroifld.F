      SUBROUTINE MPGATHEROIFLD(ISEND,ITAG,FIELD) 

!****  *MPGATHEROIFLD* - GATHER SCALAR GRID DATA FIELD
!****                    ONTO ALL PROCESSES 

!     J. BIDLOT    ECMWF   JUNE 1996  MESSAGE PASSING

!     PURPOSE.
!     --------

!     BROADCAST SCALAR GRID DATA CONTAINED IN ARRAY FIELD STORED ON
!     PROCESS ISEND  ONTO ALL PROCESSES

!*    INTERFACE.
!     ----------

!     CALL *MPGATHEROIFLD*(ISEND,ITAG,FIELD) 

!     *ISEND*     PROCESS RANK BROADCATING THE GRID FIELD
!     *ITAG*      TAG ASSOCIATED WITH A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS
!                 THE USE OF A SEPARATE TAG FOR EACH SPECIFIC CALL
!                 ENSURES THAT NO CONFLICT ARISES BETWEEN DIFFERENT
!                 MESSAGES
!     *FIELD*     INPUT/OUTPUT ARRAY CONTAINING THE FIRST FIELD

!     METHOD.
!     -------
!     MPL BROADCAST OF ARRAY FIELD TO ALL PROCESSES

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_BROADCAST
!         MPL_RECV
!         MPL_ABORT 

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE MPL_MODULE

!----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL,DIMENSION(NGX,NGY) :: FIELD
!----------------------------------------------------------------------

      REAL,DIMENSION(NGX*NGY):: ZCOMBUF

      KCOUNT=NGX*NGY 

      IF(NPROC.EQ.1) THEN
        RETURN

      ELSE

        IF(IRANK.EQ.ISEND) THEN
          IJ=0
          DO J=1,NGY
            DO I=1,NGX
              IJ=IJ+1
              ZCOMBUF(IJ)=FIELD(I,J)
            ENDDO
          ENDDO
        ENDIF

        CALL MPL_BROADCAST(ZCOMBUF(1:KCOUNT),KROOT=ISEND,KTAG=ITAG,
     &   CDSTRING='MPGATHEROIFLD:')

        IF(IRANK.NE.ISEND) THEN
          IJ=0
          DO J=1,NGY
            DO I=1,NGX
              IJ=IJ+1
              FIELD(I,J)=ZCOMBUF(IJ)
            ENDDO
          ENDDO
        ENDIF

      ENDIF

      RETURN
      END SUBROUTINE MPGATHEROIFLD
