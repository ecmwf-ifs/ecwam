      LOGICAL FUNCTION LLALLOC_GRDFLD(INDEX)

! ----------------------------------------------------------------------

!     J. BIDLOT     ECMWF    MARCH 2000

!*    PURPOSE.
!     --------

!     RETURNS WHETHER OR NOT AN ARRAY USED IN THE GRIDDED OUTPUT
!     HAS T0 BE ALLOCATED. IT DOES NOT CHECK WHETHER THE ARRAY
!     IS ALREADY ALLOCATED.

!**   INTERFACE.
!     ----------
!     *LLALLOC_GRDFLD(INDEX)*

!     *INDEX*   INTEGER  INDEX OF THE OUTPUT FIELD IN THE ARRAYS
!                        CONTROLLING THE DIFFERENT FORM OF OUTPUTS.

! ----------------------------------------------------------------------

      USE YOWCOUT  , ONLY : FFLAG    ,PFLAG    ,GFLAG    ,UFLAG     ,
     &            LFDB     ,IPFGTBL
      USE YOWSTAT  , ONLY : IASSI 
      USE YOWMPP   , ONLY : IRANK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: INDEX

      LLALLOC_GRDFLD=.FALSE.

      IF((IASSI.EQ.1 .AND. UFLAG(INDEX)) .OR.
     &   FFLAG(INDEX) .OR.
     &   PFLAG(INDEX) .OR.
     &   (.NOT.LFDB .AND. GFLAG(INDEX)) ) THEN
        LLALLOC_GRDFLD=.TRUE.
      ELSE
        IF(IRANK.EQ.IPFGTBL(INDEX)) LLALLOC_GRDFLD=.TRUE.
      ENDIF

      RETURN
      END FUNCTION LLALLOC_GRDFLD
