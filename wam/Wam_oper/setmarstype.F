      SUBROUTINE SETMARSTYPE 

! ----------------------------------------------------------------------

!**** *SETMARSTYPE* -

!     J. BIDLOT     ECMWF   JANUARY 1998  SETS MARSTYPE    

!*    PURPOSE.
!     --------
!      SETS VARIABLE MARSTYPE  

!**   INTERFACE.
!     ----------
!     *CALL SETMARSTYPE*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------
      USE YOWSTAT  , ONLY : MARSTYPE ,CDATEA   ,CDATEF   ,CDTPRO   ,
     &            IASSI    ,NENSFNB  ,NTOTENS  ,LANAONLY ,L4VTYPE

! ----------------------------------------------------------------------

      IF (IASSI .EQ. 0 ) THEN
        MARSTYPE = 'an'
      ELSE
        IF(L4VTYPE) THEN
          MARSTYPE = '4v'
        ELSE
          MARSTYPE = 'fg'
        ENDIF
      ENDIF

      IF (CDTPRO.GE.CDATEF.OR.CDATEA.EQ.CDATEF) THEN

        IF(NTOTENS.EQ.0) THEN
          MARSTYPE = 'fc'
        ELSE
          IF(NENSFNB.EQ.0) THEN
            MARSTYPE = 'cf'
          ELSE
            MARSTYPE = 'pf'
          ENDIF
        ENDIF

      ENDIF

      IF(LANAONLY) THEN
        IF (IASSI .EQ. 0 ) THEN
          MARSTYPE = 'an'
        ELSE
          MARSTYPE = 'fg'
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE SETMARSTYPE
