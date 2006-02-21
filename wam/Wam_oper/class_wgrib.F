      INTEGER FUNCTION CLASS_WGRIB(CDCLASS,IU06)

!**** *CLASS_WGRIB*  - CLASS code to GRIB code
!     March 1999:
!     This is an adapted version of the ifs function
!     CLASS_GRIB(CDCLASS) such that WAM can be built without
!     the use of the IFS library.	
!     Please update when CLASS_GRIB(CDCLASS) is modified.

!     Purpose.  To provide GRIB header code for given CLASS code.
!     --------


!**   Interface.
!     ----------
!       IGRIB=CLASS_WGRIB(CDCLASS)

!       Explicit arguments : CDCLASS - CLASS code - 2 characters (input)
!       --------------------           can be either upper or lower case
!                            CLASS_WGRIB returns value is the integer
!                                        GRIB code
!                                for ECMWF local GRIB usage of KSEC1(38)
!                            IU06 - PRINTER OUTPUT UNIT.


!        Implicit arguments :  
!        --------------------

!     Method.
!     -------
!        Table look up     

!     Externals.   
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        David Dent  *ECMWF*

!     Modifications.
!     --------------
!       Original : 98-08-10
!       D.Salmond : 03-10-02  Change to include 'test'
!     ------------------------------------------------------------------
!
      PARAMETER (JPCLASS_LIST=28)
      CHARACTER*2 CDCLASS
      CHARACTER*2 CLASSES(56)
      DATA CLASSES/'OD','RD','ER','CS','E4','AT','BE','HR','DK','FI',
     &             'FR','DE','GR','HU','IS','IE','IT','NL','NO','PT',
     &             'SI','ES','SE','CH','TR','UK','TE','EI',
     &             'od','rd','er','cs','e4','at','be','hr','dk','fi',
     &             'fr','de','gr','hu','is','ie','it','nl','no','pt',
     &             'si','es','se','ch','tr','uk','te','ei'/
      INTEGER GRIBCODES(27)
      DATA GRIBCODES/ 1,   2,   3,   4,   5, 100, 101, 102, 103, 104,
     &              105, 106, 107, 108, 109, 110, 111, 112, 113, 114,
     &              115, 116, 117, 118, 119, 120, 99, 14/

!     ------------------------------------------------------------------

!*       1.    SEARCH CLASS LIST.
!              -----------------

      ICODE=0
      DO JL=1,JPCLASS_LIST

        IF(CDCLASS.EQ.CLASSES(JL).OR.
     &     CDCLASS.EQ.CLASSES(JL+JPCLASS_LIST)) THEN
          ICODE=GRIBCODES(JL)
          EXIT
        ENDIF
      ENDDO
      IF(ICODE.NE.0) THEN
        CLASS_WGRIB=ICODE
      ELSE
        WRITE(IU06,*) ' -------------------------'
        WRITE(IU06,*) ' CLASS_WGRIB: ' 
        WRITE(IU06,*) ' ERROR DETERMINING CLASS  '
        WRITE(IU06,*) ' THE CODE FOUND IS : ', CDCLASS 
        WRITE(IU06,*) ' -------------------------'
        CALL ABORT1
      ENDIF

!     ------------------------------------------------------------------

      RETURN
      END FUNCTION CLASS_WGRIB
