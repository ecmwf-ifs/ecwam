       SUBROUTINE MERGESARCOR(IJSG, IJLG, FL1, NW1D)

!---------------------------------------------------------------------

!     PURPOSE
!     -------

!     MERGE ANALYSED SPECTRA INTO BLOCKED GRID SPECTRA.
!     MERGE ANALYSED WINDS INTO WIND BLOCKED ARRAYS.


!     AUTHOR: RENATE BROKOPF, MPM HAMBURG, GERMANY. MAY 1994.
!             JEAN BIDLOT ECMWF MAY 1999.
!     

!     EXTERNALS
!     ---------
       
!***********************************************************************

!     MODULES:
!     --------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, ICHNKFROMIJ, IPRMFROMIJ
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSARAS , ONLY : SPECW    ,NSPECW

!-----------------------------------------------------------------------

!     INTERFACE
!     ---------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IJSG, IJLG
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
      INTEGER(KIND=JWIM), INTENT(IN) :: NW1D(IJSG:IJLG)

!     VARIABLES
!     ---------

      INTEGER(KIND=JWIM) :: JP, JFRE, JANG, IP, IK

!-----------------------------------------------------------------------

!      1. UPDATING SPECTRA 
!         ----------------

      IF (NSPECW > 0) THEN

            DO JFRE=1,NFRE
              DO JANG=1,NANG
                DO JP = IJSG, IJLG
                  IF (NW1D(JP) /= 0) THEN
                    IK = ICHNKFROMIJ(JP)
                    IP = IPRMFROMIJ(JP)
                    FL1(IP, JANG, JFRE, IK) = SPECW(NW1D(JP), JANG, JFRE)
                  ENDIF
                ENDDO
              ENDDO
            ENDDO

      ENDIF


!     2. UPDATING WINDS
!        --------------

!     no longer used.

      END SUBROUTINE MERGESARCOR
