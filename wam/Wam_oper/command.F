! ======================================================================

      SUBROUTINE command (kuso, cddir, cdstep, cdroot, cdlevel, cdtype,
     &                    kerr)

! ======================================================================

!**** command
!****           crack command line options for cpfdb

!        B. HANSEN   *ECMWF*   13/03/98

!     Interface.
!     ----------
!         INPUT.
!         ------
!           *command line*   -d ccssssdxxxxyyyymmdd -t tthh  -s step
!                            -D ccssssdxxxxyyyymmdd -T tthh  -S step

!                            -p parameter -l level

!            where ccssssdxxxxyyyymmdd may be rdwavesat7k19980224
!            and   tthh                may be an21
!                  step                may be 000001
!                       ie. cc = class
!                       ie. ssss = stream
!                       ie. d = domain
!                       ie. xxxx = experiment id
!                       ie. yyyy = year
!                       ie. mm = month
!                       ie. dd = day
!                       ie. tt = type
!                       ie. hh = time
!             -d from_directory  -t from_type [ -r from_root ]
!             -D to_directory    -t to_type   [ -R to_root   ]
!             [ -x decode ] decode = decode

!         OUTPUT.
!         -------
!           *iu_grs*    output unit opened for grs file.

!     Method.
!     -------
!       NONE.

!     Modification.
!     -------------
!       NONE.


! ----------------------------------------------------------------------

      IMPLICIT NONE
!.....!!!!!!!!!!!!!

! ----------------------------------------------------------------------

      INTEGER GETCLO, GETCLA, OPTVAL, kuso, kerr, MORARG

      CHARACTER (LEN=  *)               :: CDLEVEL
      CHARACTER (LEN=  *), DIMENSION(2) :: CDTYPE
      CHARACTER (LEN=  *), DIMENSION(2) :: CDSTEP
      CHARACTER (LEN=  *), DIMENSION(2) :: CDDIR, CDROOT

      CHARACTER (LEN= 20), PARAMETER :: clopts = 'd;t;D;T;s;S;p;l;r;R;'
      CHARACTER (LEN=128)            :: clarg
      CHARACTER                      :: cloptlet

      CDLEVEL = ' '
      CDTYPE = ' '
      CDSTEP = ' '
      CDDIR = ' '
      CDROOT = ' '

!*    WRITE(kuso,'(" command Mon Mar 16 19:46:13 GMT 1998 ")')

      COMLINE : DO
        OPTVAL = GETCLO (clopts,clarg)
        cloptlet=CHAR(OPTVAL)
!       WRITE(kuso,'(" cloptlet=", a2, $ )') cloptlet
        IF ( OPTVAL .GT. 0 ) THEN      
          MORARG = GETCLA (clarg)
          CALL U2l1CR ( clarg )
          IF(MORARG.NE.0) THEN
!           WRITE(kuso,'(" clarg ", a30)') clarg(1:30)
            IF     ( cloptlet .EQ. 'd' ) THEN
              CDDIR(1) = clarg
            ELSEIF ( cloptlet .EQ. 'D' ) THEN
              CDDIR(2) = clarg
            ELSEIF ( cloptlet .EQ. 't' ) THEN
              CDTYPE(1) = clarg
            ELSEIF ( cloptlet .EQ. 'T' ) THEN
              CDTYPE(2) = clarg
            ELSEIF ( cloptlet .EQ. 's' ) THEN
              CDSTEP(1) = clarg
            ELSEIF ( cloptlet .EQ. 'S' ) THEN
              CDSTEP(2) = clarg
            ELSEIF ( cloptlet .EQ. 'r' ) THEN
              CDROOT(1) = clarg
            ELSEIF ( cloptlet .EQ. 'R' ) THEN
              CDROOT(2) = clarg
            ELSEIF ( cloptlet .EQ. 'l' ) THEN
              CDLEVEL = clarg
            ELSE
              WRITE(kuso, '("  option not supported ")')
              KERR = -1
            ENDIF 
          ENDIF
        ELSE
          EXIT COMLINE
        ENDIF
      ENDDO COMLINE
      IF ( CDDIR(1) .EQ. ' ' .OR. CDDIR(2) .EQ. ' ' .OR.
     &     CDTYPE(1) .EQ. ' ' .OR. CDTYPE(2) .EQ. ' ' .OR.
     &     CDSTEP(1) .EQ. ' ' .OR. CDSTEP(2) .EQ. ' ' .OR.
     &     CDLEVEL .EQ. ' ' ) THEN
        WRITE (kuso, '("  arguments given not sufficient ")')
        KERR = -1
      ENDIF

      RETURN
      END SUBROUTINE command


