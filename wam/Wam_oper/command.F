C ======================================================================
C
      SUBROUTINE command (kuso, cddir, cdstep, cdroot, cdlevel, cdtype,
     .                    kerr)
C
C ======================================================================
C
C**** command
C****           crack command line options for cpfdb
C
C        B. HANSEN   *ECMWF*   13/03/98
C
C     Interface.
C     ----------
C         INPUT.
C         ------
C           *command line*   -d ccssssdxxxxyyyymmdd -t tthh  -s step
C                            -D ccssssdxxxxyyyymmdd -T tthh  -S step
C
C                            -p parameter -l level
C
C            where ccssssdxxxxyyyymmdd may be rdwavesat7k19980224
C            and   tthh                may be an21
C                  step                may be 000001
C                       ie. cc = class
C                       ie. ssss = stream
C                       ie. d = domain
C                       ie. xxxx = experiment id
C                       ie. yyyy = year
C                       ie. mm = month
C                       ie. dd = day
C                       ie. tt = type
C                       ie. hh = time
C             -d from_directory  -t from_type [ -r from_root ]
C             -D to_directory    -t to_type   [ -R to_root   ]
C             [ -x decode ] decode = decode
C
C         OUTPUT.
C         -------
C           *iu_grs*    output unit opened for grs file.
C
C     Method.
C     -------
C       NONE.
C
C     Modification.
C     -------------
C       NONE.
C
C
C ----------------------------------------------------------------------
C
      IMPLICIT NONE
C.....!!!!!!!!!!!!!
C
C ----------------------------------------------------------------------
C
      INTEGER GETCLO, GETCLA, OPTVAL, kuso, kerr, MORARG
C
      CHARACTER (LEN=  *)               :: CDLEVEL
      CHARACTER (LEN=  *), DIMENSION(2) :: CDTYPE
      CHARACTER (LEN=  *), DIMENSION(2) :: CDSTEP
      CHARACTER (LEN=  *), DIMENSION(2) :: CDDIR, CDROOT
C
      CHARACTER (LEN= 20), PARAMETER :: clopts = 'd;t;D;T;s;S;p;l;r;R;'
      CHARACTER (LEN=128)            :: clarg
      CHARACTER                      :: cloptlet
C
      CDLEVEL = ' '
      CDTYPE = ' '
      CDSTEP = ' '
      CDDIR = ' '
      CDROOT = ' '
C
C*    WRITE(kuso,'(" command Mon Mar 16 19:46:13 GMT 1998 ")')
C
      COMLINE : DO
        OPTVAL = GETCLO (clopts,clarg)
        cloptlet=CHAR(OPTVAL)
C       WRITE(kuso,'(" cloptlet=", a2, $ )') cloptlet
        IF ( OPTVAL .GT. 0 ) THEN      
          MORARG = GETCLA (clarg)
          CALL U2l1CR ( clarg )
          IF(MORARG.NE.0) THEN
C           WRITE(kuso,'(" clarg ", a30)') clarg(1:30)
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
     .     CDTYPE(1) .EQ. ' ' .OR. CDTYPE(2) .EQ. ' ' .OR.
     .     CDSTEP(1) .EQ. ' ' .OR. CDSTEP(2) .EQ. ' ' .OR.
     .     CDLEVEL .EQ. ' ' ) THEN
        WRITE (kuso, '("  arguments given not sufficient ")')
        KERR = -1
      ENDIF
C
      RETURN
      END
