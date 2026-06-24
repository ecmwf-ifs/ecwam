      INTEGER FUNCTION wam_getcla(yaargument)

      CHARACTER(LEN=*) :: yaargument

      CHARACTER(LEN=1) :: yolastarg
      CHARACTER(LEN=120) :: arg

      INTEGER :: here, imorearg, ivarg
      DATA here, imorearg, ivarg, arg / 1, 0, 0, "  " /
      DATA yolastarg / " " /

      wam_getcla = 1
      CALL getarg(here,arg)
      IF ( arg (1:1) .NE. '-' ) THEN
        here = here + 1
        yaargument=arg
      ELSE
        IF (ivarg.EQ.1) THEN
          WRITE(*,*)' refused to take ', arg (1:2) ,' as argument for', &
     &    ' the option -',yolastarg
          wam_getcla = -1
        ELSE
          wam_getcla = 0
        ENDIF
      ENDIF
      ivarg=0

      END FUNCTION wam_getcla
