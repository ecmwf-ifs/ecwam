      INTEGER FUNCTION wam_getclo(yaoptions, yaargument)

      CHARACTER(LEN=1) :: yolastarg
      CHARACTER(LEN=*) :: yaoptions, yaargument
      CHARACTER(LEN=120) :: arg

      INTEGER :: here, imorearg, ivarg
      DATA here, imorearg, ivarg, arg / 1, 0, 0, "  " /
      DATA yolastarg / " " /
 
      arg=' '
      CALL getarg(here,arg)
 
      iol=len_trim(arg)
      IF (iol .EQ. 2 .AND. arg(1:1) .EQ. '-' .AND. ivarg .EQ. 0 ) THEN
        iol = len_trim(yaoptions)
        DO jl=1,iol
          wam_getclo = 0
          IF ( yaoptions(jl:jl) .EQ. arg(2:2) ) THEN
            wam_getclo = ichar(arg(2:2))
            IF (yaoptions(jl+1:jl+1) .EQ. ':' ) THEN
              yolastarg=yaoptions(jl:jl)
              ivarg=1
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ELSEIF ( ivarg .EQ. 1 ) THEN
         WRITE(*,*) ' option -', yolastarg, ' requires arguments'
         wam_getclo=-1
      ELSEIF (iol .EQ. 0) THEN
        wam_getclo=0
      ELSE
         WRITE(*,*) 'illegal option: ',arg(1:iol)
         wam_getclo=-1
      ENDIF
      here = here + 1

      END FUNCTION wam_getclo
