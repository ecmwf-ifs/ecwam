      FUNCTION IECF_len(y_char)
c
c  Calculate the 'length' of a character string.
c
c  The 'length' is the character-position of the last character in the
c  string which is neither 'BLANK', nor 'NULL'.
c
      CHARACTER*(*) y_char
      CHARACTER*1   null
c
      null = char(0)
      lgth = len(y_char)
c
      DO 10 i = lgth , 1 , -1
      IF ((y_char(i:i) .NE. ' ') .AND. (y_char(i:i) .NE. null)) goto 20
   10 CONTINUE
      i = 0
c
   20 CONTINUE
      IECF_len = i
c
      RETURN
      END
