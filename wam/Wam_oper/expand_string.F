      SUBROUTINE EXPAND_STRING(
     $     myproc,              ! %p
     $     nproc,               ! %n
     $     timestep,            ! %t
     $     max_timestep,
     $     s,                   ! %s
     $     n)
C
C    S. SAARINEM   ECMWF   MAY 1996   MESSAGE PASSING
C
C*    PURPOSE.
C     --------
C
C     	EXPAND STRING BY ADDING AN EXTENSION WHICH IS FUNCTION OF
C       THE PROCESS NUMBER, AND/OR THE TOTAL NUMBER OF PROCESSES,
C       AND/OR THE TIME STEP NUMBER, AND/OR THE INDEX OF THE STRING
C
C**   INTERFACE.
C     ----------
C     CALL EXPAND_STRING(MYPROC,NPROC,TIMESTEP,MAX_TIMESTEP,S,N)
C
C     *MYPROC*       INTEGER  PROCESS NUMBER (PE)
C     *NPROC*        INTEGER  TOTAL NUMBER OF PE'S
C     *TIMESTEP*     INTEGER  TIME STEP NUMBER (OR ANY OTHER INDEX)
C     *MAX_TIMESTEP* INTEGER  MAXIMUM TIME STEP NUMBER (OR ANY OTHER INDEX)
C     *S*            CHARACTER OR CHARACTER ARRAY THAT WILL BE EXPANDED
C     *N*            INTEGER  DIMENSION OF S 
C
C
C     METHOD.
C     -------
C     THE STRING S MUST CONTAIN THE FOLLOWING SUB-STRING THAT WILL THAN BE
C     SUBSTITUTED BY THE CORRESPONDING EXPANDED SUB STRING, THAT STARTS WITH
C     % AND p, n, t, or s
C
C     %p   -->  myproc in character mode, left justified
C     %n   -->  nproc  in character mode, left justified
C     %t   -->  timestep in character mode, left justified
C     %s   -->  the array index of S, left justified
C
C     e.g.   s='name_%n.%p'
C            call(myproc,nproc,0,0,s,1)
C
C     will produce, if nproc=11
C
C     on PE1,  s = name_11.1 
C     on PE2,  s = name_11.2 
C     on PE3,  s = name_11.3 
C
C     ...
C
C     on PE10,  s = name_11.10 
C     on PE11,  s = name_11.11
C
C
C     EXTERNALS.
C     ---------
C

C     REFERENCE.
C     ----------
C
C       NONE

C-----------------------------------------------------------------------------
      implicit none
      integer, intent(in)          :: myproc, nproc
      integer, intent(in)          :: timestep, max_timestep
      integer, intent(in)          :: n
      character*(*), intent(inout) :: s(n)
! === END OF INTERFACE BLOCK ===
      character(len=len(s))   t
      character(len=2*len(s)) tt
      integer i, j, jj, loc_p, len_t
      integer ndigs(4), num(4)
      character*6 fmt(4)

      if (n < 1) return

c*    Setup output formats
      num(1) = myproc
      num(2) = max(nproc,myproc)
      num(3) = n
      num(4) = max(max_timestep,timestep)

c*    Count number of digits in each integer
      do j=1,4
         ndigs(j) = 1
         if (num(j) /= 0) then
            ndigs(j) = 1 + log10(dble(abs(num(j))))
            if (num(j) < 0) ndigs(j) = ndigs(j) + 1 ! Room for minus sign
         endif
         ndigs(j) = min(ndigs(j),9)   ! Max 9 digits supported; i.e. '999 999 999'
         write(fmt(j),'("(i",i1,")")') ndigs(j)
      enddo

c
c*    Expand fields '%s', '%p', '%n' and '%t' with their values
c

c*    A special treatment with the sequence numbering
      if (n>1) then
         loc_p = index(s(1),'%s')
         if (loc_p > 0) then
            s(2:) = s(1)
         endif
      endif

      do i=1,n
         t = adjustl(s(i))//' '
         loc_p = index(t,'%')

         if (loc_p > 0) then
            len_t = len_trim(t)
            j = loc_p
            tt(:j-1) = t(:j-1)
            tt(j:) = ' '
            jj = j-1

            do while (j <= len_t)
               if (t(j:j) == '%') then
                  j = j + 1
                  if (j <= len_t) then
                     select case ( t(j:j) ) 
                       case ( 'p' )   ! myproc
                          write(tt(jj+1:jj+ndigs(1)),fmt(1)) myproc
                          jj = jj + ndigs(1)
                       case ( 'n' )   ! nproc
                          write(tt(jj+1:jj+ndigs(2)),fmt(2)) nproc
                          jj = jj + ndigs(2)
                       case ( 's' )   ! sequence number i=[1..n]
                          write(tt(jj+1:jj+ndigs(3)),fmt(3)) i
                          jj = jj + ndigs(3)
                       case ( 't' )   ! timestep
                          write(tt(jj+1:jj+ndigs(4)),fmt(4)) timestep
                          jj = jj + ndigs(4)
                       case default
                          tt(jj+1:jj+2) = '%'//t(j:j)
                          jj = jj + 2
                     end select
                  else
                     tt(jj+1:jj+1) = '%'
                     jj = jj + 1
                  endif
               else
                  tt(jj+1:jj+1) = t(j:j)
                  jj = jj + 1
               endif
               j = j + 1
            enddo

            t = adjustl(tt)

c*   Get also rid of any blanks in the middle of the string

            len_t = len_trim(t)
            j = 1
            do while (j < len_t)
               if (t(j:j) == ' ') then
                  t(j:) = t(j+1:)
                  len_t = len_trim(t)
               else
                  j = j + 1
               endif
            enddo

         endif

         s(i) = t
      enddo

      END SUBROUTINE expand_string
