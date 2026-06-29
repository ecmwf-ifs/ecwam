module yowcmdlo

  implicit none

  private
  public :: wam_getclo, wam_getcla, reset_wam_getclo

  integer, save :: here  = 1
  integer, save :: ivarg = 0

  character(len=1), save :: yolastarg = ' '

contains

  integer function wam_getclo(yaoptions, yaargument) result(iret)

    implicit none

    character(len=*), intent(in) :: yaoptions

    ! Retained only for interface compatibility with the old routine.
    ! It is not used directly by wam_getclo.
    character(len=*), intent(inout), optional :: yaargument

    character(len=120) :: arg

    integer :: iol
    integer :: nopt
    integer :: jl
    logical :: found

    iret = 0
    arg  = ' '

    if (here > command_argument_count()) then
      iret = 0
      return
    end if

    call get_command_argument(here, arg)

    iol = len_trim(arg)

    if (iol == 2 .and. arg(1:1) == '-' .and. ivarg == 0) then

      nopt  = len_trim(yaoptions)
      found = .false.

      do jl = 1, nopt

        if (yaoptions(jl:jl) == arg(2:2)) then

          found = .true.
          iret = ichar(arg(2:2))

          if (jl < nopt) then
            if (yaoptions(jl+1:jl+1) == ':') then
              yolastarg = yaoptions(jl:jl)
              ivarg = 1
            end if
          end if

          exit

        end if

      end do

      if (.not. found) then
        write(*,*) 'illegal option: ', trim(arg)
        iret = -1
      end if

    else if (ivarg == 1) then

      write(*,*) ' option -', yolastarg, ' requires arguments'
      iret = -1

    else if (iol == 0) then

      iret = 0

    else

      write(*,*) 'illegal option: ', arg(1:iol)
      iret = -1

    end if

    here = here + 1

  end function wam_getclo


  integer function wam_getcla(yaargument) result(iret)

    implicit none

    character(len=*), intent(out) :: yaargument

    character(len=120) :: arg

    integer :: iol

    iret = 1
    arg  = ' '
    yaargument = ' '

    if (here > command_argument_count()) then

      if (ivarg == 1) then
        write(*,*) ' option -', yolastarg, ' requires arguments'
        iret = -1
      else
        iret = 0
      end if

      ivarg = 0
      return

    end if

    call get_command_argument(here, arg)

    iol = len_trim(arg)

    if (iol == 0) then

      if (ivarg == 1) then
        write(*,*) ' option -', yolastarg, ' requires arguments'
        iret = -1
      else
        iret = 0
      end if

    else if (arg(1:1) /= '-') then

      here = here + 1
      yaargument = arg

    else

      if (ivarg == 1) then

        write(*,*) ' refused to take ', arg(1:min(2,iol)), &
                   ' as argument for the option -', yolastarg
        iret = -1

      else

        iret = 0

      end if

    end if

    ivarg = 0

  end function wam_getcla


  subroutine reset_wam_getclo()

    implicit none

    here      = 1
    ivarg     = 0
    yolastarg = ' '

  end subroutine reset_wam_getclo

end module yowcmdlo
