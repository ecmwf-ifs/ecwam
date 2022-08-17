module yowabort
implicit none
private

public :: wam_abort

interface wam_abort
  module procedure wam_abort_noargs
  module procedure wam_abort_msg
  module procedure wam_abort_msg_loc
  module procedure wam_abort_loc
end interface

contains

subroutine flush_open_units()
use parkind_wave, only : jwim
implicit none
integer(jwim) :: iu
logical :: lopened
do iu=1,99
  inquire(unit=iu,opened=lopened)
  if(lopened) then
    call flush(iu)
  end if
end do
end subroutine

subroutine wam_abort_noargs()
implicit none
#include "abor1.intfb.h"
call flush_open_units()
call abor1("!!!! *** WAVE MODEL HAS ABORTED *** !!!!")
end subroutine

subroutine wam_abort_msg(message)
implicit none
#include "abor1.intfb.h"
character(len=*), intent(in) :: message
call flush_open_units()
call abor1("!!!! *** WAVE MODEL HAS ABORTED *** !!!!"//new_line('A')//message)
end subroutine

subroutine wam_abort_msg_loc(message,file,line)
use parkind_wave, only : jwim
implicit none
#include "abor1.intfb.h"
character(len=*), intent(in) :: message
character(len=*), intent(in) :: file
integer(jwim),    intent(in) :: line
call flush_open_units()
call abor1fl(file,line,"!!!! *** WAVE MODEL HAS ABORTED *** !!!!"//new_line('A')//message)
end subroutine

subroutine wam_abort_loc(file,line)
use parkind_wave, only : jwim
implicit none
#include "abor1.intfb.h"
character(len=*), intent(in) :: file
integer(jwim),    intent(in) :: line
call flush_open_units()
call abor1fl(file,line,"!!!! *** WAVE MODEL HAS ABORTED *** !!!!")
end subroutine

end module
