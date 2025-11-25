! (C) Copyright 2013 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

#if WAM_HAVE_ATLAS
#define WAM_HAVE_ATLAS_TRACE 1
        ! Set to 0 to disable ATLAS_TRACE
#else
#define WAM_HAVE_ATLAS_TRACE 0
#endif

module wam_trace_mod

#if WAM_HAVE_ATLAS_TRACE
use atlas_Trace_module, only : atlas_Trace
#endif

implicit none

private
!-----------------------------
! atlas_Trace                !
!-----------------------------

type, public :: WAM_TRACE
#if WAM_HAVE_ATLAS_TRACE
  TYPE(atlas_Trace), POINTER :: trace => null()
#endif
contains
! Public methods

  procedure, private :: init_loc
  procedure, private :: init_labels_1

  generic, public :: init =>      &
    & init_loc,      &
    & init_labels_1

  procedure, public :: running
  procedure, public :: start
  procedure, public :: stop
  procedure, public :: pause
  procedure, public :: resume
  procedure, public :: elapsed

  procedure, public :: final
end type

!========================================================
contains
!========================================================

subroutine init_loc(this,file,line,title)
  class(WAM_TRACE) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
#if WAM_HAVE_ATLAS_TRACE
  if (associated(this%trace)) then
    call this%trace%final()
  else
    allocate(this%trace)
  endif
  call this%trace%init(file,line,title)
#endif
end subroutine

subroutine init_labels_1(this,file,line,title,label)
  class(WAM_TRACE) :: this
  character(len=*) , intent(in) :: file
  integer          , intent(in) :: line
  character(len=*) , intent(in) :: title
  character(len=*) , intent(in) :: label
#if WAM_HAVE_ATLAS_TRACE
  if (associated(this%trace)) then
    call this%trace%final()
  else
    allocate(this%trace)
  endif
  call this%trace%init(file,line,title,label)
#endif
end subroutine

!-------------------------------------------------------------------------------

function running( this )
  logical :: running
  class(WAM_TRACE) :: this
  running = .False.
#if WAM_HAVE_ATLAS_TRACE
  if (associated(this%trace)) then
    running = this%trace%running()
  endif
#endif
end function

!-------------------------------------------------------------------------------

function elapsed( this )
  use, intrinsic :: iso_c_binding, only : c_double
  real(c_double) :: elapsed
  class(WAM_TRACE) :: this
  elapsed = 0._c_double
#if WAM_HAVE_ATLAS_TRACE
  if (associated(this%trace)) then
    elapsed = this%trace%elapsed()
  endif
#endif
end function

!-------------------------------------------------------------------------------

subroutine start( this )
  class(WAM_TRACE) :: this
#if WAM_HAVE_ATLAS_TRACE
  if (associated(this%trace)) then
    call this%trace%start()
  endif
#endif
end subroutine

!-------------------------------------------------------------------------------

subroutine stop( this )
  class(WAM_TRACE) :: this
#if WAM_HAVE_ATLAS_TRACE
  if (associated(this%trace)) then
    call this%trace%stop()
  endif
#endif
end subroutine

!-------------------------------------------------------------------------------

subroutine pause( this )
  class(WAM_TRACE) :: this
#if WAM_HAVE_ATLAS_TRACE
  if (associated(this%trace)) then
    call this%trace%pause()
  endif
#endif
end subroutine

!-------------------------------------------------------------------------------

subroutine resume( this )
  class(WAM_TRACE) :: this
#if WAM_HAVE_ATLAS_TRACE
  if (associated(this%trace)) then
    call this%trace%resume()
  endif
#endif
end subroutine

!-------------------------------------------------------------------------------

subroutine final(this)
  class(WAM_TRACE), intent(inout) :: this
#if WAM_HAVE_ATLAS_TRACE
  if (associated(this%trace)) then
    call this%trace%final()
    deallocate(this%trace)
  endif
  nullify(this%trace)
#endif
end subroutine

!-------------------------------------------------------------------------------

end module

