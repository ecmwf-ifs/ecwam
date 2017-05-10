!> \file yowpd.F90
!> \brief This is an interface module. Simply call use yowpd in your fortran code to use pdlib
!> \author Thomas Huxhorn
!> \date 2013

module yowpd
  use yowError
  use yowDatapool
  use yowpdlibMain
  use yowNodepool
  use yowElementpool
  use yowSidepool
  use yowExchangeModule
  use yowRankModule
!  comment out because it cause cyclic dependency
!   use yowMpiModule
  use yowChecksModule
  implicit none
  public
end module yowpd