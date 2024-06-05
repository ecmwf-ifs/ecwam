! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

!     INLINE FUNCTION.
!     ----------------

!     Simple empirical fit to model drag coefficient
      REAL(KIND=JWRB) :: CDM, U

      CDM(U) = MAX(MIN(0.0006_JWRB+0.00008_JWRB*U, 0.001_JWRB+0.0018_JWRB*EXP(-0.05_JWRB*(U-33._JWRB))),0.001_JWRB)