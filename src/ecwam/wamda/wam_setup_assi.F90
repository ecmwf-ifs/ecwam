SUBROUTINE WAM_SETUP_ASSI
USE YOWASSI, ONLY : WAMASSI_HANDLER, GETODBRALT_HANDLER, RFL2ODB_HANDLER, WAM2ODB_HANDLER, WAM_ODB_CLOSE_HANDLER, WAM_ODB_OPEN_HANDLER
! ----------------------------------------------------
#include "wamassi.intfb.h"

#ifdef WITH_ODB
#include "getodnralt.intfb.h"
#include "rfl2odb.intfb.h"
#include "wam2odb.intfb.h"
#include "wam_odb_open.intfb.h"
#include "wam_odb_close.intfb.h"
#endif
! ----------------------------------------------------
WAMASSI_HANDLER       => WAMASSI

#ifdef WITH_ODB
GETODBRALT_HANDLER    => GETODBRALT
RFL2ODB_HANDLER       => RFL2ODB
WAM2ODB_HANDLER       => WAM2ODB
WAM_ODB_OPEN_HANDLER  => WAM_ODB_OPEN
WAM_ODB_CLOSE_HANDLER => WAM_ODB_CLOSE
#endif

END SUBROUTINE
