      SUBROUTINE CHECKCFL (MIJS, MIJL, DTC,
     &                     CFLEA,CFLWE,CFLNO,CFLSO,CFLNO2,CFLSO2, 
     &                     CFLTP,CFLTM,CFLOP,CFLOM)

! ----------------------------------------------------------------------

!**** *CHECKCFL* - 

!     J. BIDLOT    ECMWF   2005

!*    PURPOSE.
!     --------

!     CHECK THAT THE CFL CRITERIA ARE NOT VIOLATED ANYWHERE.

!**   INTERFACE.
!     ----------

!       *CALL* *CHECKCFL(MIJS, MIJL, DTC
!    &                   CFLEA,CFLWE,CFLNO,CFLSO,CFLNO2,CFLSO2, 
!    &                   CFLTP,CFLTM,CFLOP,CFLOM)*
!          *MIJS* - INDEX OF FIRST POINT
!          *MIJL* - INDEX OF LAST POINT
!          AND ALL COEFFICIENTS OF THE UPWIND SCHEME.

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!       NONE.

!     REFERENCE.
!     ----------
!       NONE.

! ----------------------------------------------------------------------

      USE YOWMPP   , ONLY : NINF
      USE YOWSHAL  , ONLY : INDEP    ,DEPTH
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON     ,KLAT2

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: MIJS, MIJL, IJ, IC, ICP
      REAL,DIMENSION(MIJS:MIJL) :: CFLEA,CFLWE,CFLNO,CFLSO,CFLNO2 
      REAL,DIMENSION(MIJS:MIJL) :: CFLSO2,CFLTP,CFLTM,CFLOP,CFLOM
      REAL,DIMENSION(MIJS:MIJL) :: DTC
      LOGICAL :: LLABORTCFL

! ----------------------------------------------------------------------

      LLABORTCFL=.FALSE.

      DO IJ = MIJS,MIJL

        IF(CFLEA(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE EAST  *'
          WRITE(IU06,*) ' * CFLEA = ',CFLEA(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(CFLWE(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *          *                   *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE WEST  *'
          WRITE(IU06,*) ' * CFLWE = ',CFLWE(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(CFLNO(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE NORTH *'
          WRITE(IU06,*) ' * CFLNO = ',CFLNO(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(CFLSO(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE SOUTH *'
          WRITE(IU06,*) ' * CFLSO = ',CFLSO(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(CFLNO2(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE NORTH *'
          WRITE(IU06,*) ' * CFLNO2 = ',CFLNO2(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(CFLSO2(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR PROPAGATION TO THE SOUTH *'
          WRITE(IU06,*) ' * CFLSO2 = ',CFLSO2(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(CFLTP(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR POSITIVE DIRECTIONS      *'
          WRITE(IU06,*) ' * CFLTP = ',CFLTP(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' * IJ = ',IJ, INDEP(IJ),DEPTH(IJ,1)
          DO IC=1,2
            ICP=KLAT(IJ,IC)
            IF(ICP.NE.NINF-1)
     &      WRITE(IU06,*)' * KLAT_',IC, INDEP(ICP),DEPTH(ICP,1)
          ENDDO
          DO IC=1,2
            ICP=KLAT2(IJ,IC)
            IF(ICP.NE.NINF-1)
     &      WRITE(IU06,*)' * KLAT2',IC, INDEP(ICP),DEPTH(ICP,1)
          ENDDO
          DO IC=1,2
            ICP=KLON(IJ,IC)
            IF(ICP.NE.NINF-1)
     &      WRITE(IU06,*)' * KLON_',IC, INDEP(ICP),DEPTH(ICP,1)
          ENDDO
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(CFLTM(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR NEGATIVE DIRECTIONS      *'
          WRITE(IU06,*) ' * CFLTM = ',CFLTM(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' * IJ = ',IJ, INDEP(IJ),DEPTH(IJ,1)
          DO IC=1,2
            ICP=KLAT(IJ,IC)
            IF(ICP.NE.NINF-1)
     &      WRITE(IU06,*)' * KLAT_',IC, INDEP(ICP),DEPTH(ICP,1)
          ENDDO
          DO IC=1,2
            ICP=KLAT2(IJ,IC)
            IF(ICP.NE.NINF-1)
     &      WRITE(IU06,*)' * KLAT2',IC, INDEP(ICP),DEPTH(ICP,1)
          ENDDO
          DO IC=1,2
            ICP=KLON(IJ,IC)
            IF(ICP.NE.NINF-1)
     &      WRITE(IU06,*)' * KLON_',IC, INDEP(ICP),DEPTH(ICP,1)
          ENDDO
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(CFLOP(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' * ============================ *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR FREQUENCY TO THE RIGHT   *'
          WRITE(IU06,*) ' * CFLOP = ',CFLOP(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(CFLOM(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ********************************'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL *'
          WRITE(IU06,*) ' *  ==========================  *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED  *'
          WRITE(IU06,*) ' * FOR FREQUENCY TO THE LEFT    *'
          WRITE(IU06,*) ' * CFLOM = ',CFLOM(IJ)
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY   *'
          WRITE(IU06,*) ' *                              *'
          WRITE(IU06,*) ' ********************************'
          LLABORTCFL=.TRUE.
        ENDIF

        IF(DTC(IJ).GT.1.) THEN
          WRITE(IU06,*) ' ************************************'
          WRITE(IU06,*) ' *                                  *'
          WRITE(IU06,*) ' * FATAL ERROR IN SUB. CHECKCFL     *'
          WRITE(IU06,*) ' * ============================     *'
          WRITE(IU06,*) ' * STABILITY CONDITION NOT FULFILLED*'
          WRITE(IU06,*) ' * DTC = ',DTC(IJ)
          WRITE(IU06,*) ' *                                  *'
          WRITE(IU06,*) ' * REDUCE IDELPRO ACCORDINGLY       *'
          WRITE(IU06,*) ' ************************************'
          LLABORTCFL=.TRUE.
        ENDIF
        IF(LLABORTCFL) THEN
          WRITE(IU06,*) ' ************************************'
          WRITE(IU06,*) ' * ABORTING IN WAVE MODEL CHECKFL !!*'
          WRITE(IU06,*) ' ************************************'
          CALL FLUSH(IU06)
          CALL ABORT1
        ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CHECKCFL
