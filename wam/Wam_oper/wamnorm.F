      SUBROUTINE WAMNORM(FIELD,NX,NY,ZMISS,ZMIN,FAVG,FMIN,FMAX,NAVG, 
     &                   LMESSPASS,IRANK,ISEND,IRECV,ITAG)

! ----------------------------------------------------------------------

!**** *WAMNORM* - COMPUTES THE AVERAGE, MINIMUM AND MAXIMUM OF FIELD 

!     J. BIDLOT       ECMWF   MAY 1998 

!*    PURPOSE.
!     --------

!     COMPUTES THE AVERAGE, MINIMUM AND MAXIMUM OF FIELD FOR ALL VALUES
!     DIFFERENT THAN ZMISS 

!**   INTERFACE.
!     ----------
!     *CALL*  WAMNORM(FIELD,NX,NY,ZMISS,FAVG,FMIN,FMAX) 
!          *FIELD*     REAL    INPUT FIELD
!          *NX*        INTEGER FIRST DIMENSION OF FIELD
!          *NY*        INTEGER SECOND DIMENSION OF FIELD
!          *ZMISS*     REAL    VALUE EXCLUDED IN NORM COMPUTATION 
!          *ZMIN*      REAL    ONLY VALUES GREATER THAN ZMIN WILL
!                              CONSIDERED IN THE NORM COMPUTATION 
!          *FAVG*      REAL    AVERAGE VALUE OF FIELD
!          *FMIN*      REAL    MINIMUM VALUE OF FIELD
!          *FMAX*      REAL    MAXIMUM VALUE OF FIELD
!          *NAVG       INTEGER NUMBER OF POINTS USED TO COMPUTE FAVG
!          *LMESSPASS* LOGICAL MESSAGE PASSING SWITCH
!          *IRANK*     INTEGER CURRENT PE 
!          *ISEND*     INTEGER PE WHERE FIELD IS DEFINED
!          *IRECV*     INTEGER PE WHERE NORM SHOULD BE WRITTEN OUT 
!          *ITAG*      INTEGER MESSAGE PASSING TAG

!     METHOD.
!     -------

!     COMPUTES THE AVERAGE, MINIMUM AND MAXIMUM OF FIELD FOR ALL VALUES
!     DIFFERENT THAN ZMISS. FIELD IS SUPPOSED TO BE ON PE ISEND, THE 
!     NORM WILL BE ON PE IRECV. 

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      REAL :: FIELD(NX,NY), ZCOMBUF(4)
      LOGICAL :: LMESSPASS

! ----------------------------------------------------------------------
!*    2. MAKE GRIDDED FIELD
!     ---------------------

      NTOT=0
      FAVG=0.
      FMIN=9999999.
      FMAX=-9999999.
      IF(LMESSPASS) THEN
        IF(IRANK.EQ.ISEND) THEN
          IF(NX.EQ.1.AND.NY.EQ.1) THEN
            FAVG=FIELD(1,1)
            FMIN=FIELD(1,1)
            FMAX=FIELD(1,1)
            NAVG=1
          ELSE
            DO J=1,NY
              DO I=1,NX
                IF(FIELD(I,J).NE.ZMISS) THEN
                  NTOT=NTOT+1
                  FAVG=FAVG+FIELD(I,J)
                  IF(FIELD(I,J).GT.ZMIN) FMIN=MIN(FMIN,FIELD(I,J))
                  FMAX=MAX(FMAX,FIELD(I,J))
                ENDIF
              ENDDO
            ENDDO
            IF(NTOT.GT.0) FAVG=FAVG/NTOT
            NAVG=NTOT
          ENDIF
        ENDIF
        IF(ISEND.NE.IRECV) THEN
          IF(IRANK.EQ.ISEND) THEN
            ZCOMBUF(1)=FAVG
            ZCOMBUF(2)=FMIN
            ZCOMBUF(3)=FMAX
            ZCOMBUF(4)=REAL(NAVG)
            CALL MPE_SEND(ZCOMBUF,4,2,IRECV,ITAG,0,0,0,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     &                    ('MPE_SEND ERROR in WAMNORM' )
          ELSEIF(IRANK.EQ.IRECV) THEN
            CALL MPE_RECV(ZCOMBUF,4,2,ISEND,ITAG,0,0,0,
     &                    KRCOUNT,KRFROM,KRTAG,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     &                       ('MPE_RECV ERROR in WAMNORM' )
            IF(KRCOUNT.NE.4) CALL MPEI_ABORT
     &      ('MPE_RECV ERROR in WAMNORM:MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     &      ('MPE_RECV ERROR in WAMNORM MISMATCHED TAGS' )

            FAVG=ZCOMBUF(1)
            FMIN=ZCOMBUF(2)
            FMAX=ZCOMBUF(3)
            NAVG=INT(ZCOMBUF(4))
          ENDIF
        ENDIF

      ELSE
        IF(NX.EQ.1.AND.NY.EQ.1) THEN
          FAVG=FIELD(1,1)
          FMIN=FIELD(1,1)
          FMAX=FIELD(1,1)
          NAVG=1
        ELSE
          DO J=1,NY
            DO I=1,NX
              IF(FIELD(I,J).NE.ZMISS) THEN
                NTOT=NTOT+1
                FAVG=FAVG+FIELD(I,J)
                IF(FIELD(I,J).GT.ZMIN) FMIN=MIN(FMIN,FIELD(I,J))
                FMAX=MAX(FMAX,FIELD(I,J))
              ENDIF
            ENDDO
          ENDDO
          IF(NTOT.GT.0) FAVG=FAVG/NTOT
          NAVG=NTOT
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE WAMNORM
