      SUBROUTINE WAMNORM(FIELD,NX,NY,ZMISS,ZMIN,FAVG,FMIN,FMAX,NAVG, 
     &                   LMESSPASS,IRANK,ISEND,IRECV,ITAG)
C
C ----------------------------------------------------------------------
C
C**** *WAMNORM* - COMPUTES THE AVERAGE, MINIMUM AND MAXIMUM OF FIELD 
C
C     J. BIDLOT       ECMWF   MAY 1998 
C
C*    PURPOSE.
C     --------
C
C     COMPUTES THE AVERAGE, MINIMUM AND MAXIMUM OF FIELD FOR ALL VALUES
C     DIFFERENT THAN ZMISS 
C
C**   INTERFACE.
C     ----------
C     *CALL*  WAMNORM(FIELD,NX,NY,ZMISS,FAVG,FMIN,FMAX) 
C          *FIELD*     REAL    INPUT FIELD
C          *NX*        INTEGER FIRST DIMENSION OF FIELD
C          *NY*        INTEGER SECOND DIMENSION OF FIELD
C          *ZMISS*     REAL    VALUE EXCLUDED IN NORM COMPUTATION 
C          *ZMIN*      REAL    ONLY VALUES GREATER THAN ZMIN WILL
C                              CONSIDERED IN THE NORM COMPUTATION 
C          *FAVG*      REAL    AVERAGE VALUE OF FIELD
C          *FMIN*      REAL    MINIMUM VALUE OF FIELD
C          *FMAX*      REAL    MAXIMUM VALUE OF FIELD
C          *NAVG       INTEGER NUMBER OF POINTS USED TO COMPUTE FAVG
C          *LMESSPASS* LOGICAL MESSAGE PASSING SWITCH
C          *IRANK*     INTEGER CURRENT PE 
C          *ISEND*     INTEGER PE WHERE FIELD IS DEFINED
C          *IRECV*     INTEGER PE WHERE NORM SHOULD BE WRITTEN OUT 
C          *ITAG*      INTEGER MESSAGE PASSING TAG
C
C     METHOD.
C     -------
C
C     COMPUTES THE AVERAGE, MINIMUM AND MAXIMUM OF FIELD FOR ALL VALUES
C     DIFFERENT THAN ZMISS. FIELD IS SUPPOSED TO BE ON PE ISEND, THE NORM
C     WILL BE PE IRECV. 
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      REAL :: FIELD(NX,NY), ZCOMBUF(4)
      LOGICAL :: LMESSPASS
C
C ----------------------------------------------------------------------
C*    2. MAKE GRIDDED FIELD
C     ---------------------
C
      NTOT=0
      FAVG=0.
      FMIN=9999999.
      FMAX=-9999999.
      IF(LMESSPASS) THEN
        IF(IRANK.EQ.ISEND) THEN
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
        IF(ISEND.NE.IRECV) THEN
          IF(IRANK.EQ.ISEND) THEN
            ZCOMBUF(1)=FAVG
            ZCOMBUF(2)=FMIN
            ZCOMBUF(3)=FMAX
            ZCOMBUF(4)=FLOAT(NAVG)
            CALL MPE_SEND(ZCOMBUF,4,2,IRECV,ITAG,0,0,0,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     &                    ('MPE_SEND ERROR in WAMNORM' )
          ELSE IF(IRANK.EQ.IRECV) THEN
            CALL MPE_RECV(ZCOMBUF,4,2,ISEND,ITAG,0,0,0,
     &                    KRCOUNT,KRFROM,KRTAG,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     &                       ('MPE_RECV ERROR in WAMNORM' )
            IF(KRCOUNT.NE.4) CALL MPEI_ABORT
     &      ('MPE_RECV ERROR in WAMNORM:MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     &      ('MPE_RECV ERROR in WAMNORM MISMATCHED TAGS' )
C
            FAVG=ZCOMBUF(1)
            FMIN=ZCOMBUF(2)
            FMAX=ZCOMBUF(3)
            NAVG=INT(ZCOMBUF(4))
          ENDIF
        ENDIF
C
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
C
      RETURN
      END
