      LOGICAL FUNCTION SET_WFLAGS(FLAG,N,NTRAIN,IFLAG)

! ----------------------------------------------------------------------

!     J. BIDLOT     ECMWF 

!*    PURPOSE.
!     --------

!     SETS A GLOBAl FLAG FROM AN ARRAY OF LOGICAL FLAGS (FLAG)
!     DEPENDING ON THE LOGIC REQUIRED BY IFLAG


!**   INTERFACE.
!     ----------
!     *SET_WFLAGS(FLAG,N,IFLAG)*

!     *FLAG* LOGICAL  INPUT FLAGS 
!     *N*    INTEGER  SIZE OF FLAG
!     *NTRAIN* INTEGER NUMBER OF WAVE TRAIN 
!     *IFLAG*INTEGER  DETERMINES THE LOGIC TO BE APPLIED TO FLAG 

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: N, IFLAG, IFL, NTRAIN
      LOGICAL :: FLAG(N)

      SET_WFLAGS=.FALSE.

      IF(IFLAG.EQ.20) THEN
        SET_WFLAGS = FLAG(1)
        IFL=2
        DO WHILE (.NOT. SET_WFLAGS .AND. IFL.LE.10)
          SET_WFLAGS = FLAG(IFL)
          IFL=IFL+1
        ENDDO
        IFL=22
        DO WHILE (.NOT. SET_WFLAGS .AND. IFL.LE.27)
          SET_WFLAGS = FLAG(IFL)
          IFL=IFL+1
        ENDDO
        IFL=34
        DO WHILE (.NOT. SET_WFLAGS .AND. IFL.LE.N)
          SET_WFLAGS = FLAG(IFL)
          IFL=IFL+1
        ENDDO

      ELSE IF(IFLAG.EQ.21) THEN
        SET_WFLAGS = FLAG(11)
        IFL=12
        DO WHILE (.NOT. SET_WFLAGS .AND. IFL.LE.16)
          SET_WFLAGS = FLAG(IFL)
          IFL=IFL+1
        ENDDO
        IFL=28
        DO WHILE (.NOT. SET_WFLAGS .AND. IFL.LE.33)
          SET_WFLAGS = FLAG(IFL)
          IFL=IFL+1
        ENDDO
        IFL=47
        DO WHILE (.NOT. SET_WFLAGS .AND. IFL.LE.(46+3*NTRAIN))
          SET_WFLAGS = FLAG(IFL)
          IFL=IFL+1
        ENDDO

      ELSE IF(IFLAG.EQ.44) THEN
        SET_WFLAGS = FLAG(44)
        IFL=45
        DO WHILE (.NOT. SET_WFLAGS .AND. IFL.LE.46)
          SET_WFLAGS = FLAG(IFL)
          IFL=IFL+1
        ENDDO

      ELSE
        WRITE(*,*) 'Nothing was done in SET_WFLAG!'
        WRITE(*,*) 'Is it normal ?'
      ENDIF

      RETURN
      END FUNCTION SET_WFLAGS
