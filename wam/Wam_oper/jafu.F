      INTEGER FUNCTION JAFU (CL, J, IAN)

C ----------------------------------------------------------------------
C
C**** *JAFU* - FUNCTION TO COMPUTE THE INDEX ARRAY FOR THE
C              ANGLES OF THE INTERACTING WAVENUMBERS.
C
C     S. HASSELMANN        MPIFM        01/12/1985.
C
C*    PURPOSE.
C     --------
C
C       INDICES DEFINING BINS IN FREQUENCY AND DIRECTION PLANE INTO
C       WHICH NONLINEAR ENERGY TRANSFER INCREMENTS ARE STORED. NEEDED
C       FOR COMPUTATION OF THE NONLINEAR ENERGY TRANSFER.
C
C**   INTERFACE.
C     ----------
C
C       *FUNCTION* *JAFU (CL, J, IAN)*
C          *CL*  - WEIGHTS.
C          *J*   - INDEX IN ANGULAR ARRAY.
C          *IAN* - NUMBER OF ANGLES IN ARRAY.
C
C     METHOD.
C     -------
C
C       SEE REFERENCE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C        S. HASSELMANN AND K. HASSELMANN,JPO, 1985 B.
C
C ----------------------------------------------------------------------
C
      IDPH = CL
      JA = J+IDPH
      IF (JA.LE.0)   JA = IAN+JA-1
      IF (JA.GE.IAN) JA = JA-IAN+1
      JAFU = JA

      RETURN
      END
