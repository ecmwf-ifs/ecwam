      subroutine MYINDEX(ndim, input, index )
C ----------------------------------------------------------------------
C
C**** *MYINDEX* -
C
C     B. Hansen  ECMWF March 1996
C
C*    PURPOSE.
C     --------
C
C     order the element of integer array input and output solution
C     in array index 
C
C     NOTE: myindex is a simple, not so efficient sort routine that
C           had to be provided to replace calls to the more general
C           routine orders which is Cray specific
C
C**   INTERFACE.
C     ----------
C
C     *CALL* *myindex (ndim, input, index )
C         ndim: dimension of array input and index
C         input: input array to be ordered
C         index: output array containing the ordered element of input
C
C     METHOD.
C     -------
C
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
      integer input(ndim), index(ndim)
      
      DO j11=1,ndim
        index(j11) = j11
      ENDDO
      
      IF ( ndim .EQ. 1 ) RETURN
      
      ihlf = ndim /2 + 1
      ir = ndim
      
      DO WHILE (.TRUE.)
        IF ( ihlf .GT. 1 ) THEN
          ihlf = ihlf - 1
          indext = index(ihlf)
          inum = input(indext)
        ELSE
          indext = index(ir)
          inum = input(indext)
          index(ir) = index(1)
          ir = ir -1 
          IF (ir .EQ. 1 ) THEN
            index(1) = indext
            RETURN
          ENDIF
        ENDIF
        ihp = ihlf
        ibg = ihlf + ihlf
        DO WHILE ( ibg .LE. ir )
          IF (ibg .LT. ir ) THEN
            IF (input(index(ibg)) .LT. input(index(ibg+1))) ibg = ibg +1
          ENDIF
          IF (inum .LT. input(index(ibg))) THEN
            index(ihp)=index(ibg)
            ihp = ibg
            ibg = ibg + ibg
          ELSE
            ibg = ir + 1
          ENDIF     
        ENDDO
        index(ihp) = indext
      ENDDO
      RETURN
      END
