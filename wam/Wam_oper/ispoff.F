      FUNCTION ispoff (klonrgg, kgy, kang, kfre, ky, kx)
c
c  Calculate the offset for a given kxlt and ixlg so as
c  to place a 2d spectrum in a grid.
c
c

      INTEGER klonrgg(kgy)
      INTEGER kgy, kang, kfre, ky, kx
c
      io = 0
      DO ji = 1, ky
        io = io + klonrgg(ji)
      ENDDO
      ispoff = (io + kx-1) * (kang*kfre)
c
      RETURN
      END
