      FUNCTION ispoff (klonrgg, kgy, kang, kfre, ky, kx)

!  Calculate the offset for a given kxlt and ixlg so as
!  to place a 2d spectrum in a grid.



      INTEGER klonrgg(kgy)
      INTEGER kgy, kang, kfre, ky, kx

      io = 0
      DO ji = 1, ky
        io = io + klonrgg(ji)
      ENDDO
      ispoff = (io + kx-1) * (kang*kfre)

      RETURN
      END FUNCTION ispoff




