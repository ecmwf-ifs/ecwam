      SUBROUTINE toto(flmin,nfre,nang,nspec,a)
      real,dimension(nspec,nang,nfre) :: a,b


      do m=1,nfre
        write(*,*) a(ij,k,m)
        do k=1,nang
          do ij = 1,nspec
!!!          if(a(ij,k,m).gt.0) then
             b(ij,k,m)= exp(a(ij,k,m))
!!!          endif
          enddo
        enddo
      enddo
      return
      END subroutine toto
