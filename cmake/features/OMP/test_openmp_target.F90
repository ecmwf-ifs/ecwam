program test_openmp_target
   
    implicit none

    integer :: a(32), i

    a = 1
    !$omp target map(tofrom:a)
    do i=1,32
        a(i) = a(i) + 2
    enddo
    !$omp end target

    if (.not. all(a == 3)) error stop

end program test_openmp_target
