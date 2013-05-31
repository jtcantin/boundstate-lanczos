        program tester
         integer i
         character(len=32) arg

         do i = 1, argc()
            call getarg(i, arg)
            write(*,*) arg
         end do
        end program
