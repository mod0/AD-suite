program test_mat_read
  real, dimension(2,2,2) :: a
  open(11, file = 'a.dat', form = 'unformatted')
  read(11) a
  write(*,*) a
  close(11)
end program test_mat_read
