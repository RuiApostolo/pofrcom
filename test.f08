program test
integer :: row,k
real,dimension(5) :: array
character :: header*60

row = 5

open(unit=72, file='processed.lammpstrj', action='read', status='old')

do j=1,15000

  write(6,*) 'header'
  do i=1,9
    read(72,*) 
    ! write(6,*)
  end do
  
  write(6,*) 'atoms from iteration ',j
  do i=1,(737-8*21)
    read(72,*) (array(k),k=1,Row)
    write(6,*) (array(k),k=1,Row)
    if (abs(array(2) - real(5)) .LE. 0.01) then
      read(72,*) (array(k),k=1,Row)
      write(6,*) (array(k),k=1,Row)
      do l=2,21
        read(72,*) (array(k),k=1,Row)
        write(6,*) (array(k),k=1,Row)
      end do
    end if
end do

end do


end program test
