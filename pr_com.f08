
program RDFcom
! From the group of Prof. Philip Camp
! subroutines by Rui Apóstolo in 2017 & 2018
! University of Edinburgh
use, intrinsic :: iso_fortran_env
implicit none
include 'variables.inc'



! notes:
! pbc condition can be improved by reading box dimensions each timestep and changing Lz into an array L[3]


debug = .false.
if (debug .eqv. .true.) open(unit=93, file='debug.log', action='write', status='replace')

! New clock start
! going to ignore the array, it gives user and system time in seconds, in this order.
call ETIME(tarray, beg_cpu_time)
cpu_time_last = 0.0_dp

open(unit=70, file='params.in', action='read', status='old')
read(70,*)
read(70,'(a)',iostat=ierror)  Inputfile
  if (ierror /= 0) call systemexit("Input file")
read(70,'(I2)',iostat=ierror)      Row
  if (ierror /= 0) call systemexit("Row")
read(70,'(I6)',iostat=ierror)      StepMax
  if (ierror /= 0) call systemexit("StepMax")
read(70,'(I6)',iostat=ierror)      IgnoreFirst
  if (ierror /= 0) call systemexit("IgnoreFirst")
read(70,'(I6)',iostat=ierror)      MolSize
  if (ierror /= 0) call systemexit("MolSize")
read(70,'(I6)',iostat=ierror)      NumFGs
  if (ierror /= 0) call systemexit("NumFGs")
read(70,'(I6)',iostat=ierror)      NumTypesFG
  if (ierror /= 0) call systemexit("NumTypesFG")
read(70,'(I6)',iostat=ierror)      SizeFG
  if (ierror /= 0) call systemexit("SizeFG")
read(70,'(I6)',iostat=ierror)      CTFG
  if (ierror /= 0) call systemexit("CTFG")
read(70,'(F10.6)',iostat=ierror)      BoxLength
  if (ierror /= 0) call systemexit("Box length")
read(70,'(I6)',iostat=ierror) ! solvent_num
  if (ierror /= 0) call systemexit("solvent_num")
read(70,'(F10.6)',iostat=ierror)   SizeBins
  if (ierror /= 0) call systemexit("SizeBins")
read(70,'(F10.6)',iostat=ierror)   LowBin
  if (ierror /= 0) call systemexit("LowBin")
read(70,'(I6)',iostat=ierror)      NumBins
  if (ierror /= 0) call systemexit("NumBins")
read(70,'(L6)',iostat=ierror)      calc2dpr
  if (ierror /= 0) call systemexit("Calculate 2D p(r)")
close(70)


! Input sanity checks
if (Inputfile == "") then
  write(6,*)      "No name for input file given in params.in"
  ! see /usr/include/sysexits.h for error codes
  call EXIT(65)
end if
if (Row < 5) then
  write(6,*)      "'Row' line in params.in must be an integer bigger than 4. X, Y and Z should be on rows 3, 4 and 5, respectively."
  call EXIT(65)
end if
if (StepMax < 0) then
  write(6,*)      "'StepMax' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (StepMax < 1001) then
  write(6,*) StepMax," steps is an unusually low number for StepMax, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if
if (StepMax-IgnoreFirst > 19999) then
  write(6,*) StepMax-IgnoreFirst," steps is an unusually high number of steps, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if
if (StepMax-IgnoreFirst < 1) then
  write(6,*)      "StepMax - IgnoreFirst  must be an integer bigger than 0."
  call EXIT(65)
else if (StepMax-IgnoreFirst < 4999) then
  write(6,*) StepMax-IgnoreFirst," steps is an unusually low number of steps, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
else if (IgnoreFirst < 0) then
  write(6,*)      "'IgnoreFirst' line in params.in must be an integer bigger than or equal to 0."
  call EXIT(65)
end if
if (MolSize < 1) then
  write(6,*)      "'MolSize' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (MolSize < 101) then
  write(6,*) MolSize," atoms is an unusually low number for MolSize, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
else if (MolSize > 1999) then
  write(6,*) MolSize," atoms is an unusually high number for MolSize, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if
if (NumFGs < 1) then
  write(6,*)      "'NumFGs' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (NumFGs > 20) then
  write(6,*) NumFGs," FGs is an unusually high number for NumFGs, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if
if (NumTypesFG < 1) then
  write(6,*)      "'NumTypesFG' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (NumTypesFG > 20) then
  write(6,*) NumTypesFG," Types of atoms in the FGs is an unusually high number for NumTypesFG, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if
if (SizeFG < 1) then
  write(6,*)      "'SizeFG' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (SizeFG < NumTypesFG) then
  write(6,*)      "'SizeFG' line in params.in must be an integer bigger than or equal to 'NumTypesFG'."
  call EXIT(65)
else if (SizeFG > 40) then
  write(6,*) SizeFG," atoms in the FGs is an unusually high number for SizeFG, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if
if (CTFG < 1) then
  write(6,*)      "'CTFG' line in params.in must be an integer bigger than 0."
  call EXIT(65)
end if
if (BoxLength <= 0.0) then
  write(6,*)      "'BoxLength' line in params.in must be a real number bigger than 0."
  call EXIT(65)
else if (BoxLength < 50.0) then
  write(6,*) BoxLength," angstrom is an unusually small size for BoxLength, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
else if (BoxLength > 200.0) then
  write(6,*) BoxLength," angstrom is an unusually large size for BoxLength, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if
! if (solvent_num < 1) then
  ! write(6,*)      "'solvent_num' line in params.in must be an integer bigger than 0"
  ! call EXIT(65)
! end if 
if (SizeBins <= 0.0) then
  write(6,*)      "'SizeBins' line in params.in must be a decimal bigger than 0.0."
  call EXIT(65)
else if (SizeBins > 1.0) then
  write(6,*) SizeBins," as a bin size might make you loose resolution, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if
if (LowBin < 0.0) then
  write(6,*)      "'LowBin' line in params.in must be a decimal bigger than or equal to 0.0."
  call EXIT(65)
else if (LowBin > 0.0) then
  write(6,*) LowBin," as the lowest bin might make you lose the closest FGcoms, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if
if (NumBins < 1) then
  write(6,*)      "'NumBins' line in params.in must be an integer bigger than 0."
  call EXIT(65)
else if (NumBins < 100) then
  write(6,*) NumBins," bins is an unusually low number for NumBins, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
else if (NumBins > 2999) then
  write(6,*) NumBins," bins is an unusually high number for NumBins, make sure it's correct."
  write(6,*) "Do you wish to continue anyway? (y/N)"
  read(5,*)  checkrow
  call userexit(checkrow)
end if









write(6,*)
write(6,*)      "Input file chosen:                  ",Inputfile
write(6,*)
write(6,*)      "Number of Rows in input file:     ",Row
write(6,*)      "Last step number:                 ",StepMax
write(6,*)      "Ignoring first:                   ",IgnoreFirst," timesteps"
write(6,*)      "Molecule Size:                    ",MolSize
write(6,*)      "Number of FGs:                    ",NumFGs
write(6,*)      "Number of different atoms in FGs: ",NumTypesFG
write(6,*)      "Number of atoms in FG:            ",SizeFG
write(6,*)      "Type number for first atom in FG: ",CTFG
write(6,*)      "Size of simulation box",BoxLength
! write(6,*) ! solvent_num
write(6,*)      "Size of histogram Bin:                    ",SizeBins
write(6,*)      "Lowest histogram bin:                     ",LowBin
write(6,*)      "Total number of bins:             ",NumBins
write(6,*)
write(6,*)
if (calc2dpr .eqv. .true.) then
  write(6,*)      "Calculating both 2D and 3D p(r)."
  allocate(results2d(NumBins,2))
else if (calc2dpr .eqv. .false.) then
  write(6,*)      "Calculating only 3D p(r)."
else
  write(6,*)      "Calculate 2D p(r) in params.in must be .true. or .false."
  ! see /usr/include/sysexits.h for error codes
  call EXIT(65)
end if

write(6,*)
write(6,*)
write(6, "(' Program started ',A)") adjustl(ctime(time()))
write(6,*)
write(6,*)

FGcount = 0_sp
allocate(Masses(1:100))
allocate(array(Row))
array = 0.0_dp
allocate(FGcoms(3,NumFGs))
allocate(IniBins(NumBins))


open(unit=71, file='masses.in', action='read', status='old')
read(71,*)
do i=0,NumTypesFG-1
  read(71,*) j, Masses(j)
end do
close(71)

open(unit=72, file=inputfile, action='read', status='old')
open(unit=81, file='3d_pr_fgcom.out', action='write', status='new')
if (calc2dpr .eqv. .true.) open(unit=82, file='2d_pr_fgcom.out', action='write', status='new')


! initialise NumBins
do i=1,NumBins
  IniBins(i) = real(LowBin) + real(i)*SizeBins
end do
write(81,*) IniBins
write(82,*) IniBins


write(6, *) "Ignoring  ",IgnoreFirst,"  timesteps"
write(6,*)
write(6,*)

! Loop over ignored input
do i=1,IgnoreFirst
  call readheader
  do j=1,MolSize
    read(72,*)
  end do
end do



call ETIME(tarray,cpu_time_now)
call texttime(int(cpu_time_now-cpu_time_last),tsl)
write(6, *) "Ignored   ",IgnoreFirst,"  timesteps, took ",tsl
write(6,*)
write(6,*)
cpu_time_last = cpu_time_now

! Loop over useful timesteps
outerloop: do i=IgnoreFirst+1,StepMax
  ! variable reset
  FGcount = 0_sp
  FGcoms = 0.0_dp
  results2d = 0.0_dp

  ! call clock and output to screen
  if (modulo(i,100) == 0) then
    ! call cpu_time(cpu_time_now) ! old clock
    call ETIME(tarray,cpu_time_now) ! new clock
    call texttime(int(cpu_time_now-cpu_time_last),tsl)
    call texttime(int(cpu_time_now-beg_cpu_time),tss)
    call texttime(int((cpu_time_now-cpu_time_last)*(((stepmax-i)/100)+1)),eta)
    write(6,*) "Timestep: ", i, &
               "  This iteration: ", tsl, &
               "Running for: ", tss, &
               "ETA: ", eta
    cpu_time_last = cpu_time_now
  end if

  call readheader

  ! Loop over this timestep's atoms
  innerloop: do j=1,(MolSize-NumFGs*(SizeFG-1))
    read(72,*) (array(k),k=1,Row)
    ! Find start of FG
    if (abs(array(2) - real(CTFG,dp)) .LE. 0.01_dp ) then
      TotalMass = 0.0_dp
      FGcount = FGcount + 1
      ! assign atommass*x, atommass*y and atommass*z to FGcom(1:3,currentFGcount)
      FGcoms(1:3,FGcount) = FGcoms(1:3,FGcount) + Masses(int(array(2)))*array(3:5)
      ! Add atommass to totalmass
      TotalMass = TotalMass + Masses(int(array(2)))
      if (debug .eqv. .true.) write(93,*) 'type ',array(2),'mass ',Masses(int(array(2))), 'fgcoms: ',FGcoms(1:3,FGcount)
      ! Loop over rest of FG
      fgloop: do l=2,SizeFG
        read(72,*) (array(k),k=1,Row)
        ! add this atom's atommass*x, atommass*y and atommass*z to FGcom(1:3,currentFGcount)
        FGcoms(1:3,FGcount) = FGcoms(1:3,FGcount) + Masses(int(array(2)))*array(3:5)
        ! Add atommass to totalmass
        TotalMass = TotalMass + Masses(int(array(2)))
        if (debug .eqv. .true.) write(93,*) 'type ',array(2),'mass ',Masses(int(array(2))), ' fgcoms: ',FGcoms(1:3,FGcount)
      end do fgloop
      FGcoms(1:3,FGcount) = FGcoms(1:3,FGcount) / TotalMass
      if (debug .eqv. .true.) write(93,*) 'after fgloop fgcoms: ',FGcoms(1:3,FGcount), ' tmass: ',TotalMass
    end if
  end do innerloop      
        if (debug .eqv. .true.) write(93,*) "before resultscalc, FGcoms+ ",FGcoms," NumBins= ",NumBins," SizeBins= ",SizeBins
      if (calc2dpr .eqv. .true.) then
        ! call prcalc2
        results2d = prcalc2d(FGcoms,NumBins,SizeBins,BoxLength)
        ! write 3d results
        write(81,*) (results2d(k,1),k=1,NumBins)
        ! write 2d results
        write(82,*) (results2d(k,2),k=1,NumBins)
      else
        ! write only 3d results
        write(81,*) prcalc3d(FGcoms,NumBins,SizeBins,BoxLength)
      end if
end do outerloop



! End of program clock call
call ETIME(tarray,end_cpu_time)
call texttime(int(end_cpu_time - beg_cpu_time), tempw)
write(6,*) "Total Time: ", ADJUSTL(tempw)



contains

!********************************
! Read vmd header - Rui Apóstolo
!********************************

subroutine readheader
implicit none
integer :: i
do i=1,9
  read(72,*)
end do
end subroutine readheader


!**************************************
! Calculate only 3D P(r) - Rui Apóstolo
!**************************************
! prcalc3d(FGcoms,Numbins,Sizebins,BoxLength)
! FGcom is in the form (1:3,SizeFG), 1:3 being xyz coords
function prcalc3d(innerarray,NB,SB,Lx)
implicit none
real(dp),intent(in),dimension(:,:) :: innerarray
real(dp),intent(in) :: SB,Lx
integer,intent(in) :: NB
real(dp),dimension(NB) :: prcalc3d
real(dp),dimension(3) :: xyz
integer :: i,j,bin
real(dp) :: r,add

prcalc3d = 0.0_dp
! 1 / ((N * N-1)/2)
add = 1.0_dp/((real(size(innerarray,2),dp)*(real(size(innerarray,2),dp)-1))/2)


do i=1,size(innerarray,2)-1
  do j = i+1,size(innerarray,2)
      ! xyz(1:3) = innerarray(1:3,i) - innerarray(1:3,j)
      xyz(1:3) = innerarray(1:3,i) - Lx*anint(innerarray(1:3,j)/Lx,dp)
      r = sqrt(xyz(1)**2 + xyz(2)**2 + xyz(3)**2)
      bin = int(r/SB)+1
      if (bin <= NB) then
        prcalc3d(bin) = prcalc3d(bin) + add
      else
        write(6,*) "Distance between FGcoms bigger than NumBins * SizeBins."
      end if
  end do
end do
end function prcalc3d


!****************************************
! Calculate 2D and 3D P(r) - Rui Apóstolo
!****************************************
! prcalc2d(FGcoms,Numbins,Sizebins)
! FGcom is in the form (1:3,SizeFG), 1:3 being xyz coords
function prcalc2d(innerarray,NB,SB,Lx)
implicit none
real(dp),intent(in),dimension(:,:) :: innerarray
real(dp),intent(in) :: SB,Lx
integer,intent(in) :: NB
real(dp),dimension(NB,2) :: prcalc2d
real(dp),dimension(3) :: xyz
integer :: i,j,bin,bin2d
real(dp) :: r,r2d,add

prcalc2d = 0.0_dp
! 1 / ((N * N-1)/2)
add = 1.0_dp/((real(size(innerarray,2),dp)*(real(size(innerarray,2),dp)-1))/2)


do i=1,size(innerarray,2)-1
  do j = i+1,size(innerarray,2)
      ! ri - rj
      ! xyz(1:3) = innerarray(1:3,i) - innerarray(1:3,j
      xyz(1:3) = innerarray(1:3,i) - Lx*anint(innerarray(1:3,j)/Lx,dp)
      ! calculating r
      r = sqrt(xyz(1)**2 + xyz(2)**2 + xyz(3)**2)
      r2d = sqrt(xyz(1)**2 + xyz(2)**2)
      bin = int(r/SB)+1
      bin2d = int(r2d/SB)+1
      ! assign 3D bin
      if (bin <= NB) then
        prcalc2d(bin,1) = prcalc2d(bin,1) + add
      else
        write(6,*) "Distance between FGcoms bigger than NumBins * SizeBins."
      end if
      ! assign 2D bin
      if (bin2d <= NB) then
        prcalc2d(bin2d,2) = prcalc2d(bin2d,2) + add
      else
        write(6,*) "Distance between FGcoms bigger than NumBins * SizeBins."
      end if
  end do
end do
end function prcalc2d


!************************************************
! Time conversion - Written by Rui Apóstolo
! Converts seconds to hours, minutes and seconds
!************************************************

subroutine texttime(seconds,ttime)
  integer(kind=4)::seconds
  integer::s,m,h,r
  character(LEN=20)::ttime,ttemp
  ttime = ""
  ttemp = ""
  s = 0
  r = 0
  m = 0
  h = 0
  ! 1234 h 56 min 78 sec
  s = modulo(seconds,60)
  r = (seconds - s)/60
  m = modulo(r,60)
  h = (r - m)/60
  if (h == 0) then
    if (m == 0) then
      write(ttemp, '(I2.2, A)') s, " seconds"
    else
      write(ttemp, '(I2.2, A, I2.2, A)') m, " min ", s, " sec"
    end if
  else
    write(ttemp, '(I4.1, A, I2.2, A, I2.2, A)') h, " h ", m, " min ", s, " sec"
  end if
  write(ttime, '(A)') ADJUSTL(ttemp)
end subroutine texttime


!********************************************
! Exits program by user choice - Rui Apóstolo
!********************************************

subroutine userexit(answer)
  character,intent(in) :: answer*1
  select case(answer)
    case("y")
      write(6,*) "Program continuing by your choice."
    case("Y")
      write(6,*) "Program continuing by your choice."
    case default
      write(6,*) "Program exit by your choice."
      call EXIT(0)
  end select
end subroutine userexit

!*****************************************************
! Exits program due to data input error - Rui Apóstolo
!*****************************************************

subroutine systemexit(error)
  character(len=*),intent(in) :: error
    write(6,*) "Program exit due to error on data input in variable: ",error
    ! see /usr/include/sysexits.h for error codes
    call EXIT(65)
end subroutine systemexit


end program RDFcom