! real kind parametrs
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)
! integer, parameter :: qp = selected_real_kind(33, 4931)


! integers
integer(sp) :: CTFG,NumFGs,SizeFG,NumTypesFG,Row,StepMax,IgnoreFirst,MolSize
integer(sp) :: FGcount
integer(sp) :: i,j,k,l
integer(sp) :: ierror


! reals
real(dp) :: TotalMass,BoxLength

! real arrays
real(dp), allocatable :: Masses(:),array(:)
real(dp), allocatable :: FGcoms(:,:),results2d(:,:)


! characters
! character :: Hash*1
character :: Inputfile*32, checkrow*1

! Logical
LOGICAL :: calc2dpr,debug

! Bin variables
integer(sp) :: NumBins
real(dp) :: SizeBins,LowBin
real(dp),allocatable :: IniBins(:)

! clock variables
! real(dp) ::  beg_cpu_time, end_cpu_time, cpu_time_now, cpu_time_last
! character(LEN=20)::tempw, tsl, tss, eta

! New clock variables
! going to ignore the array, it gives user and system time in seconds, in this order.
real, dimension(2) :: tarray 
real(sp) ::  beg_cpu_time, end_cpu_time, cpu_time_now, cpu_time_last
character(LEN=20)::tempw, tsl, tss, eta