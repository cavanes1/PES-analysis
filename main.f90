!=============================================================================
program main
  use hddata
  use geomtrans
  implicit none
  real*8, parameter :: au2cm=219474.63067d0
  real*8, parameter :: eshift=284.717676772199d0
  real*8, parameter :: bohr2angs=0.529177210903d0

  real*8 :: Rot(3,3)

  real*8, dimension(:), allocatable :: e, et
  real*8, dimension(:,:), allocatable :: hd, ckl
  real*8, dimension(:,:,:), allocatable :: geom
  real*8, dimension(:,:,:), allocatable :: cg, dcg

  real*8, allocatable :: energy(:,:)
  real*8, allocatable :: x(:,:)
  real*8, allocatable :: Bmat(:,:),intc(:)
  real*8, allocatable :: grd(:)

  real*8 :: dx
  character(3) :: sym
  real*8 :: anums
  integer :: i,j,k,ios,total,info

  call initpes

  allocate(e(nstates),et(nstates))
  allocate(hd(nstates,nstates),ckl(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  total=9999
  allocate(geom(3,natoms,total))
  allocate(energy(nstates,total))
  allocate(x(3,natoms),Bmat(ncoord,3*natoms),intc(ncoord))

  open(unit=100,file='geom.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  open(unit=101,file='energy.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  do i=1,9999
    read(101,*,iostat=ios) energy(:,i)
    do j=1,natoms
      read(100,*,iostat=ios) sym,anums,geom(:,j,i)
    end do
    if(ios.ne.0) then
      total=i-1
      exit
    end if
  end do

  print*,total, 'valid points.'
  energy=(energy+eshift)*au2cm

  open(200,file='energy.dat')
  open(201,file='quasi.dat')
  open(202,file='xyz')
  do i=1,total
    x=geom(:,:,i)
    call WilsonBMatrixAndInternalCoordinate(x, Bmat, intc, 3*natoms, ncoord)
    call Evaluate(x,e,cg,hd,dcg,ckl)
    e=e*au2cm
    hd=hd*au2cm
    write(200,"(34e15.6)") energy(:,i),e,hd(1,1),hd(2,2),hd(1,2),intc

  end do

  stop
end
!=============================================================================
