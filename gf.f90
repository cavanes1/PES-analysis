!=============================================================================
program main
  use hddata
  use geomtrans
  implicit none
  real*8, parameter :: au2cm=219474.63067d0
  real*8, parameter :: eshift=284.717676772199d0
  real*8, parameter :: bohr2angs=0.529177210903d0
  real*8, parameter :: amu2au=1822.888486d0
  real*8, parameter :: freq2cm=1378999.78d0
  real*8, parameter :: pi=dacos(-1.d0)

  real*8, allocatable :: x(:,:)
  real*8, allocatable :: Bmat(:,:),intc(:)
  real*8 :: E0
  real*8, allocatable :: dE(:), d2E(:,:), hess(:,:)
  real*8, allocatable :: mass(:), freq(:), intmode(:,:), Linv(:,:),cartmode(:,:)

  character(3) :: sym
  real*8 :: anums
  real*8 :: Rot(3,3), h,t,tol
  integer :: i,j,k,ios,total,info

  call initpes

  allocate(x(3,natoms),Bmat(ncoord,3*natoms),intc(ncoord))
  allocate(dE(ncoord),d2E(ncoord,ncoord),hess(ncoord,ncoord))
  allocate(mass(natoms),freq(ncoord),intmode(ncoord,ncoord),Linv(ncoord,ncoord))
  allocate(cartmode(3*natoms,ncoord))

  open(unit=100,file='energy.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  !reference energy
  read(100,*,iostat=ios) E0
  !ei displacement
  do i=1,ncoord
    read(100,*,iostat=ios) dE(i)
  end do
  !ei+ej displacement
  do i=1,ncoord
    do j=i,ncoord
      read(100,*,iostat=ios) d2E(i,j)
    end do
  end do

  h=1.d-3
  !finite difference hessian
  do i=1,ncoord
    do j=i,ncoord
      hess(i,j)=(d2E(i,j)-dE(i)-dE(j)+E0)/h**2
      if(i.ne.j) hess(j,i)=hess(i,j)
    end do
  end do

  open(unit=101,file='geom.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  do i=1,natoms
    read(101,*,iostat=ios) sym,anums,x(:,i),mass(i)
  end do

  call WilsonBMatrixAndInternalCoordinate(x, Bmat, intc, 3*natoms, ncoord)

  mass=mass*amu2au
  call WilsonGFMethod(Hess,Bmat,mass,freq,intmode,Linv,cartmode,ncoord,natoms)
  freq=freq*freq2cm/(2.d0*pi)

  write(*,"(f15.6)"),freq

  open(200,file='hess')
  do i=1,ncoord
    write(200,"(<ncoord>e18.8)") hess(i,:)
  end do

  stop
end
!=============================================================================
