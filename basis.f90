!========================================================================================
program main
  use hddata
  use geomtrans
  use NadVibSInterface
  implicit none
  real*8, parameter :: amu2au=1822.888486d0
  real*8, allocatable :: qPrecursor(:), qResidual(:), freqPrecursor(:), freqResidual(:)
  real*8, allocatable :: LinvPrecursor(:,:), LinvResidual(:,:), intmodeResidual(:,:)
  real*8, allocatable :: intmode(:,:),cartmode(:,:)
  real*8, allocatable :: geom(:,:),hess(:,:),Bmat(:,:),mass(:)

  character(3) :: sym
  real*8 :: anums
  integer :: ios,i,fid

  namelist /basisest/ NVS_contour 

  call FLUnit(fid)
  open(fid,file='basisest.in',delim='APOSTROPHE')
  read(unit=fid,nml=basisest)
  close(fid)

  call initpes

  allocate(mass(natoms))
  allocate(qPrecursor(ncoord),qResidual(ncoord),freqPrecursor(ncoord),freqResidual(ncoord))
  allocate(LinvPrecursor(ncoord,ncoord),LinvResidual(ncoord,ncoord))
  allocate(intmodeResidual(ncoord,ncoord))
  allocate(intmode(ncoord,ncoord),cartmode(3*natoms,ncoord))
  allocate(geom(3,natoms),hess(ncoord,ncoord),Bmat(ncoord,3*natoms))

  !Precursor
  open(unit=100,file='precursor.xyz',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  open(unit=101,file='precursor.hess',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  do i=1,natoms
    read(100,*,iostat=ios) sym,anums,geom(:,i),mass(i)
  end do
  do i=1,ncoord
    read(101,*,iostat=ios) hess(i,:)
  end do
  mass=mass*amu2au

  call WilsonBMatrixAndInternalCoordinate(geom, Bmat, qPrecursor, 3*natoms, ncoord)
  call WilsonGFMethod(Hess,Bmat,mass,freqPrecursor,intmode,&
                      LinvPrecursor,cartmode,ncoord,natoms)

  !Residual
  open(unit=200,file='residual.xyz',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  open(unit=201,file='residual.hess',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  do i=1,natoms
    read(200,*,iostat=ios) sym,anums,geom(:,i),mass(i)
  end do
  do i=1,ncoord
    read(201,*,iostat=ios) hess(i,:)
  end do
  mass=mass*amu2au

  call WilsonBMatrixAndInternalCoordinate(geom, Bmat, qResidual, 3*natoms, ncoord)
  call WilsonGFMethod(Hess,Bmat,mass,freqResidual,intmodeResidual,&
                      LinvResidual,cartmode,ncoord,natoms)

  call BasisEstimation(qPrecursor,freqPrecursor,LinvPrecursor,qResidual,freqResidual,&
                           LinvResidual,intmodeResidual,ncoord)

  stop
end
!========================================================================================
