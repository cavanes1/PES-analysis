!=============================================================================
program main
  use hddata
  use geomtrans
  implicit none
  real*8, parameter :: au2cm=219474.63067d0
  !real*8, parameter :: eshift=256.787847331183d0
  real*8, parameter :: bohr2angs=0.529177210903d0
  !General options
  DOUBLE PRECISION :: eshift    ! uniform shift on ab initio energies

  real*8 :: Rot(3,3)

  real*8, dimension(:), allocatable :: e, et
  real*8, dimension(:,:), allocatable :: hd, ckl
  real*8, dimension(:,:,:), allocatable :: geom
  real*8, dimension(:,:,:), allocatable :: cg, dcg

  character(7), allocatable :: ids(:)
  real*8, allocatable :: energy(:,:)
  real*8, allocatable :: x(:,:)
  real*8, allocatable :: Bmat(:,:),intc(:)
  real*8, allocatable :: grd(:)

  real*8 :: dx
  character(3) :: sym
  character(255) :: filename, npoints, enfDiab, epmax, w_energy, w_grad, w_fij, &
  gradcutoff, cpcutoff, deggrdbinding, deg_cap, lambda, &
  nrmediff, ediffcutoff, fixref ! dummy except for filename
  DOUBLE PRECISION,dimension(10) :: energyT, highEScale ! dummy
  real*8 :: anums,bsnum
  integer :: i,j,k,l,m,ios,total,info

  namelist /fitting/ npoints, enfDiab, epmax, w_energy, w_grad, w_fij, &
  gradcutoff, cpcutoff, deggrdbinding, deg_cap, lambda, eshift, energyT, &
  highEScale, nrmediff, ediffcutoff, fixref, natoms, nstates
  open(103,file='fit.in',delim='APOSTROPHE')
  read(unit=103,nml=fitting)
  close(103)

  call initpes

  allocate(e(nstates),et(nstates))
  allocate(hd(nstates,nstates),ckl(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  total=9999
  allocate(ids(total))
  allocate(geom(3,natoms,total))
  allocate(energy(nstates,total))
  allocate(x(3,natoms),Bmat(ncoord,3*natoms),intc(ncoord))

  open(unit=100,file='geom.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  open(unit=101,file='energy.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  open(unit=102,file='names.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  do i=1,9999
    read(101,*,iostat=ios) energy(:,i)
    read(102,*,iostat=ios) bsnum, ids(i)
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
  open(201,file='allenergies.csv')
  !open(202,file='xyz')
  open(203,file='fitener.dat')
  open(204,file='hd.dat')
  open(205,file='intc.dat')
  open(207,file='errener.dat')

  !resetting gradient output files
  do j=1, nstates
    write (filename, '("cartgrd.drt1.state", i0, ".dat")' ) j
    open (208,file=filename)
    write(208," ",advance="no")
    close (208,status="delete")
  end do

  !resetting NAC output files
  do j=1, nstates-1
    do k=j+1, nstates
      write (filename, '("cartgrd.nad.drt1.state", i0, ".drt1.state", i0, ".dat")' ) j,k
      open (209,file=filename)
      write(209," ",advance="no")
      close (209,status="delete")
    end do
  end do

  do i=1,total
    x=geom(:,:,i)
    call WilsonBMatrixAndInternalCoordinate(x, Bmat, intc, 3*natoms, ncoord)
    call Evaluate(x,e,cg,hd,dcg,ckl)
    write(201,"(A7,',',F18.12,',',F18.12,',',F18.12,',',F18.12,',',F18.12)") ids(i), e-eshift
    e=e*au2cm
    hd=hd*au2cm
    write(200,"(5es22.12)") energy(:,i)
    write(203,"(5es22.12)") e
    do j=1,nstates
      do k=j,nstates
        write(204,"(3es22.12)",advance="no") hd(j,k)
      end do
    end do
    write(204,*)
    write(205,"(27es15.6)") intc
    !full sig figs is 5f17.5
    write(207,"(5f12.0)",advance="no") e-energy(:,i)
    write(207,*) ids(i)
    
    !surface gradient printout
    do j=1, nstates
      write (filename, '("cartgrd.drt1.state", i0, ".dat")' ) j
      open (208,file=filename,position='append')
      do l=0,natoms-1
        do m=1,3
          write(208,"(33es15.6)",advance="no") cg(3*l+m,j,j)
        end do
        write(208,*)
      end do
      close (208)
    end do

    !surface NAC printout
    do j=1, nstates-1
      do k=j+1, nstates
        write (filename, '("cartgrd.nad.drt1.state", i0, ".drt1.state", i0, ".dat")' ) j,k
        open (209,file=filename,position='append')
        do l=0,natoms-1
          do m=1,3
            write(209,"(33es15.6)",advance="no") cg(3*l+m,j,k)
          end do
          write(209,*)
        end do
        close (209)
      end do
    end do

  end do

  stop
end
!=============================================================================
