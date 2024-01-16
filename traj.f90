!=============================================================================
program main
  use hddata
  implicit none
  real*8, parameter :: au2cm=219474.63067d0
  real*8, parameter :: bohr2angs=0.52917721067d0
  real*8, parameter :: ev2cm=8065.51d0
  character(3) :: sym
  real*8 :: anums
  real*8 :: Rot(3,3)
  integer :: i,j,k,l,ios,total,info,id
  integer :: num,idx,nt
  real*8 :: tx,vx
  real*8, allocatable :: geom(:,:,:),x(:,:),t(:),v(:)

  real*8, dimension(:), allocatable :: e
  real*8, dimension(:,:), allocatable :: hd,ckl
  real*8, dimension(:,:,:), allocatable :: cg, dcg

  total=99999
  call initpes
  allocate(geom(3,natoms,total),x(3,natoms),t(total),v(total))

  allocate(e(nstates))
  allocate(hd(nstates,nstates),ckl(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  read(*,*) id
  if(id.lt.0) stop 'id.lt.0!'
  print*, 'Analyze trajectory ', id

  open(unit=100,file='traj.all',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  total=0
  ios=0
  do while(ios.eq.0)
    read(100,*,iostat=ios) num
    read(100,*,iostat=ios) idx,nt,tx,vx
    read(100,*,iostat=ios) sym, x(:,1)
    read(100,*,iostat=ios) sym, x(:,2)
    read(100,*,iostat=ios) sym, x(:,3)
    read(100,*,iostat=ios) sym, x(:,4)
    read(100,*,iostat=ios) sym, x(:,5)
    read(100,*,iostat=ios) sym, x(:,6)
    read(100,*,iostat=ios) sym, x(:,7)
    read(100,*,iostat=ios) sym, x(:,8)
    read(100,*,iostat=ios) sym, x(:,9)
    read(100,*,iostat=ios) sym, x(:,10)
    read(100,*,iostat=ios) sym, x(:,11)
    if(ios.ne.0) exit
    if(idx.lt.id) cycle
    if(idx.eq.id) then
      total=total+1
      geom(:,:,total)=x
      t(total)=tx
      v(total)=vx
    end if
    if(idx.gt.id) exit
  end do
  print*,total, 'valid points in trajectory.'

  open(200,file='traj.ene.out')
  open(201,file='traj.xyz.out')
  do i=1,total
    x=geom(:,:,i)/bohr2angs
    call Evaluate(x,e,cg,hd,dcg,ckl)
    write(200,"(i5,4e14.5)") i,t(i),v(i)*ev2cm,e*au2cm

    call orientation(x,Rot,info)
    write(201,"(i4)") 11
    write(201,"(4e14.5)")  t(i),v(i)*ev2cm,e*au2cm
    write(201,"('C',2x,3f15.9)") x(:,1)*bohr2angs
    write(201,"('C',2x,3f15.9)") x(:,2)*bohr2angs
    write(201,"('C',2x,3f15.9)") x(:,3)*bohr2angs
    write(201,"('C',2x,3f15.9)") x(:,4)*bohr2angs
    write(201,"('C',2x,3f15.9)") x(:,5)*bohr2angs
    write(201,"('N',2x,3f15.9)") x(:,6)*bohr2angs
    write(201,"('H',2x,3f15.9)") x(:,7)*bohr2angs
    write(201,"('H',2x,3f15.9)") x(:,8)*bohr2angs
    write(201,"('H',2x,3f15.9)") x(:,9)*bohr2angs
    write(201,"('H',2x,3f15.9)") x(:,10)*bohr2angs
    write(201,"('H',2x,3f15.9)") x(:,11)*bohr2angs

  end do

  stop
end
!=============================================================================
