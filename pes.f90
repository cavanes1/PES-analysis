!====================================================================================
subroutine initpes
  use hddata
  use GeomTrans
  use DiabaticHamiltonian
  implicit none
  real*8, allocatable :: geom0(:),Bmat(:,:)
  integer :: nvibs
  character(3) :: sym
  real*8 :: anums
  integer :: i,j,fid,ios

  !natoms=6
  print*,natoms
  nvibs=3*natoms-6
  print*,nvibs
  !nstates=5

  !read coorninate definition
  ncoord=DefineInternalCoordinate()
  if(ncoord.ne.nvibs) stop 'ncoord.ne.nvibs!'
  print*,'System has ',ncoord, 'internal degrees of freedom.'

  !print coordinate definition
  do i=1,ncoord
    do j=1,GeometryTransformation_IntCoordDef(i)%NMotions
      write(*,"(2i5,2x,A10,2x,e12.5,4i5)") i,j,&
                GeometryTransformation_IntCoordDef(i)%motion(j)%type,&
                GeometryTransformation_IntCoordDef(i)%motion(j)%coeff,&
                GeometryTransformation_IntCoordDef(i)%motion(j)%atom
    end do
  end do

  !for coordinate transformation
  ncart=3*natoms
  nintc=ncoord
  allocate(trgt(nintc))

  !Initialize Diabatic Hamiltonian 
  call InitializeDiabaticHamiltonian(nstates,ncoord)
  print*,'The expansion of each element of DPEM has ',NHdExpansionBasis,' terms.'

  !Read HdExpansion Coefficients
  call ReadHdExpansionCoefficients(Hd_HdEC)

  !reference geometry
  allocate(refint(ncoord),geom0(3*natoms),Bmat(ncoord,3*natoms))
  call FLUnit(fid)
  open(unit=fid,file='refgeom',access='sequential',form='formatted',&
       STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  do i=1,natoms
    read(fid,*,iostat=ios) sym,anums,geom0(3*i-2:3*i)
  end do
  call WilsonBMatrixAndInternalCoordinate(geom0, Bmat, refint, 3*natoms, ncoord) 
  close(fid)

  return
end subroutine initpes
!======================================================================================
subroutine getinfo(nat,nst)
  use hddata, only: natoms,nstates
  implicit none
  integer, intent(out) :: nat,nst
  nat=natoms
  nst=nstates
  return
end
!======================================================================================
subroutine Evaluate(geom,e,cg,h,dcg,ckl)
  use hddata
  use geomtrans
  implicit none
  real*8, intent(in) :: geom(3,natoms)
  real*8, intent(out) :: e(nstates)
  real*8, intent(out) :: h(nstates,nstates)
  real*8, intent(out) :: cg(3*natoms,nstates,nstates)
  real*8, intent(out) :: dcg(3*natoms,nstates,nstates)
  real*8, intent(out) :: ckl(nstates,nstates) !AtD transformation
  real*8, allocatable :: igeom(:),dh(:,:,:),Bmat(:,:),work(:)
  real*8 :: de
  integer :: i,j,info,lwork

  lwork=999
  allocate(igeom(ncoord),dh(ncoord,nstates,nstates),Bmat(ncoord,3*natoms))
  allocate(work(lwork))

  call WilsonBMatrixAndInternalCoordinate(geom, Bmat, igeom, 3*natoms, ncoord)
  igeom=igeom-refint
  call EvaluateHd(igeom,h,dh)

  do i=1,nstates
     do j=i,nstates
       call dgemv('T',ncoord,3*natoms,1.d0,Bmat,ncoord,dh(:,i,j),1,0.d0,dcg(:,i,j),1)
       if(j.ne.i) dcg(:,j,i)=dcg(:,i,j)
     end do
  end do

  ckl=h
  call dsyev('V','U',nstates,ckl,nstates,e,work,lwork,info)
  if(info.ne.0) stop 'Failed to solve eigenvectors in Evaluate!'

  do i=1,3*natoms
    cg(i,:,:)=matmul(transpose(ckl),matmul(dcg(i,:,:),ckl))
  end do

  !derivative coupling
  do i=1,nstates-1
    do j=i+1,nstates
      de=e(j)-e(i)
      if(abs(de).lt.1d-30) de=1d-30
      cg(:,i,j)=cg(:,i,j)/de
      cg(:,j,i)=-cg(:,i,j)
    end do
  end do

  return
end
!=================================================================================
subroutine prepot
  implicit none
  call initpes
  return
end
!====FOR ANT POT==================================================================
subroutine pot(igrad,x,uu,guu,vv,gvv,dvec,ccph,repflag)
  use hddata
  implicit none
  integer, intent(in) :: igrad,repflag
  real*8, intent(in) :: x(3,natoms)
  real*8, intent(out) :: uu(nstates,nstates),guu(3,natoms,nstates,nstates)
  real*8, intent(out) :: vv(nstates),gvv(3,natoms,nstates)
  real*8, intent(out) :: dvec(3,natoms,nstates,nstates)
  real*8, intent(out) :: ccph(nstates,nstates)
  real*8, allocatable :: cg(:,:,:),dcg(:,:,:)
  integer :: i,j,k

  allocate(cg(3*natoms,nstates,nstates),dcg(3*natoms,nstates,nstates))

  call Evaluate(x,vv,cg,uu,dcg,ccph)

  do i=1,natoms
    do j=1,nstates
      gvv(:,i,j)=cg(3*i-2:3*i,j,j)
      do k=1,nstates
        guu(:,i,j,k)=dcg(3*i-2:3*i,j,k)
        dvec(:,i,j,k)=cg(3*i-2:3*i,j,k)
      end do
    end do
  end do

  return
end
!==================================================================================
