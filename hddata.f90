!=================================================================================
MODULE HdDATA
  IMPLICIT NONE
  integer :: natoms
  INTEGER :: ncoord   !total number of internal coords
  INTEGER :: nstates  !number of electronic states
  real*8, allocatable :: refint(:)
CONTAINS
!---------------------------------------------------------------------------------
SUBROUTINE EvaluateHd(igeom,hmat,dhmat)
  use DiabaticHamiltonian
  IMPLICIT NONE
  real*8, intent(in) :: igeom(ncoord)
  DOUBLE PRECISION,DIMENSION(nstates,nstates),INTENT(OUT) :: hmat
  DOUBLE PRECISION,DIMENSION(ncoord,nstates,nstates),INTENT(OUT) :: dhmat
  real*8, allocatable :: q(:),t(:),dt(:,:)
  integer :: i,j

  allocate(q(Hd_intdim),t(NHdExpansionBasis))
  allocate(dt(Hd_intdim,NHdExpansionBasis))

  q(1:Hd_intdim)=igeom(1:Hd_intdim)

  do i=1,NHdExpansionBasis
    t(i)=ExpansionBasis(q,i)
    dt(:,i)=ExpansionBasisGradient(q,i)
  end do

  hmat=0.d0
  dhmat=0.d0
  do i=1,nstates
    do j=i,nstates
      hmat(i,j)=dot_product(t,Hd_HdEC(i,j)%Array)
      dhmat(1:Hd_intdim,i,j)=matmul(dt,Hd_HdEC(i,j)%Array)
      !fill up
      if(i.ne.j) then
        hmat(j,i)=hmat(i,j)
        dhmat(:,j,i)=dhmat(:,i,j)
      end if
    end do
  end do

  return
END SUBROUTINE EvaluateHd
!---------------------------------------------------------------------------------
END MODULE HdDATA
!=================================================================================
