!Generate nadvibs.in for NadVibS
!The main procedure is:
!    1, read in the precursor information:
!           geometry through precursor.xyz
!           Hessian
!       then compute precursor vibration
!    2, select a geometry as origin and define the normal mode
!           when executing after Analyze-min, use that minimum and corresponding normal mode
!           when executing after Analyze-mex, use that mex and the normal mode of mean field
!    3, shift from old reference to the selected origin, then transform from
!       the internal coordinate adopted in fitting to normal mode
!An estimation of number of basis to use in NadVibS is also provided
module NadVibSInterface
    use DiabaticHamiltonian
    use GeomTrans
    use hddata, only: natoms,nstates,refint
    implicit none
    integer :: InternalDimension,CartesianDimension

!Parameter
    real*8::NVS_contour=10d0!Controls how large the eclipse to cover, see BasisEstimation for details

!Derived type
    !Example: type(NVS_HdExpansionCoefficient),allocatable,dimension(:,:)::NVS_HdEC
    !         NVS_HdEC(jstate,istate).Order(iorder).Array(i) is the i-th expansion coefficient
    !         in iorder-th order terms for Hd(jstate,istate)
    type NVS_HdExpansionCoefficient!Store Hd expansion coefficient in NadVibS format
        type(d2PArray),allocatable,dimension(:)::Order
    end type NVS_HdExpansionCoefficient

    type NVS_ExpansionBasisNumberingRule
        type(i2PArray),allocatable,dimension(:)::Number
    end type NVS_ExpansionBasisNumberingRule

!NadVibSInterface module only variable
    integer::NVS_NOrder
    integer,allocatable,dimension(:)::NVS_NumberOfEachOrderTerms
    type(NVS_ExpansionBasisNumberingRule),allocatable,dimension(:)::NVS_EBNR
    type(NVS_HdExpansionCoefficient),allocatable,dimension(:,:)::NVS_HdEC

contains
subroutine GenerateNadVibSInput()
    !Vibration information of precursor and residual
    real*8,dimension(InternalDimension)::qPrecursor,qResidual,freqPrecursor,freqResidual
    real*8,dimension(CartesianDimension)::rPrecursor,rResidual
    real*8,dimension(InternalDimension,InternalDimension)::HPrecursor,intmodePrecursor,LinvPrecursor,HResidual,intmodeResidual,LinvResidual
    real*8,dimension(InternalDimension,CartesianDimension)::BPrecursor,BResidual
    real*8,dimension(CartesianDimension,InternalDimension)::cartmodePrecursor,cartmodeResidual
    !Origin shift in nadvibs.in
    real*8,dimension(InternalDimension)::dshift
    real*8,dimension(InternalDimension,InternalDimension)::Tshift
    !Work space
    character*2::chartemp
    integer::i,j,k
    integer,dimension(2)::indice
    character(3) :: sym
    real*8 :: anums
    real*8::dbletemp1,dbletemp2

    real*8, allocatable :: mass(:)
    real*8, parameter :: amu2au=1822.888486d0

    allocate(mass(natoms))

    !Precursor geometries, Columbus format
    open(unit=99,file='precursor.xyz',status='old')
    do i=1,NAtoms
      read(99,*) sym,anums,rPrecursor(3*i-2:3*i),mass(i)
    end do
    close(99)
    call WilsonBMatrixAndInternalCoordinate(rPrecursor,BPrecursor,qPrecursor,CartesianDimension,InternalDimension)
    !read in Precursor Hessian in a.u.
    open(unit=99,file='precursor.hess',status='old')
    do i=1,InternalDimension
      read(99,*) HPrecursor(i,:)
    end do
    close(99)
    mass=mass*amu2au
    call WilsonGFMethod(HPrecursor,BPrecursor,mass,freqPrecursor,intmodePrecursor,LinvPrecursor,cartmodePrecursor,InternalDimension,NAtoms)
    if(minval(freqPrecursor)<0d0) stop 'Program abort: imaginary frequency found for precursor'
    !Residual vibration
    open(unit=99,file='residual.xyz',status='old')
    do i=1,NAtoms
      read(99,*) sym,anums,rResidual(3*i-2:3*i),mass(i)
    end do
    close(99)
    call WilsonBMatrixAndInternalCoordinate(rResidual,BResidual,qResidual,CartesianDimension,InternalDimension)
    !get hessian
    open(unit=99,file='residual.hess',status='old')
    do i=1,InternalDimension
      read(99,*) HResidual(i,:)
    end do
    close(99)
    mass=mass*amu2au
    call WilsonGFMethod(HResidual,BResidual,mass,freqResidual,intmodeResidual,LinvResidual,cartmodeResidual,InternalDimension,NAtoms)
    if(minval(freqResidual)<0d0) stop 'Program abort: imaginary frequency found for residual'

    !Prepare nadvibs.in
    !Shift Hd origin to normal coordinate origin
    call OriginShift(qResidual(1:InternalDimension)-refint(1:InternalDimension))
    call HdEC_Hd2NVS(intmodeResidual)!Reformat Hd expansion coefficient into NadVibS format
    do i=1,InternalDimension!Subtract the harmonic oscillator potential term
        indice=i; j=NVS_WhichExpansionBasis(2,indice)
        forall(k=1:NStates)
            NVS_HdEC(k,k).Order(2).Array(j)=NVS_HdEC(k,k).Order(2).Array(j)-0.5d0*freqResidual(i)*freqResidual(i)
        end forall
    end do
    !Definition of dshift and Tshift see Schuurman & Yarkony 2008 JCP 128 eq. (12)
    dshift=matmul(LinvPrecursor,qResidual-qPrecursor)
    Tshift=matmul(LinvPrecursor,intmodeResidual)
    open(unit=99,file='nadvibs.in',status='replace')
        write(99,'(A54)')'Angular frequency of each vibrational basis: (In a.u.)'
        write(99,*)freqResidual
        do i=1,NStates
            do j=i,NStates
                do k=0,NVS_NOrder
                    write(99,'(A2,I2,I2,A15,I2)')'Hd',j,i,'Expansion order',k
                    write(99,*)NVS_HdEC(j,i).Order(k).Array
                end do
            end do
        end do
        write(99,'(A58)')'Angular frequency of each precursor normal mode: (In a.u.)'
        write(99,*)freqPrecursor
        write(99,'(A13)')'Shift Vector:'
        write(99,*)dshift
        write(99,'(A22)')'Transformation Matrix:'
        write(99,*)Tshift
    close(99)
end subroutine GenerateNadVibSInput
!------------------------------------------------------------------------------------------
!A photoelectron spectrum usually has the brightest line near the vertical excitation,
!    which emphasizes the overlap between the basis and precursor vibrational state
!Our goal is to find the minimum basis producing a satisfying overlap (empirically > 90%),
!    and this is an integer nonlinear optimization subject to inequality constraint
!A rigorous evaluation solution will be too expensive, so we make 2 approximations:
!1. In most cases the precursor is at vibrational ground state, which is a gaussian
!   with covariance matrix diagonal in precursor normal coordinate.
!   We approximate the overlap constraint by 'the basis covers some sigma eclipse',
!   by analogy to '2 sigma covers 95%' in single variate case
!2. Define the coverage of a 1 dimensional basis function as its standard deviation, then the
!   total coverage is a cuboid with each edge equals to the widest basis along this direction
!Certainly, the smallest cuboid must be tangential to the eclipse, so we only have to solve
!    the lower and upper bound of the eclipse along each residual normal coordinate
!which can be solved by a common real nonlinear optimization subject to equality constraint:
!    Minimize and maximize the component along each residual normal coordinate,
!    subject to staying on the 2 sigma eclipse
!------------------------------------------------------------------------------------------
subroutine BasisEstimation(qPrecursor,freqPrecursor,LinvPrecursor,qResidual,freqResidual,&
                           LinvResidual,intmodeResidual,intdim)
    use NonlinearOptimization
    implicit none
    integer,intent(in)::intdim
    real*8,dimension(intdim),intent(in)::qPrecursor,qResidual,freqPrecursor,freqResidual
    real*8,dimension(intdim,intdim),intent(in)::LinvPrecursor,LinvResidual,intmodeResidual
    real*8 :: contoursq, sgn
    real*8,dimension(intdim) :: DPrecursor, basis, LowerBound, UpperBound, bound, q
    integer :: i,j
    integer*8 :: NTotalBasis
    integer,dimension(intdim) :: NBasis

    write(*,*)'The basis will cover',NVS_contour,' times of sigma eclipse'
    !Since the ellipsoid function has squared times at right hand side
    contoursq=NVS_contour**2
    !Variance of mass weighted coordinate standard deviation of precursor ground state
    DPrecursor=0.5d0/freqPrecursor

    !Calculate each lower and upper bound, then decide how much basis is required
    do i=1,intdim
      !Lower bound
      q=0d0
      sgn=1d0
      call AugmentedLagrangian(f,fd,c,cd,q,intdim,1,Precision=1d-6/dble(intdim),f_fd=f_fd,fdd=fdd,cdd=cdd)
      call f(LowerBound(i),q,intdim)
      !Upper bound
      q=0d0
      sgn=-1d0
      call AugmentedLagrangian(f,fd,c,cd,q,intdim,1,Precision=1d-6/dble(intdim),f_fd=f_fd,fdd=fdd,cdd=cdd)
      call f(UpperBound(i),q,intdim)
      UpperBound(i)=-UpperBound(i)
      bound(i)=max(dAbs(LowerBound(i)),dAbs(UpperBound(i)))
      basis(i)=freqResidual(i)*bound(i)*bound(i)+0.5d0
    end do

    write(*,*)'Refer to NormalCoverage.txt and basis.txt for suggestion on NadVibS basis'
    open(unit=99,file='NormalCoverage.txt',status='replace')
    NTotalBasis=1
    do i=1,intdim
      NBasis(i)=ceiling(basis(i))
      NTotalBasis=NTotalBasis*NBasis(i)
    end do
    write(99,'(A23)',advance='no')'Total number of basis ='
    write(99,*) NTotalBasis
    write(99,'(A4,A1,A11,A1,A11,A1,A5,A1,A9)')  &
    'mode',char(9),'lower bound',char(9),'upper bound',char(9),'basis',char(9),'raw basis'
    do i=1,intdim
      write(99,'(I4,A1,F11.5,A1,F11.5,A1,I5,A1,F9.5)')  &
      i,char(9),LowerBound(i),char(9),UpperBound(i),char(9),NBasis(i),char(9),basis(i)
    end do
    close(99)

    open(unit=99,file='basis.txt',status='replace')
    do i=1,intdim
      write(99,'(2I5)') i,NBasis(i)
    end do
    close(99)

    write(*,*)'Refer to InternalCoverage.txt to validate the basis'
    open(unit=99,file='InternalCoverage.txt',status='replace')
      write(99,'(A55)') 'This is a report on the coverage in internal coordinate'
      write(99,'(A98)') 'Please check your internal coordinate definition to make sure angles are within well-defined range'
      write(99,'(A4,A1,A8)')'   q',char(9),'coverage'
      !Check the coverage in internal coordinate, since angles should not exceed some pi
      do i=1,intdim
        sgn=0d0
        do j=1,intdim
          sgn=sgn+dAbs(intmodeResidual(i,j)*bound(j))
        end do
        write(99,'(I4,A1,F8.4)') i,char(9),sgn
      end do
    close(99)

    return
    !The merit function and constraint
    !    i and sign controls the behavior of f routines:
    !    i-th residual normal coordinate will be searched
    !    sign > 0, search for lower bound; sign < 0, search for upper bound
    contains
    !-----------------------------------------------------------------------
    subroutine f(fq,q,intdim)
      integer,intent(in)::intdim
      real*8,dimension(intdim),intent(in)::q
      real*8,intent(out)::fq
      fq=sgn*dot_product(LinvResidual(i,:),q-qResidual)
    end subroutine f
    !-----------------------------------------------------------------------
    subroutine fd(fdq,q,intdim)
      integer,intent(in)::intdim
      real*8,dimension(intdim),intent(in)::q
      real*8,dimension(intdim),intent(out)::fdq
      fdq=sgn*LinvResidual(i,:)
    end subroutine fd
    !-----------------------------------------------------------------------
    integer function f_fd(fq,fdq,q,intdim)
      integer,intent(in)::intdim
      real*8,dimension(intdim),intent(in)::q
      real*8,intent(out)::fq
      real*8,dimension(intdim),intent(out)::fdq
      fq=sgn*dot_product(LinvResidual(i,:),q-qResidual)
      fdq=sgn*LinvResidual(i,:)
      f_fd=0!Return 0
    end function f_fd
    !-----------------------------------------------------------------------
    integer function fdd(fddq,q,intdim)
      integer,intent(in)::intdim
      real*8,dimension(intdim),intent(in)::q
      real*8,dimension(intdim,intdim),intent(out)::fddq
      fddq=0d0
      fdd=0!Return 0
    end function fdd
    !-----------------------------------------------------------------------
    subroutine c(cq,q,M,intdim)
      integer,intent(in)::M,intdim
      real*8,dimension(intdim),intent(in)::q
      real*8,dimension(M),intent(out)::cq
      real*8,dimension(intdim)::temp
      temp=matmul(LinvPrecursor,q-qPrecursor)
      cq(1)=dot_product(temp/DPrecursor,temp)-contoursq
    end subroutine c
    !-----------------------------------------------------------------------
    subroutine cd(cdq,q,M,intdim)
      integer,intent(in)::M,intdim
      real*8,dimension(intdim),intent(in)::q
      real*8,dimension(intdim,M),intent(out)::cdq
      integer::i
      real*8,dimension(intdim,intdim)::temp
      temp=transpose(LinvPrecursor)
      forall(i=1:intdim)
        temp(:,i)=temp(:,i)/DPrecursor(i)
      end forall
      cdq(:,1)=2d0*matmul(matmul(temp,LinvPrecursor),q-qPrecursor)
    end subroutine cd
    !-----------------------------------------------------------------------
    integer function cdd(cddq,q,M,intdim)
      integer,intent(in)::M,intdim
      real*8,dimension(intdim),intent(in)::q
      real*8,dimension(intdim,intdim,M),intent(out)::cddq
      integer::i
      real*8,dimension(intdim,intdim)::temp
      temp=transpose(LinvPrecursor)
      forall(i=1:intdim)
        temp(:,i)=temp(:,i)/DPrecursor(i)
      end forall
      cddq(:,:,1)=2d0*matmul(temp,LinvPrecursor)
      cdd=0!Return 0
    end function cdd
    !-----------------------------------------------------------------------
end subroutine BasisEstimation
!------------------------------------------------------------------------------------------
subroutine InitializeNadVibSInterface()
    integer::i,j,iorder
    NVS_NOrder=Hd_EBNR(1).order!Determine the highest order used
    allocate(NVS_NumberOfEachOrderTerms(0:NVS_NOrder))
    allocate(NVS_EBNR(0:NVS_NOrder))
    do iorder=0,NVS_NOrder
        NVS_NumberOfEachOrderTerms(iorder)=binom(InternalDimension+iorder-1,iorder)
        allocate(NVS_EBNR(iorder).Number(NVS_NumberOfEachOrderTerms(iorder)))
        !Generate expansion basis numbering mapping for iorder-th order (NVS_EBNR(iorder))
        !My preference is to use pseudo InternalDimension+1 counter satisfying former digit >= latter digit,
        !corresponding to the direct sum of an NVS_NOrder-th order tensor's 1st dimension vector
        allocate(NVS_EBNR(iorder).Number(1).Array(iorder))
        NVS_EBNR(iorder).Number(1).Array=1
        do j=2,NVS_NumberOfEachOrderTerms(iorder)
            allocate(NVS_EBNR(iorder).Number(j).Array(iorder))
            NVS_EBNR(iorder).Number(j).Array=NVS_EBNR(iorder).Number(j-1).Array
            NVS_EBNR(iorder).Number(j).Array(1)=NVS_EBNR(iorder).Number(j).Array(1)+1!Add 1 to the 1st digit
            do i=1,iorder-1!Carry to latter digits
                if(NVS_EBNR(iorder).Number(j).Array(i)>InternalDimension) then
                    NVS_EBNR(iorder).Number(j).Array(i)=1
                    NVS_EBNR(iorder).Number(j).Array(i+1)=NVS_EBNR(iorder).Number(j).Array(i+1)+1
                end if
            end do
            do i=iorder-1,1,-1!Modify to satisfy former digit >= latter digit
                if(NVS_EBNR(iorder).Number(j).Array(i)<NVS_EBNR(iorder).Number(j).Array(i+1)) &
                    NVS_EBNR(iorder).Number(j).Array(i)=NVS_EBNR(iorder).Number(j).Array(i+1)
            end do
        end do
    end do
    allocate(NVS_HdEC(NStates,NStates))!Allocate storage space of NVS_HdEC
    do j=1,NStates
        do i=j,NStates
            allocate(NVS_HdEC(i,j).Order(0:NVS_NOrder))
            do iorder=0,NVS_NOrder
                allocate(NVS_HdEC(i,j).Order(iorder).Array(NVS_NumberOfEachOrderTerms(iorder)))
                NVS_HdEC(i,j).Order(iorder).Array=0d0
            end do
        end do
    end do
end subroutine InitializeNadVibSInterface

subroutine HdEC_Hd2NVS(L)!Transform from internal coordinate Hd_HdEC to normal mode NVS_HdEC
    !This is done by:
    !    1, select an Hd_HdEC term
    !    2, under rotation, the terms making up the multiplication become the linear combination of every mode,
    !       i.e. q = L . Q, so we go through all combination and add the contribution to NVS_HdEC
    !    3, go to 1 until all Hd_HdEC are done
    real*8,dimension(InternalDimension,InternalDimension),intent(in)::L
    integer::order,location,n,i,j
    integer,dimension(NVS_NOrder)::indice,Hdindice
    real*8::coeff
    do n=1,NHdExpansionBasis!Main loop
        order=Hd_EBNR(n).order
        if(order==0) then!Const term will not change under any transformation
            forall(i=1:NStates,j=1:NStates,i>=j)
                NVS_HdEC(i,j).Order(0).Array(1)=NVS_HdEC(i,j).Order(0).Array(1)+Hd_HdEC(i,j).Array(n)
            end forall
        else!Go through all combination in a InternalDimension+1 counter manner
            Hdindice(1:order)=Hd_EBNR(n).indice
            indice(1:order)=1!indice(i)=j means using j-th mode at i-th position in multiplication
            do while(indice(order)<=InternalDimension)!Done when the counter overflows
                coeff=1d0
                do i=1,order
                    coeff=coeff*L(Hdindice(i),indice(i))
                end do
                if(coeff/=0d0) then
                    location=NVS_WhichExpansionBasis(order,indice(1:order))
                    forall(i=1:NStates,j=1:NStates,i>=j)
                        NVS_HdEC(i,j).Order(order).Array(location)=NVS_HdEC(i,j).Order(order).Array(location)&
                            +coeff*Hd_HdEC(i,j).Array(n)
                    end forall
                end if
                indice(1)=indice(1)+1!Add 1 to the InternalDimension+1 counter
                do i=1,order-1
                    if(indice(i)>InternalDimension) then!Carry
                        indice(i)=1
                        indice(i+1)=indice(i+1)+1
                    end if
                end do
            end do
        end if
    end do
end subroutine HdEC_Hd2NVS

integer function NVS_WhichExpansionBasis(order,indiceinput)
    integer,intent(in)::order
    integer,dimension(order),intent(in)::indiceinput
    integer,dimension(order)::indice,temp
    indice=indiceinput
    call bsort(order,indice,-1)
    call bisect(1,NVS_NumberOfEachOrderTerms(order))
    contains
    recursive subroutine bisect(low,up)
        integer,intent(in)::low,up
        integer::i,bisection
        if(low==up) then
            NVS_WhichExpansionBasis=low
        else if(up-low==1) then
            do i=1,order
                if(indice(i)/=NVS_EBNR(order).Number(low).Array(i)) exit
            end do
            if(i>order) then
                NVS_WhichExpansionBasis=low
            else
                NVS_WhichExpansionBasis=up
            end if
        else
            bisection=(low+up)/2
            !Binary search according to my pseudo InternalDimension+1 counter expansion basis numbering rule
            do i=order,1,-1
                if(indice(i)/=NVS_EBNR(order).Number(bisection).Array(i)) exit
            end do
            if(i<1) then
                NVS_WhichExpansionBasis=bisection
            else
                if(indice(i)>NVS_EBNR(order).Number(bisection).Array(i)) then
                    call bisect(bisection,up)
                else
                    call bisect(low,bisection)
                end if
            end if
        end if
    end subroutine bisect
end function NVS_WhichExpansionBasis

!------------------------------------------------------------------------------------
function binom ( n, k )
!  Reference:
!    Bill Buckles, Matthew Lybanon,
!    Algorithm 515: Generation of a Vector from the Lexicographical Index,
!    ACM Transactions on Mathematical Software,
!    Volume 3, Number 2, June 1977, pages 180-182.
!  Parameters:
!    Input, integer N, K, the parameters for the binomial coefficient.
!    Output, integer BINOM, the binomial coefficient.
  implicit none
  integer binom
  integer i
  integer k
  integer k1
  integer n
  integer n1
  integer p
  integer r

  k1 = k
  p = n - k1

  if ( k1 .lt. p ) then
    p = k1
    k1 = n - p
  end if

  if ( p .eq. 0 ) then
    r = 1
  else
    r = k1 + 1
  end if

  do i = 2, p
    r = ( r * ( k1 + i ) ) / i
  end do

  binom = r

  return
end
!-----------------------------------------------------------------------------------
end module NadVibSInterface
