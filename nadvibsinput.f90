!=============================================================================
program main
  use hddata
  use geomtrans
  use NadVibSInterface
  implicit none

  character(255) :: filename, npoints, enfDiab, epmax, w_energy, w_grad, w_fij, &
  gradcutoff, cpcutoff, deggrdbinding, deg_cap, lambda, &
  nrmediff, ediffcutoff, fixref ! dummy except for filename
  DOUBLE PRECISION,dimension(10) :: energyT, highEScale ! dummy
  DOUBLE PRECISION :: eshift

  namelist /fitting/ npoints, enfDiab, epmax, w_energy, w_grad, w_fij, &
  gradcutoff, cpcutoff, deggrdbinding, deg_cap, lambda, eshift, energyT, &
  highEScale, nrmediff, ediffcutoff, fixref, natoms, nstates
  open(103,file='fit.in',delim='APOSTROPHE')
  read(unit=103,nml=fitting)
  close(103)

  call initpes
  InternalDimension=ncoord
  CartesianDimension=3*natoms
  call InitializeNadVibSInterface
  call GenerateNadVibSInput

  stop
end
!=============================================================================
