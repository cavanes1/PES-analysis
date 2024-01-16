!=============================================================================
program main
  use hddata
  use geomtrans
  use NadVibSInterface
  implicit none

  call initpes
  InternalDimension=ncoord
  CartesianDimension=3*natoms
  call InitializeNadVibSInterface
  call GenerateNadVibSInput

  stop
end
!=============================================================================
