module mod_common
  use mpi_f08        , only: MPI_REAL_RP => MPI_DOUBLE_PRECISION
  implicit none
  integer, parameter :: rp = selected_real_kind(15,307)
  integer :: ierr
  integer :: nh = 1
end module mod_common
