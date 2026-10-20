module mod_common
  use mpi_f08, only: MPI_FINALIZE
  implicit none
  !
  ! the fields are read and written in single or double precision (`sp`, `dp`), and all
  ! computations are performed in double precision (`rp`)
  !
  integer, parameter :: sp = selected_real_kind(6 , 37), &
                        dp = selected_real_kind(15,307), &
                        rp = dp
  integer, parameter :: i8 = selected_int_kind(18)
  integer(i8), parameter :: CHUNK_BYTES = 64*2_i8**20 ! size of the conversion buffer for single-precision output
  integer :: ierr
  integer :: myid = 0
contains
  subroutine abort_run(msg)
    !
    ! stops the run after an input error; must be called by all tasks
    !
    implicit none
    character(len=*), intent(in) :: msg
    if(myid == 0) print*, 'ERROR: ', msg
    if(myid == 0) print*, 'Aborting...'
    call MPI_FINALIZE(ierr)
    error stop
  end subroutine abort_run
end module mod_common
