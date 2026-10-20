module mod_decomp
  !
  ! decomposition of a grid of ng(1) x ng(2) x ng(3) points into the blocks of a dims(1) x dims(2) x dims(3)
  ! process grid (a pencil decomposition has a 1 in dims), as in 2DECOMP&FFT: along each direction, the
  ! first mod(ng,dims) tasks get one more point
  !
  use mod_common, only: abort_run
  implicit none
  private
  public distribute_grid,task_coords
contains
  subroutine distribute_grid(ng,dims,coords,n,lo,hi)
    !
    ! block lo:hi (n points per direction) of the task with coordinates `coords` (from 0) in the process grid
    !
    implicit none
    integer, intent(in ), dimension(3) :: ng,dims,coords
    integer, intent(out), dimension(3) :: n,lo,hi
    n(:)  = ng(:)/dims(:)
    lo(:) = 1 + coords(:)*n(:) + min(coords(:),mod(ng(:),dims(:)))
    where(coords(:) < mod(ng(:),dims(:))) n(:) = n(:) + 1
    hi(:) = lo(:) + n(:) - 1
  end subroutine distribute_grid
  !
  subroutine task_coords(irank,dims,order,coords)
    !
    ! coordinates (from 0) in the process grid of task irank, with ranks in 'column-major' order (the first
    ! coordinate varies fastest, irank = c1 + dims(1)*(c2 + dims(2)*c3)) or 'row-major' order (the last
    ! coordinate varies fastest, as for an MPI Cartesian topology, irank = c3 + dims(3)*(c2 + dims(2)*c1))
    !
    implicit none
    integer         , intent(in ) :: irank
    integer         , intent(in ), dimension(3) :: dims
    character(len=*), intent(in ) :: order
    integer         , intent(out), dimension(3) :: coords
    select case(order)
    case('column-major')
      coords(:) = [mod(irank,dims(1)),mod(irank/dims(1),dims(2)),irank/(dims(1)*dims(2))]
    case('row-major')
      coords(:) = [irank/(dims(2)*dims(3)),mod(irank/dims(3),dims(2)),mod(irank,dims(3))]
    case default
      call abort_run('unknown task order `'//order//'` (''column-major'' or ''row-major'').')
    end select
  end subroutine task_coords
end module mod_decomp
