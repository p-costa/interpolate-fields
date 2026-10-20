module mod_grid
  !
  ! one-dimensional grids, as in CaNS: xf(i) is the upper face of cell i (xf(0) = 0, xf(n) = l) and xc(i)
  ! its centre; the ghost points i = 0 and i = n+1 follow from periodicity, or from the mirror image about
  ! the boundary otherwise; dxf(i) = xf(i)-xf(i-1) are the cell widths and dxc(i) = xc(i+1)-xc(i) the
  ! spacings between cell centres
  !
  use mpi_f08
  use mod_common, only: rp,sp,dp,i8,myid,ierr,abort_run
  implicit none
  private
  public grid_1d,init_grid,is_same_length
  type grid_1d
    integer  :: n
    real(rp) :: l
    logical  :: is_periodic
    integer  :: iprec ! precision of the grid file, in bytes (0 for a uniform grid)
    real(rp), allocatable, dimension(:) :: xc,xf,dxc,dxf
  end type grid_1d
contains
  subroutine init_grid(n,l,fname,is_periodic,g)
    !
    ! sets a uniform grid with n cells over a length l if `fname` is blank; otherwise, reads the grid from
    ! a CaNS grid file (`grid.bin`, or `grid_x.bin`, `grid_y.bin` and `grid_z.bin` in CaNS-EIGEN) with
    ! [dzc(1:n),dzf(1:n),zc(1:n),zf(1:n)] in single or double precision (detected from the file size);
    ! the length is then that of the file, and l, if positive, is checked against it
    !
    implicit none
    integer         , intent(in ) :: n
    real(rp)        , intent(in ) :: l
    character(len=*), intent(in ) :: fname
    logical         , intent(in ) :: is_periodic
    type(grid_1d)   , intent(out) :: g
    real(sp), allocatable, dimension(:,:) :: buf_sp
    real(dp), allocatable, dimension(:,:) :: buf_dp
    integer(i8) :: fsize
    integer :: i,iunit
    character(len=16) :: cn
    !
    if(n < 1) call abort_run('the number of grid cells `ng` must be prescribed and positive.')
    g%n = n
    g%is_periodic = is_periodic
    allocate(g%xc(0:n+1),g%xf(0:n+1),g%dxc(0:n+1),g%dxf(0:n+1))
    if(len_trim(fname) == 0) then
      if(l <= 0._rp) call abort_run('the domain length `l` must be prescribed and positive for a uniform grid.')
      g%iprec = 0
      g%l = l
      do i=0,n
        g%xf(i) = (i*l)/n
      end do
      do i=1,n
        g%xc(i) = ((i-0.5_rp)*l)/n
      end do
    else
      fsize = -1
      if(myid == 0) inquire(file=trim(fname),size=fsize)
      call MPI_BCAST(fsize,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
      write(cn,'(i0)') n
      if(fsize < 0) call abort_run('could not open grid file `'//trim(fname)//'`.')
      if(     fsize == 4*int(n,i8)*storage_size(1._sp)/8) then
        g%iprec = 4
      else if(fsize == 4*int(n,i8)*storage_size(1._dp)/8) then
        g%iprec = 8
      else
        call abort_run('the size of grid file `'//trim(fname)//'` is not that of a grid with '//trim(cn)// &
                       ' cells (four arrays of n values, in single or double precision).')
      end if
      allocate(buf_dp(n,4))
      if(myid == 0) then
        open(newunit=iunit,file=trim(fname),access='stream',form='unformatted',status='old',action='read')
        if(g%iprec == 4) then
          allocate(buf_sp(n,4))
          read(iunit) buf_sp(:,:)
          buf_dp(:,:) = buf_sp(:,:)
        else
          read(iunit) buf_dp(:,:)
        end if
        close(iunit)
      end if
      call MPI_BCAST(buf_dp,4*n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      g%xc(1:n) = buf_dp(:,3)
      g%xf(1:n) = buf_dp(:,4)
      g%xf(0)   = 0._rp
      g%l       = g%xf(n)
      if(any(g%xf(1:n) <= g%xf(0:n-1)) .or. any(g%xc(1:n) <= g%xf(0:n-1)) .or. any(g%xc(1:n) >= g%xf(1:n))) &
        call abort_run('the coordinates in grid file `'//trim(fname)//'` are not increasing, '// &
                       'or not those of a grid with '//trim(cn)//' cells.')
      if(l > 0._rp .and. .not.is_same_length(l,g%l,g%iprec == 4)) &
        call abort_run('the domain length `l` is inconsistent with grid file `'//trim(fname)//'`.')
    end if
    !
    ! ghost points
    !
    if(is_periodic) then
      g%xc(0  ) = g%xc(n) - g%l
      g%xc(n+1) = g%xc(1) + g%l
      g%xf(n+1) = g%xf(1) + g%l
    else
      g%xc(0  ) = -g%xc(1)
      g%xc(n+1) = 2*g%l - g%xc(n)
      g%xf(n+1) = 2*g%l - g%xf(max(n-1,0))
    end if
    g%dxf(1:n+1) = g%xf(1:n+1) - g%xf(0:n)
    g%dxf(0)     = g%dxf(1)
    if(is_periodic) g%dxf(0) = g%dxf(n)
    g%dxc(0:n)   = g%xc(1:n+1) - g%xc(0:n)
    g%dxc(n+1)   = g%dxc(n)
    if(is_periodic) g%dxc(n+1) = g%dxc(1)
  end subroutine init_grid
  !
  logical function is_same_length(l1,l2,is_single)
    !
    ! whether two domain lengths agree within round-off; a single-precision grid file holds
    ! the sums of single-precision cell widths
    !
    implicit none
    real(rp), intent(in) :: l1,l2
    logical , intent(in) :: is_single
    real(rp) :: tol
    tol = 1.e-10_rp
    if(is_single) tol = 1.e-4_rp
    is_same_length = abs(l1-l2) <= tol*max(l1,l2)
  end function is_same_length
end module mod_grid
