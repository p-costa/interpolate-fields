  program interpolate_fields
    !
    ! interpolates CaNS and CaNS-EIGEN fields from one grid to another, by trilinear interpolation (see the README)
    !
    ! each task owns a block of the output grid and, for each variable, reads the block of the input field
    ! that brackets it, imposes the boundary conditions, interpolates and writes; the variables are processed
    ! one at a time, so that the memory footprint is that of one input and one output block per task
    !
    use mpi_f08
    use mod_common, only: rp,dp,i8,ierr,myid,abort_run
    use mod_param , only: read_input,vars,dims_in => dims,input_format,input_prefix,input_file,output_prefix, &
                          output_precision,is_reset_time,ng_i,l_i,grid_files_i,ng_o,l_o,grid_files_o,is_periodic
    use mod_grid  , only: grid_1d,init_grid,is_same_length
    use mod_bound , only: set_ghosts
    use mod_interp, only: get_weights,interp_fld
    use mod_io    , only: get_checkpoint_layout,read_fldinfo,read_block,write_field,CHUNK_BYTES
    implicit none
    !
    ! brackets and weights of the output points of a block along one direction (see `get_weights`)
    !
    type brackets
      integer , allocatable, dimension(:) :: m
      real(rp), allocatable, dimension(:) :: w
    end type brackets
    character(len=*), parameter :: cdir = 'xyz'
    character(len=256) :: fname_nml,fname_i,fname_o
    type(grid_1d)  , dimension(3) :: gi,go
    type(brackets) , dimension(3) :: br
    integer :: nproc
    integer , dimension(3) :: dims,coords,n_o,lo_o,hi_o,lo_i,hi_i,own_lo,own_hi
    type(MPI_COMM) :: comm_cart
    real(rp), allocatable, dimension(:,:,:) :: fldi,fldo
    real(rp) :: time,time_1,time_o,tol
    integer  :: istep,istep_1,istep_o
    integer(MPI_OFFSET_KIND) :: disp,disp_info,ntot_i
    integer :: iv,idir,iprec_i,iprec_o,nfields,mmax
    character(len=3) :: cloc
    real(dp) :: twi,tw(3)
    integer(i8) :: mem,mem_max
    !
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    !
    ! read the input parameters from the namelist file given as argument (default: `interp.nml`)
    !
    fname_nml = 'interp.nml'
    if(command_argument_count() > 0) call get_command_argument(1,fname_nml)
    call read_input(trim(fname_nml))
    !
    ! input and output grids
    !
    do idir=1,3
      call init_grid(ng_i(idir),l_i(idir),trim(grid_files_i(idir)),is_periodic(idir),gi(idir))
      call init_grid(ng_o(idir),l_o(idir),trim(grid_files_o(idir)),is_periodic(idir),go(idir))
      if(.not.is_same_length(gi(idir)%l,go(idir)%l,gi(idir)%iprec == 4 .or. go(idir)%iprec == 4)) &
        call abort_run('the input and output domain lengths along '//cdir(idir:idir)//' differ.')
    end do
    if(myid == 0) then
      print*, '*** interpolate-fields ***'
      do idir=1,3
        call print_grid('Input  grid, '//cdir(idir:idir)//':',gi(idir),grid_files_i(idir))
      end do
      do idir=1,3
        call print_grid('Output grid, '//cdir(idir:idir)//':',go(idir),grid_files_o(idir))
      end do
    end if
    ntot_i = product(int(ng_i(:),MPI_OFFSET_KIND))
    !
    ! process grid; each task owns a block of the output grid
    !
    dims(:) = dims_in(:)
    if(mod(nproc,product(dims(:),mask=dims(:) > 0)) /= 0 .or. &
       (all(dims(:) > 0) .and. product(dims(:)) /= nproc)) &
      call abort_run('`dims` is inconsistent with the number of MPI tasks.')
    call MPI_DIMS_CREATE(nproc,3,dims,ierr)
    if(any(dims(:) > ng_o(:))) call abort_run('`dims` exceeds the number of output grid cells.')
    call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,is_periodic,.false.,comm_cart,ierr)
    call MPI_CART_COORDS(comm_cart,myid,3,coords,ierr)
    call distribute_grid(ng_o,dims,coords,n_o,lo_o,hi_o)
    if(myid == 0) print*, 'Process grid: ', dims(:), ' (', nproc, ' tasks)'
    !
    tw(:) = 0._dp
    mem_max = 0
    do iv=1,size(vars)
      associate(var => vars(iv))
      !
      ! input file, precision, time and time step
      !
      select case(trim(input_format))
      case('legacy')
        fname_i = input_file
        call get_checkpoint_layout(trim(fname_i),ng_i,.true.,iprec_i,nfields)
        if(var%ifield > nfields) &
          call abort_run('variable `'//trim(var%name)//'` is not in the legacy checkpoint `'//trim(fname_i)//'`.')
        disp      = (var%ifield-1)*ntot_i*iprec_i
        disp_info = nfields*ntot_i*iprec_i
      case('per-variable')
        fname_i = trim(input_prefix)//'_'//trim(var%name)//'.bin'
        call get_checkpoint_layout(trim(fname_i),ng_i,.false.,iprec_i,nfields)
        disp      = 0
        disp_info = ntot_i*iprec_i
      end select
      call read_fldinfo(trim(fname_i),disp_info,iprec_i,time,istep)
      if(iv == 1) then
        time_1  = time
        istep_1 = istep
      else if(time /= time_1 .or. istep /= istep_1) then
        if(myid == 0) print*, 'WARNING: the time or time step of `', trim(var%name), &
                              '` differs from that of the first variable: ', time, istep
      end if
      select case(trim(output_precision))
      case('single')
        iprec_o = 4
      case('double')
        iprec_o = 8
      case default
        iprec_o = iprec_i
      end select
      do idir=1,3
        cloc(idir:idir) = merge('s','c',var%is_staggered(idir))
      end do
      if(myid == 0) print*, 'Variable ', trim(var%name), ' (', cloc, ') from `', trim(fname_i), &
                            '` (', merge('single','double',iprec_i == 4), ' precision), at time = ', time, &
                            'step = ', istep
      !
      ! brackets of the output points along each direction, from the input coordinates of the variable's
      ! location (faces if staggered, centres otherwise); the input points are distributed among the tasks
      ! (own_lo:own_hi, for the statistics) by the brackets of the first output point of each block, and each
      ! task reads the input block lo_i:hi_i that holds the brackets of its output block and its own points
      !
      do idir=1,3
        allocate(br(idir)%m(lo_o(idir):hi_o(idir)+1),br(idir)%w(lo_o(idir):hi_o(idir)+1))
        tol = 10*epsilon(1._rp)*max(gi(idir)%l,go(idir)%l)
        if(var%is_staggered(idir)) then
          mmax = gi(idir)%n-1
          call get_weights(gi(idir)%xf,mmax,tol,go(idir)%xf(lo_o(idir):hi_o(idir)+1),br(idir)%m,br(idir)%w)
        else
          mmax = gi(idir)%n
          call get_weights(gi(idir)%xc,mmax,tol,go(idir)%xc(lo_o(idir):hi_o(idir)+1),br(idir)%m,br(idir)%w)
        end if
        own_lo(idir) = 1
        own_hi(idir) = gi(idir)%n
        if(coords(idir) > 0          ) own_lo(idir) = br(idir)%m(lo_o(idir))+1
        if(coords(idir) < dims(idir)-1) own_hi(idir) = br(idir)%m(hi_o(idir)+1)
        lo_i(idir) = min(br(idir)%m(lo_o(idir))  ,own_lo(idir))
        hi_i(idir) = max(br(idir)%m(hi_o(idir))+1,own_hi(idir))
      end do
      !
      ! read the input block and impose the boundary conditions (the block is zeroed first, as its edges
      ! and corners along bounded directions are only set by the last boundary condition imposed there)
      !
      allocate(fldi(lo_i(1):hi_i(1),lo_i(2):hi_i(2),lo_i(3):hi_i(3)))
      fldi(:,:,:) = 0._rp
      twi = MPI_WTIME()
      call read_block(trim(fname_i),disp,ng_i,iprec_i,is_periodic,lo_i,hi_i,fldi)
      tw(1) = tw(1) + (MPI_WTIME()-twi)
      call print_stats(trim(var%name)//' (input) ',own_lo,own_hi,var%is_staggered,gi,lo_i,fldi)
      call set_ghosts(var%cbc,var%bc,var%is_staggered,gi,lo_i,hi_i,fldi)
      !
      ! interpolate
      !
      allocate(fldo(lo_o(1):hi_o(1),lo_o(2):hi_o(2),lo_o(3):hi_o(3)))
      twi = MPI_WTIME()
      call interp_fld(lo_i,lo_o,hi_o,br(1)%m,br(2)%m,br(3)%m,br(1)%w,br(2)%w,br(3)%w,fldi,fldo)
      tw(2) = tw(2) + (MPI_WTIME()-twi)
      !
      ! peak memory of the field arrays: the input block, with the read buffer of `read_block` or with
      ! the output block and the conversion buffer of `write_field`
      !
      mem = 8*size(fldo,kind=i8)
      if(iprec_o == 4) mem = mem + min(4*size(fldo,kind=i8),CHUNK_BYTES)
      mem = 8*size(fldi,kind=i8) + max(iprec_i*size(fldi,kind=i8),mem)
      mem_max = max(mem_max,mem)
      deallocate(fldi)
      call print_stats(trim(var%name)//' (output)',lo_o,hi_o,var%is_staggered,go,lo_o,fldo)
      !
      ! write the output
      !
      time_o  = time_1
      istep_o = istep_1
      if(is_reset_time) then
        time_o  = 0._rp
        istep_o = 0
      end if
      fname_o = trim(output_prefix)//'_'//trim(var%name)//'.bin'
      twi = MPI_WTIME()
      call write_field(trim(fname_o),ng_o,lo_o,hi_o,iprec_o,fldo,time_o,istep_o)
      tw(3) = tw(3) + (MPI_WTIME()-twi)
      if(myid == 0) print*, '  written to `', trim(fname_o), '` (', merge('single','double',iprec_o == 4), &
                            ' precision), at time = ', time_o, 'step = ', istep_o
      deallocate(fldo)
      do idir=1,3
        deallocate(br(idir)%m,br(idir)%w)
      end do
      end associate
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,tw     ,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,mem_max,1,MPI_INTEGER8        ,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(myid == 0) then
      print*, 'Elapsed time for reading, interpolation and writing (s, maximum over tasks): ', tw(:)
      print*, 'Memory of the field arrays (MiB, maximum over tasks and variables): ', mem_max/2._dp**20
    end if
    call MPI_FINALIZE(ierr)
  contains
    subroutine distribute_grid(ng,dims,coords,n,lo,hi)
      !
      ! block of the task with coordinates `coords` in the process grid; along each direction,
      ! the first mod(ng,dims) tasks get one more point
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
    subroutine print_grid(label,g,fname)
      implicit none
      character(len=*), intent(in) :: label,fname
      type(grid_1d)   , intent(in) :: g
      if(g%iprec == 0) then
        print*, label, g%n, ' cells, length ', g%l, ', uniform'
      else
        print*, label, g%n, ' cells, length ', g%l, ', from `', trim(fname), '` (', &
                merge('single','double',g%iprec == 4), ' precision)'
      end if
    end subroutine print_grid
    !
    subroutine print_stats(label,lo,hi,is_staggered,g,lo_p,p)
      !
      ! prints the global minimum, maximum and volume-weighted mean of p over the points lo:hi of each task,
      ! with the cell widths dxf along cell-centred directions, and the spacings dxc between cell centres
      ! along staggered directions
      !
      implicit none
      character(len=*), intent(in) :: label
      integer         , intent(in), dimension(3) :: lo,hi
      logical         , intent(in), dimension(3) :: is_staggered
      type(grid_1d)   , intent(in), dimension(3) :: g
      integer         , intent(in), dimension(3) :: lo_p
      real(rp)        , intent(in), dimension(lo_p(1):,lo_p(2):,lo_p(3):) :: p
      real(dp), allocatable, dimension(:) :: wx,wy,wz
      real(dp) :: pmin,pmax,psum(2)
      integer :: j,k
      !
      allocate(wx(lo(1):hi(1)),wy(lo(2):hi(2)),wz(lo(3):hi(3)))
      wx(:) = merge(g(1)%dxc(lo(1):hi(1)),g(1)%dxf(lo(1):hi(1)),is_staggered(1))
      wy(:) = merge(g(2)%dxc(lo(2):hi(2)),g(2)%dxf(lo(2):hi(2)),is_staggered(2))
      wz(:) = merge(g(3)%dxc(lo(3):hi(3)),g(3)%dxf(lo(3):hi(3)),is_staggered(3))
      pmin =  huge(1._dp)
      pmax = -huge(1._dp)
      psum(:) = 0._dp
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          pmin = min(pmin,minval(p(lo(1):hi(1),j,k)))
          pmax = max(pmax,maxval(p(lo(1):hi(1),j,k)))
          psum(1) = psum(1) + wz(k)*wy(j)*sum(wx(:)*p(lo(1):hi(1),j,k))
          psum(2) = psum(2) + wz(k)*wy(j)*sum(wx(:))
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,pmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,pmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,psum,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) print*, '  ', label, ': min = ', pmin, '| max = ', pmax, '| mean = ', psum(1)/psum(2)
    end subroutine print_stats
  end program interpolate_fields
