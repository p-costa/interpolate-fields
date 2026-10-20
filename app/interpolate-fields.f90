  program interpolate_fields
    !
    ! interpolates CaNS and CaNS-EIGEN fields from one grid to another, by trilinear interpolation (see the README)
    !
    ! each task owns a block of the output grid and, for each variable, reads the block of the input field
    ! that brackets it, imposes the boundary conditions, interpolates and writes; the variables are processed
    ! one at a time, so that the memory footprint is that of one input and one output block per task
    !
    use mpi_f08
    use mod_common, only: rp,dp,i8,ierr,myid,abort_run,CHUNK_BYTES
    use mod_param , only: read_input,vars,dims_in => dims,input_format,output_format,output_precision, &
                          is_compress_output,is_reset_time,ng_i,l_i,grid_files_i,ng_o,l_o,grid_files_o, &
                          is_periodic,input_file_name,output_file_name,ntiles,ntiles_tot,tile_inputs, &
                          tile_shifts,MAX_TILES
    use mod_grid  , only: grid_1d,init_grid,tile_grid,is_same_length
    use mod_bound , only: set_ghosts
    use mod_interp, only: get_weights,interp_fld
    use mod_io    , only: get_checkpoint_layout,read_fldinfo,read_block,write_field
    implicit none
    !
    ! brackets and weights of the output points of a block along one direction (see `get_weights`)
    !
    type brackets
      integer , allocatable, dimension(:) :: m
      real(rp), allocatable, dimension(:) :: w
    end type brackets
    character(len=*), parameter :: cdir = 'xyz'
    character(len=256) :: fname_nml,fname_o
    character(len=256), dimension(MAX_TILES) :: fnames_i
    type(grid_1d)  , dimension(3) :: gi,go
    type(brackets) , dimension(3) :: br
    integer :: nproc
    integer , dimension(3) :: dims,coords,n_o,lo_o,hi_o,lo_i,hi_i,own_lo,own_hi
    type(MPI_COMM) :: comm_cart
    real(rp), allocatable, dimension(:,:,:) :: fldi,fldo
    real(rp) :: time,time_1,time_o,tol,pmin,pmax,pmean
    integer  :: istep,istep_1,istep_o
    integer(MPI_OFFSET_KIND) :: disp_info,ntot_i
    integer(MPI_OFFSET_KIND), dimension(MAX_TILES) :: disps
    real(rp), dimension(MAX_TILES) :: means
    integer :: iv,idir,it,iprec_i,iprec_t,iprec_o,nfields,mmax
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
    end do
    if(myid == 0) then
      print*, '*** interpolate-fields ***'
      do idir=1,3
        call print_grid('Input  grid, '//cdir(idir:idir)//':',gi(idir),grid_files_i(idir))
      end do
    end if
    !
    ! periodic tiling: the input domain is made of ntiles(1) x ntiles(2) x ntiles(3) copies of the input grid
    !
    if(ntiles_tot > 1) then
      do idir=1,3
        call tile_grid(ntiles(idir),gi(idir))
      end do
      if(myid == 0) then
        print*, 'Tiling: ', ntiles(:), ' tiles (x fastest), with their inputs and shifts (in input cells):'
        do it=1,ntiles_tot
          print*, '  tile ', it, ': `', trim(tile_inputs(it)), '`, shift ', tile_shifts(:,it)
        end do
        do idir=1,3
          call print_grid('Tiled  grid, '//cdir(idir:idir)//':',gi(idir),grid_files_i(idir))
        end do
      end if
    end if
    do idir=1,3
      if(.not.is_same_length(gi(idir)%l,go(idir)%l,gi(idir)%iprec == 4 .or. go(idir)%iprec == 4)) &
        call abort_run('the input and output domain lengths along '//cdir(idir:idir)//' differ.')
    end do
    if(myid == 0) then
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
      ! input files (one per tile, all in the same precision), precision, and time and time step (of the
      ! first tile)
      !
      do it=ntiles_tot,1,-1 ! (backwards, so that iprec_i and disp_info end as those of the first tile)
        fnames_i(it) = input_file_name(var%name,it)
        call get_checkpoint_layout(trim(fnames_i(it)),trim(input_format),trim(var%name),ng_i,iprec_t,nfields)
        if(it < ntiles_tot .and. iprec_t /= iprec_i) &
          call abort_run('the input files of the tiles must be in the same precision.')
        iprec_i = iprec_t
        if(trim(input_format) == 'legacy') then
          if(var%ifield > nfields) &
            call abort_run('variable `'//trim(var%name)//'` is not in the legacy checkpoint `'// &
                           trim(fnames_i(it))//'`.')
          disps(it) = (var%ifield-1)*ntot_i*iprec_i
          disp_info = nfields*ntot_i*iprec_i
        else
          disps(it) = 0
          disp_info = ntot_i*iprec_i
        end if
      end do
      call read_fldinfo(trim(fnames_i(1)),trim(input_format),disp_info,iprec_i,time,istep)
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
      if(myid == 0) print*, 'Variable ', trim(var%name), ' (', cloc, ') from `', trim(fnames_i(1)), &
                            '` (', merge('single','double',iprec_i == 4), ' precision), at time = ', time, &
                            'step = ', istep
      if(myid == 0 .and. ntiles_tot > 1) print*, '  (time and step of the first of the ', ntiles_tot, ' tiles)'
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
      call read_block(fnames_i(1:ntiles_tot),trim(input_format),trim(var%name),disps(1:ntiles_tot),ng_i,ntiles, &
                      tile_shifts(:,1:ntiles_tot),iprec_i,is_periodic,lo_i,hi_i,fldi)
      tw(1) = tw(1) + (MPI_WTIME()-twi)
      call field_stats(own_lo,own_hi,var%is_staggered,gi,lo_i,fldi,pmin,pmax,pmean)
      call print_stats(trim(var%name)//' (input) ',pmin,pmax,pmean)
      call set_ghosts(var%cbc,var%bc,var%is_staggered,gi,lo_i,hi_i,fldi)
      !
      ! field operations on the input: subtract the volume-weighted mean (of each tile) or a constant offset,
      ! in double precision; the ghost points are shifted too, i.e. as the field with Dirichlet values
      ! shifted alike
      !
      if(var%is_remove_mean) then
        call tile_means(own_lo,own_hi,var%is_staggered,gi,lo_i,fldi,means(1:ntiles_tot))
        call subtract_tile_values(lo_i,hi_i,means(1:ntiles_tot),fldi)
        if(myid == 0) then
          if(ntiles_tot == 1) then
            print*, '  ', trim(var%name), ': mean of the input subtracted: ', means(1)
          else
            do it=1,ntiles_tot
              print*, '  ', trim(var%name), ': mean of tile ', it, ' subtracted: ', means(it)
            end do
          end if
        end if
      else if(var%offset /= 0._rp) then
        fldi(:,:,:) = fldi(:,:,:) - var%offset
        if(myid == 0) print*, '  ', trim(var%name), ': offset subtracted: ', var%offset
      end if
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
      !
      ! field operations on the output: subtract the residual mean (on non-uniform grids, interpolation does
      ! not preserve the mean exactly), and add a constant (the frame velocity, for u, v and w)
      !
      if(var%is_remove_mean) then
        call field_stats(lo_o,hi_o,var%is_staggered,go,lo_o,fldo,pmin,pmax,pmean)
        fldo(:,:,:) = fldo(:,:,:) - pmean
        if(myid == 0) print*, '  ', trim(var%name), ': residual mean of the output subtracted: ', pmean
      end if
      if(var%added /= 0._rp) then
        fldo(:,:,:) = fldo(:,:,:) + var%added
        if(myid == 0) print*, '  ', trim(var%name), ': frame velocity added: ', var%added
      end if
      call field_stats(lo_o,hi_o,var%is_staggered,go,lo_o,fldo,pmin,pmax,pmean)
      call print_stats(trim(var%name)//' (output)',pmin,pmax,pmean)
      !
      ! write the output
      !
      time_o  = time_1
      istep_o = istep_1
      if(is_reset_time) then
        time_o  = 0._rp
        istep_o = 0
      end if
      fname_o = output_file_name(var%name)
      twi = MPI_WTIME()
      call write_field(trim(fname_o),trim(output_format),trim(var%name),ng_o,lo_o,hi_o,iprec_o,fldo,time_o,istep_o, &
                       merge(go(1)%xf,go(1)%xc,var%is_staggered(1)), &
                       merge(go(2)%xf,go(2)%xc,var%is_staggered(2)), &
                       merge(go(3)%xf,go(3)%xc,var%is_staggered(3)),is_compress_output)
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
    subroutine field_stats(lo,hi,is_staggered,g,lo_p,p,pmin,pmax,pmean)
      !
      ! global minimum, maximum and volume-weighted mean of p, over the points lo:hi of each task (within
      ! 1:n), with the widths of the control volumes of the variable's location (`wc` or `ws` of the grids)
      !
      implicit none
      integer         , intent(in), dimension(3) :: lo,hi
      logical         , intent(in), dimension(3) :: is_staggered
      type(grid_1d)   , intent(in), dimension(3) :: g
      integer         , intent(in), dimension(3) :: lo_p
      real(rp)        , intent(in), dimension(lo_p(1):,lo_p(2):,lo_p(3):) :: p
      real(rp)        , intent(out) :: pmin,pmax,pmean
      real(dp), allocatable, dimension(:) :: wx,wy,wz
      real(dp) :: psum(2)
      integer :: j,k
      !
      allocate(wx(lo(1):hi(1)),wy(lo(2):hi(2)),wz(lo(3):hi(3)))
      wx(:) = merge(g(1)%ws(lo(1):hi(1)),g(1)%wc(lo(1):hi(1)),is_staggered(1))
      wy(:) = merge(g(2)%ws(lo(2):hi(2)),g(2)%wc(lo(2):hi(2)),is_staggered(2))
      wz(:) = merge(g(3)%ws(lo(3):hi(3)),g(3)%wc(lo(3):hi(3)),is_staggered(3))
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
      pmean = psum(1)/psum(2)
    end subroutine field_stats
    !
    subroutine tile_coords(idir,lo,hi,tc)
      !
      ! coordinate (0 to ntiles-1) along direction idir of the tile of the input points lo:hi
      !
      implicit none
      integer, intent(in ) :: idir,lo,hi
      integer, intent(out), dimension(lo:) :: tc
      integer :: i
      do i=lo,hi
        tc(i) = modulo(i-1,ntiles(idir)*ng_i(idir))/ng_i(idir)
      end do
    end subroutine tile_coords
    !
    subroutine tile_means(lo,hi,is_staggered,g,lo_p,p,means)
      !
      ! volume-weighted mean of p over each tile, from the points lo:hi of each task (within the domain),
      ! summed as in `field_stats`
      !
      implicit none
      integer         , intent(in), dimension(3) :: lo,hi
      logical         , intent(in), dimension(3) :: is_staggered
      type(grid_1d)   , intent(in), dimension(3) :: g
      integer         , intent(in), dimension(3) :: lo_p
      real(rp)        , intent(in), dimension(lo_p(1):,lo_p(2):,lo_p(3):) :: p
      real(rp)        , intent(out), dimension(:) :: means
      real(dp), allocatable, dimension(:) :: wx,wy,wz
      integer , allocatable, dimension(:) :: tcy,tcz
      real(dp), dimension(2,size(means)) :: psum
      integer :: j,k,tx,it,i1,i2,n1
      !
      allocate(wx(lo(1):hi(1)),wy(lo(2):hi(2)),wz(lo(3):hi(3)),tcy(lo(2):hi(2)),tcz(lo(3):hi(3)))
      wx(:) = merge(g(1)%ws(lo(1):hi(1)),g(1)%wc(lo(1):hi(1)),is_staggered(1))
      wy(:) = merge(g(2)%ws(lo(2):hi(2)),g(2)%wc(lo(2):hi(2)),is_staggered(2))
      wz(:) = merge(g(3)%ws(lo(3):hi(3)),g(3)%wc(lo(3):hi(3)),is_staggered(3))
      call tile_coords(2,lo(2),hi(2),tcy)
      call tile_coords(3,lo(3),hi(3),tcz)
      n1 = ng_i(1)
      psum(:,:) = 0._dp
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do tx=0,ntiles(1)-1
            i1 = max(lo(1),tx*n1+1)
            i2 = min(hi(1),(tx+1)*n1)
            if(i2 < i1) cycle
            it = 1 + tx + ntiles(1)*(tcy(j) + ntiles(2)*tcz(k))
            psum(1,it) = psum(1,it) + wz(k)*wy(j)*sum(wx(i1:i2)*p(i1:i2,j,k))
            psum(2,it) = psum(2,it) + wz(k)*wy(j)*sum(wx(i1:i2))
          end do
        end do
      end do
      call MPI_ALLREDUCE(MPI_IN_PLACE,psum,size(psum),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      means(:) = psum(1,:)/psum(2,:)
    end subroutine tile_means
    !
    subroutine subtract_tile_values(lo,hi,values,p)
      !
      ! subtracts from each point of the block p(lo(1):hi(1),...) the value of its tile
      !
      implicit none
      integer , intent(in), dimension(3) :: lo,hi
      real(rp), intent(in), dimension(:) :: values
      real(rp), intent(inout), dimension(lo(1):,lo(2):,lo(3):) :: p
      integer , allocatable, dimension(:) :: tcx,tcy,tcz
      integer :: i,j,k
      !
      allocate(tcx(lo(1):hi(1)),tcy(lo(2):hi(2)),tcz(lo(3):hi(3)))
      call tile_coords(1,lo(1),hi(1),tcx)
      call tile_coords(2,lo(2),hi(2),tcy)
      call tile_coords(3,lo(3),hi(3),tcz)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            p(i,j,k) = p(i,j,k) - values(1 + tcx(i) + ntiles(1)*(tcy(j) + ntiles(2)*tcz(k)))
          end do
        end do
      end do
    end subroutine subtract_tile_values
    !
    subroutine print_stats(label,pmin,pmax,pmean)
      implicit none
      character(len=*), intent(in) :: label
      real(rp)        , intent(in) :: pmin,pmax,pmean
      if(myid == 0) print*, '  ', label, ': min = ', pmin, '| max = ', pmax, '| mean = ', pmean
    end subroutine print_stats
  end program interpolate_fields
