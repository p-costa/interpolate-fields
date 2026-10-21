module mod_io
  !
  ! I/O of CaNS checkpoints, in single or double precision, in the formats (`fmt`):
  !   'per-variable': MPI-IO files with the ng(1)*ng(2)*ng(3) values of one field (x fastest, without halos),
  !                   followed by the time and the time step
  !   'legacy'      : single MPI-IO files with u, v, w, p and any scalars, one after the other, followed by the
  !                   time and the time step (input only)
  !   'hdf5'        : HDF5 files with one field each (see `mod_io_hdf5`)
  !
  use mpi_f08
  use mod_common , only: rp,sp,dp,i8,myid,ierr,abort_run,CHUNK_BYTES,runs_1d
  use mod_io_hdf5, only: MSG_NO_HDF5
#if defined(_USE_HDF5)
  use mod_io_hdf5, only: h5_get_field_layout,h5_read_fldinfo,h5_read_block,h5_open_output,h5_write_output_block, &
                         h5_close_output
#endif
  implicit none
  private
  public get_checkpoint_layout,read_fldinfo,read_block,write_field,out_file,open_output,write_output_block, &
         close_output
  !
  ! an output checkpoint file, open for writing a field by blocks (`open_output`, `write_output_block`,
  ! `close_output`)
  !
  type out_file
    character(len=:), allocatable :: fname,fmt,varname
    integer :: iprec
    integer(MPI_OFFSET_KIND) :: disp,disp_info ! MPI-IO: byte offsets of the field and of the time and time step
    type(MPI_FILE) :: fh
    integer(i8), dimension(2) :: h5_ids        ! HDF5: file and dataset identifiers
  end type out_file
contains
  subroutine get_checkpoint_layout(fname,fmt,varname,ng,iprec,nfields)
    !
    ! precision of a checkpoint file (iprec = 4 or 8 bytes) and, for a legacy single-file checkpoint, its
    ! number of fields (nfields >= 4; nfields = 1 otherwise): for MPI-IO files, from the file size,
    ! (ng(1)*ng(2)*ng(3)*nfields + 2)*iprec bytes; for HDF5 files, from the type of the field dataset
    !
    implicit none
    character(len=*), intent(in ) :: fname,fmt,varname
    integer         , intent(in ), dimension(3) :: ng
    integer         , intent(out) :: iprec,nfields
    type(MPI_FILE) :: fh
    integer(MPI_OFFSET_KIND) :: fsize,ntot,nvals
    integer :: ip
    logical :: is_legacy
    character(len=32) :: csize
    !
    if(fmt == 'hdf5') then
#if defined(_USE_HDF5)
      call h5_get_field_layout(fname,varname,ng,iprec)
#else
      associate(dummy_varname => varname) ! silences warnings on arguments only used with HDF5
      end associate
      call abort_run(MSG_NO_HDF5)
#endif
      nfields = 1
      return
    end if
    is_legacy = fmt == 'legacy'
    call MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    if(ierr /= MPI_SUCCESS) call abort_run('could not open checkpoint file `'//fname//'`.')
    call MPI_FILE_GET_SIZE(fh,fsize,ierr)
    call MPI_FILE_CLOSE(fh,ierr)
    ntot = product(int(ng(:),MPI_OFFSET_KIND))
    iprec   = 0
    nfields = 0
    do ip=8,4,-4
      if(mod(fsize,int(ip,MPI_OFFSET_KIND)) /= 0) cycle
      nvals = fsize/ip - 2
      if(nvals < ntot .or. mod(nvals,ntot) /= 0) cycle
      if(     is_legacy .and. nvals <  4*ntot) cycle
      if(.not.is_legacy .and. nvals /=   ntot) cycle
      iprec   = ip
      nfields = int(nvals/ntot)
      exit
    end do
    if(iprec == 0) then
      write(csize,'(i0)') fsize
      if(is_legacy) then
        call abort_run('the size of checkpoint file `'//fname//'` ('//trim(csize)//' bytes) is not that of '// &
                       'a legacy checkpoint with `ng` of the input grid (u, v, w, p and any scalars, '// &
                       'then time and step, in single or double precision).')
      else
        call abort_run('the size of checkpoint file `'//fname//'` ('//trim(csize)//' bytes) is not that of '// &
                       'a field with `ng` of the input grid (then time and step, in single or double precision).')
      end if
    end if
  end subroutine get_checkpoint_layout
  !
  subroutine read_fldinfo(fname,fmt,disp,iprec,time,istep)
    !
    ! reads the time and the time step of a checkpoint file, stored at byte offset disp of an MPI-IO file
    !
    implicit none
    character(len=*), intent(in ) :: fname,fmt
    integer(MPI_OFFSET_KIND), intent(in) :: disp
    integer         , intent(in ) :: iprec
    real(rp)        , intent(out) :: time
    integer         , intent(out) :: istep
    type(MPI_FILE) :: fh
    real(sp), dimension(2) :: fldinfo_sp
    real(dp), dimension(2) :: fldinfo
    !
    if(fmt == 'hdf5') then
#if defined(_USE_HDF5)
      call h5_read_fldinfo(fname,time,istep)
#else
      call abort_run(MSG_NO_HDF5)
#endif
      return
    end if
    call MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    if(ierr /= MPI_SUCCESS) call abort_run('could not open checkpoint file `'//fname//'`.')
    if(myid == 0) then
      if(iprec == 4) then
        call MPI_FILE_READ_AT(fh,disp,fldinfo_sp,2,MPI_REAL,MPI_STATUS_IGNORE,ierr)
        fldinfo(:) = fldinfo_sp(:)
      else
        call MPI_FILE_READ_AT(fh,disp,fldinfo   ,2,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      end if
    end if
    call MPI_FILE_CLOSE(fh,ierr)
    call MPI_BCAST(fldinfo,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    time  =      fldinfo(1)
    istep = nint(fldinfo(2))
  end subroutine read_fldinfo
  !
  subroutine read_block(fnames,fmt,varname,disps,ng,nt,shifts,iprec,is_periodic,lo,hi,p)
    !
    ! reads the points lo:hi of the input field into p, indexed as the global grid; the input domain is
    ! made of nt(1) x nt(2) x nt(3) tiles (x fastest), each the field of one file, fnames(it), of ng points
    ! (stored at byte offset disps(it) of an MPI-IO file, or in dataset `fields/<varname>` of an HDF5 file),
    ! rolled by shifts(:,it) cells (point i of the tile holds point i+shift of the file); along periodic
    ! directions, the indices outside the domain wrap around, and along the other directions they are ghost
    ! points, which are not read; for each tile, the points to read form runs of indices per direction,
    ! in increasing order, selected by the file view (MPI-IO) or by a union of hyperslabs (HDF5), and all
    ! tasks read their part of the tile (possibly none) with one collective call
    !
    implicit none
    character(len=*), intent(in), dimension(:) :: fnames
    character(len=*), intent(in) :: fmt,varname
    integer(MPI_OFFSET_KIND), intent(in), dimension(:) :: disps
    integer , intent(in), dimension(3)   :: ng,nt
    integer , intent(in), dimension(:,:) :: shifts
    integer , intent(in) :: iprec
    logical , intent(in), dimension(3) :: is_periodic
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(inout), dimension(lo(1):,lo(2):,lo(3):) :: p
    type(runs_1d), dimension(3) :: runs
    integer , dimension(3) :: nr,tc,m
    integer , allocatable, dimension(:) :: posx,posy,posz,ix,iy,iz,px,py,pz
    real(sp), allocatable, dimension(:,:,:) :: buf_sp
    real(dp), allocatable, dimension(:,:,:) :: buf_dp
    integer :: it,i,j,k
    !
    allocate(posx(lo(1):hi(1)),posy(lo(2):hi(2)),posz(lo(3):hi(3)))
    do it=1,size(fnames)
      tc(:) = [mod(it-1,nt(1)),mod((it-1)/nt(1),nt(2)),(it-1)/(nt(1)*nt(2))]
      call get_tile_runs(ng(1),nt(1),tc(1),shifts(1,it),is_periodic(1),lo(1),hi(1),runs(1),posx)
      call get_tile_runs(ng(2),nt(2),tc(2),shifts(2,it),is_periodic(2),lo(2),hi(2),runs(2),posy)
      call get_tile_runs(ng(3),nt(3),tc(3),shifts(3,it),is_periodic(3),lo(3),hi(3),runs(3),posz)
      nr(:) = [sum(runs(1)%length),sum(runs(2)%length),sum(runs(3)%length)]
      if(iprec == 4) then
        allocate(buf_sp(nr(1),nr(2),nr(3)))
      else
        allocate(buf_dp(nr(1),nr(2),nr(3)))
      end if
      if(fmt == 'hdf5') then
#if defined(_USE_HDF5)
        if(iprec == 4) then
          call h5_read_block(fnames(it),varname,runs,buf_sp=buf_sp)
        else
          call h5_read_block(fnames(it),varname,runs,buf_dp=buf_dp)
        end if
#else
        associate(dummy_varname => varname) ! silences warnings on arguments only used with HDF5
        end associate
        call abort_run(MSG_NO_HDF5)
#endif
      else
        call read_runs_mpiio(fnames(it),disps(it))
      end if
      !
      ! copy the points read to the block, in double precision
      !
      ix = pack([(i,i=lo(1),hi(1))],posx(:) > 0); px = pack(posx(:),posx(:) > 0)
      iy = pack([(j,j=lo(2),hi(2))],posy(:) > 0); py = pack(posy(:),posy(:) > 0)
      iz = pack([(k,k=lo(3),hi(3))],posz(:) > 0); pz = pack(posz(:),posz(:) > 0)
      m(:) = [size(ix),size(iy),size(iz)]
      if(iprec == 4) then
        do k=1,m(3)
          do j=1,m(2)
            do i=1,m(1)
              p(ix(i),iy(j),iz(k)) = buf_sp(px(i),py(j),pz(k))
            end do
          end do
        end do
        deallocate(buf_sp)
      else
        do k=1,m(3)
          do j=1,m(2)
            do i=1,m(1)
              p(ix(i),iy(j),iz(k)) = buf_dp(px(i),py(j),pz(k))
            end do
          end do
        end do
        deallocate(buf_dp)
      end if
    end do
  contains
    subroutine read_runs_mpiio(fname,disp)
      implicit none
      character(len=*), intent(in) :: fname
      integer(MPI_OFFSET_KIND), intent(in) :: disp
      type(MPI_DATATYPE) :: etype,type_x,type_xr,type_y,type_yr,type_glob,type_loc
      type(MPI_FILE) :: fh
      integer(MPI_ADDRESS_KIND) :: esize
      !
      if(iprec == 4) then
        etype = MPI_REAL
      else
        etype = MPI_DOUBLE_PRECISION
      end if
      call MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
      if(ierr /= MPI_SUCCESS) call abort_run('could not open checkpoint file `'//fname//'`.')
      if(any(nr(:) == 0)) then
        !
        ! no point of this tile: an empty collective read
        !
        call MPI_FILE_SET_VIEW(fh,disp,etype,etype,'native',MPI_INFO_NULL,ierr)
        if(iprec == 4) then
          call MPI_FILE_READ_ALL(fh,buf_sp,0,etype,MPI_STATUS_IGNORE,ierr)
        else
          call MPI_FILE_READ_ALL(fh,buf_dp,0,etype,MPI_STATUS_IGNORE,ierr)
        end if
        call MPI_FILE_CLOSE(fh,ierr)
        return
      end if
      !
      ! file view: rows (x), planes (x-y) and the whole field (x-y-z), each resized to its extent in the file
      !
      esize = iprec
      call MPI_TYPE_INDEXED(runs(1)%n,runs(1)%length,runs(1)%start-1,etype  ,type_x,ierr)
      call MPI_TYPE_CREATE_RESIZED(type_x,0_MPI_ADDRESS_KIND,esize*ng(1)      ,type_xr,ierr)
      call MPI_TYPE_INDEXED(runs(2)%n,runs(2)%length,runs(2)%start-1,type_xr,type_y,ierr)
      call MPI_TYPE_CREATE_RESIZED(type_y,0_MPI_ADDRESS_KIND,esize*ng(1)*ng(2),type_yr,ierr)
      call MPI_TYPE_INDEXED(runs(3)%n,runs(3)%length,runs(3)%start-1,type_yr,type_glob,ierr)
      call MPI_TYPE_COMMIT(type_glob,ierr)
      call MPI_TYPE_CONTIGUOUS(nr(1)*nr(2),etype,type_loc,ierr)
      call MPI_TYPE_COMMIT(type_loc,ierr)
      call MPI_FILE_SET_VIEW(fh,disp,etype,type_glob,'native',MPI_INFO_NULL,ierr)
      if(iprec == 4) then
        call MPI_FILE_READ_ALL(fh,buf_sp,nr(3),type_loc,MPI_STATUS_IGNORE,ierr)
      else
        call MPI_FILE_READ_ALL(fh,buf_dp,nr(3),type_loc,MPI_STATUS_IGNORE,ierr)
      end if
      call MPI_FILE_CLOSE(fh,ierr)
      call MPI_TYPE_FREE(type_x   ,ierr)
      call MPI_TYPE_FREE(type_xr  ,ierr)
      call MPI_TYPE_FREE(type_y   ,ierr)
      call MPI_TYPE_FREE(type_yr  ,ierr)
      call MPI_TYPE_FREE(type_glob,ierr)
      call MPI_TYPE_FREE(type_loc ,ierr)
    end subroutine read_runs_mpiio
  end subroutine read_block
  !
  subroutine get_tile_runs(n,nt,tc,shift,is_periodic,lo,hi,runs,pos)
    !
    ! maps the indices lo:hi along a direction of the input domain, made of nt tiles of n points, to the
    ! indices of the file of tile tc (0 to nt-1), rolled by `shift` (along periodic directions, the indices
    ! wrap around the domain; along the other directions, nt = 1 and the indices outside 1:n are ghost
    ! points, not read): the file indices to read form `runs`, in increasing order, and pos(i) is the
    ! position of the file index of i among the points read (0 if i is not in the tile, or a ghost point)
    !
    implicit none
    integer      , intent(in ) :: n,nt,tc,shift
    logical      , intent(in ) :: is_periodic
    integer      , intent(in ) :: lo,hi
    type(runs_1d), intent(out) :: runs
    integer      , intent(out), dimension(lo:) :: pos
    integer, dimension(lo:hi) :: ifile
    logical, dimension(n) :: is_read
    integer, dimension(n) :: rank
    integer :: i,ig,nruns
    logical :: is_prev
    !
    ifile(:) = 0
    is_read(:) = .false.
    do i=lo,hi
      if(is_periodic) then
        ig = modulo(i-1,nt*n)
        if(ig/n /= tc) cycle
        ifile(i) = modulo(ig-tc*n+shift,n) + 1
      else
        if(i < 1 .or. i > n) cycle
        ifile(i) = i
      end if
      is_read(ifile(i)) = .true.
    end do
    nruns = count(is_read(1:1)) + count(is_read(2:n) .and. .not.is_read(1:n-1))
    runs%n = nruns
    allocate(runs%start(nruns),runs%length(nruns))
    nruns = 0
    is_prev = .false.
    do i=1,n
      if(is_read(i)) then
        if(.not.is_prev) then
          nruns = nruns + 1
          runs%start(nruns)  = i
          runs%length(nruns) = 0
        end if
        runs%length(nruns) = runs%length(nruns) + 1
      end if
      is_prev = is_read(i)
    end do
    rank(:) = 0
    nruns = 0
    do i=1,n
      if(is_read(i)) then
        nruns = nruns + 1
        rank(i) = nruns
      end if
    end do
    pos(:) = 0
    where(ifile(:) > 0) pos(:) = rank(max(ifile(:),1))
  end subroutine get_tile_runs
  !
  subroutine write_field(fname,fmt,varname,ng,lo,hi,iprec,p,time,istep,x,y,z,is_compress)
    !
    ! writes a per-variable checkpoint file (the field, then the time and the time step; MPI-IO or HDF5) in
    ! single or double precision (iprec = 4 or 8 bytes), from the block p(lo(1):hi(1),...) of each task;
    ! HDF5 files also hold the coordinates x, y, z of the variable's location, and can be compressed, in
    ! chunks based on the smallest block of the tasks
    !
    implicit none
    character(len=*), intent(in) :: fname,fmt,varname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: iprec
    real(rp), intent(in), contiguous, dimension(lo(1):,lo(2):,lo(3):) :: p
    real(rp), intent(in) :: time
    integer , intent(in) :: istep
    real(rp), intent(in), dimension(0:) :: x,y,z
    logical , intent(in) :: is_compress
    type(out_file) :: f
    integer, dimension(3) :: nchunk
    !
    call MPI_ALLREDUCE(hi-lo+1,nchunk,3,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call open_output(fname,fmt,varname,ng,iprec,f,is_compress=is_compress,nchunk=nchunk)
    call write_output_block(f,ng,lo,hi,p)
    call close_output(f,time,istep,x,y,z)
  end subroutine write_field
  !
  subroutine open_output(fname,fmt,varname,ng,iprec,f,is_compress,nchunk,disp,disp_info,is_existing)
    !
    ! opens an output checkpoint file f for writing a field of ng points by blocks, in single or double
    ! precision (iprec = 4 or 8 bytes): an MPI-IO file with the field at byte offset disp (default: 0) and
    ! the time and the time step at byte offset disp_info (default: after the field), created unless
    ! is_existing; or an HDF5 file with dataset `fields/<varname>`, compressed if is_compress, in chunks
    ! based on blocks of nchunk points
    !
    implicit none
    character(len=*), intent(in) :: fname,fmt,varname
    integer , intent(in), dimension(3) :: ng
    integer , intent(in) :: iprec
    type(out_file), intent(out) :: f
    logical , intent(in), optional :: is_compress
    integer , intent(in), dimension(3), optional :: nchunk
    integer(MPI_OFFSET_KIND), intent(in), optional :: disp,disp_info
    logical , intent(in), optional :: is_existing
    logical :: is_existing_
    !
    f%fname   = fname
    f%fmt     = fmt
    f%varname = varname
    f%iprec   = iprec
    if(fmt == 'hdf5') then
#if defined(_USE_HDF5)
      call h5_open_output(fname,varname,ng,iprec,is_compress,nchunk,f%h5_ids)
#else
      if(present(is_compress) .or. present(nchunk)) continue ! (arguments only used with HDF5)
      call abort_run(MSG_NO_HDF5)
#endif
      return
    end if
    f%disp = 0
    if(present(disp)) f%disp = disp
    f%disp_info = f%disp + product(int(ng(:),MPI_OFFSET_KIND))*iprec
    if(present(disp_info)) f%disp_info = disp_info
    is_existing_ = .false.
    if(present(is_existing)) is_existing_ = is_existing
    if(is_existing_) then
      call MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY,MPI_INFO_NULL,f%fh,ierr)
      if(ierr /= MPI_SUCCESS) call abort_run('could not open output file `'//fname//'`.')
    else
      call MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,f%fh,ierr)
      if(ierr /= MPI_SUCCESS) call abort_run('could not create output file `'//fname//'`.')
      call MPI_FILE_SET_SIZE(f%fh,0_MPI_OFFSET_KIND,ierr)
    end if
  end subroutine open_output
  !
  subroutine write_output_block(f,ng,lo,hi,p)
    !
    ! writes the block p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) of each task (none if hi < lo) to an output
    ! file, converted to single precision by chunks of planes (`CHUNK_BYTES`) if needed; collective
    !
    implicit none
    type(out_file), intent(inout) :: f
    integer , intent(in), dimension(3) :: ng,lo,hi
    real(rp), intent(in), contiguous, dimension(lo(1):,lo(2):,lo(3):) :: p
    real(sp), allocatable, dimension(:,:,:) :: buf_sp
    type(MPI_DATATYPE) :: etype,type_glob,type_loc
    integer , dimension(3) :: n
    integer :: nzc,nchunks,ic,k0,nk
    logical :: is_empty
    !
    if(f%fmt == 'hdf5') then
#if defined(_USE_HDF5)
      call h5_write_output_block(f%h5_ids,f%iprec,lo,hi,p)
#endif
      return
    end if
    n(:) = max(hi(:)-lo(:)+1,0)
    is_empty = any(n(:) == 0)
    if(f%iprec == 4) then
      etype = MPI_REAL
    else
      etype = MPI_DOUBLE_PRECISION
    end if
    if(is_empty) then
      call MPI_FILE_SET_VIEW(f%fh,f%disp,etype,etype,'native',MPI_INFO_NULL,ierr)
      call MPI_TYPE_CONTIGUOUS(1,etype,type_loc,ierr)
    else
      call MPI_TYPE_CREATE_SUBARRAY(3,ng,n,lo-1,MPI_ORDER_FORTRAN,etype,type_glob,ierr)
      call MPI_TYPE_COMMIT(type_glob,ierr)
      call MPI_FILE_SET_VIEW(f%fh,f%disp,etype,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_TYPE_FREE(type_glob,ierr)
      call MPI_TYPE_CONTIGUOUS(n(1)*n(2),etype,type_loc,ierr)
    end if
    call MPI_TYPE_COMMIT(type_loc,ierr)
    if(f%iprec == 4) then
      !
      ! convert and write in chunks of nzc planes, so that the single-precision copy stays small;
      ! the file pointer advances within the view, and all tasks make the same number of
      ! collective calls (with no data once their block is written)
      !
      nzc = 1
      if(.not.is_empty) nzc = int(max(1_i8,min(int(n(3),i8),CHUNK_BYTES/(4*int(n(1),i8)*n(2)))))
      nchunks = (n(3)+nzc-1)/nzc
      call MPI_ALLREDUCE(MPI_IN_PLACE,nchunks,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      allocate(buf_sp(n(1),n(2),nzc))
      do ic=1,nchunks
        k0 = lo(3) + (ic-1)*nzc
        nk = max(0,min(nzc,hi(3)-k0+1))
        if(is_empty) nk = 0
        if(nk > 0) buf_sp(:,:,1:nk) = real(p(:,:,k0:k0+nk-1),sp)
        call MPI_FILE_WRITE_ALL(f%fh,buf_sp,nk,type_loc,MPI_STATUS_IGNORE,ierr)
      end do
      deallocate(buf_sp)
    else
      call MPI_FILE_WRITE_ALL(f%fh,p,n(3),type_loc,MPI_STATUS_IGNORE,ierr)
    end if
    call MPI_TYPE_FREE(type_loc,ierr)
  end subroutine write_output_block
  !
  subroutine close_output(f,time,istep,x,y,z)
    !
    ! writes the time and the time step (and, for HDF5 files, the coordinates x, y, z of the variable's
    ! location, if present), and closes an output file
    !
    implicit none
    type(out_file), intent(inout) :: f
    real(rp), intent(in) :: time
    integer , intent(in) :: istep
    real(rp), intent(in), dimension(0:), optional :: x,y,z
    real(sp), dimension(2) :: fldinfo_sp
    real(dp), dimension(2) :: fldinfo
    type(MPI_DATATYPE) :: etype
    integer :: nreals_myid
    !
    if(f%fmt == 'hdf5') then
#if defined(_USE_HDF5)
      call h5_close_output(f%h5_ids,f%fname,f%iprec,time,istep,x,y,z)
#endif
      return
    end if
    if(present(x) .or. present(y) .or. present(z)) continue ! (coordinates only written to HDF5 files)
    if(f%iprec == 4) then
      etype = MPI_REAL
    else
      etype = MPI_DOUBLE_PRECISION
    end if
    call MPI_FILE_SET_VIEW(f%fh,f%disp_info,etype,etype,'native',MPI_INFO_NULL,ierr)
    fldinfo(:) = [time,1._rp*istep]
    fldinfo_sp(:) = real(fldinfo(:),sp)
    nreals_myid = 0
    if(myid == 0) nreals_myid = 2
    if(f%iprec == 4) then
      call MPI_FILE_WRITE(f%fh,fldinfo_sp,nreals_myid,etype,MPI_STATUS_IGNORE,ierr)
    else
      call MPI_FILE_WRITE(f%fh,fldinfo   ,nreals_myid,etype,MPI_STATUS_IGNORE,ierr)
    end if
    call MPI_FILE_CLOSE(f%fh,ierr)
  end subroutine close_output
end module mod_io
