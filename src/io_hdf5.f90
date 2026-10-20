module mod_io_hdf5
  !
  ! HDF5 I/O in the format of CaNS and CaNS-EIGEN (`io_field_hdf5`, with `io_backend = 'hdf5'`): one file per
  ! variable, with the field in dataset `fields/<var>` (ng(1) x ng(2) x ng(3), x fastest), the time and the time
  ! step in `meta/time` and `meta/istep`, and the coordinates of the variable's location, ghost points included,
  ! in `grid/x`, `grid/y` and `grid/z`; and grid files `grid_<dir>.h5`, with the cell centres in `rc` and the
  ! cell faces in `rf`
  !
  ! HDF5 support is optional (`USE_HDF5=1` at build time, `_USE_HDF5` in the code); without it, only
  ! `is_hdf5_enabled`, `is_hdf5_file` and the error message `MSG_NO_HDF5` are defined
  !
#if defined(_USE_HDF5)
  use hdf5
#endif
  use mpi_f08
  use mod_common, only: rp,sp,dp,i8,myid,ierr,abort_run,CHUNK_BYTES
  implicit none
  private
  public is_hdf5_enabled,is_hdf5_file,MSG_NO_HDF5
#if defined(_USE_HDF5)
  public h5_get_field_layout,h5_read_fldinfo,h5_read_block,h5_write_field,h5_read_grid
  logical, parameter :: is_hdf5_enabled = .true.
#else
  logical, parameter :: is_hdf5_enabled = .false.
#endif
  character(len=*), parameter :: MSG_NO_HDF5 = 'HDF5 support is not compiled in: rebuild with `USE_HDF5=1`.'
contains
  logical function is_hdf5_file(fname)
    !
    ! HDF5 files are recognized by their extension, as in CaNS
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer :: i
    is_hdf5_file = .false.
    i = index(fname,'.',back=.true.)
    if(i == 0) return
    select case(fname(i:))
    case('.h5','.hdf5','.hdf')
      is_hdf5_file = .true.
    end select
  end function is_hdf5_file
#if defined(_USE_HDF5)
  !
  subroutine h5_get_field_layout(fname,varname,ng,iprec)
    !
    ! checks that dataset `fields/<varname>` of an HDF5 checkpoint has the dimensions ng, and gets its precision
    ! (iprec = 4 or 8 bytes); task 0 checks the file, and all tasks stop on an error
    !
    implicit none
    character(len=*), intent(in ) :: fname,varname
    integer         , intent(in ), dimension(3) :: ng
    integer         , intent(out) :: iprec
    integer(HID_T) :: file_id,dset,dtype,dspace
    integer(HSIZE_T), dimension(3) :: dims,maxdims
    integer(SIZE_T) :: tsize
    integer :: hdferr,tclass,rank,istat
    logical :: is_exist
    character(len=64) :: cdims
    !
    ! istat: 0 (valid), 1 (no file), 2 (no dataset), 3 (wrong dimensions), 4 (not single or double precision)
    !
    istat = 0
    iprec = 0
    dims(:) = 0
    if(myid == 0) then
      inquire(file=fname,exist=is_exist)
      if(.not.is_exist) istat = 1
      if(istat == 0) then
        call h5open_f(hdferr)
        call h5fopen_f(fname,H5F_ACC_RDONLY_F,file_id,hdferr)
        if(hdferr < 0) istat = 1
      end if
      if(istat == 0) then
        call h5lexists_f(file_id,'fields',is_exist,hdferr)
        if(is_exist) call h5lexists_f(file_id,'fields/'//varname,is_exist,hdferr)
        if(.not.is_exist) istat = 2
        if(istat == 0) then
          call h5dopen_f(file_id,'fields/'//varname,dset,hdferr)
          call h5dget_space_f(dset,dspace,hdferr)
          call h5sget_simple_extent_ndims_f(dspace,rank,hdferr)
          if(rank == 3) then
            call h5sget_simple_extent_dims_f(dspace,dims,maxdims,hdferr)
            if(any(dims(:) /= int(ng(:),HSIZE_T))) istat = 3
          else
            istat = 3
          end if
          call h5sclose_f(dspace,hdferr)
          call h5dget_type_f(dset,dtype,hdferr)
          call h5tget_class_f(dtype,tclass,hdferr)
          call h5tget_size_f(dtype,tsize,hdferr)
          call h5tclose_f(dtype,hdferr)
          call h5dclose_f(dset,hdferr)
          if(tclass == H5T_FLOAT_F .and. (tsize == 4 .or. tsize == 8)) then
            iprec = int(tsize)
          else if(istat == 0) then
            istat = 4
          end if
        end if
        call h5fclose_f(file_id,hdferr)
        call h5close_f(hdferr)
      end if
    end if
    call MPI_BCAST(istat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(iprec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    select case(istat)
    case(1)
      call abort_run('could not open HDF5 checkpoint file `'//fname//'`.')
    case(2)
      call abort_run('checkpoint file `'//fname//'` has no dataset `fields/'//varname//'`.')
    case(3)
      write(cdims,'(3(i0,1x))') ng(:)
      call abort_run('dataset `fields/'//varname//'` of `'//fname//'` is not of dimensions `ng` = '// &
                     trim(cdims)//' of the input grid.')
    case(4)
      call abort_run('dataset `fields/'//varname//'` of `'//fname//'` is not in single or double precision.')
    end select
  end subroutine h5_get_field_layout
  !
  subroutine h5_read_fldinfo(fname,time,istep)
    !
    ! reads the time and the time step of an HDF5 checkpoint (task 0, then broadcast)
    !
    implicit none
    character(len=*), intent(in ) :: fname
    real(rp)        , intent(out) :: time
    integer         , intent(out) :: istep
    integer(HID_T) :: file_id,dset
    real(dp), dimension(1) :: meta_time
    integer , dimension(1) :: meta_istep
    integer :: hdferr,istat
    logical :: is_exist
    !
    istat = 0
    meta_time(1)  = 0._dp
    meta_istep(1) = 0
    if(myid == 0) then
      call h5open_f(hdferr)
      call h5fopen_f(fname,H5F_ACC_RDONLY_F,file_id,hdferr)
      call h5lexists_f(file_id,'meta',is_exist,hdferr)
      if(is_exist) call h5lexists_f(file_id,'meta/time' ,is_exist,hdferr)
      if(is_exist) call h5lexists_f(file_id,'meta/istep',is_exist,hdferr)
      if(is_exist) then
        call h5dopen_f(file_id,'meta/time',dset,hdferr)
        call h5dread_f(dset,H5T_NATIVE_DOUBLE,meta_time,[1_HSIZE_T],hdferr)
        call h5dclose_f(dset,hdferr)
        call h5dopen_f(file_id,'meta/istep',dset,hdferr)
        call h5dread_f(dset,H5T_NATIVE_INTEGER,meta_istep,[1_HSIZE_T],hdferr)
        call h5dclose_f(dset,hdferr)
      else
        istat = 1
      end if
      call h5fclose_f(file_id,hdferr)
      call h5close_f(hdferr)
    end if
    call MPI_BCAST(istat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if(istat /= 0) call abort_run('checkpoint file `'//fname//'` has no `meta/time` and `meta/istep`.')
    call MPI_BCAST(meta_time ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(meta_istep,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
    time  = meta_time(1)
    istep = meta_istep(1)
  end subroutine h5_read_fldinfo
  !
  subroutine h5_read_block(fname,varname,nb,bstart,blen,buf_sp,buf_dp)
    !
    ! reads the points of dataset `fields/<varname>` with indices in the runs bstart(1:nb(idir),idir) to
    ! bstart+blen-1 along each direction (a union of hyperslabs, in increasing order) into buf_sp or buf_dp
    ! (whichever is present, in the precision of the file); all tasks read with one collective call
    !
    implicit none
    character(len=*), intent(in) :: fname,varname
    integer , intent(in), dimension(3)   :: nb
    integer , intent(in), dimension(2,3) :: bstart,blen
    real(sp), intent(out), dimension(:,:,:), optional :: buf_sp
    real(dp), intent(out), dimension(:,:,:), optional :: buf_dp
    integer(HID_T) :: file_id,dset,filespace,memspace,plist_id
    integer(HSIZE_T), dimension(3) :: offset,count,nr
    integer :: hdferr,istat,ib,jb,kb,iop
    !
    call h5open_f(hdferr)
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdferr)
    call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD%MPI_VAL,MPI_INFO_NULL%MPI_VAL,hdferr)
    call h5fopen_f(fname,H5F_ACC_RDONLY_F,file_id,istat,access_prp=plist_id)
    call h5pclose_f(plist_id,hdferr)
    if(istat < 0) call abort_run('could not open checkpoint file `'//fname//'`.')
    call h5dopen_f(file_id,'fields/'//varname,dset,hdferr)
    call h5dget_space_f(dset,filespace,hdferr)
    iop = H5S_SELECT_SET_F
    do kb=1,nb(3)
      do jb=1,nb(2)
        do ib=1,nb(1)
          offset(:) = [bstart(ib,1),bstart(jb,2),bstart(kb,3)] - 1
          count(:)  = [blen(ib,1)  ,blen(jb,2)  ,blen(kb,3)  ]
          call h5sselect_hyperslab_f(filespace,iop,offset,count,hdferr)
          iop = H5S_SELECT_OR_F
        end do
      end do
    end do
    if(present(buf_sp)) then
      nr(:) = shape(buf_sp)
    else
      nr(:) = shape(buf_dp)
    end if
    call h5screate_simple_f(3,nr,memspace,hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdferr)
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdferr)
    if(present(buf_sp)) then
      call h5dread_f(dset,H5T_NATIVE_REAL  ,buf_sp,nr,hdferr,mem_space_id=memspace,file_space_id=filespace, &
                     xfer_prp=plist_id)
    else
      call h5dread_f(dset,H5T_NATIVE_DOUBLE,buf_dp,nr,hdferr,mem_space_id=memspace,file_space_id=filespace, &
                     xfer_prp=plist_id)
    end if
    if(hdferr < 0) call abort_run('could not read dataset `fields/'//varname//'` of `'//fname//'`.')
    call h5pclose_f(plist_id,hdferr)
    call h5sclose_f(memspace,hdferr)
    call h5sclose_f(filespace,hdferr)
    call h5dclose_f(dset,hdferr)
    call h5fclose_f(file_id,hdferr)
    call h5close_f(hdferr)
  end subroutine h5_read_block
  !
  subroutine h5_write_field(fname,varname,ng,lo,hi,iprec,p,time,istep,x,y,z,is_compress)
    !
    ! writes an HDF5 checkpoint in single or double precision (iprec = 4 or 8 bytes): each task writes its
    ! block p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), converted to single precision by chunks of planes if needed;
    ! then task 0 writes the time, the time step and the coordinates x(0:ng(1)+1), y(0:ng(2)+1), z(0:ng(3)+1)
    ! of the variable's location; with compression, the dataset is chunked and gzip-compressed (with the
    ! shuffle filter), as in CaNS (`hdf5_checkpoint_chunk_dims`)
    !
    implicit none
    character(len=*), intent(in) :: fname,varname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: iprec
    real(rp), intent(in), contiguous, dimension(lo(1):,lo(2):,lo(3):) :: p
    real(rp), intent(in) :: time
    integer , intent(in) :: istep
    real(rp), intent(in), dimension(0:) :: x,y,z
    logical , intent(in) :: is_compress
    integer, parameter :: COMPRESSION_LEVEL = 4
    integer(i8), parameter :: CHUNK_MAX = 2*2_i8**20 ! 2 MiB, as in CaNS
    real(sp), allocatable, dimension(:,:,:) :: buf_sp
    integer(HID_T) :: file_id,group_id,dset,filespace,memspace,plist_id,dcpl_id,dtype
    integer(HSIZE_T), dimension(3) :: dims,offset,count,chunk_dims
    integer , dimension(3) :: n,nn
    integer :: hdferr,istat,idir,nzc,nchunks,ic,k0,nk
    logical :: is_reduced
    !
    n(:) = hi(:)-lo(:)+1
    call h5open_f(hdferr)
    dtype = H5T_NATIVE_DOUBLE ! (the HDF5 predefined types are set by `h5open_f`)
    if(iprec == 4) dtype = H5T_NATIVE_REAL
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdferr)
    call h5pset_fapl_mpio_f(plist_id,MPI_COMM_WORLD%MPI_VAL,MPI_INFO_NULL%MPI_VAL,hdferr)
    call h5fcreate_f(fname,H5F_ACC_TRUNC_F,file_id,istat,access_prp=plist_id)
    call h5pclose_f(plist_id,hdferr)
    if(istat < 0) call abort_run('could not create output file `'//fname//'`.')
    call h5gcreate_f(file_id,'fields',group_id,hdferr)
    dims(:) = ng(:)
    call h5screate_simple_f(3,dims,filespace,hdferr)
    call h5pcreate_f(H5P_DATASET_CREATE_F,dcpl_id,hdferr)
    if(is_compress) then
      !
      ! chunks: the smallest block over the tasks, halved along each direction in turn down to `CHUNK_MAX`
      !
      call MPI_ALLREDUCE(n,nn,3,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
      do while(product(int(nn,i8))*iprec > CHUNK_MAX)
        is_reduced = .false.
        do idir=1,3
          if(nn(idir) > 1) then
            nn(idir) = max(1,nn(idir)/2)
            is_reduced = .true.
            if(product(int(nn,i8))*iprec <= CHUNK_MAX) exit
          end if
        end do
        if(.not.is_reduced) exit
      end do
      chunk_dims(:) = nn(:)
      call h5pset_chunk_f(dcpl_id,3,chunk_dims,hdferr)
      call h5pset_shuffle_f(dcpl_id,hdferr)
      call h5pset_deflate_f(dcpl_id,COMPRESSION_LEVEL,hdferr)
    end if
    call h5dcreate_f(group_id,varname,dtype,filespace,dset,hdferr,dcpl_id=dcpl_id)
    call h5pclose_f(dcpl_id,hdferr)
    call h5gclose_f(group_id,hdferr)
    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdferr)
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdferr)
    if(iprec == 4) then
      !
      ! convert and write in chunks of nzc planes; all tasks make the same number of collective calls
      ! (with empty selections once their block is written)
      !
      nzc = int(max(1_i8,min(int(n(3),i8),CHUNK_BYTES/(4*int(n(1),i8)*n(2)))))
      nchunks = (n(3)+nzc-1)/nzc
      call MPI_ALLREDUCE(MPI_IN_PLACE,nchunks,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      allocate(buf_sp(n(1),n(2),nzc))
      count(:) = [n(1),n(2),nzc]
      call h5screate_simple_f(3,count,memspace,hdferr)
      do ic=1,nchunks
        k0 = lo(3) + (ic-1)*nzc
        nk = max(0,min(nzc,hi(3)-k0+1))
        if(nk > 0) then
          buf_sp(:,:,1:nk) = real(p(:,:,k0:k0+nk-1),sp)
          offset(:) = [lo(1),lo(2),k0] - 1
          count(:)  = [n(1),n(2),nk]
          call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,count,hdferr)
          call h5sselect_hyperslab_f(memspace ,H5S_SELECT_SET_F,[0_HSIZE_T,0_HSIZE_T,0_HSIZE_T],count,hdferr)
        else
          call h5sselect_none_f(filespace,hdferr)
          call h5sselect_none_f(memspace ,hdferr)
        end if
        call h5dwrite_f(dset,dtype,buf_sp,count,hdferr,mem_space_id=memspace,file_space_id=filespace, &
                        xfer_prp=plist_id)
      end do
      deallocate(buf_sp)
    else
      offset(:) = lo(:) - 1
      count(:)  = n(:)
      call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offset,count,hdferr)
      call h5screate_simple_f(3,count,memspace,hdferr)
      call h5dwrite_f(dset,dtype,p,count,hdferr,mem_space_id=memspace,file_space_id=filespace,xfer_prp=plist_id)
    end if
    if(hdferr < 0) call abort_run('could not write dataset `fields/'//varname//'` of `'//fname//'`.')
    call h5pclose_f(plist_id,hdferr)
    call h5sclose_f(memspace,hdferr)
    call h5sclose_f(filespace,hdferr)
    call h5dclose_f(dset,hdferr)
    call h5fclose_f(file_id,hdferr)
    !
    ! metadata
    !
    if(myid == 0) then
      call h5fopen_f(fname,H5F_ACC_RDWR_F,file_id,hdferr)
      call h5gcreate_f(file_id,'meta',group_id,hdferr)
      call write_1d(group_id,'time',dtype,[time])
      call write_1d_int(group_id,'istep',[istep])
      call h5gclose_f(group_id,hdferr)
      call h5gcreate_f(file_id,'grid',group_id,hdferr)
      call write_1d(group_id,'x',dtype,x)
      call write_1d(group_id,'y',dtype,y)
      call write_1d(group_id,'z',dtype,z)
      call h5gclose_f(group_id,hdferr)
      call h5fclose_f(file_id,hdferr)
    end if
    call h5close_f(hdferr)
  contains
    subroutine write_1d(loc_id,name,dtype,a)
      implicit none
      integer(HID_T)  , intent(in) :: loc_id,dtype
      character(len=*), intent(in) :: name
      real(rp)        , intent(in), dimension(:) :: a
      integer(HID_T) :: space,dset
      integer(HSIZE_T), dimension(1) :: dims
      integer :: hdferr
      dims(1) = size(a)
      call h5screate_simple_f(1,dims,space,hdferr)
      call h5dcreate_f(loc_id,name,dtype,space,dset,hdferr)
      if(dtype == H5T_NATIVE_REAL) then
        call h5dwrite_f(dset,dtype,real(a,sp),dims,hdferr)
      else
        call h5dwrite_f(dset,dtype,real(a,dp),dims,hdferr)
      end if
      call h5dclose_f(dset,hdferr)
      call h5sclose_f(space,hdferr)
    end subroutine write_1d
    subroutine write_1d_int(loc_id,name,a)
      implicit none
      integer(HID_T)  , intent(in) :: loc_id
      character(len=*), intent(in) :: name
      integer         , intent(in), dimension(:) :: a
      integer(HID_T) :: space,dset
      integer(HSIZE_T), dimension(1) :: dims
      integer :: hdferr
      dims(1) = size(a)
      call h5screate_simple_f(1,dims,space,hdferr)
      call h5dcreate_f(loc_id,name,H5T_NATIVE_INTEGER,space,dset,hdferr)
      call h5dwrite_f(dset,H5T_NATIVE_INTEGER,a,dims,hdferr)
      call h5dclose_f(dset,hdferr)
      call h5sclose_f(space,hdferr)
    end subroutine write_1d_int
  end subroutine h5_write_field
  !
  subroutine h5_read_grid(fname,n,istat,iprec,xc,xf)
    !
    ! reads the cell centres `rc` and faces `rf` (n values each) of a grid file `grid_<dir>.h5` written by CaNS;
    ! called by one task: istat = 0 (valid), 1 (no file), 2 (missing datasets or wrong size), 3 (not real)
    !
    implicit none
    character(len=*), intent(in ) :: fname
    integer         , intent(in ) :: n
    integer         , intent(out) :: istat,iprec
    real(dp)        , intent(out), dimension(n) :: xc,xf
    integer(HID_T) :: file_id
    integer :: hdferr,iprec_c,iprec_f
    logical :: is_exist
    !
    istat = 0
    iprec = 0
    inquire(file=fname,exist=is_exist)
    if(.not.is_exist) then
      istat = 1
      return
    end if
    call h5open_f(hdferr)
    call h5fopen_f(fname,H5F_ACC_RDONLY_F,file_id,hdferr)
    if(hdferr < 0) then
      istat = 1
    else
      call read_dset('rc',xc,iprec_c)
      if(istat == 0) call read_dset('rf',xf,iprec_f)
      if(istat == 0) iprec = min(iprec_c,iprec_f)
      call h5fclose_f(file_id,hdferr)
    end if
    call h5close_f(hdferr)
  contains
    subroutine read_dset(name,a,ip)
      implicit none
      character(len=*), intent(in ) :: name
      real(dp)        , intent(out), dimension(:) :: a
      integer         , intent(out) :: ip
      integer(HID_T) :: dset,dspace,dtype
      integer(HSIZE_T), dimension(1) :: dims,maxdims
      integer(SIZE_T) :: tsize
      integer :: hdferr,rank,tclass
      logical :: is_exist
      ip = 0
      call h5lexists_f(file_id,name,is_exist,hdferr)
      if(.not.is_exist) then
        istat = 2
        return
      end if
      call h5dopen_f(file_id,name,dset,hdferr)
      call h5dget_space_f(dset,dspace,hdferr)
      call h5sget_simple_extent_ndims_f(dspace,rank,hdferr)
      dims(1) = 0
      if(rank == 1) call h5sget_simple_extent_dims_f(dspace,dims,maxdims,hdferr)
      call h5sclose_f(dspace,hdferr)
      call h5dget_type_f(dset,dtype,hdferr)
      call h5tget_class_f(dtype,tclass,hdferr)
      call h5tget_size_f(dtype,tsize,hdferr)
      call h5tclose_f(dtype,hdferr)
      if(rank /= 1 .or. dims(1) /= int(size(a),HSIZE_T)) then
        istat = 2
      else if(tclass /= H5T_FLOAT_F) then
        istat = 3
      else
        ip = int(tsize)
        dims(1) = size(a)
        call h5dread_f(dset,H5T_NATIVE_DOUBLE,a,dims,hdferr)
      end if
      call h5dclose_f(dset,hdferr)
    end subroutine read_dset
  end subroutine h5_read_grid
#endif
end module mod_io_hdf5
