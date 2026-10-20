module mod_io
  !
  ! MPI-IO of CaNS checkpoints: per-variable files (the ng(1)*ng(2)*ng(3) values of one field, x fastest,
  ! without halos, followed by the time and the time step) and legacy single files (u, v, w, p and any
  ! scalars, one after the other, followed by the time and the time step), in single or double precision
  !
  use mpi_f08
  use mod_common, only: rp,sp,dp,i8,myid,ierr,abort_run
  implicit none
  private
  public get_checkpoint_layout,read_fldinfo,read_block,write_field,CHUNK_BYTES
  integer(i8), parameter :: CHUNK_BYTES = 64*2_i8**20 ! size of the buffer for single-precision output
contains
  subroutine get_checkpoint_layout(fname,ng,is_legacy,iprec,nfields)
    !
    ! detects from the size of a checkpoint file, (ng(1)*ng(2)*ng(3)*nfields + 2)*iprec bytes, its precision
    ! (iprec = 4 or 8 bytes) and, for a legacy single-file checkpoint, its number of fields (nfields >= 4);
    ! nfields = 1 for a per-variable file
    !
    implicit none
    character(len=*), intent(in ) :: fname
    integer         , intent(in ), dimension(3) :: ng
    logical         , intent(in ) :: is_legacy
    integer         , intent(out) :: iprec,nfields
    type(MPI_FILE) :: fh
    integer(MPI_OFFSET_KIND) :: fsize,ntot,nvals
    integer :: ip
    character(len=32) :: csize
    !
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
  subroutine read_fldinfo(fname,disp,iprec,time,istep)
    !
    ! reads the time and the time step, stored at byte offset disp of a checkpoint file
    !
    implicit none
    character(len=*), intent(in ) :: fname
    integer(MPI_OFFSET_KIND), intent(in) :: disp
    integer         , intent(in ) :: iprec
    real(rp)        , intent(out) :: time
    integer         , intent(out) :: istep
    type(MPI_FILE) :: fh
    real(sp), dimension(2) :: fldinfo_sp
    real(dp), dimension(2) :: fldinfo
    !
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
  subroutine read_block(fname,disp,ng,iprec,is_periodic,lo,hi,p)
    !
    ! reads the points lo:hi of a field stored at byte offset disp of a file into p, indexed as the global
    ! grid; along periodic directions, the indices outside 1:ng wrap around, and along the other directions
    ! they are ghost points, which are not read; the file view selects the points to read (at most two runs
    ! of indices per direction, in increasing order), and all tasks read their blocks with one collective call
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer(MPI_OFFSET_KIND), intent(in) :: disp
    integer , intent(in), dimension(3) :: ng
    integer , intent(in) :: iprec
    logical , intent(in), dimension(3) :: is_periodic
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(inout), dimension(lo(1):,lo(2):,lo(3):) :: p
    integer , dimension(2,3) :: bstart,blen
    integer , dimension(3) :: nb,nr,ia,ib
    integer , allocatable, dimension(:) :: posx,posy,posz
    real(sp), allocatable, dimension(:,:,:) :: buf_sp
    real(dp), allocatable, dimension(:,:,:) :: buf_dp
    type(MPI_DATATYPE) :: etype,type_x,type_xr,type_y,type_yr,type_glob,type_loc
    type(MPI_FILE) :: fh
    integer(MPI_ADDRESS_KIND) :: esize
    integer :: i,j,k,jj,kk
    !
    allocate(posx(lo(1):hi(1)),posy(lo(2):hi(2)),posz(lo(3):hi(3)))
    call get_file_runs(ng(1),is_periodic(1),lo(1),hi(1),nb(1),bstart(:,1),blen(:,1),ia(1),ib(1),posx)
    call get_file_runs(ng(2),is_periodic(2),lo(2),hi(2),nb(2),bstart(:,2),blen(:,2),ia(2),ib(2),posy)
    call get_file_runs(ng(3),is_periodic(3),lo(3),hi(3),nb(3),bstart(:,3),blen(:,3),ia(3),ib(3),posz)
    nr(:) = [sum(blen(1:nb(1),1)),sum(blen(1:nb(2),2)),sum(blen(1:nb(3),3))]
    !
    ! file view: rows (x), planes (x-y) and the whole field (x-y-z), each resized to its extent in the file
    !
    if(iprec == 4) then
      etype = MPI_REAL
    else
      etype = MPI_DOUBLE_PRECISION
    end if
    esize = iprec
    call MPI_TYPE_INDEXED(nb(1),blen(1:nb(1),1),bstart(1:nb(1),1)-1,etype  ,type_x,ierr)
    call MPI_TYPE_CREATE_RESIZED(type_x,0_MPI_ADDRESS_KIND,esize*ng(1)      ,type_xr,ierr)
    call MPI_TYPE_INDEXED(nb(2),blen(1:nb(2),2),bstart(1:nb(2),2)-1,type_xr,type_y,ierr)
    call MPI_TYPE_CREATE_RESIZED(type_y,0_MPI_ADDRESS_KIND,esize*ng(1)*ng(2),type_yr,ierr)
    call MPI_TYPE_INDEXED(nb(3),blen(1:nb(3),3),bstart(1:nb(3),3)-1,type_yr,type_glob,ierr)
    call MPI_TYPE_COMMIT(type_glob,ierr)
    call MPI_TYPE_CONTIGUOUS(nr(1)*nr(2),etype,type_loc,ierr)
    call MPI_TYPE_COMMIT(type_loc,ierr)
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
    if(ierr /= MPI_SUCCESS) call abort_run('could not open checkpoint file `'//fname//'`.')
    call MPI_FILE_SET_VIEW(fh,disp,etype,type_glob,'native',MPI_INFO_NULL,ierr)
    if(iprec == 4) then
      allocate(buf_sp(nr(1),nr(2),nr(3)))
      call MPI_FILE_READ_ALL(fh,buf_sp,nr(3),type_loc,MPI_STATUS_IGNORE,ierr)
    else
      allocate(buf_dp(nr(1),nr(2),nr(3)))
      call MPI_FILE_READ_ALL(fh,buf_dp,nr(3),type_loc,MPI_STATUS_IGNORE,ierr)
    end if
    call MPI_FILE_CLOSE(fh,ierr)
    call MPI_TYPE_FREE(type_x   ,ierr)
    call MPI_TYPE_FREE(type_xr  ,ierr)
    call MPI_TYPE_FREE(type_y   ,ierr)
    call MPI_TYPE_FREE(type_yr  ,ierr)
    call MPI_TYPE_FREE(type_glob,ierr)
    call MPI_TYPE_FREE(type_loc ,ierr)
    !
    ! copy the points read to the block, in double precision
    !
    if(iprec == 4) then
      do k=ia(3),ib(3)
        kk = posz(k)
        do j=ia(2),ib(2)
          jj = posy(j)
          do i=ia(1),ib(1)
            p(i,j,k) = buf_sp(posx(i),jj,kk)
          end do
        end do
      end do
    else
      do k=ia(3),ib(3)
        kk = posz(k)
        do j=ia(2),ib(2)
          jj = posy(j)
          do i=ia(1),ib(1)
            p(i,j,k) = buf_dp(posx(i),jj,kk)
          end do
        end do
      end do
    end if
  end subroutine read_block
  !
  subroutine get_file_runs(n,is_periodic,lo,hi,nb,bstart,blen,ia,ib,pos)
    !
    ! maps the indices lo:hi along a direction with n points to the file indices 1:n (modulo n if periodic;
    ! otherwise, the indices outside 1:n are ghost points, not read): the file indices to read form nb runs
    ! (one, or two if the range wraps around), starting at bstart(1:nb) and of lengths blen(1:nb), in
    ! increasing order; ia:ib are the indices read, and pos(i) the position of the file index of i among
    ! the points read
    !
    implicit none
    integer, intent(in ) :: n
    logical, intent(in ) :: is_periodic
    integer, intent(in ) :: lo,hi
    integer, intent(out) :: nb
    integer, intent(out), dimension(2) :: bstart,blen
    integer, intent(out) :: ia,ib
    integer, intent(out), dimension(lo:) :: pos
    integer :: i,ig,ig0,nn
    !
    bstart(:) = 1
    blen(:)   = 0
    pos(:)    = 0
    if(.not.is_periodic) then
      ia = max(lo,1)
      ib = min(hi,n)
      nb = 1
      bstart(1) = ia
      blen(1)   = ib-ia+1
      do i=ia,ib
        pos(i) = i-ia+1
      end do
    else
      ia = lo
      ib = hi
      nn  = hi-lo+1
      ig0 = modulo(lo-1,n)+1
      if(nn >= n) then
        nb = 1
        blen(1) = n
      else if(ig0+nn-1 <= n) then
        nb = 1
        bstart(1) = ig0
        blen(1)   = nn
      else
        nb = 2
        bstart(:) = [1,ig0]
        blen(:)   = [ig0+nn-1-n,n-ig0+1]
      end if
      do i=lo,hi
        ig = modulo(i-1,n)+1
        if(ig >= bstart(nb)) then
          pos(i) = sum(blen(1:nb-1)) + ig-bstart(nb)+1
        else
          pos(i) = ig-bstart(1)+1
        end if
      end do
    end if
  end subroutine get_file_runs
  !
  subroutine write_field(fname,ng,lo,hi,iprec,p,time,istep)
    !
    ! writes a per-variable checkpoint file (the field, then the time and the time step) in single or double
    ! precision (iprec = 4 or 8 bytes); each task writes its block p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),
    ! converted to single precision by chunks of planes (`CHUNK_BYTES`) if needed
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: iprec
    real(rp), intent(in), contiguous, dimension(lo(1):,lo(2):,lo(3):) :: p
    real(rp), intent(in) :: time
    integer , intent(in) :: istep
    real(sp), allocatable, dimension(:,:,:) :: buf_sp
    real(sp), dimension(2) :: fldinfo_sp
    real(dp), dimension(2) :: fldinfo
    type(MPI_DATATYPE) :: etype,type_glob,type_loc
    type(MPI_FILE) :: fh
    integer(MPI_OFFSET_KIND) :: disp
    integer , dimension(3) :: n
    integer :: nreals_myid,nzc,nchunks,ic,k0,nk
    !
    n(:) = hi(:)-lo(:)+1
    if(iprec == 4) then
      etype = MPI_REAL
    else
      etype = MPI_DOUBLE_PRECISION
    end if
    call MPI_TYPE_CREATE_SUBARRAY(3,ng,n,lo-1,MPI_ORDER_FORTRAN,etype,type_glob,ierr)
    call MPI_TYPE_COMMIT(type_glob,ierr)
    call MPI_TYPE_CONTIGUOUS(n(1)*n(2),etype,type_loc,ierr)
    call MPI_TYPE_COMMIT(type_loc,ierr)
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD,fname,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    if(ierr /= MPI_SUCCESS) call abort_run('could not create output file `'//fname//'`.')
    call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,ierr)
    call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,etype,type_glob,'native',MPI_INFO_NULL,ierr)
    if(iprec == 4) then
      !
      ! convert and write in chunks of nzc planes, so that the single-precision copy stays small;
      ! the file pointer advances within the view, and all tasks make the same number of
      ! collective calls (with no data once their block is written)
      !
      nzc = int(max(1_i8,min(int(n(3),i8),CHUNK_BYTES/(4*int(n(1),i8)*n(2)))))
      nchunks = (n(3)+nzc-1)/nzc
      call MPI_ALLREDUCE(MPI_IN_PLACE,nchunks,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      allocate(buf_sp(n(1),n(2),nzc))
      do ic=1,nchunks
        k0 = lo(3) + (ic-1)*nzc
        nk = max(0,min(nzc,hi(3)-k0+1))
        if(nk > 0) buf_sp(:,:,1:nk) = real(p(:,:,k0:k0+nk-1),sp)
        call MPI_FILE_WRITE_ALL(fh,buf_sp,nk,type_loc,MPI_STATUS_IGNORE,ierr)
      end do
      deallocate(buf_sp)
    else
      call MPI_FILE_WRITE_ALL(fh,p     ,n(3),type_loc,MPI_STATUS_IGNORE,ierr)
    end if
    disp = product(int(ng(:),MPI_OFFSET_KIND))*iprec
    call MPI_FILE_SET_VIEW(fh,disp,etype,etype,'native',MPI_INFO_NULL,ierr)
    fldinfo(:) = [time,1._rp*istep]
    fldinfo_sp(:) = real(fldinfo(:),sp)
    nreals_myid = 0
    if(myid == 0) nreals_myid = 2
    if(iprec == 4) then
      call MPI_FILE_WRITE(fh,fldinfo_sp,nreals_myid,etype,MPI_STATUS_IGNORE,ierr)
    else
      call MPI_FILE_WRITE(fh,fldinfo   ,nreals_myid,etype,MPI_STATUS_IGNORE,ierr)
    end if
    call MPI_FILE_CLOSE(fh,ierr)
    call MPI_TYPE_FREE(type_glob,ierr)
    call MPI_TYPE_FREE(type_loc ,ierr)
  end subroutine write_field
end module mod_io
