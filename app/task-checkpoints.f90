  program task_checkpoints
    !
    ! converts between a checkpoint saved with one file per task and a single checkpoint (see the README):
    !   mode = 'merge': per-task files -> one checkpoint (per-variable `.bin` files, a legacy single file, or HDF5)
    !   mode = 'split': one checkpoint -> per-task files
    ! a per-task file holds the variables `task_variables`, one after the other, over the block of its task
    ! (without halos) in a task_dims(1) x task_dims(2) x task_dims(3) decomposition of the grid (a pencil
    ! decomposition has a 1 in task_dims), followed by the time and the time step, in single or double
    ! precision (detected from the file size); the tasks of this program process the files of ranks
    ! myid, myid+nproc, ..., one round at a time, and need not match task_dims
    !
    use mpi_f08
    use mod_common , only: rp,sp,dp,i8,ierr,myid,abort_run
    use mod_decomp , only: distribute_grid,task_coords
    use mod_io     , only: get_checkpoint_layout,read_fldinfo,read_block,out_file,open_output,write_output_block, &
                           close_output
    use mod_io_hdf5, only: is_hdf5_enabled
    implicit none
    integer, parameter :: MAX_VARS = 128, LEN_NAME = 32, LEN_FILE = 256
    !
    ! input parameters, namelist `&tasks`
    !
    character(len=8)       :: mode = ''
    integer, dimension(3)  :: ng = 0,task_dims = 0
    character(len=16)      :: task_order = 'row-major'
    character(len=LEN_FILE):: task_prefix = '',task_suffix = '.bin'
    integer                :: task_digits = 5
    character(len=LEN_NAME), dimension(MAX_VARS) :: task_variables = '',variables = ''
    character(len=16)      :: task_precision = 'input',checkpoint_format = 'per-variable',checkpoint_precision = 'input'
    character(len=LEN_FILE):: checkpoint_prefix = '',checkpoint_file = ''
    real(rp), dimension(MAX_VARS) :: offset = 0._rp
    logical                :: is_compress = .false.,is_reset_time = .false.
    namelist /tasks/ mode,ng,task_dims,task_order,task_prefix,task_suffix,task_digits,task_variables, &
                     task_precision,checkpoint_format,checkpoint_prefix,checkpoint_file,checkpoint_precision, &
                     variables,offset,is_compress,is_reset_time
    !
    character(len=LEN_FILE) :: fname_nml
    character(len=256) :: errmsg
    integer :: nproc,nvars,nf,nrounds,iunit,istat,iv
    integer(i8) :: ntot
    !
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    !
    ! read and check the input parameters, from the namelist file given as argument (default: `tasks.nml`)
    !
    fname_nml = 'tasks.nml'
    if(command_argument_count() > 0) call get_command_argument(1,fname_nml)
    open(newunit=iunit,file=trim(fname_nml),status='old',action='read',iostat=istat,iomsg=errmsg)
    if(istat /= 0) call abort_run('could not open the namelist file `'//trim(fname_nml)//'`: '//trim(errmsg))
    read(iunit,nml=tasks,iostat=istat,iomsg=errmsg)
    if(istat /= 0) call abort_run('could not read namelist `&tasks` from `'//trim(fname_nml)//'`: '//trim(errmsg))
    close(iunit)
    if(trim(mode) /= 'merge' .and. trim(mode) /= 'split') call abort_run('`mode` must be ''merge'' or ''split''.')
    if(any(ng(:) < 1)) call abort_run('`ng` must be prescribed and positive.')
    if(any(task_dims(:) < 1) .or. any(task_dims(:) > ng(:))) &
      call abort_run('`task_dims` must be prescribed, positive and at most `ng`.')
    if(trim(task_order) /= 'row-major' .and. trim(task_order) /= 'column-major') &
      call abort_run('`task_order` must be ''row-major'' or ''column-major''.')
    if(len_trim(task_prefix) == 0) call abort_run('`task_prefix` must be prescribed.')
    if(task_digits < 1 .or. task_digits > 9) call abort_run('`task_digits` must be between 1 and 9.')
    nvars = count(task_variables(:) /= '')
    if(nvars == 0) then
      nvars = 4
      task_variables(1:4) = ['u','v','w','p']
    end if
    if(any(task_variables(1:nvars) == '')) call abort_run('`task_variables` must not have blank entries.')
    do iv=1,nvars
      if(any(task_variables(1:iv-1) == task_variables(iv))) call abort_run('`task_variables` has repeated names.')
    end do
    if(any(offset(nvars+1:) /= 0._rp)) call abort_run('`offset` set beyond the number of `task_variables`.')
    select case(trim(checkpoint_format))
    case('per-variable','hdf5')
      if(len_trim(checkpoint_prefix) == 0) call abort_run('`checkpoint_prefix` must be prescribed.')
    case('legacy')
      if(len_trim(checkpoint_file) == 0) call abort_run('`checkpoint_file` must be prescribed.')
    case default
      call abort_run('`checkpoint_format` must be ''per-variable'', ''legacy'' or ''hdf5''.')
    end select
    if(trim(checkpoint_format) == 'hdf5' .and. .not.is_hdf5_enabled) &
      call abort_run('HDF5 checkpoint requested, but HDF5 support is not compiled in: rebuild with `USE_HDF5=1`.')
    if(all(trim(task_precision)       /= ['input ','single','double']) .or. &
       all(trim(checkpoint_precision) /= ['input ','single','double'])) &
      call abort_run('`task_precision` and `checkpoint_precision` must be ''input'', ''single'' or ''double''.')
    nf = product(task_dims(:))
    nrounds = (nf+nproc-1)/nproc
    ntot = product(int(ng(:),i8))
    if(myid == 0) then
      print*, '*** task-checkpoints: ', trim(mode), ' ***'
      print*, 'Grid: ', ng(:), '| per-task files: ', task_dims(:), '(', nf, ' files, ', trim(task_order), &
              ' rank order), processed in ', nrounds, ' rounds by ', nproc, ' tasks'
      print*, 'Variables of the per-task files: ', (trim(task_variables(iv))//' ',iv=1,nvars)
    end if
    select case(trim(mode))
    case('merge')
      call merge_files()
    case('split')
      call split_files()
    end select
    call MPI_FINALIZE(ierr)
  contains
    function task_file(irank) result(fname)
      !
      ! per-task file of rank irank: <task_prefix><irank, with at least task_digits digits><task_suffix>
      !
      implicit none
      integer, intent(in) :: irank
      character(len=:), allocatable :: fname
      character(len=16) :: crank,cfmt
      write(cfmt,'(a,i0,a)') '(i0.',task_digits,')'
      write(crank,cfmt) irank
      fname = trim(task_prefix)//trim(crank)//trim(task_suffix)
    end function task_file
    !
    function checkpoint_name(name) result(fname)
      implicit none
      character(len=*), intent(in) :: name
      character(len=:), allocatable :: fname
      select case(trim(checkpoint_format))
      case('legacy')
        fname = trim(checkpoint_file)
      case('hdf5')
        fname = trim(checkpoint_prefix)//'_'//trim(name)//'.h5'
      case default
        fname = trim(checkpoint_prefix)//'_'//trim(name)//'.bin'
      end select
    end function checkpoint_name
    !
    subroutine task_block(irank,lo,hi,nloc)
      !
      ! block lo:hi (nloc points) of the per-task file of rank irank; an empty block (lo > hi) if irank >= nf
      !
      implicit none
      integer    , intent(in ) :: irank
      integer    , intent(out), dimension(3) :: lo,hi
      integer(i8), intent(out) :: nloc
      integer, dimension(3) :: coords,n
      if(irank >= nf) then
        lo(:) = 1
        hi(:) = 0
        nloc  = 0
        return
      end if
      call task_coords(irank,task_dims,trim(task_order),coords)
      call distribute_grid(ng,task_dims,coords,n,lo,hi)
      nloc = product(int(n(:),i8))
    end subroutine task_block
    !
    subroutine stop_on_errors(is_error,msg)
      !
      ! stops the run if a task found an error (in the files it processes), which it prints
      !
      implicit none
      logical         , intent(in) :: is_error
      character(len=*), intent(in) :: msg
      logical :: is_error_any
      call MPI_ALLREDUCE(is_error,is_error_any,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
      if(is_error) print*, 'ERROR (task ', myid, '): ', trim(msg)
      if(is_error_any) call abort_run('errors in the per-task files (see above).')
    end subroutine stop_on_errors
    !
    subroutine merge_files()
      !
      ! per-task files -> one checkpoint, with the variables `variables` (default: all `task_variables`); with
      ! a subset, only these are written (into an existing, complete file for a legacy checkpoint), e.g. to
      ! merge one field again
      !
      implicit none
      integer , dimension(3) :: lo,hi
      integer , dimension(nrounds) :: iprec_t
      integer(i8) :: nloc,fsize,ipos
      real(sp), allocatable, dimension(:,:,:) :: buf_sp
      real(rp), allocatable, dimension(:,:,:) :: p
      real(sp), dimension(2) :: info_sp
      real(dp), dimension(2) :: info,info0,tmin,tmax
      integer , dimension(MAX_VARS) :: iwrite
      integer :: k,irank,iprec_min,iprec_max,iprec_o,nwrite,jv,iv,istep_o
      integer(MPI_OFFSET_KIND) :: disp,disp_info
      real(rp) :: time_o
      type(out_file) :: f
      logical :: is_error,is_subset,is_existing
      character(len=LEN_FILE) :: fname
      character(len=512) :: msg
      !
      ! variables to write
      !
      nwrite = count(variables(:) /= '')
      if(nwrite == 0) then
        nwrite = nvars
        variables(1:nvars) = task_variables(1:nvars)
      end if
      do jv=1,nwrite
        iwrite(jv) = findloc(task_variables(1:nvars),variables(jv),dim=1)
        if(iwrite(jv) == 0) call abort_run('variable `'//trim(variables(jv))//'` is not in `task_variables`.')
      end do
      is_subset = nwrite < nvars
      !
      ! the per-task files of this task: precision (from the size), time and time step
      !
      is_error = .false.
      msg = ''
      iprec_t(:) = 0
      tmin(:) =  huge(1._dp)
      tmax(:) = -huge(1._dp)
      info0(:) = 0._dp
      do k=1,nrounds
        irank = myid + (k-1)*nproc
        if(irank >= nf) exit
        call task_block(irank,lo,hi,nloc)
        fname = task_file(irank)
        fsize = -1
        inquire(file=trim(fname),size=fsize)
        if(fsize < 0) then
          is_error = .true.
          msg = 'could not open `'//trim(fname)//'`.'
          exit
        else if(fsize == (nvars*nloc+2)*4) then
          iprec_t(k) = 4
        else if(fsize == (nvars*nloc+2)*8) then
          iprec_t(k) = 8
        else
          is_error = .true.
          write(msg,'(a,i0,a,i0,a,i0,a)') 'file `'//trim(fname)//'` (', fsize, ' bytes) is not that of rank ', &
                                          irank, ' with ', nvars, ' variables, in single or double precision.'
          exit
        end if
        open(newunit=iunit,file=trim(fname),access='stream',form='unformatted',status='old',action='read')
        ipos = 1 + nvars*nloc*iprec_t(k)
        if(iprec_t(k) == 4) then
          read(iunit,pos=ipos) info_sp(:)
          info(:) = info_sp(:)
        else
          read(iunit,pos=ipos) info(:)
        end if
        close(iunit)
        tmin(:) = min(tmin(:),info(:))
        tmax(:) = max(tmax(:),info(:))
        if(irank == 0) info0(:) = info(:)
      end do
      call stop_on_errors(is_error,msg)
      call MPI_ALLREDUCE(MPI_IN_PLACE,tmin,2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,tmax,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(info0,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(minval(iprec_t,mask=iprec_t > 0),iprec_min,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(maxval(iprec_t,mask=iprec_t > 0),iprec_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
      !
      ! the time step must be the same in all files; the time may not (it is read per task at each restart
      ! of a run saving one file per task), so only a warning is issued, and that of rank 0 is kept
      !
      if(tmin(2) /= tmax(2)) call abort_run('inconsistent time step across the per-task files.')
      if(tmin(1) /= tmax(1) .and. myid == 0) &
        print*, 'WARNING: inconsistent time across the per-task files (min, max): ', tmin(1), tmax(1), &
                '; the time of rank 0 is kept.'
      select case(trim(checkpoint_precision))
      case('single')
        iprec_o = 4
      case('double')
        iprec_o = 8
      case default
        if(iprec_min /= iprec_max) &
          call abort_run('per-task files in single and double precision: set `checkpoint_precision`.')
        iprec_o = iprec_min
      end select
      time_o  = info0(1)
      istep_o = nint(info0(2))
      if(is_reset_time) then
        time_o  = 0._rp
        istep_o = 0
      end if
      if(myid == 0) print*, 'Per-task files in ', merge('single','double',iprec_min == 4), &
                            trim(merge(' and double','           ',iprec_min /= iprec_max)), ' precision, at time = ', &
                            info0(1), 'step = ', nint(info0(2))
      if(trim(checkpoint_format) == 'legacy' .and. is_subset) then
        fsize = -1
        if(myid == 0) inquire(file=trim(checkpoint_file),size=fsize)
        call MPI_BCAST(fsize,1,MPI_INTEGER8,0,MPI_COMM_WORLD,ierr)
        if(fsize /= (nvars*ntot+2)*iprec_o) &
          call abort_run('merging some variables into a legacy checkpoint requires an existing, complete file.')
      end if
      !
      ! write the variables, one at a time, each task writing the block of one file per round
      !
      do jv=1,nwrite
        iv = iwrite(jv)
        fname = checkpoint_name(task_variables(iv))
        if(trim(checkpoint_format) == 'legacy') then
          disp        = (iv-1)*ntot*iprec_o
          disp_info   = nvars*ntot*iprec_o
          is_existing = is_subset .or. jv > 1
          call open_output(trim(fname),trim(checkpoint_format),trim(task_variables(iv)),ng,iprec_o,f, &
                           disp=disp,disp_info=disp_info,is_existing=is_existing)
        else
          call open_output(trim(fname),trim(checkpoint_format),trim(task_variables(iv)),ng,iprec_o,f, &
                           is_compress=is_compress,nchunk=ng/task_dims)
        end if
        do k=1,nrounds
          irank = myid + (k-1)*nproc
          call task_block(irank,lo,hi,nloc)
          allocate(p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
          if(irank < nf) then
            open(newunit=iunit,file=task_file(irank),access='stream',form='unformatted',status='old',action='read')
            ipos = 1 + (iv-1)*nloc*iprec_t(k)
            if(iprec_t(k) == 4) then
              allocate(buf_sp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
              read(iunit,pos=ipos) buf_sp
              p(:,:,:) = buf_sp(:,:,:)
              deallocate(buf_sp)
            else
              read(iunit,pos=ipos) p
            end if
            close(iunit)
            if(offset(iv) /= 0._rp) p(:,:,:) = p(:,:,:) - offset(iv)
          end if
          call write_output_block(f,ng,lo,hi,p)
          deallocate(p)
        end do
        call close_output(f,time_o,istep_o)
        if(myid == 0) print*, '  ', trim(task_variables(iv)), ' written to `', trim(fname), '` (', &
                              merge('single','double',iprec_o == 4), ' precision)', &
                              trim(merge(' with offset subtracted','                       ',offset(iv) /= 0._rp))
      end do
      if(myid == 0) print*, 'Merged at time = ', time_o, 'step = ', istep_o
    end subroutine merge_files
    !
    subroutine split_files()
      !
      ! one checkpoint -> per-task files, with all `task_variables`
      !
      implicit none
      character(len=LEN_FILE), dimension(MAX_VARS) :: fnames
      integer(MPI_OFFSET_KIND), dimension(MAX_VARS) :: disps
      integer, dimension(MAX_VARS) :: iprec_c
      integer(MPI_OFFSET_KIND) :: disp_info
      integer , dimension(3) :: lo,hi
      integer(i8) :: nloc
      real(rp), allocatable, dimension(:,:,:) :: p
      real(rp) :: time
      integer  :: istep,k,irank,nfields,iprec_o,iunit_t
      integer , dimension(3,1) :: no_shift
      logical :: is_error
      character(len=512) :: msg
      !
      ! the checkpoint: precision of each variable, and time and time step (of the first variable)
      !
      do iv=nvars,1,-1 ! (backwards, so that disp_info ends as that of the first variable)
        fnames(iv) = checkpoint_name(task_variables(iv))
        call get_checkpoint_layout(trim(fnames(iv)),trim(checkpoint_format),trim(task_variables(iv)),ng, &
                                   iprec_c(iv),nfields)
        if(trim(checkpoint_format) == 'legacy') then
          if(nfields < nvars) call abort_run('the legacy checkpoint has fewer fields than `task_variables`.')
          disps(iv) = (iv-1)*ntot*iprec_c(iv)
          disp_info = nfields*ntot*iprec_c(iv)
        else
          disps(iv) = 0
          disp_info = ntot*iprec_c(iv)
        end if
      end do
      call read_fldinfo(trim(fnames(1)),trim(checkpoint_format),disp_info,iprec_c(1),time,istep)
      select case(trim(task_precision))
      case('single')
        iprec_o = 4
      case('double')
        iprec_o = 8
      case default
        if(any(iprec_c(1:nvars) /= iprec_c(1))) &
          call abort_run('checkpoint files in single and double precision: set `task_precision`.')
        iprec_o = iprec_c(1)
      end select
      if(is_reset_time) then
        time  = 0._rp
        istep = 0
      end if
      if(myid == 0) print*, 'Checkpoint `', trim(fnames(1)), '`, ... at time = ', time, 'step = ', istep, &
                            '; per-task files in ', merge('single','double',iprec_o == 4), ' precision'
      !
      ! write the per-task files, one per task and round, reading their blocks from the checkpoint
      !
      no_shift(:,:) = 0
      do k=1,nrounds
        irank = myid + (k-1)*nproc
        call task_block(irank,lo,hi,nloc)
        is_error = .false.
        msg = ''
        if(irank < nf) then
          open(newunit=iunit_t,file=task_file(irank),access='stream',form='unformatted',status='replace', &
               action='write',iostat=istat,iomsg=errmsg)
          is_error = istat /= 0
          if(is_error) msg = 'could not create `'//task_file(irank)//'`: '//trim(errmsg)
        end if
        call stop_on_errors(is_error,msg)
        do iv=1,nvars
          allocate(p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
          call read_block(fnames(iv:iv),trim(checkpoint_format),trim(task_variables(iv)),disps(iv:iv),ng, &
                          [1,1,1],no_shift,iprec_c(iv),[.false.,.false.,.false.],lo,hi,p)
          if(offset(iv) /= 0._rp) p(:,:,:) = p(:,:,:) - offset(iv)
          if(irank < nf) then
            if(iprec_o == 4) then
              write(iunit_t) real(p,sp)
            else
              write(iunit_t) real(p,dp)
            end if
          end if
          deallocate(p)
        end do
        if(irank < nf) then
          if(iprec_o == 4) then
            write(iunit_t) real([time,1._rp*istep],sp)
          else
            write(iunit_t) real([time,1._rp*istep],dp)
          end if
          close(iunit_t)
        end if
        if(myid == 0) print*, '  round ', k, ' of ', nrounds, ' written'
      end do
      if(myid == 0) print*, 'Split into ', nf, ' files `', task_file(0), '`, ...'
    end subroutine split_files
  end program task_checkpoints
