module mod_param
  !
  ! input parameters, read from a namelist file (see the README)
  !
  use mod_common , only: rp,i8,myid,abort_run
  use mod_io_hdf5, only: is_hdf5_enabled,is_hdf5_file
  implicit none
  private
  public read_input,var_t,input_file_name,output_file_name,MAX_TILES
  public dims,input_format,input_prefix,input_file,output_format,output_prefix,output_precision, &
         is_compress_output,is_reset_time,ng_i,l_i,grid_files_i,ng_o,l_o,grid_files_o,is_periodic,vars, &
         ntiles,ntiles_tot,tile_inputs,tile_shifts
  integer, parameter :: MAX_VARS = 128, & ! maximum number of variables processed
                        MAX_SCAL = 999, & ! maximum number of scalars, `s_001` to `s_999`, as in CaNS
                        MAX_TILES = 256   ! maximum number of tiles
  integer, parameter :: LEN_NAME = 32, LEN_FILE = 256
  !
  ! a variable to process: name, location (staggered, i.e. on the cell faces, or cell-centred, per
  ! direction), boundary conditions, position in a legacy single-file checkpoint (u, v, w, p and
  ! the scalars, in this order; 0 for other variables), and field operations: removal of the
  ! volume-weighted mean or of a constant offset, and a constant added to the output (frame velocity)
  !
  type var_t
    character(len=LEN_NAME)            :: name
    logical         , dimension(3)     :: is_staggered
    character(len=1), dimension(0:1,3) :: cbc
    real(rp)        , dimension(0:1,3) :: bc
    integer                            :: ifield
    logical                            :: is_remove_mean = .false.
    real(rp)                           :: offset = 0._rp, added = 0._rp
  end type var_t
  !
  ! &interp
  !
  integer                , protected, dimension(3) :: dims
  character(len=16)      , protected :: input_format,output_format,output_precision
  character(len=LEN_FILE), protected :: input_prefix,input_file,output_prefix
  logical                , protected :: is_compress_output,is_reset_time
  !
  ! &grid_in and &grid_out
  !
  integer                , protected, dimension(3) :: ng_i,ng_o
  real(rp)               , protected, dimension(3) :: l_i,l_o
  character(len=LEN_FILE), protected, dimension(3) :: grid_files_i,grid_files_o
  !
  ! derived: periodic directions, and the list of variables to process
  !
  logical                , protected, dimension(3) :: is_periodic
  type(var_t), allocatable, protected, dimension(:) :: vars
  !
  ! &tiling: number of tiles per direction, and the input (prefix, or legacy file) and the shift (in input
  ! cells) of each tile, x fastest
  !
  integer                , protected, dimension(3) :: ntiles
  integer                , protected :: ntiles_tot
  character(len=LEN_FILE), protected, dimension(MAX_TILES) :: tile_inputs
  integer                , protected, dimension(3,MAX_TILES) :: tile_shifts
contains
  subroutine read_input(fname)
    implicit none
    character(len=*), intent(in) :: fname
    character(len=LEN_NAME), dimension(MAX_VARS) :: variables
    character(len=1), dimension(0:1,3,3)        :: cbcvel
    real(rp)        , dimension(0:1,3,3)        ::  bcvel
    character(len=1), dimension(0:1,3)          :: cbcpre
    real(rp)        , dimension(0:1,3)          ::  bcpre
    character(len=1), dimension(0:1,3,MAX_SCAL) :: cbcscal
    real(rp)        , dimension(0:1,3,MAX_SCAL) ::  bcscal
    type(var_t)     , dimension(MAX_VARS)       :: vars_nml
    logical         , dimension(MAX_VARS)       :: is_remove_mean
    real(rp)        , dimension(MAX_VARS)       :: offset
    real(rp)        , dimension(3)              :: vel_frame
    integer :: iunit,istat,iv,jv,nvars,nvars_nml,ivel,iscal
    character(len=256) :: errmsg
    character(len=LEN_NAME) :: name
    !
    open(newunit=iunit,file=fname,status='old',action='read',iostat=istat,iomsg=errmsg)
    if(istat /= 0) call abort_run('could not open the namelist file `'//fname//'`: '//trim(errmsg))
    call read_nml_interp(iunit,fname,variables)
    call read_nml_grid(iunit,fname,'grid_in' ,ng_i,l_i,grid_files_i)
    call read_nml_grid(iunit,fname,'grid_out',ng_o,l_o,grid_files_o)
    call read_nml_bc(iunit,fname,cbcvel,bcvel,cbcpre,bcpre,cbcscal,bcscal)
    call read_nml_variables(iunit,fname,nvars_nml,vars_nml)
    call read_nml_operations(iunit,fname,is_remove_mean,offset,vel_frame)
    call read_nml_tiling(iunit,fname)
    close(iunit)
    !
    ! &interp
    !
    if(any(dims(:) < 0)) call abort_run('`dims` must not be negative.')
    select case(trim(input_format))
    case('per-variable','hdf5')
      if(len_trim(input_prefix) == 0) &
        call abort_run('`input_prefix` must be prescribed for `input_format = '''//trim(input_format)//'''`.')
    case('legacy')
      if(len_trim(input_file) == 0) &
        call abort_run('`input_file` must be prescribed for `input_format = ''legacy''`.')
    case default
      call abort_run('`input_format` must be ''per-variable'', ''legacy'' or ''hdf5''.')
    end select
    select case(trim(output_format))
    case('per-variable','hdf5')
    case default
      call abort_run('`output_format` must be ''per-variable'' or ''hdf5''.')
    end select
    if(len_trim(output_prefix) == 0) call abort_run('`output_prefix` must be prescribed.')
    if(.not.is_hdf5_enabled .and. (trim(input_format) == 'hdf5' .or. trim(output_format) == 'hdf5' .or. &
                                   any([(is_hdf5_file(trim(grid_files_i(iv))),iv=1,3)]) .or. &
                                   any([(is_hdf5_file(trim(grid_files_o(iv))),iv=1,3)]))) &
      call abort_run('HDF5 input or output requested, but HDF5 support is not compiled in: '// &
                     'rebuild with `USE_HDF5=1`.')
    if(is_compress_output .and. trim(output_format) /= 'hdf5' .and. myid == 0) &
      print*, 'WARNING: `is_compress_output` is ignored for `output_format = ''per-variable''`.'
    select case(trim(output_precision))
    case('input','single','double')
    case default
      call abort_run('`output_precision` must be ''input'', ''single'' or ''double''.')
    end select
    !
    ! list of variables, with their locations and boundary conditions
    !
    nvars = count(variables(:) /= '')
    if(nvars == 0) then
      nvars = 4
      variables(1:4) = ['u','v','w','p']
    end if
    if(any(variables(1:nvars) == '')) call abort_run('`variables` must not have blank entries.')
    allocate(vars(nvars))
    do iv=1,nvars
      name = adjustl(variables(iv))
      if(any(variables(1:iv-1) == variables(iv))) call abort_run('variable `'//trim(name)//'` is listed twice.')
      vars(iv)%name = name
      vars(iv)%is_staggered(:) = .false.
      select case(trim(name))
      case('u','v','w')
        ivel = index('uvw',trim(name))
        vars(iv)%is_staggered(ivel) = .true.
        vars(iv)%cbc(:,:) = cbcvel(:,:,ivel)
        vars(iv)%bc(:,:)  =  bcvel(:,:,ivel)
        vars(iv)%ifield   = ivel
      case('p')
        vars(iv)%cbc(:,:) = cbcpre(:,:)
        vars(iv)%bc(:,:)  =  bcpre(:,:)
        vars(iv)%ifield   = 4
      case default
        if(is_scalar_name(name,iscal)) then
          vars(iv)%cbc(:,:) = cbcscal(:,:,iscal)
          vars(iv)%bc(:,:)  =  bcscal(:,:,iscal)
          vars(iv)%ifield   = 4+iscal
        else
          jv = findloc(vars_nml(1:nvars_nml)%name,name,dim=1)
          if(jv == 0) &
            call abort_run('variable `'//trim(name)//'` is not u, v, w, p or a scalar `s_<nnn>`: '// &
                           'its location and boundary conditions must be prescribed in a `&variable` group.')
          vars(iv) = vars_nml(jv)
          vars(iv)%ifield = 0
        end if
      end select
      if(any(vars(iv)%cbc(:,:) /= 'P' .and. vars(iv)%cbc(:,:) /= 'D' .and. vars(iv)%cbc(:,:) /= 'N')) &
        call abort_run('the boundary conditions of variable `'//trim(name)//'` are missing or invalid '// &
                       '(`cbcvel`, `cbcpre` or `cbcscal` in `&bc`, or `cbc` in `&variable`; '// &
                       '''P'', ''D'' or ''N'').')
      if(any((vars(iv)%cbc(0,:) == 'P') .neqv. (vars(iv)%cbc(1,:) == 'P'))) &
        call abort_run('variable `'//trim(name)//'`: a periodic boundary condition must be set on both sides.')
      if(iv == 1) is_periodic(:) = vars(iv)%cbc(0,:) == 'P'
      if(any((vars(iv)%cbc(0,:) == 'P') .neqv. is_periodic(:))) &
        call abort_run('variable `'//trim(name)//'`: the periodic directions must be the same for all variables.')
      if(trim(input_format) == 'legacy' .and. vars(iv)%ifield == 0) &
        call abort_run('variable `'//trim(name)//'` cannot be read from a legacy single-file checkpoint '// &
                       '(u, v, w, p and scalars only).')
    end do
    !
    ! field operations, per variable in the order of `variables`
    !
    if(any(is_remove_mean(nvars+1:)) .or. any(offset(nvars+1:) /= 0._rp)) &
      call abort_run('`is_remove_mean` or `offset` set beyond the number of `variables`.')
    do iv=1,nvars
      if(is_remove_mean(iv) .and. offset(iv) /= 0._rp) &
        call abort_run('variable `'//trim(vars(iv)%name)//'`: `is_remove_mean` and `offset` are exclusive.')
      vars(iv)%is_remove_mean = is_remove_mean(iv)
      vars(iv)%offset         = offset(iv)
      ivel = index('uvw',trim(vars(iv)%name))
      if(ivel > 0 .and. len_trim(vars(iv)%name) == 1) vars(iv)%added = vel_frame(ivel)
    end do
    do ivel=1,3
      if(vel_frame(ivel) /= 0._rp .and. .not.any(vars(:)%name == 'uvw'(ivel:ivel)) .and. myid == 0) &
        print*, 'WARNING: `vel_frame(', ivel, ')` ignored (', 'uvw'(ivel:ivel), ' not in `variables`).'
    end do
    do jv=1,nvars_nml
      if(myid == 0 .and. .not.any(vars(:)%name == vars_nml(jv)%name)) &
        print*, 'WARNING: `&variable` group of `', trim(vars_nml(jv)%name), '` ignored (not in `variables`).'
    end do
    !
    ! tiling, along periodic directions only; shifts along uniform directions only, as whole-cell shifts
    ! of a non-uniform grid would move the values to cells of other sizes
    !
    if(any(ntiles(:) < 1)) call abort_run('`ntiles` must be positive.')
    if(any(ntiles(:) > 1 .and. .not.is_periodic(:))) &
      call abort_run('tiling (`ntiles` > 1) is only possible along periodic directions.')
    if(product(int(ntiles(:),i8)) > MAX_TILES) call abort_run('too many tiles.')
    ntiles_tot = product(ntiles(:))
    if(any(tile_inputs(ntiles_tot+1:) /= '') .or. any(tile_shifts(:,ntiles_tot+1:) /= 0)) &
      call abort_run('`tile_inputs` or `tile_shifts` set beyond the number of tiles.')
    do iv=1,3
      if(any(tile_shifts(iv,1:ntiles_tot) /= 0) .and. .not.is_periodic(iv)) &
        call abort_run('tile shifts are only possible along periodic directions.')
      if(any(tile_shifts(iv,1:ntiles_tot) /= 0) .and. len_trim(grid_files_i(iv)) > 0) &
        call abort_run('tile shifts are only possible along directions with a uniform input grid.')
    end do
    where(tile_inputs(1:ntiles_tot) == '')
      tile_inputs(1:ntiles_tot) = merge(input_file,input_prefix,trim(input_format) == 'legacy')
    end where
    !
    ! the output must not overwrite the input
    !
    do iv=1,nvars
      do jv=1,ntiles_tot
        if(input_file_name(vars(iv)%name,jv) == output_file_name(vars(iv)%name)) &
          call abort_run('the output file of variable `'//trim(vars(iv)%name)//'` would overwrite its input.')
      end do
    end do
  end subroutine read_input
  !
  function input_file_name(name,itile) result(fname)
    !
    ! input file of variable `name` for tile `itile`: `<prefix>_<name>.bin` or `.h5`, or a legacy file,
    ! with the prefix or file of the tile (by default, `input_prefix` or `input_file`)
    !
    implicit none
    character(len=*), intent(in) :: name
    integer         , intent(in) :: itile
    character(len=:), allocatable :: fname
    select case(trim(input_format))
    case('legacy')
      fname = trim(tile_inputs(itile))
    case('hdf5')
      fname = trim(tile_inputs(itile))//'_'//trim(name)//'.h5'
    case default
      fname = trim(tile_inputs(itile))//'_'//trim(name)//'.bin'
    end select
  end function input_file_name
  !
  function output_file_name(name) result(fname)
    !
    ! output file of variable `name`: `<output_prefix>_<name>.bin` or `.h5`
    !
    implicit none
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: fname
    if(trim(output_format) == 'hdf5') then
      fname = trim(output_prefix)//'_'//trim(name)//'.h5'
    else
      fname = trim(output_prefix)//'_'//trim(name)//'.bin'
    end if
  end function output_file_name
  !
  subroutine read_nml_interp(iunit,fname,variables)
    implicit none
    integer         , intent(in ) :: iunit
    character(len=*), intent(in ) :: fname
    character(len=LEN_NAME), intent(out), dimension(MAX_VARS) :: variables
    integer :: istat
    character(len=256) :: errmsg
    namelist /interp/ dims,variables,input_format,input_prefix,input_file, &
                      output_format,output_prefix,output_precision,is_compress_output,is_reset_time
    !
    dims(:)          = [1,0,0]
    variables(:)     = ''
    input_format     = 'per-variable'
    input_prefix     = ''
    input_file       = ''
    output_format    = 'per-variable'
    output_prefix    = ''
    output_precision = 'input'
    is_compress_output = .false.
    is_reset_time    = .false.
    rewind(iunit)
    read(iunit,nml=interp,iostat=istat,iomsg=errmsg)
    if(istat /= 0) call abort_run('could not read namelist `&interp` from `'//fname//'`: '//trim(errmsg))
  end subroutine read_nml_interp
  !
  subroutine read_nml_grid(iunit,fname,group,ng_g,l_g,grid_files_g)
    !
    ! reads `&grid_in` or `&grid_out`
    !
    implicit none
    integer         , intent(in ) :: iunit
    character(len=*), intent(in ) :: fname,group
    integer         , intent(out), dimension(3) :: ng_g
    real(rp)        , intent(out), dimension(3) :: l_g
    character(len=*), intent(out), dimension(3) :: grid_files_g
    integer :: istat
    character(len=256) :: errmsg
    integer                 , dimension(3) :: ng
    real(rp)                , dimension(3) :: l
    character(len=LEN_FILE) , dimension(3) :: grid_files
    namelist /grid_in/  ng,l,grid_files
    namelist /grid_out/ ng,l,grid_files
    !
    ng(:) = 0
    l(:)  = 0._rp
    grid_files(:) = ''
    rewind(iunit)
    if(group == 'grid_in') then
      read(iunit,nml=grid_in ,iostat=istat,iomsg=errmsg)
    else
      read(iunit,nml=grid_out,iostat=istat,iomsg=errmsg)
    end if
    if(istat /= 0) call abort_run('could not read namelist `&'//group//'` from `'//fname//'`: '//trim(errmsg))
    ng_g(:) = ng(:)
    l_g(:)  = l(:)
    grid_files_g(:) = grid_files(:)
  end subroutine read_nml_grid
  !
  subroutine read_nml_bc(iunit,fname,cbcvel,bcvel,cbcpre,bcpre,cbcscal,bcscal)
    !
    ! reads `&bc` (optional, if no variable needs it): boundary conditions of u, v, w, p and the scalars, as in CaNS
    !
    implicit none
    integer         , intent(in ) :: iunit
    character(len=*), intent(in ) :: fname
    character(len=1), intent(out), dimension(0:1,3,3) :: cbcvel
    real(rp)        , intent(out), dimension(0:1,3,3) ::  bcvel
    character(len=1), intent(out), dimension(0:1,3)   :: cbcpre
    real(rp)        , intent(out), dimension(0:1,3)   ::  bcpre
    character(len=1), intent(out), dimension(0:1,3,MAX_SCAL) :: cbcscal
    real(rp)        , intent(out), dimension(0:1,3,MAX_SCAL) ::  bcscal
    integer :: istat
    character(len=256) :: errmsg
    namelist /bc/ cbcvel,bcvel,cbcpre,bcpre,cbcscal,bcscal
    !
    cbcvel(:,:,:)  = ' '
    bcvel(:,:,:)   = 0._rp
    cbcpre(:,:)    = ' '
    bcpre(:,:)     = 0._rp
    cbcscal(:,:,:) = ' '
    bcscal(:,:,:)  = 0._rp
    rewind(iunit)
    read(iunit,nml=bc,iostat=istat,iomsg=errmsg)
    if(istat > 0) call abort_run('could not read namelist `&bc` from `'//fname//'`: '//trim(errmsg))
  end subroutine read_nml_bc
  !
  subroutine read_nml_variables(iunit,fname,nvars_nml,vars_nml)
    !
    ! reads the `&variable` groups (optional, one per variable other than u, v, w, p and the scalars)
    !
    implicit none
    integer         , intent(in ) :: iunit
    character(len=*), intent(in ) :: fname
    integer         , intent(out) :: nvars_nml
    type(var_t)     , intent(out), dimension(:) :: vars_nml
    integer :: istat,iscal
    character(len=256) :: errmsg
    character(len=LEN_NAME) :: name
    logical         , dimension(3)     :: is_staggered
    character(len=1), dimension(0:1,3) :: cbc
    real(rp)        , dimension(0:1,3) :: bc
    namelist /variable/ name,is_staggered,cbc,bc
    !
    nvars_nml = 0
    rewind(iunit)
    do
      name = ''
      is_staggered(:) = .false.
      cbc(:,:) = ' '
      bc(:,:)  = 0._rp
      read(iunit,nml=variable,iostat=istat,iomsg=errmsg)
      if(istat < 0) exit
      if(istat > 0) call abort_run('could not read namelist `&variable` from `'//fname//'`: '//trim(errmsg))
      name = adjustl(name)
      if(name == '') call abort_run('a `&variable` group has no `name`.')
      if(any(name == ['u','v','w','p']) .or. is_scalar_name(name,iscal)) &
        call abort_run('`&variable` group of `'//trim(name)//'`: u, v, w, p and the scalars take their '// &
                       'locations from their names and their boundary conditions from `&bc`.')
      if(any(vars_nml(1:nvars_nml)%name == name)) &
        call abort_run('`&variable` group of `'//trim(name)//'` given twice.')
      if(nvars_nml == size(vars_nml)) call abort_run('too many `&variable` groups.')
      nvars_nml = nvars_nml + 1
      vars_nml(nvars_nml)%name            = name
      vars_nml(nvars_nml)%is_staggered(:) = is_staggered(:)
      vars_nml(nvars_nml)%cbc(:,:)        = cbc(:,:)
      vars_nml(nvars_nml)%bc(:,:)         = bc(:,:)
      vars_nml(nvars_nml)%ifield          = 0
    end do
  end subroutine read_nml_variables
  !
  subroutine read_nml_operations(iunit,fname,is_remove_mean,offset,vel_frame)
    !
    ! reads `&operations` (optional): field operations, per variable in the order of `variables`
    !
    implicit none
    integer         , intent(in ) :: iunit
    character(len=*), intent(in ) :: fname
    logical         , intent(out), dimension(MAX_VARS) :: is_remove_mean
    real(rp)        , intent(out), dimension(MAX_VARS) :: offset
    real(rp)        , intent(out), dimension(3)        :: vel_frame
    integer :: istat
    character(len=256) :: errmsg
    namelist /operations/ is_remove_mean,offset,vel_frame
    !
    is_remove_mean(:) = .false.
    offset(:)         = 0._rp
    vel_frame(:)      = 0._rp
    rewind(iunit)
    read(iunit,nml=operations,iostat=istat,iomsg=errmsg)
    if(istat > 0) call abort_run('could not read namelist `&operations` from `'//fname//'`: '//trim(errmsg))
  end subroutine read_nml_operations
  !
  subroutine read_nml_tiling(iunit,fname)
    !
    ! reads `&tiling` (optional)
    !
    implicit none
    integer         , intent(in) :: iunit
    character(len=*), intent(in) :: fname
    integer :: istat
    character(len=256) :: errmsg
    namelist /tiling/ ntiles,tile_inputs,tile_shifts
    !
    ntiles(:)        = 1
    tile_inputs(:)   = ''
    tile_shifts(:,:) = 0
    rewind(iunit)
    read(iunit,nml=tiling,iostat=istat,iomsg=errmsg)
    if(istat > 0) call abort_run('could not read namelist `&tiling` from `'//fname//'`: '//trim(errmsg))
  end subroutine read_nml_tiling
  !
  logical function is_scalar_name(name,iscal)
    !
    ! whether name is that of a scalar, `s_<nnn>` as in CaNS, with number iscal
    !
    implicit none
    character(len=*), intent(in ) :: name
    integer         , intent(out) :: iscal
    iscal = 0
    is_scalar_name = len_trim(name) == 5 .and. name(1:2) == 's_' .and. verify(name(3:5),'0123456789') == 0
    if(is_scalar_name) read(name(3:5),'(i3)') iscal
    is_scalar_name = is_scalar_name .and. iscal >= 1
  end function is_scalar_name
end module mod_param
