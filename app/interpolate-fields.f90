  program interpolate_fields
    use mpi_f08
    use mod_common, only: rp,ierr
    use mod_bound , only: makehalo,updthalo,set_bc
    use mod_io    , only: load
    implicit none
    !
    ! input domain parameters
    !
    real(rp), parameter,       dimension(3) :: l     = [1.5_rp,3._rp,1._rp]
    integer , parameter,       dimension(3) :: ni    = [ 64, 64, 64]
    integer , parameter,       dimension(3) :: no    = [128,128,128]
    real(rp), parameter,       dimension(3) :: dlo   = l(:)/no(:)
    real(rp), parameter,       dimension(3) :: dli   = l(:)/ni(:)
    !
    ! boundary conditions
    !
    ! velocity
    character(len=1), parameter, dimension(0:1,3,3) :: cbcvel = &
      reshape(['P','P','P','P','P','P',  & ! u lower,upper bound in x,y,z
               'P','P','P','P','P','P',  & ! v lower,upper bound in x,y,z
               'P','P','P','P','P','P'], & ! w lower,upper bound in x,y,z
              shape(cbcvel))
    real(rp)        , parameter, dimension(0:1,3,3) ::  bcvel = &
        reshape([0._rp,0._rp,0._rp,0._rp,0._rp,0._rp,   &
                 0._rp,0._rp,0._rp,0._rp,0._rp,0._rp,   &
                 0._rp,0._rp,0._rp,0._rp,0._rp,0._rp],  &
                shape(bcvel))
    !
    ! pressure
    character(len=1), parameter, dimension(0:1,3) :: cbcpre = &
      reshape(['P','P','P','P','P','P'],shape(cbcpre))
    real(rp)        , parameter, dimension(0:1,3) ::  bcpre = &
      reshape([0._rp,0._rp,0._rp,0._rp,0._rp,0._rp],shape(bcpre))
    !
    ! default BC, in case it is needed for another field; commented for now:
    !!character(len=1), parameter, dimension(0:1,3) :: cbc = cbcpre
    !!real(rp)        , parameter, dimension(0:1,3) :: bc  =  bcpre
    !
    ! file names
    !
    character(len=*), parameter             :: input_file  = 'data/fld_i.bin', &
                                               output_file = 'data/fld_o.bin'
    !
    ! local problem sizes
    !
    integer, dimension(3)                   :: nni,nno,lo_i,hi_i,lo_o,hi_o
    !
    ! MPI stuff
    !
    integer                                 :: myid,nproc,dims(3),coords(3)
    type(MPI_DATATYPE)                      :: halo(3)
    type(MPI_COMM)                          :: comm_cart
    logical, dimension(3)                   :: periods
    integer, dimension(0:1,3)               :: nb
    logical, dimension(0:1,3)               :: is_bound
    !
    ! computational variables
    !
    real(rp), allocatable, dimension(:,:,:) :: ui,vi,wi,pi
    real(rp), allocatable, dimension(:,:,:) :: uo,vo,wo,po
    real(rp)                                :: time
    integer                                 :: istep
    !
    ! other variables
    !
    integer :: idir
    !
    ! initialize MPI
    !
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    !
    ! create processor grid
    !
    dims(:) = [0,0,1]
    call MPI_DIMS_CREATE(nproc,2,dims(1:2),ierr)
    !
    periods(:) = .false.; where(cbcpre(0,:)//cbcpre(1,:) == 'PP') periods(:) = .true.
    call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,.true.,comm_cart)
    call MPI_CART_COORDS(comm_cart,myid,3,coords)
    !
    ! decompose the domain
    !
    call distribute_grid(ni,dims,coords,[1,1,1],nni,lo_i,hi_i)
    call distribute_grid(no,dims,coords,[1,1,1],nno,lo_o,hi_o)
    !
    ! allocate input and output arrays
    !
    allocate(ui(0:nni(1)+1,0:nni(2)+1,0:nni(3)+1), &
             vi(0:nni(1)+1,0:nni(2)+1,0:nni(3)+1), &
             wi(0:nni(1)+1,0:nni(2)+1,0:nni(3)+1), &
             pi(0:nni(1)+1,0:nni(2)+1,0:nni(3)+1), &
             uo(0:nno(1)+1,0:nno(2)+1,0:nno(3)+1), &
             vo(0:nno(1)+1,0:nno(2)+1,0:nno(3)+1), &
             wo(0:nno(1)+1,0:nno(2)+1,0:nno(3)+1), &
             po(0:nno(1)+1,0:nno(2)+1,0:nno(3)+1))
    !
    ! determine neighbors
    !
    call MPI_CART_SHIFT(comm_cart,0,1,nb(0,1),nb(1,1),ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,nb(0,2),nb(1,2),ierr)
    nb(:,3) = MPI_PROC_NULL
    is_bound(:,:) = .false.
    where(nb(:,:) == MPI_PROC_NULL) is_bound(:,:) = .true.
    !
    ! generate halo datatypes
    !
    do idir=1,3
      call makehalo(idir,1,nni,halo(idir))
    end do
    !
    ! read input data
    !
    call load('r',input_file,MPI_COMM_WORLD,myid,ni,[1,1,1],lo_i,hi_i,ui,vi,wi,pi,time,istep)
    if(myid.eq.0) print*, 'Loaded field at time = ', time, 'step = ',istep,'.'
    !
    ! impose boundary conditions
    !
    do idir = 1,3
      call updthalo(1,halo(idir),nb(:,idir),idir,ui)
      call updthalo(1,halo(idir),nb(:,idir),idir,vi)
      call updthalo(1,halo(idir),nb(:,idir),idir,wi)
      call updthalo(1,halo(idir),nb(:,idir),idir,pi)
    end do
    !
    if(is_bound(0,1)) then
      call set_bc(cbcvel(0,1,1),0,1,1,.false.,bcvel(0,1,1),dli(1),ui)
      call set_bc(cbcvel(0,1,2),0,1,1,.true. ,bcvel(0,1,2),dli(1),vi)
      call set_bc(cbcvel(0,1,3),0,1,1,.true. ,bcvel(0,1,3),dli(1),wi)
      call set_bc(cbcpre(0,1  ),0,1,1,.true. ,bcpre(0,1  ),dli(1),pi)
    end if
    if(is_bound(1,1)) then
      call set_bc(cbcvel(1,1,1),1,1,1,.false.,bcvel(1,1,1),dli(1),ui)
      call set_bc(cbcvel(1,1,2),1,1,1,.true. ,bcvel(1,1,2),dli(1),vi)
      call set_bc(cbcvel(1,1,3),1,1,1,.true. ,bcvel(1,1,3),dli(1),wi)
      call set_bc(cbcpre(1,1  ),1,1,1,.true. ,bcpre(1,1  ),dli(1),pi)
    end if
    if(is_bound(0,2)) then
      call set_bc(cbcvel(0,2,1),0,2,1,.true. ,bcvel(0,2,1),dli(2),ui)
      call set_bc(cbcvel(0,2,2),0,2,1,.false.,bcvel(0,2,2),dli(2),vi)
      call set_bc(cbcvel(0,2,3),0,2,1,.true. ,bcvel(0,2,3),dli(2),wi)
      call set_bc(cbcpre(0,2  ),0,2,1,.true. ,bcpre(0,2  ),dli(2),pi)
     end if
    if(is_bound(1,2)) then
      call set_bc(cbcvel(1,2,1),1,2,1,.true. ,bcvel(1,2,1),dli(2),ui)
      call set_bc(cbcvel(1,2,2),1,2,1,.false.,bcvel(1,2,2),dli(2),vi)
      call set_bc(cbcvel(1,2,3),1,2,1,.true. ,bcvel(1,2,3),dli(2),wi)
      call set_bc(cbcpre(1,2  ),1,2,1,.true. ,bcpre(1,2  ),dli(2),pi)
    end if
    if(is_bound(0,3)) then
      call set_bc(cbcvel(0,3,1),0,3,1,.true. ,bcvel(0,3,1),dli(3),ui)
      call set_bc(cbcvel(0,3,2),0,3,1,.true. ,bcvel(0,3,2),dli(3),vi)
      call set_bc(cbcvel(0,3,3),0,3,1,.false.,bcvel(0,3,3),dli(3),wi)
      call set_bc(cbcpre(0,3  ),0,3,1,.true. ,bcpre(0,3  ),dli(3),pi)
    end if
    if(is_bound(1,3)) then
      call set_bc(cbcvel(1,3,1),1,3,1,.true. ,bcvel(1,3,1),dli(3),ui)
      call set_bc(cbcvel(1,3,2),1,3,1,.true. ,bcvel(1,3,2),dli(3),vi)
      call set_bc(cbcvel(1,3,3),1,3,1,.false.,bcvel(1,3,3),dli(3),wi)
      call set_bc(cbcpre(1,3  ),1,3,1,.true. ,bcpre(1,3  ),dli(3),pi)
    end if
    !
    ! interpolate field from grid 'i' to mesh 'o'
    !
    call interp_fld([.true. ,.false.,.false.],lo_i,lo_o,hi_o,dli,dlo,ui,uo)
    call interp_fld([.false.,.true. ,.false.],lo_i,lo_o,hi_o,dli,dlo,vi,vo)
    call interp_fld([.false.,.false.,.true. ],lo_i,lo_o,hi_o,dli,dlo,wi,wo)
    call interp_fld([.false.,.false.,.false.],lo_i,lo_o,hi_o,dli,dlo,pi,po)
    !
    call load('w',output_file,MPI_COMM_WORLD,myid,no,[1,1,1],lo_o,hi_o,uo,vo,wo,po,time,istep)
    call MPI_FINALIZE(ierr)
  contains
    subroutine interp_fld(is_staggered,lo_i,lo_o,hi_o,dli,dlo,fldi,fldo)
      implicit none
      logical , intent(in ), dimension(3) :: is_staggered
      integer , intent(in ), dimension(3) :: lo_i,lo_o,hi_o
      real(rp), intent(in ), dimension(3) :: dli,dlo
      real(rp), intent(in ), dimension(lo_i(1)-1:,lo_i(2)-1:,lo_i(3)-1:) :: fldi
      real(rp), intent(out), dimension(lo_o(1)-1:,lo_o(2)-1:,lo_o(3)-1:) :: fldo
      real(rp), dimension(3) :: ds
      real(rp) :: deltax,deltay,deltaz
      integer  :: i,j,k,ii,ji,ki,iip,jip,kip,iim,jim,kim
      real(rp) :: xo,yo,zo,xp,yp,zp,xm,ym,zm
      real(rp) :: f000,f001,f010,f100,f011,f101,f110,f111,val
      ds(:) = 0.5_rp
      where(.not.is_staggered(:)) ds(:) = 0._rp
      !
      do k=lo_o(3),hi_o(3)
        zo = (k-ds(3))*dlo(3)
        ki = nint(zo/dli(3)+ds(3))
        kip = ki + 1
        kim = ki - 1
        if(abs((kip-ds(3))*dli(3)-zo) <= abs((kim-ds(3))*dli(3)-zo)) then
          kip = kip
          kim = ki
        else
          kip = ki
          kim = kim
        endif
        zm = (kim-ds(3))*dli(3)
        zp = (kip-ds(3))*dli(3)
        deltaz = (zo-zm)/(zp-zm)
        do j=lo_o(2),hi_o(2)
          yo  = (j-ds(2))*dlo(2)
          ji  = nint(yo/dli(2)+ds(2))
          jip = ji + 1
          jim = ji - 1
          if(abs((jip-ds(2))*dli(2)-yo) <= abs((jim-ds(2))*dli(2)-yo)) then
            jip = jip
            jim = ji
          else
            jip = ji
            jim = jim
          endif
          ym = (jim-ds(2))*dli(2)
          yp = (jip-ds(2))*dli(2)
          deltay = (yo-ym)/(yp-ym)
          do i=lo_o(1),hi_o(1)
            xo  = (i-ds(1))*dlo(1)
            ii  = nint(xo/dli(1)+ds(1))
            iip = ii + 1
            iim = ii - 1
            if(abs((iip-ds(1))*dli(1)-xo) <= abs((iim-ds(1))*dli(1)-xo)) then
              iip = iip
              iim = ii
            else
              iip = ii
              iim = iim
            endif
            xm = (iim-ds(1))*dli(1)
            xp = (iip-ds(1))*dli(1)
            deltax = (xo-xm)/(xp-xm)
            !
            f000 = fldi(iim,jim,kim)
            f001 = fldi(iim,jim,kip)
            f010 = fldi(iim,jip,kim)
            f011 = fldi(iim,jip,kip)
            f100 = fldi(iip,jim,kim)
            f101 = fldi(iip,jim,kip)
            f110 = fldi(iip,jip,kim)
            f111 = fldi(iip,jip,kip)
            val = f000*(1.-deltax)*(1.-deltay)*(1.-deltaz) + &
                  f001*(1.-deltax)*(1.-deltay)*(   deltaz) + &
                  f010*(1.-deltax)*(   deltay)*(1.-deltaz) + &
                  f011*(1.-deltax)*(   deltay)*(   deltaz) + &
                  f100*(   deltax)*(1.-deltay)*(1.-deltaz) + &
                  f101*(   deltax)*(1.-deltay)*(   deltaz) + &
                  f110*(   deltax)*(   deltay)*(1.-deltaz) + &
                  f111*(   deltax)*(   deltay)*(   deltaz)
            fldo(i,j,k) = val
          end do
        end do
      end do
    end subroutine interp_fld
    !
    subroutine distribute_grid(ng,dims,coords,lo_g,n,lo,hi)
      implicit none
      integer, intent(in ), dimension(3) :: ng,dims,coords,lo_g
      integer, intent(out), dimension(3) :: n,lo,hi
      n(:) = ng(:)/dims(:)
      where(coords(:)+1 <= mod(ng(:),dims(:))) n(:) = n(:) + 1
      lo(:) = lo_g(:)   + (coords(:)  )*n(:)
      hi(:) = lo_g(:)-1 + (coords(:)+1)*n(:)
      where(coords(:)+1 >  mod(ng(:),dims(:)))
        lo(:) = lo(:) +    mod(ng(:),dims(:))
        hi(:) = hi(:) +    mod(ng(:),dims(:))
      end where
    end subroutine distribute_grid
  end program interpolate_fields
