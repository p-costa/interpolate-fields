module mod_bound
  use mod_common, only: rp
  use mod_grid  , only: grid_1d
  implicit none
  private
  public set_ghosts
contains
  subroutine set_ghosts(cbc,bc,is_staggered,g,lo,hi,p)
    !
    ! imposes the boundary conditions along the non-periodic directions on a block p(lo(1):hi(1),...) of a
    ! field, indexed as the global grid, as in CaNS (`set_bc`): a cell-centred variable gets the ghost
    ! points i = 0 and i = n+1; a staggered variable gets the points on the boundary faces, i = 0 and
    ! i = n (where Dirichlet conditions prescribe the value, and Neumann conditions a one-sided difference)
    !
    ! the directions are processed in order, each over the whole block, so that the edges and corners
    ! are set consistently by the last direction
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    real(rp)        , intent(in), dimension(0:1,3) :: bc
    logical         , intent(in), dimension(3)     :: is_staggered
    type(grid_1d)   , intent(in), dimension(3)     :: g
    integer         , intent(in), dimension(3)     :: lo,hi
    real(rp), intent(inout), dimension(lo(1):,lo(2):,lo(3):) :: p
    integer :: idir,n
    !
    do idir=1,3
      if(g(idir)%is_periodic) cycle
      n = g(idir)%n
      if(.not.is_staggered(idir)) then
        if(lo(idir) <= 0) then
          select case(cbc(0,idir))
          case('D')
            call set_plane(idir,lo,0  ,1  ,2*bc(0,idir)               ,-1._rp,p)
          case('N')
            call set_plane(idir,lo,0  ,1  ,-g(idir)%dxc(0)*bc(0,idir) , 1._rp,p)
          end select
        end if
        if(hi(idir) >= n+1) then
          select case(cbc(1,idir))
          case('D')
            call set_plane(idir,lo,n+1,n  ,2*bc(1,idir)               ,-1._rp,p)
          case('N')
            call set_plane(idir,lo,n+1,n  , g(idir)%dxc(n)*bc(1,idir) , 1._rp,p)
          end select
        end if
      else
        if(lo(idir) <= 0) then
          select case(cbc(0,idir))
          case('D')
            call set_plane(idir,lo,0  ,0  ,bc(0,idir)                 , 0._rp,p)
          case('N')
            call set_plane(idir,lo,0  ,1  ,-g(idir)%dxf(0)*bc(0,idir) , 1._rp,p)
          end select
        end if
        if(hi(idir) >= n) then
          select case(cbc(1,idir))
          case('D')
            call set_plane(idir,lo,n  ,n  ,bc(1,idir)                 , 0._rp,p)
          case('N')
            call set_plane(idir,lo,n  ,n-1, g(idir)%dxf(n)*bc(1,idir) , 1._rp,p)
          end select
        end if
      end if
    end do
  end subroutine set_ghosts
  !
  subroutine set_plane(idir,lo,ig,isrc,factor,sgn,p)
    !
    ! sets plane ig along direction idir: p(ig) = factor + sgn*p(isrc), or p(ig) = factor if sgn = 0
    !
    implicit none
    integer , intent(in) :: idir
    integer , intent(in), dimension(3) :: lo
    integer , intent(in) :: ig,isrc
    real(rp), intent(in) :: factor,sgn
    real(rp), intent(inout), dimension(lo(1):,lo(2):,lo(3):) :: p
    if(sgn == 0._rp) then
      select case(idir)
      case(1)
        p(ig,:,:) = factor
      case(2)
        p(:,ig,:) = factor
      case(3)
        p(:,:,ig) = factor
      end select
    else
      select case(idir)
      case(1)
        p(ig,:,:) = factor + sgn*p(isrc,:,:)
      case(2)
        p(:,ig,:) = factor + sgn*p(:,isrc,:)
      case(3)
        p(:,:,ig) = factor + sgn*p(:,:,isrc)
      end select
    end if
  end subroutine set_plane
end module mod_bound
