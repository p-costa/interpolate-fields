module mod_interp
  use mod_common, only: rp
  implicit none
  private
  public get_weights,interp_fld
contains
  subroutine get_weights(xi,mmax,tol,xo,m,w)
    !
    ! for each output coordinate xo(i), finds the input points xi(m(i)) and xi(m(i)+1) that bracket it, with
    ! m(i) the largest index in [0,mmax] such that xi(m(i)) <= xo(i) (binary search), and the linear weight
    ! w(i) of xi(m(i)+1); output points that coincide with an input point (within tol) get a weight of
    ! exactly 0 or 1, so that their values are reproduced exactly
    !
    implicit none
    real(rp), intent(in ), dimension(0:) :: xi
    integer , intent(in )                :: mmax
    real(rp), intent(in )                :: tol
    real(rp), intent(in ), dimension(:)  :: xo
    integer , intent(out), dimension(:)  :: m
    real(rp), intent(out), dimension(:)  :: w
    integer :: i,ml,mu,mm
    !
    do i=1,size(xo)
      if(xo(i) >= xi(mmax)) then
        ml = mmax
      else
        ml = 0
        mu = mmax
        do while(mu-ml > 1)
          mm = (ml+mu)/2
          if(xi(mm) <= xo(i)) then
            ml = mm
          else
            mu = mm
          end if
        end do
      end if
      m(i) = ml
      if(     abs(xo(i)-xi(ml  )) <= tol) then
        w(i) = 0._rp
      else if(abs(xo(i)-xi(ml+1)) <= tol) then
        w(i) = 1._rp
      else
        w(i) = (xo(i)-xi(ml))/(xi(ml+1)-xi(ml))
      end if
    end do
  end subroutine get_weights
  !
  subroutine interp_fld(lo_i,lo_o,hi_o,mx,my,mz,wx,wy,wz,fldi,fldo)
    !
    ! trilinear interpolation of fldi onto the output block lo_o:hi_o, with the brackets m and weights w
    ! of each output point along x, y and z (see `get_weights`); both arrays are indexed as their global grid
    !
    implicit none
    integer , intent(in ), dimension(3) :: lo_i,lo_o,hi_o
    integer , intent(in ), dimension(lo_o(1):) :: mx
    integer , intent(in ), dimension(lo_o(2):) :: my
    integer , intent(in ), dimension(lo_o(3):) :: mz
    real(rp), intent(in ), dimension(lo_o(1):) :: wx
    real(rp), intent(in ), dimension(lo_o(2):) :: wy
    real(rp), intent(in ), dimension(lo_o(3):) :: wz
    real(rp), intent(in ), dimension(lo_i(1):,lo_i(2):,lo_i(3):) :: fldi
    real(rp), intent(out), dimension(lo_o(1):,lo_o(2):,lo_o(3):) :: fldo
    integer  :: i,j,k,ii,jj,kk
    real(rp) :: a,b,c
    !
    do k=lo_o(3),hi_o(3)
      kk = mz(k)
      c  = wz(k)
      do j=lo_o(2),hi_o(2)
        jj = my(j)
        b  = wy(j)
        do i=lo_o(1),hi_o(1)
          ii = mx(i)
          a  = wx(i)
          fldo(i,j,k) = (1._rp-c)*((1._rp-b)*((1._rp-a)*fldi(ii,jj  ,kk  ) + a*fldi(ii+1,jj  ,kk  ))   + &
                                          b *((1._rp-a)*fldi(ii,jj+1,kk  ) + a*fldi(ii+1,jj+1,kk  ))) + &
                                c *((1._rp-b)*((1._rp-a)*fldi(ii,jj  ,kk+1) + a*fldi(ii+1,jj  ,kk+1))   + &
                                          b *((1._rp-a)*fldi(ii,jj+1,kk+1) + a*fldi(ii+1,jj+1,kk+1)))
        end do
      end do
    end do
  end subroutine interp_fld
end module mod_interp
