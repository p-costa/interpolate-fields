module mod_bound
  use mpi_f08
  use mod_common, only: rp,MPI_REAL_RP,ierr
  private
  public makehalo,updthalo,set_bc
contains
  subroutine makehalo(idir,nh,n,halo)
    implicit none
    integer           , intent(in ) :: idir,nh
    integer           , intent(in ), dimension(3) :: n
    type(MPI_DATATYPE), intent(out) :: halo
    integer, dimension(3) :: nn
    nn(:) = n(:) + 2*nh
    select case(idir)
    case(1)
      call MPI_TYPE_VECTOR(nn(2)*nn(3),nh            ,nn(1)            ,MPI_REAL_RP,halo,ierr)
    case(2)
      call MPI_TYPE_VECTOR(      nn(3),nh*nn(1)      ,nn(1)*nn(2)      ,MPI_REAL_RP,halo,ierr)
    case(3)
      call MPI_TYPE_VECTOR(          1,nh*nn(1)*nn(2),nn(1)*nn(2)*nn(3),MPI_REAL_RP,halo,ierr)
    end select
    call MPI_TYPE_COMMIT(halo,ierr)
  end subroutine makehalo
  !
  subroutine updthalo(nh,halo,nb,idir,p)
    implicit none
    integer           , intent(in) :: nh ! number of ghost points
    type(MPI_DATATYPE), intent(in) :: halo
    integer           , intent(in), dimension(0:1) :: nb
    integer           , intent(in) :: idir
    real(rp)          , dimension(1-nh:,1-nh:,1-nh:), intent(inout) :: p
    integer           , dimension(3) :: lo,hi
    !
    !  this subroutine updates the halo that store info
    !  from the neighboring computational sub-domain
    !
    lo(:) = lbound(p)+nh
    hi(:) = ubound(p)-nh
    select case(idir)
    case(1) ! x direction
      call MPI_SENDRECV(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        p(hi(1)+1   ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh  ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    case(2) ! y direction
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
                        p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    case(3) ! z direction
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,lo(3)     ),1,halo,nb(0),0, &
                        p(lo(1)-nh,lo(2)-nh,hi(3)+1   ),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,hi(3)-nh+1),1,halo,nb(1),0, &
                        p(lo(1)-nh,lo(2)-nh,lo(3)-nh  ),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    end select
  end subroutine updthalo
  !
  subroutine set_bc(ctype,ibound,idir,nh,centered,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in) :: ibound,idir,nh
    logical , intent(in) :: centered
    real(rp), intent(in) :: rvalue,dr
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: factor,sgn
    integer  :: n,dh
    !
    n = size(p,idir) - 2*nh
    factor = rvalue
    if(ctype == 'D'.and.centered) then
      factor = 2.*factor
      sgn    = -1.
    end if
    if(ctype == 'N') then
      if(     ibound == 0) then
        factor = -dr*factor ! n.b.: only valid for nh /= 1 or factor /= 0
      else if(ibound == 1) then
        factor =  dr*factor ! n.b.: only valid for nh /= 1 or factor /= 0
      end if
      sgn    = 1.
    end if
    !
    dh = nh-1
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !$OMP WORKSHARE
        p(0  :  0-dh,:,:) = p(n:n-dh,:,:)
        p(n+1:n+1+dh,:,:) = p(1:1+dh,:,:)
        !$OMP END WORKSHARE
      case(2)
        !$OMP WORKSHARE
        p(:,0  :  0-dh,:) = p(:,n:n-dh,:)
        p(:,n+1:n+1+dh,:) = p(:,1:1+dh,:)
        !$OMP END WORKSHARE
      case(3)
        !$OMP WORKSHARE
        p(:,:,0  :  0-dh) = p(:,:,n:n-dh)
        p(:,:,n+1:n+1+dh) = p(:,:,1:1+dh)
        !$OMP END WORKSHARE
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if     (ibound == 0) then
            !$OMP WORKSHARE
            p(0  :  0-dh,:,:) = factor+sgn*p(1:1+dh,:,:)
            !$OMP END WORKSHARE
          else if(ibound == 1) then
            !$OMP WORKSHARE
            p(n+1:n+1+dh,:,:) = factor+sgn*p(n:n-dh,:,:)
            !$OMP END WORKSHARE
          end if
        case(2)
          if     (ibound == 0) then
            !$OMP WORKSHARE
            p(:,0  :  0-dh,:) = factor+sgn*p(:,1:1+dh,:)
            !$OMP END WORKSHARE
          else if(ibound == 1) then
            !$OMP WORKSHARE
            p(:,n+1:n+1+dh,:) = factor+sgn*p(:,n:n-dh,:)
            !$OMP END WORKSHARE
          end if
        case(3)
          if     (ibound == 0) then
            !$OMP WORKSHARE
            p(:,:,0  :  0-dh) = factor+sgn*p(:,:,1:1+dh)
            !$OMP END WORKSHARE
          else if(ibound == 1) then
            !$OMP WORKSHARE
            p(:,:,n+1:n+1+dh) = factor+sgn*p(:,:,n:n-dh)
            !$OMP END WORKSHARE
          end if
        end select
      else if(.not.centered.and.ctype == 'D') then
        select case(idir)
        case(1)
          if     (ibound == 0) then
            !$OMP WORKSHARE
            p(0:0-dh,:,:) = factor
            !$OMP END WORKSHARE
          else if(ibound == 1) then
            !$OMP WORKSHARE
            p(n+1   ,:,:) = p(n-1,:,:) ! unused
            p(n:n+dh,:,:) = factor
            !$OMP END WORKSHARE
          end if
        case(2)
          if     (ibound == 0) then
            !$OMP WORKSHARE
            p(:,0:0-dh,:) = factor
            !$OMP END WORKSHARE
          else if(ibound == 1) then
            !$OMP WORKSHARE
            p(:,n+1   ,:) = p(:,n-1,:) ! unused
            p(:,n:n+dh,:) = factor
            !$OMP END WORKSHARE
          end if
        case(3)
          if     (ibound == 0) then
            !$OMP WORKSHARE
            p(:,:,0:0-dh) = factor
            !$OMP END WORKSHARE
          else if(ibound == 1) then
            !$OMP WORKSHARE
            p(:,:,n+1   ) = p(:,:,n-1) ! unused
            p(:,:,n:n+dh) = factor
            !$OMP END WORKSHARE
          end if
        end select
      else if(.not.centered.and.ctype == 'N') then
        select case(idir)
        case(1)
          if     (ibound == 0) then
            !$OMP WORKSHARE
            !p(0,:,:) = 1./3.*(-2.*factor+4.*p(1  ,:,:)-p(2  ,:,:))
            p(0:0-dh,:,:) = 1.*factor + p(1  :  1+dh,:,:)
            !$OMP END WORKSHARE
          else if(ibound == 1) then
            !$OMP WORKSHARE
            !p(n,:,:) = 1./3.*(-2.*factor+4.*p(n-1,:,:)-p(n-2,:,:))
            p(n+1,:,:) = p(n,:,:) ! unused
            p(n:n+dh,:,:) = 1.*factor + p(n-1:n-1-dh,:,:)
            !$OMP END WORKSHARE
          end if
        case(2)
          if     (ibound == 0) then
            !$OMP WORKSHARE
            !p(:,0  ,:) = 1./3.*(-2.*factor+4.*p(:,1,:)-p(:,2  ,:))
            p(:,0:0-dh,:) = 1.*factor + p(:,1  :  1+dh,:)
            !$OMP END WORKSHARE
          else if(ibound == 1) then
            !$OMP WORKSHARE
            !p(:,n,:) = 1./3.*(-2.*factor+4.*p(:,n-1,:)-p(:,n-2,:))
            p(:,n+1,:) = p(:,n,:) ! unused
            p(:,n:n+dh,:) = 1.*factor + p(:,n-1:n-1-dh,:)
            !$OMP END WORKSHARE
          end if
        case(3)
          if     (ibound == 0) then
            !$OMP WORKSHARE
            !p(:,:,0) = 1./3.*(-2.*factor+4.*p(:,:,1  )-p(:,:,2  ))
            p(:,:,0:0-dh) = 1.*factor + p(:,:,1  :  1+dh)
            !$OMP END WORKSHARE
          else if(ibound == 1) then
            !$OMP WORKSHARE
            !p(:,:,n) = 1./3.*(-2.*factor+4.*p(:,:,n-1)-p(:,:,n-2))
            p(:,:,n+1) = p(:,:,n) ! unused
            p(:,:,n:n+dh) = 1.*factor + p(:,:,n-1:n-1-dh)
            !$OMP END WORKSHARE
          end if
        end select
      end if
    end select
  end subroutine set_bc
end module mod_bound
