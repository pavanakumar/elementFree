!> \file boundary.f90
!> The boundary condition information is stored 
!> in this module

!> \brief The boundary condition module where the 
!>        BC implementation is placed (also contains
!>        the size, index and normal information for 
!>        the boundary points)
module Boundary
  use connectivity
  implicit none 
  integer, dimension(:), allocatable :: wall_pts, ext_pts
  integer, dimension(:), allocatable :: in_pts, out_pts
  integer, dimension(:), allocatable :: sym_pts
  integer, dimension(:), allocatable :: bcNodeIDs
  integer :: wall_sz=0, ext_sz=0, in_sz=0, out_sz=0, sym_sz=0
  integer :: bc_sz=0
  real (kind = 8), dimension(:,:), allocatable :: wall_nxyz
  real (kind = 8), dimension(:,:), allocatable :: ext_nxyz
  real (kind = 8), dimension(:,:), allocatable :: in_nxyz
  real (kind = 8), dimension(:,:), allocatable :: out_nxyz
  real (kind = 8), dimension(:,:), allocatable :: sym_nxyz
  real (kind = 8) :: q_inf(5), M_inf
  logical :: boundaryAllocated = .false.
contains 
  !> \brief
  subroutine AllocateBoundary
    implicit none
    if( wall_sz .ne. 0 ) then 
       allocate( wall_pts(wall_sz) )
       allocate( wall_nxyz(3,wall_sz) )
    end if
    if( ext_sz .ne. 0  ) then 
       allocate( ext_pts(ext_sz) )
       allocate( ext_nxyz(3,ext_sz) )
    end if
    if( in_sz .ne. 0 ) then 
       allocate( in_pts(in_sz) )
       allocate( in_nxyz(3,in_sz) )
    end if
    if( out_sz .ne. 0 ) then 
       allocate( out_pts(out_sz) )
       allocate( out_nxyz(3,out_sz) )
    end if
    if( sym_sz .ne. 0 ) then 
       allocate( sym_pts(sym_sz) )
       allocate( sym_nxyz(3,sym_sz) )
    end if
    allocate( bcNodeIDs( wall_sz + ext_sz + in_sz + out_sz + sym_sz ) )
    boundaryAllocated = .true.
  end subroutine AllocateBoundary
  !> \brief
  subroutine DeallocateBoundary( )
    implicit none
    if( boundaryAllocated .eqv. .true. ) then
       if( wall_sz .ne. 0 ) deallocate( wall_pts, wall_nxyz )
       if( ext_sz  .ne. 0 ) deallocate( ext_pts, ext_nxyz )
       if( in_sz   .ne. 0 ) deallocate( in_pts, in_nxyz )
       if( out_sz  .ne. 0 ) deallocate( out_pts, out_nxyz )
       if( sym_sz  .ne. 0 ) deallocate( sym_pts, sym_nxyz )
       deallocate( bcNodeIDs )
       boundaryAllocated = .false.
    else
       write(*,*) 'Attempting to allocated unallocated boundary module'
       stop
    end if
  end subroutine DeallocateBoundary
  !> \brief
  subroutine ApplyWallBC( q, tot_sz )
    implicit none
    integer :: i,tot_sz
    real (kind=8), dimension(5,tot_sz) :: q
    real (kind=8) :: Un 
    do i = 1, wall_sz
      Un = dot_product( q(2:4,wall_pts(i)), wall_nxyz(1:3,i) )
      q(2:4,wall_pts(i)) = q(2:4,wall_pts(i)) - wall_nxyz(1:3,i) * Un
    end do
  end subroutine ApplyWallBC
  !> \brief
  subroutine ApplySymBC( q, tot_sz )
    implicit none
    integer :: i,tot_sz
    real (kind=8), dimension(5,tot_sz) :: q
    real (kind=8) :: Un 
    do i = 1, sym_sz
      Un = dot_product( q(2:4,sym_pts(i)), sym_nxyz(1:3,i) )
      q(2:4,sym_pts(i)) = q(2:4,sym_pts(i)) - sym_nxyz(1:3,i) * Un
    end do
  end subroutine ApplySymBC
  !> \brief
  subroutine ApplyExteriorBC( q, tot_sz )
    use Constants
    !!! _plus/_minus -> Outgoing/Incoming Characteristics
    !!! _ff   -> Far-field BC values
    !!! _i    -> Interior values
    !!! _inf  -> Free-stream values
    implicit none
    integer :: tot_sz,i,ii
    real (kind=8),dimension(5,tot_sz) :: q
    real (kind=8) :: v_plus, v_minus, q_i(5)
    real (kind=8) :: un_i, a_i, un_inf, a_inf
    real (kind=8) :: un_ff, a_ff, s
    do i = 1, out_sz
      q_i     = q(:,out_pts(i))
      un_i    = dot_product( q_i(2:4), out_nxyz(:,i) )
      a_i     = dsqrt( gama * q_i(5) / q_i(1) )
      un_inf  = dot_product( q_inf(2:4), out_nxyz(:,i) )
      a_inf   = dsqrt( gama * q_inf(5) / q_inf(1) )
      v_plus  = un_i   + ( 2.0d0 * a_i ) * invgm1
      v_minus = un_inf - ( 2.0d0 * a_inf ) * invgm1
      un_ff   = 0.50d0  * ( v_plus + v_minus )
      a_ff    = 0.250d0 * ( v_plus - v_minus ) * gm1
      if( un_ff .lt. 0.0e0 ) then
      !!! Inlet
        q_i(2:4) = q_inf(2:4) + ( un_ff - un_inf ) * out_nxyz(:,i)
        s       = ( q_inf(1) ** gama ) / q_inf(5)
      !!! Outlet
      else
        q_i(2:4) = q_i(2:4) + ( un_ff - un_i ) * out_nxyz(:,i)
        s       = ( q_i(1) ** gama ) / q_i(5)
      end if
        q_i(1)   = ( a_ff * a_ff * s * invgama ) ** invgm1
        q_i(5)   = a_ff * a_ff * q_i(1) * invgama
        q(:,out_pts(i)) = q_i
    end do
  end subroutine ApplyExteriorBC
  !> \brief
  subroutine ApplyInflowBC( q, tot_sz )
    implicit none
    integer :: tot_sz,i
    real (kind=8),dimension(5,tot_sz) :: q
    do i = 1, in_sz
      q(:,in_pts(i)) = q_inf
    end do
  end subroutine ApplyInflowBC
  !> \brief
  subroutine ApplyExtrapolationBC( q, tot_sz )
    implicit none
    integer :: tot_sz,i
    real (kind=8),dimension(5,tot_sz) :: q
    do i = 1, ext_sz
      q(:,ext_pts(i)) = q(:,adjncy(xadj(ext_pts(i))))
    end do
  end subroutine ApplyExtrapolationBC
  !> \brief
  subroutine ApplyBC( q, tot_sz )
    implicit none
    integer :: tot_sz
    real (kind=8),dimension(5,tot_sz) :: q
    if( wall_sz.ne.0 ) call ApplyWallBC( q, tot_sz )
    if(  sym_sz.ne.0 ) call ApplySymBC( q, tot_sz )
    if(  out_sz.ne.0 ) call ApplyExteriorBC( q, tot_sz )
    if(  ext_sz.ne.0 ) call ApplyExtrapolationBC( q, tot_sz )
    if( in_sz.ne.0 ) call ApplyInflowBC( q, tot_sz )
  end subroutine ApplyBC
  !> \brief Unit test functions
  subroutine printExtBCPoints
    implicit none
    integer i
    do i = 1, ext_sz
      write(*,*) xyz( :, adjncy( xadj( ext_pts(i) ) ) )
    end do
  end subroutine printExtBCPoints
  !> \brief
  subroutine setFreeStream(q0)
    implicit none
    real (kind=8) :: q0(5)
    q_inf = q0
  end subroutine setFreeStream
  !> \brief
  subroutine printBoundaryStats
    implicit none
     if( wall_sz .ne. 0 ) write(*,*) "Total wall pts = ", wall_sz
     if( ext_sz .ne. 0  ) write(*,*) "Total extrapolation pts  = ", ext_sz
     if( out_sz .ne. 0  ) write(*,*) "Total outflow pts  = ", out_sz
     if( sym_sz .ne. 0  ) write(*,*) "Total symmetry pts  = ", sym_sz
     if( in_sz .ne. 0  ) write(*,*) "Total inlet pts  = ", in_sz
     write(*,*) "Total BC pts   = ", bc_sz
  end subroutine printBoundaryStats
  !> \brief
  subroutine writeTestBC
    implicit none
    integer :: i
    !!! Wall test
    if( wall_sz .ne. 0 ) then
      open( unit=111, file="testWall.dat" )
      write(111,*) "variables=x,y,z,nx,ny,nz"
      do i = 1, wall_sz
        write(111,*) xyz(:,wall_pts(i)), wall_nxyz(:,i)
      end do
      close(111)
   end if
   !!! Inflow test
   if( in_sz .ne. 0 ) then
      open( unit=111, file="testIn.dat" )
      write(111,*) "variables=x,y,z,nx,ny,nz"
      do i = 1, in_sz
        write(111,*) xyz(:,in_pts(i)), in_nxyz(:,i)
      end do
      close(111)
   end if
  end subroutine writeTestBC
end module Boundary

