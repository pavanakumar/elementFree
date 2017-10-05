!> @file connectivity.f90
!> The connectivity information 
!> (directed graph is stored here)

!> @brief The module that stores all connectivity and cloud size inforamtion
module Connectivity
  implicit none
  integer , dimension(:) , allocatable :: xadj , adjncy
  real (kind = 8) , dimension(:,:) , allocatable :: xyz
  real (kind = 8) , dimension(:) , allocatable :: min_dist 
  integer :: size = 0, adj_size = 0 , padding = 0 , max_ngbr
  logical :: connectivityAllocated = .false.
contains
  !> \brief Allocate connectivity
  subroutine AllocateConnectivity( )
    implicit none
    if( size .eq. 0 .or. adj_size .eq. 0 ) then
       write(*,*) 'Must try to allocate a non-zero size for connectivity'
       write(*,*) 'Xadj size = ' , size + 1 , 'Adj Size = ' , adj_size
       stop
    else
       allocate( xadj( size + 1 ) , min_dist( size ) )
       allocate( xyz(3,size + padding ) )
       allocate( adjncy(adj_size) )
       connectivityAllocated = .true.
       !write(*,*) 'Allocated xadj of size ' , size + 1 , " and adjncy of size " , adj_size 
    end if
  end subroutine AllocateConnectivity
  !> \brief Deallocate connectivity 
  subroutine DeallocateConnectivity( )
    implicit none
    if( connectivityAllocated .eqv. .true. ) then
       deallocate( xadj , adjncy )
       deallocate( xyz )
       connectivityAllocated = .false.
    else
       write(*,*) 'Warning trying to deallocate unallocated connectivity module'
    end if
  end subroutine DeallocateConnectivity
  !> \brief Check connectivity
  subroutine CheckConnectivity( )
    implicit none
    integer :: i , j 
    write(*,*) "Checking connectivity ..."
    do i = 1 , size
      do j = xadj(i) , xadj(i + 1) - 1  
        if( adjncy(j) .eq. i ) then
          if( j .ne. xadj(i) ) then
            adjncy(j) = adjncy(j-1)
            write(*,*) 'Warning: Node ' , i , ' is its own neighbour fixing it by repeating node'
          else
            write(*,*) 'Error: First node in ' , i , ' repeats ... check connectivity'
            stop 
          end if
        end if
        if( adjncy(j) .eq. 0 ) then
          write(*,*) 'Error: node number "0" found in connectivity'
          stop
        end if
        if( adjncy(j) .gt. size ) then
          write(*,*) 'Error: node number > max size found in connectivity'
          stop
        end if
      end do
    end do    
    write(*,*) "Connectivity ok ..."
  end subroutine CheckConnectivity
end module Connectivity


