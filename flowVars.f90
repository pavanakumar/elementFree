module flowVars
  use Boundary
  use Connectivity
  implicit none
  real (kind = 8) , dimension(:,:) , allocatable   :: U , U_temp , dF , Limiter
  real (kind = 8) , dimension(:) , allocatable     :: local_dt
  real (kind = 8) , dimension(:,:,:) , allocatable :: gradU
  contains
  !> \brief Allocate flow variables
  subroutine allocateFlowVars
    implicit none
    allocate( U(5,size+padding) , U_temp(5,size+padding) )
    allocate( dF(5,size+padding) , Limiter(5,size+padding) )
    allocate( local_dt(size+padding) )
    allocate( gradU(3,5,size+padding) )
  end subroutine allocateFlowVars
  !> \brief Deallocate flow variables
  subroutine deallocateFlowVars
    implicit none
    deallocate( U, U_temp , dF , Limiter , local_dt , gradU )
  end subroutine deallocateFlowVars
  !> \brief Write solution at boundary only
  subroutine writeFlowSolutionBC( file_name )
    implicit none
    integer :: i
    character (len=*) :: file_name
    open( unit=333 , file=trim(file_name) )
    write(333,*) 'variables=rho,u,v,w,p'
    do i = 1 , bc_sz
      write(333,*) U(1:5,bcNodeIDs(i))
    end do
    close(333)
  end subroutine writeFlowSolutionBC
  !> \brief Write solution
  subroutine writeFlowSolution( file_name )
    implicit none
    integer :: i
    character (len=*) :: file_name
    open( unit=333 , file=trim(file_name) )
    write(333,*) 'filetype=solution'
    write(333,*) 'variables=rho,u,v,w,p'
    do i = 1 , size
      write(333,*) U(1:5,i)
    end do
    close(333)
  end subroutine writeFlowSolution
  !> \brief Read solution
  subroutine readFlowSolution( file_name )
    implicit none
    integer :: i
    character (len=*) :: file_name
    open( unit=333 , file=trim(file_name), status='OLD', action='READ' )
    read(333,*) ! Header 1
    read(333,*) ! Header 2
    do i = 1 , size
      read(333,*) U(1:5,i)
    end do
    close(333)
  end subroutine readFlowSolution
  !> \brief
  subroutine setFlowValues( q )
    implicit none
    real (kind = 8) :: q(5) 
    integer i
    do i = 1 , size
      U(:,i)     = q
    end do
    U_temp = U
  end subroutine setFlowValues
end module flowVars

