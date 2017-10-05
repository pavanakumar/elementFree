!> \file main.f90
!> The program to test the GFL fortran modules
program finitePoint
  use Connectivity
  use Boundary
  use LeastSquares
  use Flow
  use FlowVars
  implicit none
  !!! Open input files and put in various modules
  call processFiles
  !!! Init flow values
  call initFlowVars
  !!! Run the solver 
  call runSolver
  !!! Clean up
  write(*,*) "Cleaning up ..."
  call deallocateConnectivity
  call deallocateBoundary
  call deallocateLSCoeff
  call deallocateFlowVars
end program finitePoint
!> \brief Convert an integer to string
function toStr(k)
    integer, intent(in) :: k
    character(len=512) :: toStr
    write (toStr, *) k
    toStr = adjustl(toStr)
    toStr = trim(toStr)
end function toStr
!> \brief Read anandh dat files
subroutine processFiles
  use Connectivity
  use Boundary
  use LeastSquares
  implicit none
  logical :: connectivityCheck = .false., testLSGrad = .false.
  character (len=256) :: lsTypeString
  namelist /lsInput/ lsTypeString, connectivityCheck, testLSGrad
  !!! Read input namelist file
  open( unit=666, file="input.nam", status="old" )
  read(666,nml=lsInput)
  close(666) 
  !!! Read the input files
  write(*,*) "Reading connectivity and bc data ..."
  call readGrid
  if( connectivityCheck .eqv. .true. ) &
    call checkConnectivity
  !!! form the least-sqiares coeffs
  call formLS( trim(lsTypeString) )
  !!! Test the coeffs
  if( testLSGrad .eqv. .true. ) then
    call testConstantGradFunction
    call testLinearGradFunction
    call testQuadraticGradFunction
  end if
  write(*,*) "Finished reading connectivity and bc data ..."
end subroutine processFiles
!> \brief Form the least-squares coeffs
subroutine formLS( args )
  use LeastSquares
  implicit none
  character (len=*) :: args
  call allocateLSCoeff
  !!! Calculate the least-squares coefficients
  if( args .eq. "taylor" .or. args .eq. "TAYLOR") then
    write(*,*) "Forming the Taylor series LS coeff ..."
    call formLSCoefficients
  else
    write(*,*) "Forming the Polynomial LS coeff ..."
    call formLSCoefficientsPoly
  end if
end subroutine formLS
!> \brief Initilize the flow variables
subroutine initFlowVars
  use Constants
  use Boundary
  use Flow
  use FlowVars
  implicit none
  real (kind=8) :: axis(3), angle, u_inf(3)
  namelist /input/ M_inf, axis, angle
  write(*,*) "Allocating and setting up flow variables ..."
  call allocateFlowVars
  !!! Read input namelist file
  open( unit=666, file="input.nam", status="old" )
  read(666,nml=input)
  close(666)
  angle = angle * degToRad
  u_inf = (/ 1.0, 0.0, 0.0   /) !! Assume angle is measured wrt x axis
  call rotateAxisAngle( u_inf, axis, angle )
  write(*,*) 'Freestream Conditions'
  write(*,*) '====================='
  write(*,*) 'Mach Number     = ' , M_inf
  write(*,*) 'Angle of Attack = ' , angle * radToDeg
  write(*,*) 'Freestream velocity = ' , u_inf
  q_inf(1)   = 1.0d0
  q_inf(2:4) = u_inf(1:3)
  q_inf(5)   = 1.0d0 / ( gama * M_inf )
  call setFlowValues( q_inf )
end subroutine initFlowVars
!> \brief Read solver inputs from namelist file
subroutine runSolver
  use Flow
  use FlowVars
  implicit none
  integer :: numIter, numPass, i, solIter = 100, resIter = 1
  integer :: j
  real(kind = 8):: rkCoeff(10)
  character (len=512) :: solFilePrefix, toStr
  logical :: is_restart = .FALSE.
  namelist /solver/ numIter, CFL, numPass, & 
                    rkCoeff, resIter, solIter, solFilePrefix
  open( unit=888, file="input.nam", status="old" )
  read(888,nml=solver)
  close(888)
  write(*,*) "Solver Inputs"
  write(*,*) "============="
  write(*,*) "Total Iterations = ", numIter
  write(*,*) "CFL number       = ", CFL
  write(*,*) "Number of passes = ", numPass
  write(*,*) "RK Coefficients  = ", rkCoeff(1:numPass)
  write(*,*) "Residue print    = ", resIter
  write(*,*) "Solution print   = ", solIter
  write(*,*) "Solution file name = ", trim(solFilePrefix)
  !!! Explicit time step loop
  do i = 1 , numIter
    call ApplyWallBC( U , size )
    do j = 1 , numPass
      call LocalTimeStep( U , local_dt )
      call FirstOrderResidue( U , dF );
      call RKUpdate( U , U_temp , dF , local_dt , j , numPass , rkCoeff(j) )
      call ApplyBC( U , size );
    end do
    if( mod(i,resIter) .eq. 0 ) & 
      write(*,*) "Iteration = ", i, "Maximum residue = ", max_res
    if( mod(i,solIter) .eq. 0 ) then
      call writeFlowSolution( trim(solFilePrefix)//"_"//trim(toStr(i))//".dat" )
!      call writeFlowSolutionBC( trim(solFilePrefix)//"_"//trim(toStr(i))//".bc.dat" )
    end if
  end do

end subroutine runSolver
!> \brief Axis angle rotation
!>        Rodrigues' rotation formula
!>        Can be more efficiently done
subroutine rotateAxisAngle( vec, axis, angle )
  implicit none
  real(kind=8) :: vec(3),axis(3),angle
  real(kind=8) :: N(3,3), N2(3,3)
  !!! Skew matrix construction
  N(2,1) = axis(3)
  N(3,1) = -axis(2)
  N(1,2) = -axis(3)
  N(1,3) = axis(2)
  N(2,3) = -axis(1)
  N(3,2) = axis(1)
  N2 = MATMUL( N, N )
  N = N * dsin(angle) + N2 * ( 1.0d0 - dcos(angle) )
  vec = vec + MATMUL( N, vec )
end subroutine rotateAxisAngle

subroutine readGrid
  use connectivity
  use boundary
  implicit none
  character (len=*), parameter :: gridFile = "blcone.dat"
  character (len=*), parameter :: normalFile = "blcone_normal.dat"
  integer :: num_bcs, nngbr, bc_flag, iorder
  integer :: node_id, i, bc_node_id
  real (kind=8) :: tgt(3), norm(3)
  integer :: wall_count = 0, ext_count = 0, out_count = 0, sym_count = 0
  integer :: in_count = 0

  open( unit=111, file=gridFile, status='old' )
  open( unit=222, file=normalFile, status='old' )

  read(111,*) size, num_bcs
  allocate( xyz(3,size), xadj(size + 1), min_dist(size) )

  !!! Read nodes from blcone and check their bc_flag 
  !!! and increment counter to get sizes
  xadj(1) = 1
  adj_size = 0
  do i = 1, size
    read(111,*) node_id, xyz(1:3,i), nngbr, bc_flag, iorder
    read(111,*) ! ignore neighbours
    xadj(i+1) = xadj(i) + nngbr
    adj_size = adj_size + nngbr
    if( bc_flag .ne. 0 ) then
      read(222,*) bc_node_id! ignore tangents and normals
      if( bc_node_id .ne. node_id ) then
        write(*,*) 'Make sure the blcone_normal.dat is sorted'
        write(*,*) 'mv blcone_normal.dat blcone_normal.old'
        write(*,*) 'sort -n blcone_normal.old > blcone_normal.dat'
        stop
      end if
      !!! Count wall bc points
      if( bc_flag.eq.1 ) wall_sz = wall_sz + 1
      if( bc_flag.eq.2 ) out_sz  = out_sz + 1
      if( bc_flag.eq.4 ) sym_sz  = sym_sz + 1
    end if
  end do

  ! re-read file
  allocate( adjncy(adj_size) )
  call AllocateBoundary

  rewind(111)
  rewind(222)

  read(111,*) ! ignore size header
  do i = 1, size
    read(111,*) node_id, xyz(1:3,i), nngbr, bc_flag, iorder
    if( nngbr .ne. ( xadj(i+1) - xadj(i) ) ) then
      write(*,*) "Error: Number of neighbours is not consistent in Adjncy list"
      write(*,*) nngbr, xadj(i+1) - xadj(i)
      stop
    end if
    read(111,*) adjncy( xadj(i) : xadj(i+1) - 1 )
     if( bc_flag .ne. 0 ) then
      read(222,*) bc_node_id, tgt(1:3), norm(1:3)
      !!! Wall boundary
      if( bc_flag.eq.1 ) then
        wall_count = wall_count + 1
        wall_pts(wall_count) = node_id + 1
        wall_nxyz(:,wall_count) = -norm
      end if
      !!! Riemann boundary
      if( bc_flag.eq.2 ) then
        out_count = out_count + 1
        out_pts(out_count) = node_id + 1
        out_nxyz(:,out_count) = -norm
      end if
      !!! Symmetry boundary
      if( bc_flag.eq.4 ) then
        sym_count = sym_count + 1
        sym_pts(sym_count) = node_id + 1
        sym_nxyz(:,sym_count) = -norm
      end if
    end if
  end do

  do i = 1, adj_size
    adjncy(i) = adjncy(i) + 1
  end do
  connectivityAllocated = .true.
  close(111)
  close(222)
  !! Check bcs
  open( 777, file="bc_check.dat" )
  write(777,*) 'variables=x,y,z,nx,ny,nz'
  do i = 1, wall_sz
    write(777,*) xyz(:,wall_pts(i)), wall_nxyz(:,i)
  end do
  call printBoundaryStats
end subroutine readGrid

