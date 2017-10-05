!> @file leastSquares.f90
!> The least-squares derivative calculation
!> for the full stencil 

!> @brief Module that store all the least-squares coeff 
!>        abd has both poly/Taylor series based LS coeff
module LeastSquares
  use Connectivity
  implicit none
  real (kind = 8), dimension(:,:) , allocatable :: n_io, r_io , error
  real (kind = 8), dimension(:) , allocatable :: S_io, testFunction
  real (kind = 8), dimension(:) , allocatable :: praveen_deviation
  integer :: ls_pow = 1
  logical :: leastSquaresAllocated = .false. , polyLSCoeff = .false.
contains
  !> \brief 
  subroutine AllocateLSCoeff()
    implicit none
    if( connectivityAllocated .eqv. .true. ) then
       allocate( n_io(3,adj_size) )
       allocate( r_io(3,adj_size) )
       allocate( S_io(adj_size) , testFunction( size ) , error( 3 , size ) )
       allocate( praveen_deviation( size ) )
       leastSquaresAllocated = .true.
    else
       write(*,*) 'Must Allocate the Connectivity module before attempting allocation of LS module'
       stop
    end if
  end subroutine AllocateLSCoeff
  !> \brief 
  subroutine DeallocateLSCoeff()
    implicit none
    if( leastSquaresAllocated .eqv. .true. ) then
       deallocate( n_io ,r_io, S_io , testFunction , error )
    else 
       write(*,*) 'Attempting to deallocate unallocated leastSquares module'
       stop
    end if
  end subroutine DeallocateLSCoeff
  !> \brief
  subroutine FindMinDistance
    implicit none
    integer i,j,cnt
    real ( kind = 8 ) :: dist( max_ngbr ) , dxyz(3)
    do i = 1 , size
      do j = xadj(i) , xadj(i+1) - 1
        dxyz = xyz(:,adjncy(j)) - xyz(:,i) 
        dist( j - xadj(i) + 1 ) = dot_product( dxyz , dxyz )
      end do
      min_dist(i) = dsqrt( minval( dist( 1 : ( xadj(i+1) - xadj(i) ) ) ) )
    end do
  end subroutine FindMinDistance
  !> \brief
  subroutine formLSCoefficients()
    implicit none
    integer i,j,info,iwork(3)
    real (kind = 8) :: val(6),term(3),den,mag,a_norm=1.0d0
    character :: UPLO = 'U'
    real (kind = 8) :: condNum( size ), work(9)
    if( leastSquaresAllocated .eqv. .true. ) then
      polyLSCoeff = .false.
      max_ngbr = maxval( xadj( 2 : size + 1 ) - xadj( 1 : size ) )
      call FindMinDistance
!!! Calculate the Delta_x of point from neighbours
!!! Calculate the weights of neighbour points
      do i = 1 , size
        do j = xadj(i) , xadj(i + 1) - 1 
          r_io(:,j) = 0.50d0 * ( xyz(:,adjncy(j)) - xyz(:,i) ) !!! r_io = r_ij / 2
          S_io(j) = 1.0d0 / ( dot_product( r_io(:,j) , r_io(:,j) ) ** ls_pow ) !! 1/d^(2*ls_pow) weighting 
        end do
      end do
!!! Calculate the LS coeffs 
      do i = 1 , size
        val = 0.0d0
        do j = xadj(i) , xadj(i + 1) - 1 
          val(1) = val(1) + S_io(j) * r_io(1,j) ** 2
          val(2) = val(2) + S_io(j) * r_io(1,j) * r_io(2,j)
          val(3) = val(3) + S_io(j) * r_io(2,j) ** 2
          val(4) = val(4) + S_io(j) * r_io(1,j) * r_io(3,j)
          val(5) = val(5) + S_io(j) * r_io(2,j) * r_io(3,j)
          val(6) = val(6) + S_io(j) * r_io(3,j) ** 2
        end do
!!! Find the Cholesky decomposition ( LL^T x = b )
        call DPPTRF( UPLO , 3 , val , info )
        call DPPCON( UPLO , 3 , val , a_norm , condNum(i) , work , iwork , info )
        do j = xadj(i) , xadj(i + 1) - 1
          term = r_io(:,j) * S_io(j)
!!! Back Substitute to get the coeff (Lx* = b, L^Tx = x* )
          call DPPTRS( UPLO, 3, 1, val, term , 3, info )
!!! Now replace the weight with the magnitude 
!!! of the LS solution vector
          S_io(j)   = dsqrt( dot_product( term , term ) )
          n_io(:,j) = term / S_io(j)
        end do
      end do
    else
      write(*,*) 'Allocate the leastSquares module before calling its module functions'
      stop
    end if
    write(*,*) 'Max condition number = ' , maxval(condNum)
!    call formPraveenDeviation
  end subroutine formLSCoefficients

  subroutine formLSCoefficientsPoly
    implicit none
    integer i,j,info
    real (kind = 8) :: val(10),term(4),den,mag
    character :: UPLO = 'U'
 
    if( leastSquaresAllocated .eqv. .true. ) then
       polyLSCoeff = .true.
       max_ngbr = maxval( xadj(2:size+1) - xadj(1:size) )
       call FindMinDistance
!!! Calculate the Delta_x of point from neighbours
!!! Calculate the weights of neighbours of points
       do i = 1 , size
          do j = xadj(i) , xadj(i + 1) - 1 
             r_io(:,j) = 0.50d0 * ( xyz(:,adjncy(j)) - xyz(:,i) ) 
             S_io(j) = 1.0d0 / ( dot_product( r_io(:,j) , r_io(:,j) ) ** ls_pow ) 
          end do
       end do
!!! Calculate the Poly LS coeffs 
       do i = 1 , size
          val = 0.0d0
          do j = xadj(i) , xadj(i + 1) - 1
             val(1)  = val(1)  + 1.0d0 * S_io(j)
             val(2)  = val(2)  + r_io(1,j) * S_io(j)
             val(3)  = val(3)  + r_io(1,j) ** 2 * S_io(j)
             val(4)  = val(4)  + r_io(2,j) * S_io(j)
             val(5)  = val(5)  + r_io(2,j) * r_io(1,j) * S_io(j)
             val(6)  = val(6)  + r_io(2,j) ** 2 * S_io(j)
             val(7)  = val(7)  + r_io(3,j) * S_io(j)
             val(8)  = val(8)  + r_io(3,j) * r_io(1,j) * S_io(j)
             val(9)  = val(9)  + r_io(3,j) * r_io(2,j) * S_io(j)
             val(10) = val(10) + r_io(3,j) ** 2 * S_io(j)
         end do
!!! Find the Cholesky decomposition ( LL^T x = b )
          call DPPTRF( UPLO , 4 , val , info )
          do j = xadj(i) , xadj(i + 1) - 1
             term(1) = S_io(j)
             term(2:4) = r_io(1:3,j) * S_io(j)
!!! Back Substitute to get the coeff (Lx* = b, L^Tx = x* )
             call DPPTRS( UPLO, 4, 1, val, term , 4, info )
!!! Now replace the weight with the magnitude 
!!! of the LS solution vector
             S_io(j)   = dsqrt( dot_product( term(2:4) , term(2:4) ) )
             n_io(:,j) = term(2:4) / S_io(j)
          end do
       end do
    else  !!! Error handling 
       write(*,*) 'Allocate the leastSquares module before calling its module functions'
       stop
    end if
!    call formPraveenDeviation
  end subroutine formLSCoefficientsPoly

  !!! Polynomial weighted Least-Squares using SVD 
  subroutine formLSCoefficientsPolyWSVD
    implicit none

  end subroutine formLSCoefficientsPolyWSVD

  subroutine formPraveenDeviation
    implicit none
    real(kind = 8) :: deviation( max_ngbr ) , mag_r_io
    integer :: i,j
    do i = 1 , size
      do j = xadj(i) , xadj(i + 1) - 1 
        !! Calculate the dot product of least-squares coeff with
        !! the edge vector for all neighbours of i and get the maximum
        !! deviation
        mag_r_io = dsqrt( dot_product( r_io(:,j) , r_io(:,j) ) )
        deviation( j - xadj(i) + 1 ) = acos( dot_product( r_io(:,j) , n_io(:,j) ) / mag_r_io ) 
      end do
      praveen_deviation(i) = maxval( dabs( deviation( 1 : ( xadj(i+1) - xadj(i) ) ) ) )
    end do
  end subroutine formPraveenDeviation
  !> \brief Unit test
  !> Test if the LS can capture linear function exactly
  subroutine linearFunction
    implicit none
    integer i
    do i = 1 , size 
       testFunction(i) = sum( xyz(:,i) )
    end do
  end subroutine linearFunction
  !> \brief Test if the LS can capture linear function exactly
  subroutine constantFunction
    implicit none
    testFunction = 2.0d0
  end subroutine constantFunction
  !> \brief Unit test if the LS can capture quadratic function exactly
  subroutine quadraticFunction
    implicit none
    integer i
    do i = 1 , size 
       testFunction(i) = sum ( xyz(:,i) ** 2 ) 
    end do
  end subroutine quadraticFunction
  !> \brief
  subroutine testLinearGradFunction
    implicit none
    integer i,j,k
    real (kind=8) :: dfun , polyFactor = 0.5d0 , sine = -1.0d0
    if( polyLSCoeff .eqv. .true. ) then
      polyFactor = 0.50d0 ; sine = 1.0d0
    end if
    call linearFunction
    error = 0.0d0
    do i = 1 , size       
       do j = xadj(i) , xadj(i + 1) - 1
          dfun = polyFactor * ( testFunction(adjncy(j)) + sine * testFunction(i) )
          do k = 1 , 3 
             error(k,i) = error(k,i) + dfun * n_io(k,j) * S_io(j)
          end do
       end do
    end do
!    do i = 1 , size
!      error(:,i) = 1.0d0 - error(:,i)
!    end do
!    write(*,*) "Error in linear function gradient = " , maxval(dabs(error))
    open( unit=999 , file="testLinear.dat" )
    write(999,*) "variables=calc,actual,error"
    do i = 1 , size
      write(999,*) norm2(error(1:3,i)),sqrt(3.0),sqrt(3.0)-norm2(error(1:3,i))
    end do
    close(999)
  end subroutine testLinearGradFunction
  !> \brief
  subroutine testConstantGradFunction
    implicit none
    integer i,j,k
    real (kind=8) dfun
    call constantFunction
    error = 0.0d0
    do i = 1 , size       
       do j = xadj(i) , xadj(i + 1) - 1
          if( polyLSCoeff .eqv. .true. ) then
            dfun = testFunction(adjncy(j))
          else 
            dfun = testFunction(adjncy(j)) - testFunction(i)
          end if
          do k = 1 , 3 
            error(k,i) = error(k,i) + dfun * n_io(k,j) * S_io(j)
          end do
       end do
    end do
    write(*,*) "Error in constant function gradient = " , maxval(dabs(error))
  end subroutine testConstantGradFunction
  !> \brief
  subroutine testQuadraticGradFunction
    implicit none
    integer i,j,k
    real (kind=8) :: dfun , polyFactor = 0.5d0 , sine = -1.0d0 , denom(3)
    real (kind=8) :: exact(size)
    if( polyLSCoeff .eqv. .true. ) then
      polyFactor = 0.50d0 ; sine = 1.0d0
    end if
    call quadraticFunction
    error = 0.0d0
    do i = 1 , size       
       do j = xadj(i) , xadj(i + 1) - 1
          dfun = polyFactor * ( testFunction(adjncy(j)) + sine * testFunction(i) )
          do k = 1 , 3 
            error(k,i) = error(k,i) + dfun * n_io(k,j) * S_io(j)
          end do
       end do
    end do
    do i = 1 , size       
      exact(i) = norm2( 2.0 * xyz(:,i) )
    end do
    open( unit=999 , file="testQuadratic.dat" )
    write(999,*) "variables=calc,actual,error"
    do i = 1 , size
      write(999,*) norm2(error(1:3,i)) , exact(i) , &
                   exact(i) - norm2(error(1:3,i))
    end do
    close(999)
  end subroutine testQuadraticGradFunction
  !> \brief
  subroutine getPolyLSError( U , error )
    implicit none
    integer i,j,k,info
    real (kind = 8) , intent(in) :: U(5 , size + padding )
    real (kind = 8) :: error( 5 , size + padding )
    real (kind = 8) :: val(10) , term(5,4) , den , mag
    real (kind = 8) :: rloc(3) , Sloc , t_err(5) 
    character :: UPLO = 'U'
 
!! dilip telephone exchange == 9945848448
!!! Calculate the Delta_x of point from neighbours
!!! Calculate the weights of neighbours of points
     do i = 1 , size
       val = 0.0d0; term = 0.0d0
       do j = xadj(i) , xadj(i + 1) - 1 
         rloc = 0.50d0 * ( xyz(:,adjncy(j)) - xyz(:,i) ) 
         Sloc = 1.0d0 / ( dot_product( rloc , rloc ) ** ls_pow ) 
         val(1)    = val(1)  + 1.0d0 * Sloc
         val(2)    = val(2)  + rloc(1) * Sloc
         val(3)    = val(3)  + rloc(1) ** 2 * Sloc
         val(4)    = val(4)  + rloc(2) * Sloc
         val(5)    = val(5)  + rloc(2) * rloc(1) * Sloc
         val(6)    = val(6)  + rloc(2) ** 2 * Sloc
         val(7)    = val(7)  + rloc(3) * Sloc
         val(8)    = val(8)  + rloc(3) * rloc(1) * Sloc
         val(9)    = val(9)  + rloc(3) * rloc(2) * Sloc
         val(10)   = val(10) + rloc(3) ** 2 * Sloc
         do k = 1 , 5
           term(k,1)   = term(k,1) + U( k , adjncy(j) ) * Sloc
           term(k,2:4) = term(k,2:4) + U( k , adjncy(j) ) * rloc(1:3) * Sloc
         end do
      end do
      call DPPTRF( UPLO , 4 , val , info )
      k = 1 ; !! mass error
      call DPPTRS( UPLO , 4 , 1 , val , term(k,:) , 4 , info )
      error(k,i) = dabs( ( term(k,1) - U(k,i) ) ) / U(k,i) 
      !!! Momentum error
      do k = 2 , 4
        call DPPTRS( UPLO , 4 , 1 , val , term(k,:) , 4 , info )
        error(k,i) = dabs( ( term(k,1) - U(k,i) ) )
      end do
      Sloc = dsqrt( sum( U(2:4,i) ** 2 ) )
      error(2:4,i) = error(2:4,i) / Sloc
      !!! Energy error
      k = 5 ; !! mass error
      call DPPTRS( UPLO , 4 , 1 , val , term(k,:) , 4 , info )
      error(k,i) = dabs( ( term(k,1) - U(k,i) ) ) / U(k,i) 
    end do
  end subroutine getPolyLSError
end module LeastSquares

