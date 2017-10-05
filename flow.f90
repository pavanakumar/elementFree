!> @file flow.f90
!> The flow realted functions and inputs occur here
!> @brief The flow module houses all flux related operations
module Flow
  use Connectivity
  use LeastSquares
  use Constants
  implicit none
  real (kind = 8) :: max_res(5), max_xyz(3), max_l2_res(5), CFL = 0.3d0
  real (kind = 8) :: tol = 1.0d-6
contains 
  !> \brief Calculate the local dt for local time stepping
  subroutine ReciprocalTimeScale( U, speed )
    implicit none
    real (kind = 8), intent(in)  :: U( 5, size + padding )
    real (kind = 8), intent(out) :: speed(size)
    integer i,j
    speed = 0.0d0
    do i = 1, size
      do j = xadj(i), xadj(i + 1) - 1
        speed(i) = speed(i) + ( dabs( dot_product( U(2:4,i), n_io(:,j) ) ) &
                      + dsqrt( gama * U(5,i) / U(1,i) ) ) * S_io(j)
      end do
    end do
  end subroutine ReciprocalTimeScale
  !> \brief Calculate the local dt for local time stepping
  subroutine LocalTimeStep( U, local_dt )
    implicit none
    real (kind = 8), intent(in)  :: U( 5, size + padding )
    real (kind = 8), intent(out) :: local_dt(size)
    integer i,j
    call ReciprocalTimeScale( U, local_dt )
    do i = 1, size
      local_dt(i) = CFL / local_dt(i)
    end do
  end subroutine LocalTimeStep
  !> \brief Remember that this gradient is calculated based on the 
  !>        mid-point LS coeffs so actual grad is half of this value
  subroutine CalculateGradU( U, gradU )
    implicit none
    real (kind = 8), intent(out) :: gradU(3,5,size + padding)
    real (kind = 8), intent(in)  :: U(5,size + padding)
    integer :: i,j,k
    real (kind=8) :: temp_f(5), UL(5), UR(5)
    real (kind=8) :: dfun(5), polyFactor = 1.0d0, sine = -1.0d0
    if( polyLSCoeff .eqv. .true. ) then
      polyFactor = 0.50d0 ; sine = 1.0d0
    end if
    gradU = 0.0d0
    do i = 1, size
      do j = xadj(i), xadj(i + 1) - 1
        dfun = polyFactor * ( U(:,adjncy(j)) + sine * U(:,i) )
        do k = 1, 3
          gradU(k,:,i) = gradU(k,:,i) + dfun * n_io(k,j) * S_io(j)
        end do
      end do
    end do
  end subroutine CalculateGradU
  !> \brief
  subroutine CalculateLimiter( U, gradU, Limiter )
    implicit none
    real (kind = 8), intent(in)  :: gradU(3,5,size + padding), U(5,size + padding)
    real (kind = 8), intent(out) :: Limiter(5,size + padding)
    integer :: i,j,k,cnt
    real (kind=8) :: Umax(5), Umin(5), U_tilde(5)
    real (kind=8) :: tempLimiter(5,max_ngbr), polyFactor = 1.0d0 
    if( polyLSCoeff .eqv. .true. ) polyFactor = 2.0d0 
    do i = 1, size
      !!! Find the max and min state value among all neighbouring nodes
      do k = 1, 5
        Umax(k) = max( U(k,i), maxval( U(k, adjncy(xadj(i) : xadj(i + 1) - 1)) ) )
        Umin(k) = min( U(k,i), minval( U(k, adjncy(xadj(i) : xadj(i + 1) - 1)) ) )
      end do
      !!! Calculate the grad_U \cdot r_io with every neighbouring point
      !!! for constructing the Limiter function
      cnt = 1
      do j = xadj(i), xadj(i + 1) - 1
        do k = 1, 5 
          U_tilde(k) = polyFactor * dot_product( gradU(:,k,i), r_io(:,j) )
          if ( U_tilde(k) .gt. 0.0d0 ) then 
            tempLimiter(k,cnt) = min( 1.0d0, ( Umax(k) - U(k,i) ) / U_tilde(k) )
          else if( U_tilde(k) .lt. 0.0d0 ) then
            tempLimiter(k,cnt) = min( 1.0d0, ( Umin(k) - U(k,i) ) / U_tilde(k) )
          else
            tempLimiter(k,cnt) = 1.0d0
          end if
          if( tempLimiter(k,cnt) .lt. 0.0d0 ) then
            write(*,*) "Non-positive limiter value found "
          end if
        end do !! loop over k ends
        cnt = cnt + 1
      end do !! loop over j ends
      !!! Construct limiter function
      do k = 1, 5
        Limiter(k,i) = minval( tempLimiter( k, 1 : xadj(i+1) - xadj(i) ) )
      end do
    end do !! loop over i ends
  end subroutine CalculateLimiter
  !> \brief Differential limiter implementation (Venkatakrishnan)
  subroutine CalculateLimiterVenkat( U, gradU, Limiter )
    implicit none
    real (kind = 8), intent(in)  :: gradU(3,5,size + padding), U(5,size + padding)
    real (kind = 8), intent(out) :: Limiter(5,size + padding)
    real (kind=8) :: Umax(5), Umin(5), U_tilde(5), dUmax, omega, Klimit=0.01d0
    real (kind=8) :: numerator, denominator
    real (kind=8) :: tempLimiter(5,max_ngbr), dU(5,max_ngbr)
    integer :: i,j,k,cnt

    do i = 1, size    !size = npts
      omega = Klimit * ( min_dist(i)**3 )
      !!! Find the max and min state value among all neighbouring nodes
      do k = 1, 5
        Umax(k) = max( U(k,i), maxval( U(k, adjncy(xadj(i) : xadj(i + 1) - 1)) ) )
        Umin(k) = min( U(k,i), minval( U(k, adjncy(xadj(i) : xadj(i + 1) - 1)) ) )
      end do
         
      do k = 1,5
        cnt = 1
        do j = xadj(i), xadj(i + 1) - 1
          dU(k,cnt) =  U(k,adjncy(j)) - U(k,i)
          cnt = cnt + 1
        end do
      end do

      do k = 1,5
        cnt = 1
        do j = xadj(i), xadj(i + 1) - 1
          if( dU(k,cnt) .gt. 0.d0 ) then
            dUmax = Umax(k) - U(k,i)
          else
            dUmax = Umin(k)  - U(k,i)
          end if    
          numerator = dUmax**2 + 2.0d0 * dUmax * dU(k,cnt) + omega
          denominator = dUmax**2 + 2.0d0 * ( dU(k,cnt)**2 ) + dUmax * dU(k,cnt) + omega
          tempLimiter(k,cnt) = numerator / denominator
          cnt = cnt + 1
        end do
      end do

      do k = 1, 5
        Limiter(k,i) = minval( tempLimiter( k, 1 : xadj(i+1) - xadj(i) ) )
      end do

    end do !! loop over i ends
  end subroutine CalculateLimiterVenkat
  !> \brief Calculate the second order residuals for parallel version
  subroutine SecondOrderResidue( U, gradU, Limiter, dF_i )
    implicit none
    real (kind = 8), intent(in)  :: gradU(3, 5, size + padding ), U( 5, size + padding )
    real (kind = 8), intent(in)  :: Limiter(5, size + padding )
    real (kind = 8), intent(out) :: dF_i( 5, size + padding  )
    integer :: i,j,k
    real (kind=8) :: F_i(5), UL(5), UR(5)
    real (kind=8) :: F_o(5), polyFactor = 0.50d0
    if( polyLSCoeff .eqv. .true. ) polyFactor = 1.0d0
    dF_i = 0.0d0
    do i = 1, size
      do j = xadj(i), xadj(i + 1) - 1
        do k = 1, 5
          UL(k) = U(k,i) + polyFactor * Limiter(k,i) * dot_product( gradU(:,k,i), r_io(:,j) )
          UR(k) = U(k,adjncy(j)) - & 
                  polyFactor * Limiter(k,adjncy(j)) * dot_product( gradU(:,k,adjncy(j)), r_io(:,j) )
        end do
        call roe( UL, UR, n_io(:,j), F_o )
        !! Should think of a way to remove this branch statement 
        if( polyLSCoeff .eqv. .false. ) then
          F_i = normalFlux( U(:,i), n_io(:,j) )
          F_o = F_o - F_i !!! remember dx_io = (xo - xi) in ls module
        end if
        dF_i(:,i) = dF_i(:,i) + F_o * S_io(j)
      end do
    end do
  end subroutine SecondOrderResidue
  !> \brief Calculate the First Order Residual
  subroutine FirstOrderResidue( U, dF_i )
    implicit none
    real (kind = 8), intent(in)  :: U( 5, size + padding )
    real (kind = 8), intent(out) :: dF_i( 5, size + padding  )
    integer :: i,j
    real (kind=8) :: F_io(5), F_i(5)
    dF_i = 0.0d0
    !!! private (i,j,F_io,F_i,polyLSCoeff) shared(xadj,U,adjncy,dF_i)
    do i = 1, size
      do j = xadj(i), xadj(i + 1) - 1
        call roe( U(:,i), U(:,adjncy(j)), n_io(:,j), F_io )
        !! Should think of a way to remove this branch statement
        if( polyLSCoeff .eqv. .false. ) then
          F_i = normalFlux( U(:,i), n_io(:,j) )
          F_io = F_io - F_i !!! remember dx_ij = (xj - xi) in ls module
        end if
        dF_i(:,i) = dF_i(:,i) + F_io * S_io(j)
      end do
    end do
  end subroutine FirstOrderResidue
  !> \brief Remember primitive form [ rho, u, ..., p ]
  subroutine ConvertToPrimitive( U )
    implicit none
    real (kind = 8), intent(inout) :: U( 5, size + padding )
    integer i
    real (kind=8) :: u_sqr
    do i = 1, size
      !!! u,v,w = (rhou,rhov,rhow)/rho
      U(2:4,i) = U(2:4,i) / U(1,i)
      !!! p = (g - 1) (rho et - 1/2 rho u^2 )
      u_sqr = dot_product( U(2:4,i), U(2:4,i) )
      U(5,i) = gm1 * ( U(5,i) - 0.50d0 * U(1,i) * u_sqr )
    end do
  end subroutine ConvertToPrimitive
  !> \brief Remember the conservative form [ rho, rho u,.., rho e ] 
  subroutine ConvertToConservative( U )
    implicit none
    real (kind = 8), intent(inout) :: U( 5, size + padding )
    integer i
    do i = 1, size
      !!! rho et = p / ( g -1 ) + 1/2 * rho * u^2
      U(5,i) = U(5,i) / gm1 + 0.50d0 * U(1,i) * dot_product( U(2:4,i), U(2:4,i) )
      !!! (rhou,rhov,rhow) = (u,v,w)rho
      U(2:4,i) = U(2:4,i) * U(1,i)
    end do
  end subroutine ConvertToConservative
  !> \brief e = p / (rho (g -1)) + 1/2 u^2
  !> rhoe + p = p( 1+ 1/(g-1) ) + 1/2 rho u^2
  !> (rhoe + p) = rho/rho * ( g/(g-1) p + 1/2 rho u^2 )
  !> (rhoe+p) = rho ( g/(g-1) p / rho + 1/2 u^2 )
  !> (rhoe+p)un = rho un ( g/(g-1) p / rho + 1/2 u^2 )
  function normalFlux( q, normal ) result( Fn )
    implicit none
    real (kind=8) :: Fn(5), q(5), normal(3)
    integer :: i
    Fn(1) = q(1) * dot_product( q(2:4), normal(1:3) )
    Fn(2:4) = Fn(1) * q(2:4) + q(5) * normal(1:3)
    Fn(5) = Fn(1) * ( gbygm1 * q(5) / q(1) + &
         0.50d0 * dot_product( q(2:4), q(2:4) ) )
  end function normalFlux
  !> \brief Using the Van-Leer Scheme from 
  !> Nikhil Shindes thesis   
  subroutine van_leer_hifun(qL, qR, normal, F)
    implicit none
    real (kind=8), intent(in)  :: qL(5), qR(5), normal(3)
    real (kind=8), intent(out) :: F(5)
    real (kind=8) :: cL, cR, MnL, MnR, UnL,UnR, FL(5)
    cR = dsqrt( gama * qR(5) / qR(1) )
    cL = dsqrt( gama * qL(5) / qL(1) )
    UnR = dot_product( qR(2:4), normal )
    UnL = dot_product( qL(2:4), normal )
    MnR =  UnR / cR
    MnL =  UnL / cL 
    if( MnL .le. -1.0d0 ) then
      FL = 0.0d0
    else if( MnL .ge. 1.0d0 ) then
      FL = normalFlux( qL, normal )
    else
      FL(1)   = 0.250d0 * qL(1) * cL * ( 1.0d0 + MnL ) ** 2
      FL(2:4) = FL(1) * ( qL(2:4) + ( 2.0d0 - MnL ) * cL * normal(1:3) * invgama )
      FL(5)   = FL(1) * ( invgm1 * cL ** 2 + 0.50d0 * dot_product( qL(2:4), qL(2:4) ) )
    end if
    if( MnR .le. -1.0d0 ) then
      F = normalFlux( qR, normal )
    else if( MnR .ge. 1.0d0 ) then
      F = 0.0d0
    else
      F(1)   = - 0.250d0 * qR(1) * cR * ( 1.0d0 - MnR ) ** 2
      F(2:4) = F(1) * ( qR(2:4) - ( 2.0d0 + MnR ) * cR * normal(1:3) * invgama )
      F(5)   = F(1) * ( invgm1 * cR ** 2 + 0.50d0 * dot_product( qR(2:4), qR(2:4) ) )
    end if
    F = F + FL
  end subroutine van_leer_hifun
  !> \brief The KFVS formulation is not 
  !> verbatim from Nikhil Shinde 
  !> thesis (HiFUN) as it contains typos
  !> qL/R contain the primitive variables
  !> { \rho, u, v, w, p }L/R
  subroutine kfvs(qL, qR, normal, F)
    implicit none
    real (kind = 8), intent(in)  :: qL(5), qR(5), normal(3)
    real (kind = 8), intent(out) :: F(5)
    real (kind = 8) :: FL(5)
    call kfvs_splflx( qL, normal, 1.0d0, FL )
    call kfvs_splflx( qR, normal, -1.0d0, F  )
    F = F + FL
  end subroutine kfvs
  !> \brief Using Anandhanarayanan technique here
  subroutine kfvs_splflx( q, normal, sign_, F )
    implicit none
    integer i
    real (kind=8) :: F(5), q(5), normal(3), sign_
    real (kind=8) :: An, Bn, Sn, Un, sqrtBeta, term
    sqrtBeta = dsqrt( 0.50d0 * q(1) / q(5) )
    Un       = dot_product( q(2:4), normal(1:3) )
    Sn       = Un * sqrtBeta
    An       = 0.50d0 * ( 1.0d0 + sign_ * derf( Sn ) )
    Bn       = sign_ * dexp( - Sn * Sn ) / sqrtBeta * ( 0.50d0 * invsqrtpi ) 
    term     = 0.50d0 * ( gp1bygm1 * q(5) / q(1) + dot_product( q(2:4), q(2:4) ) )
    F(1)     = q(1) * ( Un * An + Bn )
    F(2:4)   = F(1) * q(2:4) + q(5) * normal(1:3) * An
    F(5)     = F(1) * term + 0.50d0 * q(5) * Un * An
  end subroutine kfvs_splflx
  !> \brief Simple Minded Roe FDS
  subroutine roe(qL, qR, normal, F)
    implicit none
    integer i
    real (kind=8) :: h0L, h0R, Rfac, rhoRoe, aRoe, &
         h0Roe, drho, dp, VnL, VnR, VnRoe, qSqRoe, &
         dV, drhoTerm
    real (kind=8) :: fL(5), fR(5), df1(5), df2(5), df3(5)
    real (kind=8) :: dvel(3), velRoe(3)
    real (kind=8) :: qL(5), qR(5), normal(3), F(5)
!!!  Left and Right Flow values
    h0L = gbygm1 * qL(5) / qL(1) + 0.50d0 * &
         dot_product( qL(2:4), qL(2:4) )
    h0R = gbygm1 * qR(5) / qR(1) + 0.50d0 * &
         dot_product( qR(2:4), qR(2:4) )
    VnL = dot_product( qL(2:4), normal(1:3) )
    VnR = dot_product( qR(2:4), normal(1:3) )
!!!  Roe Averaged States
    Rfac = dsqrt( qR(1) / qL(1) )
    rhoRoe = Rfac * qL(1)
    velRoe = ( qL(2:4) + Rfac * qR(2:4) ) / ( 1.0d0 + Rfac )
    h0Roe = (h0L + Rfac * h0R) / ( 1.0d0 + Rfac )
    qSqRoe = dot_product(velRoe, velRoe)
    aRoe = dsqrt( gm1 * (h0Roe - 0.50d0 * qSqRoe) )
    VnRoe = (VnL + Rfac * VnR) / (1.0d0 + Rfac)
!!!  Differnetial values for characteristics
    drho = qR(1)   - qL(1)
    dvel = qR(2:4) - qL(2:4)
    dp   = qR(5)   - qL(5)
    dV = dot_product( dvel, normal )
    fL = normalFlux( qL, normal )
    fR = normalFlux( qR, normal )
!!! Form Inteface roe averaged fluxes df1 - df3
!!! Calculate df1
    drhoTerm = drho - dp / aRoe **2
    df1(1)   = dabs(VnRoe) * drhoTerm
    df1(2:4) = dabs(VnRoe) * ( velRoe * drhoTerm + &
         rhoRoe * (dvel - normal * dV) )
    df1(5)   = dabs(VnRoe) * ( 0.50d0 * qSqRoe * drhoTerm + &
         rhoRoe * ( dot_product(velRoe, dvel) - VnRoe * dV ) )
!!! Calculate df2
    df2(1)   = 0.50d0 * dabs(VnRoe + aRoe) * (dp / aRoe ** 2 + rhoRoe * dV / aRoe)
    df2(2:4) = df2(1) * ( velRoe + normal * aRoe )
    df2(5) = df2(1) * ( h0Roe + VnRoe * aRoe )
!!! Calculate df3
    df3(1)   = 0.50d0 * dabs(VnRoe - aRoe) * (dp / aRoe ** 2 - rhoRoe * dV / aRoe)
    df3(2:4) = df3(1) * ( velRoe - normal * aRoe )
    df3(5)   = df3(1) * ( h0Roe - VnRoe * aRoe )
!!! Sum fluxes 0.5 * (fL + fR - df1 - df2 -df3 ) to yield flux at the interface
    F = 0.50d0 * (fL + fR - df1 - df2 - df3)
  end subroutine roe
  !> \brief
  subroutine RKUpdate( U, U_old, dF_i, local_dt, my_pass, max_pass, rk_coeff )
    implicit none
    real (kind = 8) :: U( 5, size + padding ), U_old( 5, size + padding )
    real (kind = 8) :: dF_i( 5, size + padding ), local_dt( size )
    integer :: i, iloc(1), max_pass, my_pass
    real (kind=8) :: rk_coeff
    call ConvertToConservative( U )
    if( my_pass .eq. 1 ) U_old = U
    do i = 1, size
      U(:,i)    = U_old(:,i) - rk_coeff * local_dt(i) * dF_i(:,i)
    end do
    call ConvertToPrimitive( U )
    if( my_pass .eq. max_pass ) then
      iloc = maxloc( dabs( U(1,1:size) - U_old(1,1:size) ) )
      max_xyz = xyz(:,iloc(1))
      max_res(:) = dabs( U(:,iloc(1)) - U_old(:,iloc(1)) )
    end if
  end subroutine RKUpdate
  !> \brief
  subroutine setValues( U, q )
    implicit none
    real (kind = 8) :: U( 5, size + padding )
    real (kind = 8) :: q(5)
    integer i
    do i = 1, size
      U(1:5,i)      = q(1:5)
    end do
  end subroutine setValues
end module Flow

