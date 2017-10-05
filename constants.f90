module Constants
  implicit none

  !!! Math constants
  real (kind = 8) :: pi
  real (kind = 8) :: invpi 
  real (kind = 8) :: invsqrtpi 
  real (kind = 8) :: degToRad , radToDeg
  !!! Flow constants for ideal gas
  real (kind = 8) :: gama     = 1.40d0 
  real (kind = 8) :: gm1      = 0.40d0
  real (kind = 8) :: gbygm1   = 3.50d0
  real (kind = 8) :: invgama  = 1.0d0 / 1.40d0
  real (kind = 8) :: gp1bygm1 = 6.0d0 
  real (kind = 8) :: invgm1   = 2.50d0
  parameter ( pi = 3.141592653589793238462643383279502880d0 , & 
              invpi = 0.3183098861837906715377675267450287240d0 , &
              invsqrtpi = dsqrt( invpi ) , &
              degToRad = pi / 180.0d0 , &
              radToDeg = 180.0 / pi )

end module Constants


