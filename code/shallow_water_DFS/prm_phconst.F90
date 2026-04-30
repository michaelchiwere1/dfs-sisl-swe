module prm_phconst

  real(8),parameter :: ER = 6.37122d6     !! Radius of the Earth [m]
  real(8),parameter :: ER_INV = 1.0d0/ER
  
  real(8),parameter :: PI = 3.1415926535897932d0
  real(8),parameter :: PI2 = 2.0d0*PI
  real(8),parameter :: OMG = 7.292d-5     !! Angular velocity of the earth’s rotation [/s]
  real(8),parameter :: GRAV = 9.80616d0   !! Gravity acceleration [m/s**2]
  
end module prm_phconst
