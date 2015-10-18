
module constants_and_switches_mod

!-----------------------------------------------------
!
!  Defines useful constants for Earth in mks units.
!
!-----------------------------------------------------

implicit none
private

!-----------------------------------------------------------------------
!------------ physical constants -------------------

real, public, parameter :: RADIUS = 2575.e3      !  radius of titan (meters)
real, public, parameter :: OMEGA  = 4.56068e-6   !  rotation rate of planet (1/sec)
real, public, parameter :: GRAV   = 1.35         !  acceleration due to gravity (m/s2)
real, public, parameter :: RDGAS  = 296.8        !  gas constant for dry air (J/Kg/deg)
real, public, parameter :: KAPPA  = 2./7.        !  RDGAS / CP
real, public, parameter :: cp     = RDGAS/KAPPA  !  spec heat cap of dry air (J/kg/deg)

!------------ water vapor constants -----------------

real, public, parameter :: RVGAS = 518.3         !  gas constant for methane vapor (J/Kg/deg)
real, public, parameter :: DENS_H2O = 422.62     !  density of liquid methane (Kg/m3)
real, public, parameter :: HLV = 4.9e5           !  latent heat of evaporation (J/Kg)
real, public, parameter :: HLF = 5.868e4         !  latent heat of fusion (J/Kg)
real, public, parameter :: HLS = 5.6868e5        !  latent heat of sublimation (J/Kg)
real, public, parameter :: TFREEZE = 90.6        !  temp where fresh methane freezes (deg K)
real, public, parameter :: T_triple = 90.7       !  temperature at triple point 
real, public, parameter :: p_triple = 11696.0    !  pressure at triple point

!-- constants in Bolton (1980) water vapor formulae ---

real, public, parameter :: pa_to_mb = 0.01       ! Pascals to mb conversion
real, public, parameter :: const1 = 2840.0       ! used for saturation temperature formula
real, public, parameter :: const2 = 3.5          !
real, public, parameter :: const3 = 4.805        !
real, public, parameter :: const4 = 55.0         !
real, public, parameter :: const5 = 0.81         ! used for equivalent potential temperature formula
real, public, parameter :: const6 = 3376.0       !
real, public, parameter :: const7 = 2.54         !

!------------ miscellaneous constants -----------------

real, public, parameter :: STEFAN  =  5.6734e-8  !  Stefan-Boltzmann constant (W/m2/deg4)
real, public, parameter :: VONKARM =  0.40       !  Von Karman constant
real, public, parameter :: PI      =  3.14159265358979323846 ! is it enough?
!-----------------------------------------------------------------------


! switches
! in addition to temp and surface pressure, the following fields are read in

integer, public, parameter             ::    &
       vor_and_div     = 1,                  &
       u_and_v         = 2,                  &
       all             = 3

integer, public, parameter             ::    &
       FMS             = 1,                  &
       NCEP_spectral   = 2,                  &
       CCM3            = 3,                  &
       ERA40           = 4

real, public, parameter ::  reference_sea_level_pres = 1.467e5

real, public, parameter :: gspval = -9.e-14! missing data flag

end module constants_and_switches_mod

