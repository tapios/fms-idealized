module radiation_mod

! ==================================================================================
! ==================================================================================

   use fms_mod,               only: open_file, check_nml_error, close_file, &
                                    mpp_pe, mpp_root_pe, error_mesg, FATAL

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks

   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, get_time, &
                                    operator(+), operator(-), operator(/=)

   use        transforms_mod, only: get_deg_lat

!==================================================================================
implicit none
private
!==================================================================================

! version information 
 
character(len=128) :: version='$Id: radiation.f90 $'
character(len=128) :: tag='homemade'

!==================================================================================

! public interfaces

public :: radiation_init, radiation_down, radiation_up, radiation_end
!==================================================================================

! = Orbital forcing: compute insolation based on specified paleo year. =
! Added by ian: feb 2009.
! Code modification overview: (1) In radiation_init, compute values of orbital
! parameters based on paleo year (orb_year). (2) In radiation_down, compute TOA
! downward shortwave radiation (solar) using these parameters.
! ===

! module variables
logical :: initialized =.false.
real    :: secs_per_day    = 86400.
real    :: solar_constant  = 1365. ! W/m^2
real    :: day             = 0.0
! modif omp: winter/summer hemisphere
real    :: del_sol         = 1.4 ! only used if perpetual_equinox=.true.
real    :: odp             = 1.0 ! odp = optical depth parameter, added by rwills
real    :: lw_tau0_eqtr    = 6.0
real    :: lw_tau0_pole    = 1.5
real    :: atm_abs         = 0.2
real    :: sw_diff         = 0.0
real    :: lw_linear_frac  = 0.1
real    :: albedo_value    = 0.3
real    :: lw_tau_exponent  = 4.0 
real    :: sw_tau_exponent  = 4.0 

logical :: fixed_day          = .false.
real    :: fixed_day_value    = 0.       ! time in days from spring equinox 
real    :: days_in_year       = 360 ! how many model days in solar year
logical :: perpetual_equinox  = .false. ! whether to use perpetual equinox rather than seasonally varying insolation
logical :: annual_mean        = .true. ! whether to use annual mean

! Orbital forcing: Either specify orbital_year, in which case orbital parameters are computed,
! or leave orbital_year blank and specify orbital parameters orb_ecc, orb_obl, and orb_long_perh.
! If both are left blank, orbital_year=1990 is assumed. Commented text below includes [1990 values].
! ian: Feb 2009
real    :: R_UNDEF         = 1.e36    ! undefined real
real    :: orb_year        = 1.e36    ! calendar year (AD) used to determine orbital parameters [1990]
real    :: orb_ecc         = 1.e36    ! eccentricity of Earth's orbit (unitless) [0.017]
real    :: orb_obl         = 1.e36    ! obliquity of Earth's orbit (degrees) [23.45]
real    :: orb_long_perh   = 1.e36    ! angle in Earth's orbit between equinox and perihelion (degrees) [281.4]
! end orbital forcing

real, allocatable, dimension(:,:)   :: insolation, lw_tau0, sw_tau0, h0, h0_arg
real, allocatable, dimension(:,:)   :: b_surf
real, allocatable, dimension(:,:,:) :: b, tdt_rad, tdt_solar
real, allocatable, dimension(:,:,:) :: lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux
real, allocatable, dimension(:,:,:) :: lw_dtrans
real, allocatable, dimension(:,:,:) :: lw_tau, sw_tau
real, allocatable, dimension(:,:)   :: net_lw_surf, net_lw_toa, sw_down_toa
real, allocatable, dimension(:,:)   :: albedo
real, allocatable, dimension(:,:)   :: h0_ann, h0_arg_ann, insolation_ann, insolation_annual_mean


real, save :: pi, deg_to_rad , rad_to_deg

namelist/radiation_nml/ albedo_value, lw_linear_frac, solar_constant, del_sol, odp, &
                        lw_tau0_eqtr, lw_tau0_pole, atm_abs, lw_tau_exponent, sw_tau_exponent, &
                        perpetual_equinox, annual_mean, fixed_day, fixed_day_value, days_in_year, &
                        orb_year, orb_ecc, orb_obl, orb_long_perh
                      
!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_lwup_toa, id_swdn_sfc, id_swdn_toa, id_lwdn_sfc, id_lwup_sfc, id_net_lw_surf, &
           id_tdt_rad, id_flux_rad, id_flux_lw, id_flux_sw, id_tdt_solar

character(len=10), parameter :: mod_name = 'two_stream'

real :: missing_value = -999.


contains



! ==================================================================================
! ==================================================================================


subroutine radiation_init(is, ie, js, je, num_levels, axes, Time)

!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit
!-----------------------------------------------------------------------------------------
! === Start insolation parameters (for series approximation of orbital parameters as function of paleo year) ===
! obliquity
real, dimension(5) :: obamp=(/ -2462.22,-857.32,-629.32,-414.28,-311.76 /)
real, dimension(5) :: obrate=(/ 31.61,32.62,24.17,31.98,44.83 /)
real, dimension(5) :: obphas=(/ 251.90,280.83,128.31,292.73,15.37 /)
! eccentricity
real, dimension(6) :: ecamp=(/ 0.018608,0.016275,-0.013007,0.009888,-0.003367,0.003331 /)
real, dimension(6) :: ecrate=(/ 4.20721,7.34609,17.85726,17.22055,16.84673,5.1990790 /)
real, dimension(6) :: ecphas=(/ 28.6201,193.7888,308.3070,320.1996,279.3770,87.1950 /)
! longitude of perihelion
real :: mvamp=7391.02
real :: mvrate=31.61
real :: mvphas=251.90
! other
real :: yr1950 ! years after 1950 [negative for paleo year]
real :: psecdeg = 1.0/3600.0 ! arc sec to deg conversion
real :: cossum ! for summing a cosine series
real :: sinsum ! for summing a sine series
real :: fvelp ! intermediate variable for computing orb_long_perh
real    :: time_in_ann, beta, lambda0_ann, lambda_m_ann, lambda_ann, delta_ann  ! for intermediate steps during insolation calculation
integer :: i ! index
integer :: j, k, day_ann
integer :: lat_max
logical :: calc_orb =.true. ! whether orbital parameters need to be calculated

real, allocatable, dimension(:,:)   ::   lat                 ! latitude in radians 
real, allocatable, dimension(:)     ::   deg_lat
allocate(lat(is:ie, js:je))
allocate(deg_lat(js:je))
call get_deg_lat(deg_lat)

pi = 4.0*atan(1.)

do j=js,je
   lat(:,j) = deg_lat(j)*pi/180.
enddo

! === End insolation parameters ===

! read namelist and copy to logfile

unit = open_file ('input.nml', action='read')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=radiation_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'radiation_nml')
enddo
10 call close_file (unit)

unit = open_file ('logfile.out', action='append')
if ( mpp_pe() == 0 ) then
  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
  write (unit, nml=radiation_nml)
endif
call close_file (unit)

pi    = 4.0*atan(1.)
deg_to_rad = 2.*pi/360.
rad_to_deg = 360./2./pi

                ! === START INSOLATION CALCULATION === !

!
!      Find orbital parameters === (ecc, obliq, long of perihelion)
!                            Ian, Feb 2009

!########################################################################!
! EITHER DEFINE ECCENTRICITY, OBLIQUITY, PERIHELION OR SET THEM AT GIVEN YEAR ORB_YEAR
!########################################################################!


! parse input: find out whether orbital parameters need to be calculated
if ( orb_ecc .ne. R_UNDEF .and. orb_obl .ne. R_UNDEF &
           .and. orb_long_perh .ne. R_UNDEF &
           .and. orb_year == R_UNDEF) then ! all 3 orbital parameters are specified, orb_year is not
  calc_orb=.false. ! do not calculate orbital parameters
else if ( orb_ecc == R_UNDEF .and. orb_obl == R_UNDEF &
             .and. orb_long_perh == R_UNDEF ) then ! no orbital parameters are explicitly specified
  calc_orb=.true. ! calculate orbital parameters parameters based on orb_year
  if ( orb_year == R_UNDEF ) then ! orbital year not specified
    orb_year=1990 ! use default orbital year
  endif
else ! 1-2 orbital parameters are specified, or orb_year AND some parameters are both specified
  call error_mesg('radiation', 'insolation: must specify EITHER orb_year OR&
             & orbital parameters (orb_ecc,orb_obl,orb_long_perh)', FATAL)
end if

!########################################################################!
! NO ORBITAL PARAMATERS ARE NEEDED FOR IDEALIZED EQUINOX INSOLATION
!########################################################################!

if ( perpetual_equinox ) calc_orb=.false. ! no need to calculate orbital parameters if running perpetual equinox

!########################################################################!
! COMPUTE ORBITAL PARAMETERS AT GIVEN YEAR ORB_YEAR IF ECCENTRICITY AND OBLIQUITY AND PERIHELION ARE UNDEFINED
!########################################################################!

if ( calc_orb ) then ! calculate orbital parameters
  ! = This uses highly truncated series approximation to Berger (1978) =
  ! Drawing on Berger (1978) data table from NCDC and the CCSM3 shr_orb_mod.F90

  yr1950=orb_year-1950

  ! orb_obl: obliquity in degrees
  orb_obl=23.320556
  do i = 1, 5
    orb_obl=orb_obl+obamp(i)*psecdeg*cos((obrate(i)*psecdeg*yr1950+ &
    &       obphas(i))*deg_to_rad)
  end do

  ! orb_ecc: eccentricity (dimensionless)
  cossum=0.
  sinsum=0.
  do i = 1, 6
    cossum = cossum+ecamp(i)*cos((ecrate(i)*psecdeg*yr1950+ecphas(i))*deg_to_rad)
    sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*yr1950+ecphas(i))*deg_to_rad)
  end do
  orb_ecc=sqrt(cossum**2+sinsum**2)

  ! orb_long_perh: moving vernal equinox longitude of perihelion
  ! note: this is omega+180 deg (see lambda definition, Berger 1978 appendix)
  ! fvelp: (intermediate variable in calculation) fixed vernal equinox long of perh
  fvelp=atan(sinsum/cossum)
  if ( cossum < 0. ) then
    fvelp=fvelp+pi
  else if ( cossum > 0. .and. sinsum < 0. ) then
    fvelp=fvelp+2*pi
  end if
  orb_long_perh=fvelp/deg_to_rad+ 50.439273*psecdeg*yr1950 + 3.392506+180.+ &
  &             mvamp*psecdeg*sin((mvrate*psecdeg*yr1950 + mvphas)*deg_to_rad)
end if

! Write orbital parameter values to log file.
! It would repeat this for each processor, writing it nproc times.
! Write it out only for one processor using 'if ( mpp...' .
if ( mpp_pe() == mpp_root_pe() ) then
  if ( calc_orb ) then
    write(*,'(a,f7.0)') 'INSOLATION: Calculated orbital parameters for year (AD) ',&
    &           orb_year
    write(*,'(a,f11.6)') 'INSOLATION: eccentricity = ', orb_ecc
    write(*,'(a,f11.6)') 'INSOLATION: obliquity = ', orb_obl
    write(*,'(a,f11.6)') 'INSOLATION: longitude of perihelion = ', orb_long_perh
  elseif ( perpetual_equinox ) then
    write(*,*) 'INSOLATION: running in perpetual equinox'
  else
    write(*,*)  'INSOLATION: Orbital parameters specified in input namelist'
    write(*,'(a,f11.6)') 'INSOLATION: eccentricity = ', orb_ecc
    write(*,'(a,f11.6)') 'INSOLATION: obliquity = ', orb_obl
    write(*,'(a,f11.6)') 'INSOLATION: longitude of perihelion = ', orb_long_perh
  end if
end if

!########################################################################!
! CALCULATE INSOLATION PROFILE THAT CORRESPONDS TO ORBITAL CONFIGURATIONS SPECIFIED ABOVE:
! (i)    ANNUAL-MEAN PROFILE                      IF FIXED_DAY = .FALSE.
! (ii)   PROFILE AT A FIXED DAY = FIXED_DAY_VALUE IF FIXED_DAY = .TRUE.                                              
!########################################################################!
!                            Ian, Feb 2009; Xavier, May 2012
allocate (h0_ann                    (ie-is+1, je-js+1))
allocate (h0_arg_ann                (ie-is+1, je-js+1))
allocate (insolation_ann            (ie-is+1, je-js+1))
allocate (insolation_annual_mean    (ie-is+1, je-js+1))

insolation_annual_mean = 0

do k = 1,  days_in_year

  ! === START ANNUAL INSOLATION CALCULATION === ! 
  ! Calculate annual mean insolation as function of latitude, and orbital parameters (Ian, Feb. 2009)
  
  ! The calendar is referenced to the vernal equinox (day=0)
  
  ! lambda (or solar longitude) is the angular distance along Earth's orbit measured from vernal equinox (21 March). 
  ! Estimate lambda from calendar day using an approximation from Berger 1978 section 3.

  if (fixed_day) then
     time_in_ann = fixed_day_value * secs_per_day
  else
     time_in_ann = 1.*k*secs_per_day  
  endif   

  beta=sqrt(1.-orb_ecc**2)
  lambda0_ann  = 2.*pi*time_in_ann/(days_in_year*secs_per_day)
  lambda_m_ann = lambda0_ann &
         - 2.*( (1./2.*orb_ecc+1./8.*orb_ecc**3)*(1.+beta)*sin(-orb_long_perh*deg_to_rad) &
         - 1./4.*orb_ecc**2*(1./2.+beta)*sin(-2.*orb_long_perh*deg_to_rad) &
         + 1./8.*orb_ecc**3*(1./3.+beta)*(sin(-3.*orb_long_perh*deg_to_rad)) )
  lambda_ann = lambda_m_ann &
         + (2.*orb_ecc-1./4.*orb_ecc**3)*sin(lambda_m_ann-orb_long_perh*deg_to_rad) &
         + (5./4.)*orb_ecc**2*sin(2*(lambda_m_ann-orb_long_perh*deg_to_rad)) &
         + (13./12.)*orb_ecc**3*sin(3*(lambda_m_ann-orb_long_perh*deg_to_rad))
  
  delta_ann = asin(sin(orb_obl*deg_to_rad)*sin(lambda_ann)); ! declination of the sun

  ! hour angle 
  h0_arg_ann   = -tan(lat(:,:)) * tan(delta_ann)
  where (h0_arg_ann .le. -1.)
     h0_ann = pi
  elsewhere (h0_arg_ann .ge. 1.)
     h0_ann = 0.
  elsewhere
     h0_ann = acos(h0_arg_ann)
  endwhere
  
  ! Insolation: Berger 1978 eq (10)
  insolation_ann=solar_constant/pi*(1.+orb_ecc*cos(lambda_ann-orb_long_perh*deg_to_rad))**2 / &
  &   (1.-orb_ecc**2)**2 * ( h0_ann*sin(lat(:,:))*sin(delta_ann) + &
  &   cos(lat(:,:))*cos(delta_ann)*sin(h0_ann) )

  insolation_annual_mean = insolation_annual_mean + insolation_ann / days_in_year 

enddo
 
               ! === END INSOLATION CALCULATION ===  !

initialized = .true.

allocate (b                (ie-is+1, je-js+1, num_levels))
allocate (tdt_rad          (ie-is+1, je-js+1, num_levels))
allocate (tdt_solar        (ie-is+1, je-js+1, num_levels))

allocate (sw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_flux          (ie-is+1, je-js+1, num_levels+1))

allocate (lw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (lw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (lw_flux          (ie-is+1, je-js+1, num_levels+1))

allocate (rad_flux         (ie-is+1, je-js+1, num_levels+1))

allocate (b_surf           (ie-is+1, je-js+1))
allocate (net_lw_toa       (ie-is+1, je-js+1))
allocate (net_lw_surf      (ie-is+1, je-js+1))
allocate (sw_down_toa      (ie-is+1, je-js+1))
allocate (albedo           (ie-is+1, je-js+1))

allocate (lw_dtrans        (ie-is+1, je-js+1, num_levels))
allocate (lw_tau           (ie-is+1, je-js+1, num_levels+1))
allocate (sw_tau           (ie-is+1, je-js+1, num_levels+1))

allocate (h0               (ie-is+1, je-js+1)) 
allocate (h0_arg           (ie-is+1, je-js+1)) 
allocate (insolation       (ie-is+1, je-js+1)) 

allocate (lw_tau0          (ie-is+1, je-js+1))
allocate (sw_tau0          (ie-is+1, je-js+1))


!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

    id_lwup_toa = &
    register_diag_field ( mod_name, 'lwup_toa', axes(1:2), Time, &
               'outgoing longwave radiation', &
               'watts/m2', missing_value=missing_value               )
    id_swdn_sfc = &
    register_diag_field ( mod_name, 'swdn_sfc', axes(1:2), Time, &
               'SW flux down at surface', &
               'watts/m2', missing_value=missing_value               )
    id_swdn_toa = &
    register_diag_field ( mod_name, 'swdn_toa', axes(1:2), Time, &
               'SW flux down at TOA', &
               'watts/m2', missing_value=missing_value               )
    id_lwup_sfc = &
    register_diag_field ( mod_name, 'lwup_sfc', axes(1:2), Time, &
               'LW flux up at surface', &
               'watts/m2', missing_value=missing_value               )

    id_lwdn_sfc = &
    register_diag_field ( mod_name, 'lwdn_sfc', axes(1:2), Time, &
               'LW flux down at surface', &
               'watts/m2', missing_value=missing_value               )

    id_net_lw_surf = &
    register_diag_field ( mod_name, 'net_lw_surf', axes(1:2), Time, &
               'Net upward LW flux at surface', &
               'W/m2', missing_value=missing_value               )

    id_tdt_rad = &
        register_diag_field ( mod_name, 'tdt_rad', axes(1:3), Time, &
               'Temperature tendency due to radiation', &
               'K/s', missing_value=missing_value               )

    id_tdt_solar = &
        register_diag_field ( mod_name, 'tdt_solar', axes(1:3), Time, &
               'Temperature tendency due to SW radiation', &
               'K/s', missing_value=missing_value               )

    id_flux_rad = &
        register_diag_field ( mod_name, 'flux_rad', axes(half), Time, &
               'Total radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_lw = &
        register_diag_field ( mod_name, 'flux_lw', axes(half), Time, &
               'Net longwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw = &
        register_diag_field ( mod_name, 'flux_sw', axes(half), Time, &
               'Net shortwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )

return
end subroutine radiation_init


! ==================================================================================

subroutine radiation_down (is, js, Time_diag, lat, p_half, t,         &
                           net_surf_sw_down, surf_lw_down)

! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

integer, intent(in)                 :: is, js
integer                             :: seconds, days
real                                :: time_in
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat
real, intent(out) , dimension(:,:)  :: net_surf_sw_down
real, intent(out) , dimension(:,:)  :: surf_lw_down
real, intent(in) , dimension(:,:,:) :: t, p_half

integer :: i, j, k, n

logical :: used

real    :: beta, lambda0, lambda_m, lambda, delta  ! for intermediate steps during insolation calculation

! ==================================================================================

             ! === Start TOA downwelling shortwave calculation ===  !

!########################################################################!
! IDEALIZED POLYNOMIAL REPRESENTATION FOR EQUINOX INSOLATION THAT APPROXIMATES EARTH INSOLATION IN THE ANNUAL-MEAN
!########################################################################!
if ( perpetual_equinox ) then
! Paul, 2006
  insolation = 0.25*solar_constant*(1.0 + del_sol* (1. - 3.*sin(lat(:,:))**2)/4. )

!########################################################################!
! ANNUAL-MEAN INSOLATION COMPUTED FROM BERGER ANALYTICAL FORMULATION FOR DAILY-MEAN INSOLATION
!########################################################################!
else if (annual_mean) then
! Ian, 2009; Xavier, 2012
  insolation = insolation_annual_mean(:,:) 

!########################################################################!
! INSOLATION FOR DAY = FIXED_DAY_VALUE COMPUTED FROM BERGER ANALYTICAL FORMULATION FOR DAILY-MEAN INSOLATION
!########################################################################!
else if (fixed_day) then
! Ian, 2009; Xavier, 2012
  insolation = insolation_annual_mean(:,:)

!########################################################################!
! BERGER ANALYTICAL FORMULATION FOR DAILY-MEAN INSOLATION
!########################################################################!
else ! use seasonally varying insolation based on orb_year

  call get_time(Time_diag, seconds, days)
  time_in = seconds + secs_per_day*days

! Calculate daily mean insolation as function of latitude, calendar day of year, and orbital parameters (Ian, Feb. 2009)
! The calendar is referenced to the vernal equinox (day=0)
! lambda (or solar longitude) is the angular distance along Earth's orbit measured from vernal equinox (21 March). 
! Estimate lambda from calendar day using an approximation from Berger 1978 section 3.
  beta=sqrt(1.-orb_ecc**2)
  lambda0  = 2.*pi*time_in/(days_in_year*secs_per_day)
  lambda_m = lambda0 &
         -2.*( (1./2.*orb_ecc+1./8.*orb_ecc**3)*(1.+beta)*sin(-orb_long_perh*deg_to_rad) &
         -1./4.*orb_ecc**2*(1./2.+beta)*sin(-2.*orb_long_perh*deg_to_rad) &
         +1./8.*orb_ecc**3*(1./3.+beta)*(sin(-3.*orb_long_perh*deg_to_rad)) )
  lambda = lambda_m & 
        + (2.*orb_ecc-1./4.*orb_ecc**3)*sin(lambda_m-orb_long_perh*deg_to_rad) &
        + (5./4.)*orb_ecc**2*sin(2*(lambda_m-orb_long_perh*deg_to_rad)) &
        + (13./12.)*orb_ecc**3*sin(3*(lambda_m-orb_long_perh*deg_to_rad))
  
  delta = asin(sin(orb_obl*deg_to_rad)*sin(lambda)); ! declination of the sun

  ! hour angle 
  h0_arg  = -tan(lat(:,:)) * tan(delta)
  where (h0_arg .le. -1.)
     h0 = pi
  elsewhere (h0_arg .ge. 1.)
     h0 = 0.
  elsewhere
     h0 = acos(h0_arg)
  endwhere
  
  ! Insolation: Berger 1978 eq (10)
  insolation = solar_constant/pi * (1.+orb_ecc*cos(lambda-orb_long_perh*deg_to_rad))**2 / &
  &   (1.-orb_ecc**2)**2 * ( h0*sin(lat(:,:))*sin(delta) + &
  &   cos(lat(:,:))*cos(delta)*sin(h0) )
  
end if

             ! === End TOA downwelling shortwave calculation ===  !

! ==================================================================================

lw_tau0 = lw_tau0_eqtr + (lw_tau0_pole - lw_tau0_eqtr)*sin(lat(:,:))**2
sw_tau0 = (1.0 - sw_diff*sin(lat(:,:))**2)*atm_abs
lw_tau0 = lw_tau0*odp

! set a constant albedo for testing
albedo(:,:) = albedo_value

n = size(t,3)

do k = 1, n+1
   lw_tau(:,:,k)  = lw_tau0 * (lw_linear_frac * p_half(:,:,k)/pstd_mks      &
       + (1.0 - lw_linear_frac) * (p_half(:,:,k)/pstd_mks )**lw_tau_exponent)
   sw_tau(:,:,k)  = sw_tau0 * (p_half(:,:,k)/pstd_mks )**sw_tau_exponent
end do

b = stefan*t**4

do k = 1, n
  lw_dtrans(:,:,k) = exp(-(lw_tau(:,:,k+1)-lw_tau(:,:,k)))
end do

lw_down(:,:,1) = 0.0
do k = 1,n
    lw_down(:,:,k+1) = lw_down(:,:,k)*lw_dtrans(:,:,k) + b(:,:,k)*(1.0 - lw_dtrans(:,:,k))
end do

do k = 1,n+1
   sw_down(:,:,k) = insolation(:,:)*exp(-sw_tau(:,:,k))
end do

surf_lw_down     = lw_down(:,:,n+1)
sw_down_toa      = sw_down(:,:,1)
net_surf_sw_down = sw_down(:,:,n+1)*(1. - albedo(:,:))

    
!------- downward lw flux surface -------
      if ( id_lwdn_sfc > 0 ) then
          used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
      endif
!------- incoming sw flux toa -------
      if ( id_swdn_toa > 0 ) then
          used = send_data ( id_swdn_toa, sw_down_toa, Time_diag)
      endif
!------- downward sw flux surface -------
      if ( id_swdn_sfc > 0 ) then
          used = send_data ( id_swdn_sfc, net_surf_sw_down, Time_diag)
      endif

return
end subroutine radiation_down

! ==================================================================================

subroutine radiation_up (is, js, Time_diag, lat, p_half, t_surf, t, tdt)

! Now complete the radiation calculation by computing the upward and net fluxes.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat
real, intent(in) , dimension(:,:)   :: t_surf
real, intent(in) , dimension(:,:,:) :: t, p_half
real, intent(inout), dimension(:,:,:) :: tdt


integer :: i, j, k, n

logical :: used

n = size(t,3)
b_surf = stefan*t_surf**4

! compute upward longwave flux by integrating upward
lw_up(:,:,n+1) = b_surf
do k = n,1,-1
   lw_up(:,:,k) = lw_up(:,:,k+1)*lw_dtrans(:,:,k) + b(:,:,k)*(1.0 - lw_dtrans(:,:,k))
end do

! compute shortwave flux
do k = 1,n+1
   sw_up(:,:,k) = albedo(:,:)*sw_down(:,:,n+1)
end do

! net fluxes (positive up)
lw_flux  = lw_up - lw_down
sw_flux  = sw_up - sw_down
rad_flux = lw_flux + sw_flux

do k = 1,n
  tdt_rad(:,:,k) = ( rad_flux(:,:,k+1) - rad_flux(:,:,k) )  &
             *grav/(cp_air*(p_half(:,:,k+1)-p_half(:,:,k)))
  tdt_solar(:,:,k) = ( sw_flux(:,:,k+1) - sw_flux(:,:,k))  &
             *grav/(cp_air*(p_half(:,:,k+1)-p_half(:,:,k)))
  tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
end do


net_lw_toa  = lw_up(:,:,1)
net_lw_surf = lw_flux(:, :, n+1)        

!------- outgoing lw flux toa (olr) -------
      if ( id_lwup_toa > 0 ) then
          used = send_data ( id_lwup_toa, net_lw_toa, Time_diag)
      endif
!------- upward lw flux surface -------
      if ( id_lwup_sfc > 0 ) then
          used = send_data ( id_lwup_sfc, b_surf, Time_diag)
      endif
!------- net upward lw flux surface -------
      if ( id_net_lw_surf > 0 ) then
          used = send_data ( id_net_lw_surf, net_lw_surf, Time_diag)
      endif
!------- temperature tendency due to radiation ------------
      if ( id_tdt_rad > 0 ) then
         used = send_data ( id_tdt_rad, tdt_rad, Time_diag)
      endif
!------- temperature tendency due to solar radiation ------------
      if ( id_tdt_solar > 0 ) then
         used = send_data ( id_tdt_solar, tdt_solar, Time_diag)
      endif
!------- total radiative flux (at half levels) -----------
      if ( id_flux_rad > 0 ) then
         used = send_data ( id_flux_rad, rad_flux, Time_diag)
      endif
!------- longwave radiative flux (at half levels) --------
      if ( id_flux_lw > 0 ) then 
         used = send_data ( id_flux_lw, lw_flux, Time_diag)
      endif
!------- shortwave radiative flux (at half levels) --------
      if ( id_flux_sw > 0 ) then
         used = send_data ( id_flux_sw, sw_flux, Time_diag)
      endif

return
end subroutine radiation_up

! ==================================================================================

                                                                      
subroutine radiation_end
                                                                                                      
deallocate (b, tdt_rad, tdt_solar) 
deallocate (lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux)
deallocate (b_surf, net_lw_toa, net_lw_surf, sw_down_toa, albedo)
deallocate (lw_dtrans, lw_tau, sw_tau)
deallocate (insolation, lw_tau0, sw_tau0, h0, h0_arg, h0_ann, h0_arg_ann) 

end subroutine radiation_end

! ==================================================================================

end module radiation_mod




