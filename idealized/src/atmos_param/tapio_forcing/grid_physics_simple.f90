module grid_physics

  ! GRID_PHYSICS provides routines that can be combined to form a "physics" 
  ! package in an idealized GCM. 
  !
  ! The namelist GRID_PHYS.NML contains parameters for the parameterizations. 
  !
  ! Radiative transfer is parameterized as Newtonian relaxation toward the
  ! radiative-equilibrium state of a semi-gray atmosphere with surface temperature 
  ! specified in SURFACE_TEMPERATURE_FORCED. 
  !
  ! The surface temperature of the radiative-equilibrium state depends on the
  ! parameters DELH and TSFC_SP, the pole-to-equator temperature contrast and the
  ! surface temperature at the south pole, respectively. The extent of the frictional 
  ! boundary layer is determined by SIGMA_B. For the Newtonian relaxation 
  ! parameterizations, SIGMA_B also determines the height of the low-latitudinal layer 
  ! in which the relaxation coefficient is increased compared with the interior 
  ! atmospheric relaxation coefficient.  

use     constants_mod, only: kappa, cp_air, grav, rdgas, rvgas, stefan, &
                             hlv, tfreeze

use           fms_mod, only: error_mesg, fatal, file_exist,       &
                             open_namelist_file, check_nml_error, &
                             mpp_pe, mpp_root_pe, close_file,     &
                             write_version_number, stdlog

use  time_manager_mod, only: time_type
use  diag_manager_mod, only: register_diag_field, send_data

  implicit none

  
  real, private ::             &    ! constants and defaults for parameters  
       pi            =     3.141592653, &
       day           = 86400., &    ! [s]
       reference_sea_level_press = 1.e5, &       
       delh          =   120., &    ! pole-equator surface temperature contrast
       tsfc_sp       =   320., &    ! surface temperature at south pole
       phi0          =     0., &    ! heating offset (degrees)
       lat_tropic    =    32.092, & ! latitude within/outside surface temps are flat
       t_strat       =   200., &    ! "skin" temperature of atmosphere
       scale_height_ratio = 3.5, &  ! ratio of water vapor and pressure scale heights
       ka_days       =    50., &    ! latitude-independent part of Newtonian 
       ka,                     &    ! relaxation coefficient
       ks_days       =     4., &    ! Newtonian damping coefficient in low
       ks,                     &    ! latitudes near surface
       Cdrag         = 5.e-6        ! quadratic "drag coefficient" [1/m]

  real, public ::              &
       sigma_b       = 0.85       ! extent of `PBL' (also used to compute relaxation rates)


  logical, private ::          &    ! switches
       flatxt = .false.,       &
       flatt  = .false.
         
  namelist /grid_phys_list/    &
       tsfc_sp,       delh,            t_strat,             &
       ka_days,       ks_days,         sigma_b,             &
       reference_sea_level_press,      Cdrag,               &
       scale_height_ratio,                                  &
       flatt,         flatxt,          lat_tropic,          &
       phi0

  private grid_phys_list

   character(len=128) :: version='$Id: grid_physics_simple.f90 $'
   character(len=128) :: tag='homemade'
   character(len=14)  :: mod_name = "grid_physics"

! diagnostic IDs
   integer id_tdt, id_teq
   real :: missing_value = -1.e10
contains
  
  !-------------------------------------------------------------------------------

  subroutine grid_phys_init(axes, Time)
    ! reads namelist file with parameters and echoes parameter values

    integer, intent(in) :: axes(4)
    type(time_type), intent(in) :: Time 

    integer  unit, io, ierr

      if (file_exist('input.nml')) then
         unit = open_namelist_file ( )
         ierr=1; do while (ierr /= 0)
            read  (unit, nml=grid_phys_list, iostat=io, end=100)
            ierr = check_nml_error (io, 'grid_phys_list')
         enddo
  100     call close_file (unit)
      endif

!     ----- write version info and namelist to log file -----

      call write_version_number (version,tag)
      if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=grid_phys_list)

!     ---------------write info to screen--------------------

            
      if(mpp_pe()==mpp_root_pe()) then
         write(*, *) 'Parameters in radiative equilibrium surface temperature:'
         write(*,10) 'South pole surface temperature:             ', tsfc_sp, 'K'
         write(*,10) 'Equator-pole surface temperature contrast:  ', delh, 'K'
         write(*,10) 'Offset of temperature maximum from equator: ', phi0, 'deg'

         write(*, *) 'Relaxation toward radiative equilibrium state of a'
         write(*, *) 'semi-gray atmosphere with an isothermal stratosphere'
         write(*,10) 'at temperature t_strat =', t_strat, 'K.'
         write(*, *)
         write(*,10) 'Relaxation time scale in interior atmosphere:       ', &
              ka_days, 'days'
         write(*,10) 'Relaxation time scale in low latitudes near surface:', &
              ks_days, 'days'
      endif
      ka = 1./(ka_days * day)
      ks = 1./(ks_days * day)
      
      if(mpp_pe()==mpp_root_pe()) then
      write(*, *)
      write(*,10) 'Ratio of absorber scale height over pressure scale height:',&
           scale_height_ratio
      write(*, *)
   endif

   if(mpp_pe()==mpp_root_pe()) then
      write(*, *)
      write(*, *) 'Quadratic friction in boundary layer.'
      write(*,10) 'Extent of frictional layer: sigma_b =', sigma_b
      write(*,20) 'Drag coefficient:             Cdrag =', Cdrag, 'm**(-1)'
      write(*, *)
   endif
    
10  format(1x,a,1x,f8.2,1x,a)
20  format(1x,a,1x,e8.2,1x,a)

! register fields with diagnostic manager
      id_teq = register_diag_field( mod_name, 'teq',                         &
           axes(1:3), Time, 'equilibrium temperature', 'deg_K',              &
           missing_value=missing_value)

      id_tdt = register_diag_field( mod_name, 'tdt_total_grid_physics',      &
           axes(1:3), Time, 'total temperature tendency', 'deg_K/sec',       &
           missing_value=missing_value)

  end subroutine grid_phys_init
  
  !------------------------------------------------------------------

  subroutine compute_grid_physics( is,       ie,       js,       je,      &
                   Time, delta_t, lat,    dt_ug,    dt_vg,    dt_tg,      &   
                      dt_qg,   dt_trg,       ug,       vg,       tg,      &
                         qg,      trg,   p_half,   p_full,                &
                ug_previous,        vg_previous,        tg_previous,      &
                qg_previous,       trg_previous,    p_half_previous,      &
            p_full_previous)

    integer, intent(in)                     :: is, ie, js, je
    type(time_type), intent(in)            ::  Time
    real, intent(in)                       ::  delta_t
    real, intent(in), dimension (:,:,:)     :: ug, vg, tg, qg, p_half, p_full
    real, intent(in), dimension (:,:,:,:)   :: trg
    real, intent(in), dimension (:,:)       :: lat
    
    real, optional, dimension (:,:,:)       :: ug_previous, vg_previous, & 
         tg_previous, qg_previous, p_half_previous, p_full_previous
    real, optional, dimension (:,:,:,:)     :: trg_previous
    
    real, intent(inout), dimension(:,:,:)   :: dt_ug, dt_vg, dt_tg, dt_qg
    real, intent(inout), dimension(:,:,:,:) :: dt_trg


    if(present(ug_previous)) then
          call surface_drag(ug_previous, vg_previous, p_half_previous, &  
               p_full_previous, dt_ug, dt_vg)
          call diabatic_forcing(Time, is, js, tg_previous,                 &
               p_half_previous, p_full_previous, lat, dt_tg)
    else
       call error_mesg('compute_grid_physics',                             &
            'ug_previous is not present',fatal)
    end if

  end subroutine compute_grid_physics

  !---------------------------------------------------------------------
  
  subroutine compute_grid_adjustment(tg, qg, p_half, p_full)
    
    real, intent(in), dimension (:,:,:) :: tg, qg, p_full
    real, intent(in), dimension (:,:,:) :: p_half
    
    ! dummy routine
    
  end subroutine compute_grid_adjustment
  
  !-----------------------------------------------------------------
  
  subroutine surface_drag(ug, vg, p_half, p_full, dt_ug, dt_vg)
    
    real, intent(in), dimension (:,:,:)    :: ug, vg, p_full, p_half
    real, intent(inout), dimension(:,:,:)  :: dt_ug, dt_vg
    real, dimension(size(ug,1),size(ug,2)) :: sigma, sigma_norm, sigma_max
    real, dimension(size(ug, 1), size(ug, 2), size(ug, 3)) :: vg_norm
    integer :: k, num_level
    
    ! pseudo=parameter for number of levels
    num_level = size(ug,3)
    
    vg_norm = sqrt(ug**2 + vg**2)
    do k = 1, num_level
       sigma(:, :)  = p_full(:,:,k) / p_half(:,:,num_level+1)
       sigma_norm   = (sigma - sigma_b) / (1.0 - sigma_b)
       sigma_max    = max(sigma_norm, 0.0)
       dt_ug(:,:,k) = dt_ug(:,:,k) - Cdrag * sigma_max * vg_norm(:,:,k)*ug(:,:,k)
       dt_vg(:,:,k) = dt_vg(:,:,k) - Cdrag * sigma_max * vg_norm(:,:,k)*vg(:,:,k)
    end do
    
  end subroutine surface_drag

  !-----------------------------------------------------------------

  subroutine diabatic_forcing(Time, is, js, tg, p_half, p_full, lat, dt_tg)
     
    implicit none

    integer, intent(in)                   :: is, js
    type(time_type), intent(in)           :: Time
    ! temperature and pressure fields
    real, intent(in), dimension(:, :, :)  :: tg, p_full
    real, intent(in), dimension(:, :, :)  :: p_half 
    real, intent(in), dimension(:, :)     :: lat
    
    ! time tendency in temperature field
    real, intent(inout), dimension(:,:,:) :: dt_tg 
            
    ! local variables
    real, dimension(size(tg, 1), size(tg, 2), size(tg, 3)) ::              &
         heating_rate
 
    integer :: num_level
    logical :: used


    ! pseudo-parameters
    num_level = size(tg, 3)        

    if (size(lat,1) .ne. size(tg,1) &
         .or. size(lat,2) .ne. size(tg,2)) then
       call error_mesg('diabatic_forcing',                             &
            'mismatched argument list',fatal)
       ! terminate execution
    end if

    heating_rate = relax_to_gray_equilibrium(Time, tg, lat,      &
         p_full, p_half, is, js)
    
    dt_tg        = dt_tg + heating_rate

! send data to diagnostic manager
       if(id_tdt > 0) used = send_data(id_tdt, heating_rate, Time,      &
            is, js)

  end subroutine diabatic_forcing

  !-----------------------------------------------------------------------
  
  function relax_to_gray_equilibrium(Time, temp, lat, p_full, p_half, is, js)&
       result (heating_rate)
  
    ! computes heating rate from relaxation towards radiative
    ! equilibrium state of a gray atmosphere with an isothermal 
    ! stratosphere 
    
    implicit none

    type(time_type), intent(in)           :: Time

    integer, intent(in) :: is, js

    real, intent(in), dimension(:, :) ::                                   &
         lat
    real, intent(in)                  ::                                   &
         temp(:, :, :),                                                    &
         p_full(:, :, :),                                                  &
         p_half(:, :, :)
    
    ! local variables
    real, dimension(size(temp, 1), size(temp, 2)) ::                       & 
         sin_lat_2,                                                        &
         cos_lat_2,                                                        &
         cos_lat_4,                                                        &
         cos_lat_8,                                                        & 
         sigma,                                                            &
         sigma_norm,                                                       &
         sigma_max,                                                        &
         kt,                                                               &
         t_sfc,                                                            &
         t_strat_var,                                                      &
         optical_thickness

    real, dimension(size(temp, 1), size(temp, 2), size(temp, 3)) ::        &
         p_full_ref,                                                       &
         temp_eq,                                                          &
         heating_rate

    integer :: num_level, k 
    logical used

    num_level = size(temp, 3)
    cos_lat_2 = cos(lat)*cos(lat)
    cos_lat_4 = cos_lat_2**2
    cos_lat_8 = cos_lat_2**4

    if (size(lat,1) .ne. size(p_full,1) &
         .or. size(lat,2) .ne. size(p_full,2)) then
       call error_mesg('relax_to_gray_equilibrium',                        &
            'mismatched argument list',fatal)
       ! terminate execution
    end if

    ! normalized full level pressure and functions thereof
    p_full_ref        = p_full / reference_sea_level_press

    ! surface temperature in background state    
    t_sfc             = surface_temperature_forced(lat) 

    ! latitude-dependent optical thickness
    optical_thickness = (t_sfc / t_strat)**4 - 1.

    do k = num_level, 1, -1
       ! compute current sigma-level to get height-dependence of cooling coeff
       sigma      = p_full(:, :, k) / p_half(:, :, num_level+1)
       sigma_norm = (sigma - sigma_b) / (1.0 - sigma_b)
       sigma_max  = max(sigma_norm, 0.0) ! sigma_max = 0 outside `PBL'

       ! Newtonian damping coefficient
       ! changed by CW 2/24 for new runs
       ! was ^4 for all runs before 3/17/04, changed to ^8 on 3/17
       kt         = ka + (ks - ka) * sigma_max *                              &
            cos( lat - phi0*pi/180.0 )**8
!!$       kt         = ka + (ks - ka) * sigma_max * cos_lat_8(:,:) 
!       kt         = ka + (ks - ka) * sigma_max * cos_lat_4(:,:) 
       
       ! temperature in background state is radiative equilibrium 
       ! of a gray atmosphere (cf. Goody and Yung, pp. 392 ff.)
       temp_eq(:,:,k)    =                                                 &
            t_strat * (1. + optical_thickness                          &
                     * p_full_ref(:, :, k)**scale_height_ratio)**(1./4.)

       ! relax towards background state
       heating_rate(:, :, k) = ( temp_eq(:,:,k)  - temp(:, :, k) ) * kt
    end do
    if (id_teq > 0) used = send_data ( id_teq, temp_eq,  Time, is, js)

  end function relax_to_gray_equilibrium

  !--------------------------------------------------------------------------

  function surface_temperature(temp, p_full, p_half)

    ! estimate the surface temperature from the temperature on lowest
    ! full level by assuming that the potential temperature on the surface
    ! is equal to the potential temperature on the lowest full level

    implicit none

    real, intent(in) ::                                                    &
         temp(:, :, :),    & ! temperature field
         p_full(:, :, :),  & ! full-level pressure
         p_half(:, :, :)     ! half-level pressure

    real, dimension(size(temp, 1), size(temp, 2)) ::                       &
         surface_temperature

    integer :: num_level
    num_level = size(temp, 3)

    surface_temperature = temp(:, :, num_level)                            &
         * ( p_half(:, :, num_level+1) / p_full(:, :, num_level) )**kappa
    
  end function surface_temperature

  !-----------------------------------------------------------------------

  function surface_temperature_forced(lat)

    ! returns a prescribed surface temperature

    implicit none

    real, intent(in) ::                                                    &
         lat(:, :)

    real, dimension(size(lat, 1), size(lat, 2)) ::                         &
         surface_temperature_forced,                                       &
         latitude_dependence

    integer :: i

    real match
    
    latitude_dependence = cos(lat)**2 + 2. * sin(phi0*pi/180) * ( 1 + sin(lat) )
    surface_temperature_forced = tsfc_sp + delh * latitude_dependence

    ! make surface temperature flat inside/outside lat_tropic
    match = tsfc_sp + delh * cos(lat_tropic*pi/180.0)**2
    do i=1,size(lat,2)
       if(flatt .and. (abs(lat(1,i)*180./pi) < lat_tropic) ) &
            surface_temperature_forced(:, i) = match
       if(flatxt .and. (abs(lat(1,i)*180./pi) > lat_tropic) ) &
            surface_temperature_forced(:, i) = match
    enddo
       
  end function surface_temperature_forced
  
end module grid_physics


