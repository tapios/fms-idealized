module atmosphere_mod

use                  fms_mod, only: set_domain, write_version_number, field_size, file_exist, &
                                    mpp_pe, mpp_root_pe, error_mesg, FATAL, read_data, &
                                    write_data, nullify_domain, open_namelist_file, &
                                    check_nml_error, stdlog, close_file

! grid_physics added by fridoo 01 feb 2012. Needed for tapio_forcing, used for the dry model 
use             grid_physics, only: grid_phys_init, compute_grid_physics,&
                                    surface_temperature
 

use            constants_mod, only: grav, pi, dens_h2o

use           transforms_mod, only: trans_grid_to_spherical, trans_spherical_to_grid, &
                                    get_deg_lat, get_deg_lon, get_grid_boundaries, grid_domain, &
                                    spectral_domain, get_grid_domain, get_lon_max, get_lat_max, &
                                    area_weighted_global_mean ! added fridoo Sept 2012 for water conservation scheme (tmerlis)

use         global_integral_mod, only: mass_weighted_global_integral ! added fridoo Sept 2012 for water conservation scheme (tmerlis)

use         time_manager_mod, only: time_type, set_time, get_time, &
                                    operator( + ), operator( - ), operator( < )

use     press_and_geopot_mod, only: compute_pressures_and_heights

use    spectral_dynamics_mod, only: spectral_dynamics_init, spectral_dynamics, spectral_dynamics_end, &
                                    get_num_levels, get_axis_id, spectral_diagnostics, get_initial_fields, &
                                    get_surf_geopotential, diffuse_surf_water 

use          tracer_type_mod, only: tracer_type

use           hs_forcing_mod, only: hs_forcing_init, hs_forcing

use        field_manager_mod, only: MODEL_ATMOS

use       tracer_manager_mod, only: get_number_tracers 

use         surface_flux_mod, only: surface_flux
use     vert_turb_driver_mod, only: vert_turb_driver_init,                    &
                                    vert_turb_driver, vert_turb_driver_end
use            vert_diff_mod, only: gcm_vert_diff_init, gcm_vert_diff_down, &
                                    gcm_vert_diff_up, gcm_vert_diff_end, & 
                                    surf_diff_type, gcm_vert_diff
use            radiation_mod, only: radiation_init, radiation_down, radiation_up, &
                                    radiation_end

use       dry_convection_mod, only: dry_convection_init, dry_convection

use          mixed_layer_mod, only: mixed_layer_init, mixed_layer, mixed_layer_end
use          lscale_cond_mod, only: lscale_cond_init, lscale_cond 
use  qe_moist_convection_mod, only: qe_moist_convection_init, &
                                    qe_moist_convection,   &
                                    qe_moist_convection_end
use         diag_manager_mod, only: register_diag_field, send_data


implicit none
private 
!=================================================================================================================================

character(len=128) :: version= &
'$Id: atmosphere.f90,v 13.0 2006/03/28 21:17:28 fms Exp $'
      
character(len=128) :: tagname= &
'$Name: latest $'
character(len=10), parameter :: mod_name='atmosphere'

!=================================================================================================================================

public :: atmosphere_init, atmosphere, atmosphere_end

!=================================================================================================================================

logical :: turb = .false.
logical :: ldry_convection = .false.
logical :: do_virtual = .false. ! whether virtual temp used in gcm_vert_diff
logical :: lwet_convection = .false.
logical :: two_stream = .false.
logical :: mixed_layer_bc = .false.

! add fridoo 1 feb 2012
logical :: dry_model=.false.
logical :: tapio_forcing=.false.
logical :: hs=.false.
! end fridoo

! add fridoo Sept 2012 (tmerlis)
logical :: atmos_water_correction = .false.
! end fridoo

real :: roughness_heat = 0.05
real :: roughness_moist = 0.05
real :: roughness_mom = 0.05

! added by LJJ
logical :: bucket = .true. 
real :: init_bucket_depth = 20. ! default initial bucket depth in m LJJ
real :: init_bucket_depth_land = 20. 
real :: land_left       = 90.
real :: land_right      = 270.
real :: land_bottom     = -10.
real :: land_top        = 10.
real :: max_bucket_depth_land = 1000. ! default large value
real :: robert_bucket = 0.04   ! default robert coefficient for bucket depth LJJ
real :: raw_bucket = 0.53       ! default raw coefficient for bucket depth LJJ
real :: damping_coeff_bucket = 200. ! default damping coefficient for diffusing of surface water - default is no diffusion. [degrees/year] 
! end addition by LJJ

namelist/atmosphere_nml/ turb, ldry_convection, lwet_convection, roughness_heat, two_stream, mixed_layer_bc,               &
                         roughness_moist, roughness_mom, do_virtual, tapio_forcing, dry_model, hs, atmos_water_correction, &
                         bucket, init_bucket_depth, init_bucket_depth_land, land_left, land_right, land_bottom, land_top,  &
                         robert_bucket, raw_bucket, max_bucket_depth_land, damping_coeff_bucket

!===================================================================
integer, parameter :: num_time_levels = 2
integer :: is, ie, js, je, num_levels, num_tracers, nhum


real, allocatable, dimension(:,:,:) :: p_half, p_full
real, allocatable, dimension(:,:,:) :: p_half_prev, p_full_prev
real, allocatable, dimension(:,:,:) :: z_half, z_full

type(tracer_type), allocatable, dimension(:) :: tracer_attributes
real, allocatable, dimension(:,:,:,:,:) :: grid_tracers
real, allocatable, dimension(:,:,:    ) :: psg, wg_full
real, allocatable, dimension(:,:,:,:  ) :: ug, vg, tg

real, allocatable, dimension(:,:,:)   :: bucket_depth                ! added by LJJ
real, allocatable, dimension(:,:    ) :: dt_psg, dt_bucket, bucket_diffusion, filt   ! added by LJJ
real, allocatable, dimension(:,:,:  ) :: dt_ug, dt_vg, dt_tg
real, allocatable, dimension(:,:,:,:) :: dt_tracers

real, allocatable, dimension(:)   :: deg_lat, rad_lat, deg_lon
real, allocatable, dimension(:,:) :: rad_lat_2d, rad_lon_2d

real, allocatable, dimension(:,:)   ::                       &
     z_surf,               &   ! surface height
     t_surf,               &   ! surface temperature
     q_surf,               &   ! surface moisture
     u_surf,               &   ! surface U wind
     v_surf,               &   ! surface V wind
     rough_mom,            &   ! momentum roughness length for surface_flux
     rough_heat,           &   ! heat roughness length for surface_flux
     rough_moist,          &   ! moisture roughness length for surface_flux
     depth_change_lh,      &   ! tendency in bucket depth due to latent heat transfer     ! added by LJJ
     depth_change_cond,    &   ! tendency in bucket depth due to condensation rain        ! added by LJJ
     depth_change_conv,    &   ! tendency in bucket depth due to convection rain          ! added by LJJ
     gust,                 &   ! gustiness constant
     flux_t,               &   ! surface sensible heat flux
     flux_q,               &   ! surface moisture flux
     flux_r,               &   ! surface radiation flux
     flux_u,               &   ! surface flux of zonal mom.
     flux_v,               &   ! surface flux of meridional mom.
     drag_m,               &   ! momentum drag coefficient
     drag_t,               &   ! heat drag coefficient
     drag_q,               &   ! moisture drag coefficient
     w_atm,                &   ! wind speed
     ustar,                &   ! friction velocity
     bstar,                &   ! buoyancy scale
     qstar,                &   ! moisture scale
     dhdt_surf,            &   ! d(sensible heat flux)/d(surface temp)
     dedt_surf,            &   ! d(latent heat flux)/d(surface temp)???
     dedq_surf,            &   ! d(latent heat flux)/d(surface moisture)???
     drdt_surf,            &   ! d(upward longwave)/d(surface temp)
     dhdt_atm,             &   ! d(sensible heat flux)/d(atmos.temp)
     dedq_atm,             &   ! d(latent heat flux)/d(atmospheric mixing rat.)
     dtaudv_atm,           &   ! d(stress component)/d(atmos wind)
     fracland,             &   ! fraction of land in gridbox
     rough                     ! roughness for vert_turb_driver


real, allocatable, dimension(:,:,:) ::                                        &
     diff_m,               &   ! momentum diffusion coeff.
     diff_t,               &   ! temperature diffusion coeff.
     diss_heat,            &   !  heat dissipated by vertical diffusion
     flux_tr,              &   ! surface tracer flux
     non_diff_dt_ug,       &   ! zonal wind tendency except from vertical diffusion
     non_diff_dt_vg,       &   ! merid. wind tendency except from vertical diffusion
     non_diff_dt_tg,       &   ! temperature tendency except from vertical diffusion
     non_diff_dt_qg,       &   ! moisture tendency except from vertical diffusion
     conv_dt_tg,           &   ! temperature tendency from convection
     conv_dt_qg,           &   ! moisture tendency from convection
     cond_dt_tg,           &   ! temperature tendency from condensation
     cond_dt_qg                ! moisture tendency from condensation

logical, allocatable, dimension(:,:) ::                                       &
     avail,                &   ! generate surf. flux (all true)
     land,                 &   ! land points (all false)
     coldT                     ! should precipitation be snow at this point

real, allocatable, dimension(:,:) ::                                          &
     klzbs,                &   ! stored level of zero buoyancy values         
     cape,                 &   ! convectively available potential energy      
     cin,                  &   ! convective inhibition (this and the above are before the adjustment)                                                    
     invtau_q_relaxation,  &   ! temperature relaxation time scale           
     invtau_t_relaxation,  &   ! humidity relaxation time scale               
     rain,                 &   !                                              
     snow

real, allocatable, dimension(:,:) ::                                          &
     ocean_mask                !  used for initializing bucket depth        ! added by LJJ

real, allocatable, dimension(:,:,:) ::                                        &
     t_ref,                &   ! relaxation temperature for bettsmiller scheme
     q_ref                     ! relaxation moisture for bettsmiller scheme

real, allocatable, dimension(:) ::                                            &
     sin_lat                   ! sine of latitude

integer ::                                                                    &
     id_diff_dt_ug,        &   ! zonal wind tendency from vertical diffusion
     id_diff_dt_vg,        &   ! merid. wind tendency from vertical diffusion
     id_diff_dt_tg,        &   ! temperature tendency from vertical diffusion
     id_diff_dt_qg,        &   ! moisture tendency from vertical diffusion
     id_conv_rain,         &   ! rain from convection
     id_cond_rain,         &   ! rain from condensation
     id_conv_dt_tg,        &   ! temperature tendency from convection
     id_conv_dt_qg,        &   ! temperature tendency from convection
     id_cond_dt_tg,        &   ! temperature tendency from convection
     id_cond_dt_qg,        &   ! temperature tendency from convection 
     id_drag_m,            &
     id_drag_q,            &
     id_drag_t,            &
     id_bucket_depth,      &   ! bucket depth vaiable for output added by LJJ
     id_bucket_depth_conv, &   ! bucket depth variation induced by convection, added by LJJ
     id_bucket_depth_cond, &   ! bucket depth variation induced by condensation, added by LJJ
     id_bucket_depth_lh,   &   ! bucket depth variation induced by LH, added by LJJ
     id_bucket_diffusion       ! diffused surface water depth 

integer, allocatable, dimension(:,:) ::  & convflag                  
! indicates which qe convection subroutines are used              

logical :: used
integer, dimension(4) :: axes
                                                                               
real, allocatable, dimension(:,:)   ::                                        &
     net_surf_sw_down,      &   ! net sw flux at surface
     surf_lw_down               ! downward lw flux at surface

integer :: previous, current, future
logical :: module_is_initialized =.false.
character(len=4) :: ch_tmp1, ch_tmp2

integer         :: dt_integer
real            :: dt_real
type(time_type) :: Time_step

integer, dimension(4) :: axis_id

type(surf_diff_type) :: Tri_surf ! used by gcm_vert_diff 

!=================================================================================================================================
contains
!=================================================================================================================================

subroutine atmosphere_init(Time_init, Time, Time_step_in)

type (time_type), intent(in)  :: Time_init, Time, Time_step_in

integer :: seconds, days, lon_max, lat_max, ntr, nt, i, j, ierr, io
integer, dimension(4) :: siz
real, dimension(2) :: time_pointers
character(len=64) :: file, tr_name
character(len=256) :: message

integer :: unit


if(module_is_initialized) return

call write_version_number(version, tagname)

Time_step = Time_step_in
call get_time(Time_step, seconds, days)
dt_integer   = 86400*days + seconds
dt_real      = float(dt_integer)

unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=atmosphere_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'atmosphere_nml')
enddo
10 call close_file (unit)

if ( mpp_pe() == mpp_root_pe() )   write (stdlog(), nml=atmosphere_nml)

call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
allocate (tracer_attributes(num_tracers))
call spectral_dynamics_init(Time, Time_step, tracer_attributes, dry_model, nhum)
call get_grid_domain(is, ie, js, je)
call get_num_levels(num_levels)

allocate (p_half       (is:ie, js:je, num_levels+1))
allocate (p_half_prev  (is:ie, js:je, num_levels+1))
allocate (z_half       (is:ie, js:je, num_levels+1))
allocate (p_full       (is:ie, js:je, num_levels))
allocate (p_full_prev  (is:ie, js:je, num_levels))
allocate (z_full       (is:ie, js:je, num_levels))
allocate (wg_full      (is:ie, js:je, num_levels))
allocate (psg          (is:ie, js:je, num_time_levels))
allocate (ug           (is:ie, js:je, num_levels, num_time_levels))
allocate (vg           (is:ie, js:je, num_levels, num_time_levels))
allocate (tg           (is:ie, js:je, num_levels, num_time_levels))
allocate (grid_tracers (is:ie, js:je, num_levels, num_time_levels, num_tracers ))

allocate (dt_psg     (is:ie, js:je))
allocate (dt_bucket  (is:ie, js:je))                        ! added by LJJ
allocate (bucket_diffusion(is:ie, js:je))                   ! added by LJJ
allocate (filt       (is:ie, js:je))
allocate (dt_ug      (is:ie, js:je, num_levels))
allocate (dt_vg      (is:ie, js:je, num_levels))
allocate (dt_tg      (is:ie, js:je, num_levels))
allocate (dt_tracers (is:ie, js:je, num_levels, num_tracers))

allocate (deg_lat    (       js:je))
allocate (deg_lon    (is:ie       ))
allocate (rad_lat    (       js:je))
allocate (rad_lat_2d (is:ie, js:je))
allocate (rad_lon_2d (is:ie, js:je))

allocate(z_surf      (is:ie, js:je))
allocate(t_surf      (is:ie, js:je))
allocate(bucket_depth (is:ie, js:je, num_time_levels))        ! added by LJJ
allocate(depth_change_lh(is:ie, js:je))                       ! added by LJJ
allocate(depth_change_cond(is:ie, js:je))                     ! added by LJJ
allocate(depth_change_conv(is:ie, js:je))                     ! added by LJJ
allocate(q_surf      (is:ie, js:je))
allocate(u_surf      (is:ie, js:je))
allocate(v_surf      (is:ie, js:je))
allocate(rough_mom   (is:ie, js:je))
allocate(rough_heat  (is:ie, js:je))
allocate(rough_moist (is:ie, js:je))
allocate(gust        (is:ie, js:je))
allocate(flux_t      (is:ie, js:je))
allocate(flux_q      (is:ie, js:je))
allocate(flux_r      (is:ie, js:je))
allocate(flux_u      (is:ie, js:je))
allocate(flux_v      (is:ie, js:je))
allocate(drag_m      (is:ie, js:je))
allocate(drag_t      (is:ie, js:je))
allocate(drag_q      (is:ie, js:je))
allocate(w_atm       (is:ie, js:je))
allocate(ustar       (is:ie, js:je))
allocate(bstar       (is:ie, js:je))
allocate(qstar       (is:ie, js:je))
allocate(dhdt_surf   (is:ie, js:je))
allocate(dedt_surf   (is:ie, js:je))
allocate(dedq_surf   (is:ie, js:je))
allocate(drdt_surf   (is:ie, js:je))
allocate(dhdt_atm    (is:ie, js:je))
allocate(dedq_atm    (is:ie, js:je))
allocate(dtaudv_atm  (is:ie, js:je))
allocate(land        (is:ie, js:je))
allocate(avail       (is:ie, js:je))
allocate(fracland    (is:ie, js:je))
allocate(rough       (is:ie, js:je))
allocate(sin_lat     (       js:je))
allocate(diff_t      (is:ie, js:je, num_levels))
allocate(diff_m      (is:ie, js:je, num_levels))
allocate(diss_heat   (is:ie, js:je, num_levels))
allocate(flux_tr     (is:ie, js:je,             num_tracers))

allocate(non_diff_dt_ug  (is:ie, js:je, num_levels))
allocate(non_diff_dt_vg  (is:ie, js:je, num_levels))
allocate(non_diff_dt_tg  (is:ie, js:je, num_levels))
allocate(non_diff_dt_qg  (is:ie, js:je, num_levels))


allocate(net_surf_sw_down        (is:ie, js:je))
allocate(surf_lw_down            (is:ie, js:je))
allocate(conv_dt_tg  (is:ie, js:je, num_levels))
allocate(conv_dt_qg  (is:ie, js:je, num_levels))
allocate(cond_dt_tg  (is:ie, js:je, num_levels))
allocate(cond_dt_qg  (is:ie, js:je, num_levels))

allocate(coldT        (is:ie, js:je))
allocate(klzbs        (is:ie, js:je))
allocate(cape         (is:ie, js:je))
allocate(cin          (is:ie, js:je))
allocate(invtau_q_relaxation  (is:ie, js:je))
allocate(invtau_t_relaxation  (is:ie, js:je))
allocate(rain         (is:ie, js:je))
allocate(snow         (is:ie, js:je))
allocate(convflag     (is:ie, js:je))
allocate(ocean_mask   (is:ie, js:je))              ! added by LJJ

allocate(     t_ref  (is:ie, js:je, num_levels))
allocate(     q_ref  (is:ie, js:je, num_levels))



call get_surf_geopotential(z_surf)
z_surf = z_surf/grav 
            
t_ref = 0.0; q_ref = 0.0

coldT = .false.
rain = 0.0; snow = 0.0

land = .false. 
avail = .true.
rough_mom = roughness_mom
rough_heat = roughness_heat
rough_moist = roughness_moist
gust = 1.0 
q_surf = 0.0
u_surf = 0.0
v_surf = 0.0
fracland = 0.0 ! fraction of each gridbox that is land
flux_tr = 0.0

! added by LJJ
! default to all ocean
! ocean_mask = 1.0

!if(file_exist('INPUT/ocean_mask.nc')) then
!  call read_data('INPUT/ocean_mask.nc', 'ocean_mask', ocean_mask)
!  where (ocean_mask > 0.0)
!     bucket_depth(:,:,1) = init_bucket_depth
!     bucket_depth(:,:,2) = init_bucket_depth
!  elsewhere
!     bucket_depth(:,:,1) = init_bucket_depth_land
!     bucket_depth(:,:,2) = init_bucket_depth_land
!  endwhere
!else 
!  bucket_depth = init_bucket_depth   !initializing depth of bucket
!endif 

call get_deg_lat(deg_lat)
call get_deg_lon(deg_lon)


if(bucket) then
  do i = is, ie
     do j = js, je
        if (deg_lat(j) >= land_bottom .and. deg_lat(j) <= land_top .and. deg_lon(i) >= land_left .and. deg_lon(i) <= land_right) then
           ocean_mask(i,j) = 0.0         ! land
           bucket_depth(i,j,1)  = init_bucket_depth_land
           bucket_depth(i,j,2)  = init_bucket_depth_land
       else
           ocean_mask(i,j) = 1.0
           bucket_depth(i,j,1) = init_bucket_depth
           bucket_depth(i,j,2) = init_bucket_depth
        endif
      enddo
   enddo
else
   ocean_mask = 1.0
   bucket_depth = init_bucket_depth   
endif

dt_bucket = 0.0
bucket_diffusion = 0.0
filt = 0.0
! end addition by LJJ

p_half = 0.; z_half = 0.; p_full = 0.; z_full = 0.
wg_full = 0.; psg = 0.; ug = 0.; vg = 0.; tg = 0.; grid_tracers = 0.
p_half_prev = 0.; p_full_prev = 0.
dt_psg = 0.; dt_ug  = 0.; dt_vg  = 0.; dt_tg  = 0.; dt_tracers = 0.


!--------------------------------------------------------------------------------------------------------------------------------

file = 'INPUT/atmosphere.res.nc'
if(file_exist(trim(file))) then
   call get_lon_max(lon_max)
   call get_lat_max(lat_max)
   call field_size(trim(file), 'ug', siz)
   if(lon_max /= siz(1) .or. lat_max /= siz(2)) then
      write(message,*) 'Resolution of restart data does not match resolution specified on namelist. Restart data: lon_max=', &
                     siz(1),', lat_max=',siz(2),'  Namelist: lon_max=',lon_max,', lat_max=',lat_max
      call error_mesg('atmosphere_init', message, FATAL)
   endif
   call nullify_domain()
! call read_data(trim(file), 'previous', previous)           ! No interface of read_data exists to read integer scalars
! call read_data(trim(file), 'current',  current)            ! No interface of read_data exists to read integer scalars
   call read_data(trim(file), 'time_pointers', time_pointers) ! Getaround for no interface to read integer scalars
   previous = int(time_pointers(1))                           ! Getaround for no interface to read integer scalars
   current  = int(time_pointers(2))                           ! Getaround for no interface to read integer scalars
   do nt=1,num_time_levels
         call read_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain, timelevel=nt)
         call read_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain, timelevel=nt)
         call read_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain, timelevel=nt)
         call read_data(trim(file), 'psg', psg(:,:,  nt), grid_domain, timelevel=nt)
   do ntr = 1,num_tracers
            tr_name = trim(tracer_attributes(ntr)%name)
            call read_data(trim(file), trim(tr_name), grid_tracers(:,:,:,nt,ntr), grid_domain, timelevel=nt)
   enddo ! end loop over tracers
enddo ! end loop over time levels
   call read_data(trim(file), 'wg_full', wg_full, grid_domain)
else
   previous = 1; current = 1
   call get_initial_fields(ug(:,:,:,1), vg(:,:,:,1), tg(:,:,:,1), psg(:,:,1), grid_tracers(:,:,:,1,:))
endif

   !----------------------------------------------------------------------------------------------------------------------------
if(dry_model) then 
 call compute_pressures_and_heights(tg(:,:,:,current), psg(:,:,current), z_full, z_half, p_full, p_half)
else
  call compute_pressures_and_heights( &
          tg(:,:,:,current), psg(:,:,current), z_full, z_half, p_full, p_half, grid_tracers(:,:,:,current,nhum))
     p_full_prev = p_full
     p_half_prev = p_half
endif


call get_deg_lat(deg_lat)
call get_deg_lon(deg_lon)
do j=js,je
  rad_lat_2d(:,j) = deg_lat(j)*pi/180.
enddo
do i=is,ie  
  rad_lon_2d(i,:) = deg_lon(i)*pi/180.
enddo


if(mixed_layer_bc) then
 ! need an initial condition for the mixed layer temperature
 ! may be overwritten by restart file
 ! choose an unstable initial condition to allow moisture
 ! to quickly enter the atmosphere avoiding problems with the convection scheme
 t_surf = tg(:,:,num_levels,current)+1.0
 !call mixed_layer_init(is, ie, js, je, num_levels, t_surf, get_axis_id(), Time)
 call mixed_layer_init(is, ie, js, je, num_levels, t_surf, bucket_depth, ocean_mask, get_axis_id(), Time)
endif

! fridoo modification 01 feb 2012 
if(tapio_forcing) then 
      call grid_phys_init(get_axis_id(), Time)
elseif(hs) then
      call hs_forcing_init(get_axis_id(), Time)
endif
!end fridoo
 
if(turb) then 
! need to call gcm_vert_diff_init even if using gcm_vert_diff (rather than
! gcm_vert_diff_down) because the variable sphum is not initialized
! otherwise in the vert_diff module
 call gcm_vert_diff_init (Tri_surf, ie-is+1, je-js+1, num_levels, .true., do_virtual)
end if

if(ldry_convection) call dry_convection_init(get_axis_id(), Time)


if(.not.dry_model) then ! loop added by fridoo 01 feb 2012
   call lscale_cond_init()

   axes = get_axis_id()

   id_cond_dt_qg = register_diag_field(mod_name, 'dt_qg_condensation',        &
        axes(1:3), Time, 'Moisture tendency from condensation','kg/kg/s')
   id_cond_dt_tg = register_diag_field(mod_name, 'dt_tg_condensation',        &
        axes(1:3), Time, 'Temperature tendency from condensation','K/s')
   id_cond_rain = register_diag_field(mod_name, 'condensation_rain',          &
        axes(1:2), Time, 'Rain from condensation','kg/m/m/s')

   if(bucket) then
      id_bucket_depth = register_diag_field(mod_name, 'bucket_depth',            &         ! added by LJJ
           axes(1:2), Time, 'Depth of surface reservoir', 'm')
      id_bucket_depth_conv = register_diag_field(mod_name, 'bucket_depth_conv',  &         ! added by LJJ
        axes(1:2), Time, 'Tendency of bucket depth induced by Convection', 'm/s')
      id_bucket_depth_cond = register_diag_field(mod_name, 'bucket_depth_cond',  &         ! added by LJJ
        axes(1:2), Time, 'Tendency of bucket depth induced by Condensation', 'm/s')
      id_bucket_depth_lh = register_diag_field(mod_name, 'bucket_depth_lh',      &         ! added by LJJ
        axes(1:2), Time, 'Tendency of bucket depth induced by LH', 'm/s')
      id_bucket_diffusion = register_diag_field(mod_name, 'bucket_diffusion',  &    !added by SDG
        axes(1:2), Time, 'Diffusion rate of bucket','m/s')
   endif

   id_drag_m = register_diag_field(mod_name, 'drag_coeff_mo',                  &
        axes(1:2), Time, 'Drag coefficient for momentum flux at surface','m/s')
   id_drag_q = register_diag_field(mod_name, 'drag_coeff_lh',                  &
        axes(1:2), Time, 'Drag coefficient for moisture flux at surface','m/s')
   id_drag_t = register_diag_field(mod_name, 'drag_coeff_sh',                  &
        axes(1:2), Time, 'Drag coefficient for sensible heat flux at surface','m/s')
endif

if(lwet_convection) then  
   call qe_moist_convection_init()
   id_conv_dt_qg = register_diag_field(mod_name, 'dt_qg_convection',          &
        axes(1:3), Time, 'Moisture tendency from convection','kg/kg/s')
   id_conv_dt_tg = register_diag_field(mod_name, 'dt_tg_convection',          &
        axes(1:3), Time, 'Temperature tendency from convection','K/s')
   id_conv_rain = register_diag_field(mod_name, 'convection_rain',            &
        axes(1:2), Time, 'Rain from convection','kg/m/m/s')
endif

if(two_stream) call radiation_init(is, ie, js, je, num_levels, get_axis_id(), Time)

if(turb) then
   call vert_turb_driver_init (ie-is+1,je-js+1,num_levels,get_axis_id(),Time)

   axes = get_axis_id()
   id_diff_dt_ug = register_diag_field(mod_name, 'dt_ug_diffusion',        & 
        axes(1:3), Time, 'zonal wind tendency from diffusion','m/s^2')
   id_diff_dt_vg = register_diag_field(mod_name, 'dt_vg_diffusion',        &
        axes(1:3), Time, 'meridional wind tendency from diffusion','m/s^2')
   id_diff_dt_tg = register_diag_field(mod_name, 'dt_tg_diffusion',        &
        axes(1:3), Time, 'temperature diffusion tendency','T/s')
   id_diff_dt_qg = register_diag_field(mod_name, 'dt_qg_diffusion',        &
           axes(1:3), Time, 'moisture diffusion tendency','T/s')
endif

module_is_initialized = .true.

return
end subroutine atmosphere_init

!=================================================================================================================================

subroutine atmosphere(Time)
type(time_type), intent(in) :: Time

real    :: delta_t
type(time_type) :: Time_next

real, dimension(is:ie, js:je, num_levels) :: tg_tmp, qg_tmp 
real inv_cp_air 

! added fridoo sept 2012 (tmerlis)
real :: atmos_water_previous, atmos_water_future, evap_minus_precip, water_correction_factor
! end fridoo

integer  :: i,j

if(.not.module_is_initialized) then
  call error_mesg('atmosphere','atmosphere module is not initialized',FATAL)
endif

dt_ug  = 0.0
dt_vg  = 0.0
dt_tg  = 0.0
dt_psg = 0.0
dt_bucket = 0.0                ! added by LJJ
bucket_diffusion = 0.0         ! added by LJJ
filt      = 0.0                ! added by LJJ
dt_tracers = 0.0 

conv_dt_tg  = 0.0
conv_dt_qg  = 0.0
cond_dt_tg  = 0.0
cond_dt_qg  = 0.0

evap_minus_precip = 0.0

if(current == previous) then
  delta_t = dt_real
else
  delta_t = 2*dt_real
endif

! added fridoo sept 2012 (tmerlis)

if(atmos_water_correction) then
      atmos_water_previous = mass_weighted_global_integral( &
           grid_tracers(:,:,:,previous,nhum), psg(:,:,previous))
endif

! end fridoo

if(ldry_convection) then
   grid_tracers(:,:,:,:,nhum) = 0.0
   call dry_convection(           Time,             tg(:,:,:,previous), &
                         p_full(:,:,:),                  p_half(:,:,:), &
                     conv_dt_tg(:,:,:) )
   dt_tg = dt_tg + conv_dt_tg
endif


! fridoo modification 01 feb 2012
if(tapio_forcing) then
   call compute_grid_physics(        1,                        ie-is+1, &
                                     1,                        je-js+1, &
                                  Time,                        delta_t, &
          rad_lat_2d(:,:             ),          dt_ug(:,:,:         ), &
               dt_vg(:,:,:           ),          dt_tg(:,:,:         ), &
          dt_tracers(:,:,:,nhum      ),     dt_tracers(:,:,:,:       ), &
                    ug(:,:,:,previous),             vg(:,:,:,previous), &
                    tg(:,:,:,previous),   grid_tracers(:,:,:,previous,nhum),&
          grid_tracers(:,:,:,previous,:),       p_half(:,:,:         ), &
              p_full(:,:,:           ),             ug(:,:,:,previous), &
                  vg(:,:,:,previous  ),             tg(:,:,:,previous), &
        grid_tracers(:,:,:,previous,nhum),grid_tracers(:,:,:,previous,:),&
              p_half(:,:,:         ),           p_full(:,:,:           ) )
elseif(hs) then
   call hs_forcing(1, ie-is+1, 1, je-js+1, delta_t, Time, rad_lat_2d(:,:),    &
                p_half(:,:,:         ),       p_full(:,:,:           ), &
                    ug(:,:,:,previous),           vg(:,:,:,previous  ), &
                    tg(:,:,:,previous), grid_tracers(:,:,:,previous,:), &
                    ug(:,:,:,previous),           vg(:,:,:,previous  ), &
                    tg(:,:,:,previous), grid_tracers(:,:,:,previous,:), &
                 dt_ug(:,:,:         ),        dt_vg(:,:,:           ), &
                 dt_tg(:,:,:         ),   dt_tracers(:,:,:,:))
endif
! end fridoo modification

if (lwet_convection) then
   rain = 0.0; snow = 0.0
   call qe_moist_convection ( delta_t,              tg(:,:,:,previous),      &
    grid_tracers(:,:,:,previous,nhum),                     p_full_prev,      &
                          p_half_prev,                           coldT,      &
                                 rain,                            snow,      &
                           conv_dt_tg,                      conv_dt_qg,      &
                                q_ref,                        convflag,      &
                                klzbs,                            cape,      &
                                  cin,             invtau_q_relaxation,      &
                  invtau_t_relaxation,                           t_ref)

   tg_tmp = conv_dt_tg + tg(:,:,:,previous) 
   qg_tmp = conv_dt_qg + grid_tracers(:,:,:,previous,nhum)
!  note the delta's are returned rather than the time derivatives

   conv_dt_tg = conv_dt_tg/delta_t
   conv_dt_qg = conv_dt_qg/delta_t
   depth_change_conv = rain/dens_h2o     ! added by LJJ
   rain       = rain/delta_t
                       
   if(atmos_water_correction) evap_minus_precip = evap_minus_precip - area_weighted_global_mean(rain) 
                                                    
   dt_tg = dt_tg + conv_dt_tg
   dt_tracers(:,:,:,nhum) = dt_tracers(:,:,:,nhum) + conv_dt_qg
                                                                            
   if(id_conv_dt_qg > 0) used = send_data(id_conv_dt_qg, conv_dt_qg, Time)
   if(id_conv_dt_tg > 0) used = send_data(id_conv_dt_tg, conv_dt_tg, Time)
   if(id_conv_rain > 0) used = send_data(id_conv_rain, rain, Time)
                                                                             
else 

   tg_tmp = tg(:,:,:,previous) 
   qg_tmp = grid_tracers(:,:,:,previous,nhum)

endif

rain = 0.0

if(.not.dry_model) then      !loop added by fridoo 01 feb 2012
   call lscale_cond (         tg_tmp,                          qg_tmp,        &
                         p_full_prev,                     p_half_prev,        &
                               coldT,                            rain,        &
                                snow,                      cond_dt_tg,        &
                          cond_dt_qg )
                                                                          
   cond_dt_tg = cond_dt_tg/delta_t
   cond_dt_qg = cond_dt_qg/delta_t
   depth_change_cond = rain/dens_h2o        ! added by LJJ
   rain       = rain/delta_t
                                      
   if(atmos_water_correction) evap_minus_precip = evap_minus_precip - area_weighted_global_mean(rain) 
                                       
   dt_tg = dt_tg + cond_dt_tg
   dt_tracers(:,:,:,nhum) = dt_tracers(:,:,:,nhum) + cond_dt_qg

                                                                               
   if(id_cond_dt_qg > 0) used = send_data(id_cond_dt_qg, cond_dt_qg, Time)
   if(id_cond_dt_tg > 0) used = send_data(id_cond_dt_tg, cond_dt_tg, Time)
   if(id_cond_rain > 0) used = send_data(id_cond_rain, rain, Time)
endif

! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

if(two_stream) then
   call radiation_down(is, js, Time,                   &
                       rad_lat_2d(:,:),                &
                       p_half(:,:,:),                  &
                       tg(:,:,:,previous),             &
                       net_surf_sw_down(:,:),          &
                       surf_lw_down(:,:))
end if


! fridoo modification 01 feb 2012
if(tapio_forcing) then 
   t_surf = surface_temperature(tg(:,:,:,previous), p_full, p_half)
endif
! end fridoo modification

call surface_flux(                                                          &
                  tg(:,:,num_levels,previous),                              &
   grid_tracers(:,:,num_levels,previous,nhum),                              &
                  ug(:,:,num_levels,previous),                              &
                  vg(:,:,num_levels,previous),                              &
                       p_full(:,:,num_levels),                              &
           z_full(:,:,num_levels)-z_surf(:,:),                              &
                     p_half(:,:,num_levels+1),                              &
                                  t_surf(:,:),                              &
                                  t_surf(:,:),                              &
                                  q_surf(:,:),                              &
                    bucket_depth(:,:,current),                              &     ! added by LJJ
                         depth_change_lh(:,:),                              &     ! added by LJJ
                       depth_change_conv(:,:),                              &     ! added by LJJ
                       depth_change_cond(:,:),                              &     ! added by LJJ
                                  u_surf(:,:),                              &
                                  v_surf(:,:),                              &
                               rough_mom(:,:),                              &
                              rough_heat(:,:),                              &
                             rough_moist(:,:),                              &
                                    gust(:,:),                              &
                                  flux_t(:,:),                              &
                                  flux_q(:,:),                              &
                                  flux_r(:,:),                              &
                                  flux_u(:,:),                              &
                                  flux_v(:,:),                              &
                                  drag_m(:,:),                              &
                                  drag_t(:,:),                              &
                                  drag_q(:,:),                              &
                                   w_atm(:,:),                              &
                                   ustar(:,:),                              &
                                   bstar(:,:),                              &
                                   qstar(:,:),                              &
                               dhdt_surf(:,:),                              &
                               dedt_surf(:,:),                              &
                               dedq_surf(:,:),                              &
                               drdt_surf(:,:),                              &
                                dhdt_atm(:,:),                              &
                                dedq_atm(:,:),                              &
                              dtaudv_atm(:,:),                              &
                                      delta_t,                              &
                                    land(:,:),                              &
                                   avail(:,:)  )

   if(id_drag_m > 0) used = send_data(id_drag_m, drag_m, Time)
   if(id_drag_q > 0) used = send_data(id_drag_q, drag_q, Time)
   if(id_drag_t > 0) used = send_data(id_drag_t, drag_t, Time)

! Now complete the radiation calculation by computing the upward and net fluxes.

if(two_stream) then
   call radiation_up(is, js, Time,                   &
                     rad_lat_2d(:,:),                &
                     p_half(:,:,:),                  &
                     t_surf(:,:),                    &
                     tg(:,:,:,previous),             &
                     dt_tg(:,:,:))
end if

if(turb) then

   call vert_turb_driver(            1,                              1, &
                                  Time,                 Time+Time_step, &
                               delta_t,                  fracland(:,:), &
                         p_half(:,:,:),                  p_full(:,:,:), &
                         z_half(:,:,:),                  z_full(:,:,:), &
                            ustar(:,:),                     bstar(:,:), &
                            rough(:,:),             ug(:,:,:,current ), &
                    vg(:,:,:,current ),             tg(:,:,:,current ), &
      grid_tracers(:,:,:,current,nhum),             ug(:,:,:,previous), &
                    vg(:,:,:,previous),             tg(:,:,:,previous), &
     grid_tracers(:,:,:,previous,nhum),                   dt_ug(:,:,:), &
                          dt_vg(:,:,:),                   dt_tg(:,:,:), &
                dt_tracers(:,:,:,nhum),                  diff_t(:,:,:), &
                         diff_m(:,:,:),                      gust(:,:)  )

!
!! Don't zero these derivatives as the surface flux depends implicitly
!! on the lowest level values
!! However it should be noted that these derivatives do not take into
!! account the change in the Monin-Obukhov coefficients, and so are not
!! very accurate.
!
!!$   dtaudv_atm = 0.0



! loop added by fridoo 01 feb 2012
   if(.not.mixed_layer_bc .and. .not.tapio_forcing .and. .not.hs) then
       call error_mesg('atmosphere','no diffusion implentation for non-mixed layer b.c.',FATAL) 
   endif 
!end addition

! We must use gcm_vert_diff_down and _up rather than gcm_vert_diff as the surface flux
! depends implicitly on the surface values
!
! Don't want to do time splitting for the implicit diffusion step in case
! of compensation of the tendencies
!

   non_diff_dt_ug  = dt_ug
   non_diff_dt_vg  = dt_vg
   non_diff_dt_tg  = dt_tg
   non_diff_dt_qg  = dt_tracers(:,:,:,nhum)

   call gcm_vert_diff_down (1, 1,                                          &
                            delta_t,             ug(:,:,:,previous),       &
                            vg(:,:,:,previous),  tg(:,:,:,previous),       &
                            grid_tracers(:,:,:,previous,nhum),             &
                            grid_tracers(:,:,:,previous,:), diff_m(:,:,:), &
                            diff_t(:,:,:),                  p_half(:,:,:), &
                            p_full(:,:,:),                  z_full(:,:,:), &
                            flux_u(:,:),                    flux_v(:,:),   &
                            dtaudv_atm(:,:),                               &
                            flux_tr(:,:,:),                                &
                            dt_ug(:,:,:),                    dt_vg(:,:,:), &
                            dt_tg(:,:,:),          dt_tracers(:,:,:,nhum), &
                            dt_tracers(:,:,:,:),         diss_heat(:,:,:), &
                            Tri_surf)

!
! update surface temperature
!
   if(mixed_layer_bc) then
         call mixed_layer(                                                 &
                              is,                                          &
                              ie,                                          &
                              js,                                          &
                              je,                                          &
                              Time,                                        &
                              t_surf(:,:),                                 &
                              flux_t(:,:),                                 &
                              flux_q(:,:),                                 &
                              flux_r(:,:),                                 &
                              flux_u(:,:),                                 &
                                  dt_real,                                 &
                    net_surf_sw_down(:,:),                                 &
                        surf_lw_down(:,:),                                 &
                            Tri_surf,                                      &
                           dhdt_surf(:,:),                                 &
                           dedt_surf(:,:),                                 &
                           dedq_surf(:,:),                                 &
                           drdt_surf(:,:),                                 &
                            dhdt_atm(:,:),                                 &
                            dedq_atm(:,:))

   
   endif

   if(atmos_water_correction) evap_minus_precip = evap_minus_precip + area_weighted_global_mean( flux_q )

   call gcm_vert_diff_up (1, 1, delta_t, Tri_surf,  &
                          dt_tg(:,:,:),     dt_tracers(:,:,:,nhum))

   

   if(id_diff_dt_ug > 0) used = send_data(id_diff_dt_ug, dt_ug - non_diff_dt_ug, Time)
   if(id_diff_dt_vg > 0) used = send_data(id_diff_dt_vg, dt_vg - non_diff_dt_vg, Time)
   if(id_diff_dt_tg > 0) used = send_data(id_diff_dt_tg, dt_tg - non_diff_dt_tg, Time)
   if(id_diff_dt_qg > 0) used = send_data(id_diff_dt_qg, dt_tracers(:,:,:,nhum) - non_diff_dt_qg, Time)
endif
 


Time_next = Time + Time_step

if(previous == current) then
  future = num_time_levels + 1 - current
else
  future = previous
endif

! added by LJJ
if(bucket) then
   ! bucket time tendency
   dt_bucket = depth_change_cond + depth_change_conv - depth_change_lh
   !change in bucket depth in one leapfrog timestep [m]                                 

   !diffuse_surf_water transforms dt_bucket to spherical, diffuses water, and transforms back
   call diffuse_surf_water(dt_bucket,bucket_depth(:,:,previous),delta_t,damping_coeff_bucket,bucket_diffusion)

   ! use the raw filter in leapfrog time stepping

   filt(:,:) = bucket_depth(:,:,previous) - 2.0 * bucket_depth(:,:,current)

   if(previous == current) then
      bucket_depth(:,:,future ) = bucket_depth(:,:,previous) + dt_bucket
      bucket_depth(:,:,current) = bucket_depth(:,:,current ) + robert_bucket &
        *(bucket_depth(:,:,previous) - 2.0*bucket_depth(:,:,current) + bucket_depth(:,:,future)) * raw_bucket
   else
      bucket_depth(:,:,current) = bucket_depth(:,:,current ) + robert_bucket &
        *(bucket_depth(:,:,previous) - 2.0*bucket_depth(:,:,current)) * raw_bucket 
      bucket_depth(:,:,future ) = bucket_depth(:,:,previous) + dt_bucket
      bucket_depth(:,:,current) = bucket_depth(:,:,current) + robert_bucket * bucket_depth(:,:,future) * raw_bucket
   endif

   bucket_depth(:,:,future) = bucket_depth(:,:,future) + robert_bucket * (filt(:,:) + bucket_depth(:,:, future)) &
                           * (raw_bucket - 1.0)  

   where (bucket_depth <= 0.) bucket_depth = 0.

   ! truncate surface reservoir over land points
   do i=is,ie
      do j=js,je
         if (ocean_mask(i,j) < 0.5) then
            if (bucket_depth(i,j,future) > max_bucket_depth_land) then
               bucket_depth(i,j,future) = max_bucket_depth_land
            endif
         endif
      enddo
   enddo

   if(id_bucket_depth > 0) used = send_data(id_bucket_depth, bucket_depth(:,:,future), Time)
   if(id_bucket_depth_conv > 0) used = send_data(id_bucket_depth_conv, depth_change_conv(:,:), Time)
   if(id_bucket_depth_cond > 0) used = send_data(id_bucket_depth_cond, depth_change_cond(:,:), Time)
   if(id_bucket_depth_lh > 0) used = send_data(id_bucket_depth_lh, depth_change_lh(:,:), Time)

   bucket_diffusion=bucket_diffusion/delta_t !convert to m/s for output

   if(id_bucket_diffusion > 0) used = send_data(id_bucket_diffusion, bucket_diffusion(:,:), Time)

   bucket_diffusion=bucket_diffusion*delta_t !convert back to m 

endif
! end addition by LJJ

call spectral_dynamics(Time, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
                       tg(:,:,:,future), tracer_attributes, grid_tracers(:,:,:,:,:), future, &
                       dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers, wg_full, p_full, p_half, z_full)

! call complete_robert_filter(tracer_attributes)! moved inside the spectral_dynamics: ZTAN 01/18/2012

if(dry_model) then
  call compute_pressures_and_heights(tg(:,:,:,future), psg(:,:,future), z_full, z_half, p_full, p_half)
else
  ! first do current time for p_prev vars to be used in conv/cond 
  call compute_pressures_and_heights( tg(:,:,:,current), psg(:,:,current), &
       z_full, z_half, p_full_prev, p_half_prev, grid_tracers(:,:,:,current,nhum))
  call compute_pressures_and_heights( &
     tg(:,:,:,future), psg(:,:,future), z_full, z_half, p_full, p_half, grid_tracers(:,:,:,future,nhum))
endif

! added fridoo sept 2012
!if(.not.dry_model) then
   if(atmos_water_correction) then
      atmos_water_future = mass_weighted_global_integral( &
           grid_tracers(:,:,:,future,nhum), psg(:,:,future))
      water_correction_factor = (atmos_water_previous + delta_t*evap_minus_precip)/atmos_water_future
      grid_tracers(:,:,:,future,nhum) = water_correction_factor*grid_tracers(:,:,:,future,nhum)
   endif
!endif
! end fridoo

call spectral_diagnostics(Time_next, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
                          tg(:,:,:,future), wg_full, grid_tracers(:,:,:,:,:), future)

previous = current 
current  = future

return 
end subroutine atmosphere

!=================================================================================================================================

subroutine atmosphere_end
integer :: ntr, nt
character(len=64) :: file, tr_name

if(.not.module_is_initialized) return

file='RESTART/atmosphere.res'
call nullify_domain()
!call write_data(trim(file), 'previous', previous) ! No interface exists to write a scalar
!call write_data(trim(file), 'current',  current)  ! No interface exists to write a scalar
call write_data(trim(file), 'time_pointers', (/real(previous),real(current)/)) ! getaround for no interface to write a scalar
do nt=1,num_time_levels
  call write_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'psg', psg(:,:,  nt), grid_domain)
  do ntr = 1,num_tracers
    tr_name = trim(tracer_attributes(ntr)%name)
    call write_data(trim(file), tr_name, grid_tracers(:,:,:,nt,ntr), grid_domain)
  enddo
enddo
call write_data(trim(file), 'wg_full', wg_full, grid_domain)

deallocate (p_half, z_half, p_full, z_full, wg_full, psg, ug, vg, tg, grid_tracers)
deallocate (dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers)
deallocate (deg_lat, rad_lat, rad_lat_2d)
deallocate (dt_bucket, filt)
deallocate (p_half_prev, p_full_prev)

call set_domain(grid_domain)

if(two_stream)      call radiation_end
if(lwet_convection) call qe_moist_convection_end
if(turb)            call gcm_vert_diff_end
if(mixed_layer_bc)  call mixed_layer_end(t_surf, bucket_depth)

call spectral_dynamics_end(tracer_attributes)
deallocate(tracer_attributes)

module_is_initialized = .false.

end subroutine atmosphere_end

!=================================================================================================================================

end module atmosphere_mod
