module mixed_layer_mod

!
! Implementation of mixed layer boundary condition
!

use                  fms_mod, only: set_domain, write_version_number, &
                                    mpp_pe, mpp_root_pe, error_mesg, FATAL, WARNING

use                  fms_mod, only: stdlog, check_nml_error, close_file,&
                                    open_namelist_file, stdout, file_exist, &
                                    read_data, write_data, open_file, &
                                    nullify_domain

use            constants_mod, only: HLV, PI, RHO_CP, CP_AIR, CP_OCEAN, OMEGA, RADIUS 

use         diag_manager_mod, only: register_diag_field, send_data

use       time_manager_mod,   only: time_type 

use           transforms_mod, only: get_deg_lat, grid_domain, get_lat_max

use            vert_diff_mod, only: surf_diff_type

use          mpp_domains_mod, only: mpp_global_field

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: mixed_layer.f90 $'
      
character(len=128) :: tagname= &
'$Name:  $'
character(len=128), parameter :: mod_name='mixed_layer'

!=================================================================================================================================

public :: mixed_layer_init, replicate, mrdnl_gradient, gauss_smooth_field, mixed_layer, mixed_layer_end

!=================================================================================================================================

logical :: evaporation = .true.
logical :: ekman_layer = .false.
logical :: gauss_smooth = .true.
real    :: qflux_amp = 0.0
real    :: qflux_width = 16.0  ! width of qflux region in degrees
real    :: depth = 40.0
real    :: depth_land = 1.0
logical :: load_qflux = .false.


namelist/mixed_layer_nml/ evaporation, qflux_amp, depth, depth_land, qflux_width, load_qflux, ekman_layer

!=================================================================================================================================


logical :: module_is_initialized =.false.
logical :: used

integer :: iter
integer, dimension(4) :: axes
integer ::                                                                    &
     id_t_surf,            &   ! surface temperature
     id_flux_lhe,          &   ! latent heat flux at surface
     id_flux_oceanq,       &   ! oceanic Q flux 
     id_flux_u,            &   ! surface stress
     id_flux_t                 ! sensible heat flux at surface

real, allocatable, dimension(:,:)   ::                                        &
     ocean_qflux,           &   ! Q-flux 
     rad_lat_2d                 ! latitude in radians 

real, allocatable, dimension(:)   ::     &
     deg_lat                            

real, allocatable, dimension(:,:)   ::                                        &
     gamma_t,               &   ! Used to calculate the implicit
     gamma_q,               &   ! correction to the diffusion in
     fn_t,                  &   ! the lowest layer
     fn_q,                  &   ! 
     en_t,                  &   !
     en_q,                  &   !
     alpha_t,               &   !
     alpha_q,               &   !
     alpha_lw,              &   !
     beta_t,                &   !
     beta_q,                &   !
     beta_lw,               &   !
     t_surf_dependence,     &   !
     corrected_flux,        &   !
     eff_heat_capacity,     &   ! Effective heat capacity
     delta_t_surf,          &   ! Increment in surface temperature
     depth_map                  ! 2d depth for vaiable heat capacity

real, allocatable, dimension(:,:)   ::                                        &
     rad_lat_replicate,         &
     rad_lat_diff_replicate,    &
     flux_u_replicate,          &
     flux_m_replicate,          &
     t_surf_replicate,          &
     rad_lat_global,            &
     rad_lat_diff_global,       &
     coriolis_global,           &
     flux_u_global,             &
     flux_m_global,             &
     t_surf_global,             &
     d_t_surf_d_lat_global,     &
     flux_h_integrand_global,   &
     flux_h_global,             &
     d_flux_h_cos_d_lat_global, &
     div_flux_h_global,         &
     smooth_div_flux_h_global  
      
real, allocatable, dimension(:)   ::                                          &
     flux_h_1d,        &
     flux_m_1d,        &
     flux_u_1d,        &
     t_surf_1d

real, allocatable, dimension(:,:)  ::                                         &
     rad_lat_del,      &
     rad_lat_diff                                                  

real, allocatable, dimension(:,:,:) ::                                                     &
     smooth_factor,    &
     smooth_field_X
    
real inv_cp_air


!=================================================================================================================================
contains
!=================================================================================================================================

subroutine mixed_layer_init(is, ie, js, je, num_levels, t_surf, bucket_depth, ocean_mask, axes, Time)

type(time_type), intent(in)       :: Time
real, intent(inout), dimension(:,:) :: t_surf
real, intent(inout), dimension(:,:,:) :: bucket_depth
real, intent(inout), dimension(:,:) :: ocean_mask
integer, intent(in), dimension(4) :: axes
integer, intent(in) :: is, ie, js, je, num_levels 

integer :: j
real    :: rad_qwidth
integer:: ierr, io, unit
integer :: lat_max

if(module_is_initialized) return

call write_version_number(version, tagname)

unit = open_namelist_file ()
ierr=1
do while (ierr /= 0)
  read  (unit, nml=mixed_layer_nml, iostat=io, end=10)
  ierr = check_nml_error (io, 'mixed_layer_nml')
enddo
10 call close_file (unit)

if ( mpp_pe() == mpp_root_pe() )   write (stdlog(), nml=mixed_layer_nml)

call get_lat_max(lat_max)

allocate(rad_lat_2d              (is:ie, js:je))
allocate(ocean_qflux             (is:ie, js:je))

allocate(deg_lat                 (js:je))

allocate(gamma_t                 (is:ie, js:je))
allocate(gamma_q                 (is:ie, js:je))
allocate(en_t                    (is:ie, js:je))
allocate(en_q                    (is:ie, js:je))
allocate(fn_t                    (is:ie, js:je))
allocate(fn_q                    (is:ie, js:je))
allocate(alpha_t                 (is:ie, js:je))
allocate(alpha_q                 (is:ie, js:je))
allocate(alpha_lw                (is:ie, js:je))
allocate(beta_t                  (is:ie, js:je))
allocate(beta_q                  (is:ie, js:je))
allocate(beta_lw                 (is:ie, js:je))
allocate(delta_t_surf            (is:ie, js:je))
allocate(eff_heat_capacity       (is:ie, js:je))
allocate(corrected_flux          (is:ie, js:je))
allocate(t_surf_dependence       (is:ie, js:je))
allocate(depth_map               (is:ie, js:je))

if (ekman_layer) then
allocate(flux_u_1d                  (js:je))
allocate(flux_m_1d                  (js:je))
allocate(t_surf_1d                  (js:je))
allocate(rad_lat_replicate          (is:ie, js:je))
allocate(rad_lat_diff_replicate     (is:ie, js:je))
allocate(flux_u_replicate           (is:ie, js:je))
allocate(flux_m_replicate           (is:ie, js:je))
allocate(t_surf_replicate           (is:ie, js:je))

allocate(rad_lat_global             (is:ie, 1:lat_max))
allocate(rad_lat_diff_global        (is:ie, 1:lat_max))
allocate(coriolis_global            (is:ie, 1:lat_max))
allocate(flux_u_global              (is:ie, 1:lat_max))
allocate(flux_m_global              (is:ie, 1:lat_max))
allocate(t_surf_global              (is:ie, 1:lat_max))
allocate(d_t_surf_d_lat_global      (is:ie, 1:lat_max)) 
allocate(flux_h_integrand_global    (is:ie, 1:lat_max))  
allocate(flux_h_global              (is:ie, 1:lat_max))
allocate(d_flux_h_cos_d_lat_global  (is:ie, 1:lat_max))
allocate(div_flux_h_global          (is:ie, 1:lat_max))
allocate(smooth_div_flux_h_global   (is:ie, 1:lat_max))

allocate(rad_lat_del                (is:ie, 1:lat_max))                                    
allocate(rad_lat_diff               (is:ie, 1:lat_max))                                     
allocate(smooth_factor              (is:ie, 1:lat_max, 1:lat_max))                                    
allocate(smooth_field_X             (is:ie, 1:lat_max, 1:lat_max)) 
endif

!
!see if restart file exists for the surface temperature
!
if (file_exist('INPUT/mixed_layer.res.nc')) then

   call nullify_domain()
   call read_data(trim('INPUT/mixed_layer.res'), 't_surf',   t_surf, grid_domain)
   call read_data(trim('INPUT/mixed_layer.res'), 'bucket_depth', bucket_depth, grid_domain)

else if (file_exist('INPUT/swamp.res')) then
         unit = open_file (file='INPUT/swamp.res', &
                           form='native', action='read')
         call read_data (unit, t_surf)
         call close_file (unit)
  call error_mesg('mixed_layer','mixed_layer restart file not found, using swamp restart file', WARNING)
else
  call error_mesg('mixed_layer','mixed_layer restart file not found', WARNING)
endif

id_t_surf = register_diag_field(mod_name, 't_surf',        &
                                axes(1:2), Time, 'surface temperature','K')
id_flux_t = register_diag_field(mod_name, 'flux_t',        &
                                axes(1:2), Time, 'sensible heat flux up at surface','watts/m2')
id_flux_lhe = register_diag_field(mod_name, 'flux_lhe',        &
                                 axes(1:2), Time, 'latent heat flux up at surface','watts/m2')
id_flux_oceanq = register_diag_field(mod_name, 'flux_oceanq',        &
                                 axes(1:2), Time, 'oceanic Q-flux','watts/m2')

! set up depth_map for spatially varying heat capacity  ! Added by LJJ
where (ocean_mask > 0.0)
  depth_map(:,:) = depth
elsewhere 
  depth_map(:,:) = depth_land
endwhere

! End LJJ addition

! latitude will be needed for oceanic q flux
call get_deg_lat(deg_lat)
do j=js,je
  rad_lat_2d(:,j) = deg_lat(j)*PI/180.
enddo

! calculate ocean Q flux
rad_qwidth = qflux_width*PI/180.
ocean_qflux = qflux_amp*(1-2.*rad_lat_2d**2/rad_qwidth**2) * &
        exp(- ((rad_lat_2d)**2/(rad_qwidth)**2))/cos(rad_lat_2d)


! zero Q flux over land         ! Added by LJJ
where (ocean_mask > 0.0)
     ocean_qflux(:,:) = ocean_qflux(:,:)
elsewhere
     ocean_qflux(:,:) = 0.0
endwhere

! End LJJ addition

! load Q flux 
if (load_qflux) then
  call read_data('INPUT/ocean_qflux.nc', 'ocean_qflux',  ocean_qflux)
endif

inv_cp_air = 1.0 / CP_AIR 

module_is_initialized = .true.

return
end subroutine mixed_layer_init

!=================================================================================================================================

subroutine replicate(ni, nf, field_2D, field_1D)

   integer :: i

   integer, intent(in) ::                                               &    
       ni,                                                              & 
       nf 

    real, dimension(:), intent(in) ::                                   &
       field_1D             ! 1d input field

    real, dimension(:,:), intent(inout) ::                              &
       field_2D             ! 2d output field

    field_2D = 0
    do i = 1, (nf-ni+1)
       field_2D(i,:) = field_1D(:)
    enddo

end subroutine replicate

subroutine mrdnl_gradient(field_X, rad_lat, d_field_X_d_lat)

    integer :: l, k

    real, dimension(:,:), intent(in) ::                                 &
         field_X             ! 2d input field

    real, dimension(:,:), intent(in) ::                                 &
         rad_lat             ! 2d latitude in radian

    real, dimension(:,:), intent(out) ::                                &
         d_field_X_d_lat     ! meridional gradient of field_X


    rad_lat_del  = 0
    rad_lat_diff = 0
    do l=2, (size(rad_lat,2)-1)
       rad_lat_del(:,l)   = (rad_lat(:,l+1)-rad_lat(:,l))/(rad_lat(:,l)-rad_lat(:,l-1))
       rad_lat_diff(:,l)  = (rad_lat(:,l+1)-rad_lat(:,l))
    enddo
    rad_lat_del(:,1)                = (rad_lat(:,2)-rad_lat(:,1))/(rad_lat(:,1)-rad_lat(:,size(rad_lat,2)))
    rad_lat_del(:,size(rad_lat,2))  = (rad_lat(:,1)-rad_lat(:,size(rad_lat,2)))/(rad_lat(:,size(rad_lat,2))-rad_lat(:,size(rad_lat,2)-1))
    rad_lat_diff(:,1)               = rad_lat(:,2) - rad_lat(:,1)  
    rad_lat_diff(:,size(rad_lat,2)) = rad_lat(:,1) - rad_lat(:,size(rad_lat,2))


    d_field_X_d_lat = 0
    do k=2, (size(field_X,2)-1)
       d_field_X_d_lat(:,k) = (field_X(:,k+1) - field_X(:,k)*(1-rad_lat_del(:,k)**2) - field_X(:,k-1)*(rad_lat_del(:,k)**2))/((1+rad_lat_del(:,k))*rad_lat_diff(:,k))
    enddo
end subroutine mrdnl_gradient

subroutine gauss_smooth_field(field_X, rad_lat, stand_dev, gauss_smooth_field_X)

    integer :: i, j

    real, intent(in) ::                                                    &
          stand_dev            

    real, dimension(:,:), intent(in) ::                                    &
         field_X             ! 2d input field

    real, dimension(:,:), intent(in) ::                                    &
         rad_lat             ! 2d latitude in radian

    real ::                                                                &
         coeff_max

    real, dimension(:,:), intent(out) ::                                   &
         gauss_smooth_field_X

   coeff_max             = 0
   smooth_factor         = 0
   smooth_field_X        = 0 
   gauss_smooth_field_X  = 0
    
   coeff_max = 1/(stand_dev*sqrt(2*pi))
   do j = 2, (size(field_X,2)-1)
      do i = 2, (size(field_X,2)-1)
         smooth_factor(:,i,j) = coeff_max*exp(-((rad_lat(:,i)-rad_lat(:,j))**2)/(2*(stand_dev**2))) &
                   *(rad_lat(:,i+1)-rad_lat(:,i-1))/2
         smooth_field_X(:,i,j) = field_X(:,i)*smooth_factor(:,i,j)
      enddo
      gauss_smooth_field_X(:,j) = sum(smooth_field_X(:,1:size(field_X,2),j),2)/sum(smooth_factor(:,1:size(field_X,2),j),2)
   enddo

end subroutine gauss_smooth_field


!=================================================================================================================================


!=================================================================================================================================

subroutine mixed_layer (                                               &
     is,                                                               &
     ie,                                                               &
     js,                                                               &
     je,                                                               &
     Time,                                                             &
     t_surf,                                                           &
     flux_t,                                                           &
     flux_q,                                                           &
     flux_r,                                                           &
     flux_u,                                                           &
     dt,                                                               &
     net_surf_sw_down,                                                 &
     surf_lw_down,                                                     &
     Tri_surf,                                                         &
     dhdt_surf,                                                        &
     dedt_surf,                                                        &
     dedq_surf,                                                        &
     drdt_surf,                                                        &
     dhdt_atm,                                                         &
     dedq_atm)         




! ---- arguments -----------------------------------------------------------
type(time_type), intent(in)       :: Time
real, intent(in),  dimension(:,:) :: &
     net_surf_sw_down, surf_lw_down
real, intent(in), dimension(:,:) :: &
     flux_t,    flux_q,     flux_r,    flux_u
real, intent(inout), dimension(:,:) :: t_surf
real, intent(in), dimension(:,:) :: &
   dhdt_surf, dedt_surf, dedq_surf, &
   drdt_surf, dhdt_atm, dedq_atm  
real, intent(in) :: dt
type(surf_diff_type), intent(inout) :: Tri_surf
integer, intent(in) :: is, ie, js, je
integer :: lat_max, nh_lath, sh_lath, j
real :: rad_lat_std

if(.not.module_is_initialized) then
  call error_mesg('mixed_layer','mixed_layer module is not initialized',FATAL)
endif

! Need to calculate the implicit changes to the lowest level delta_q and delta_t
! - see the discussion in vert_diff.tech.ps
                                                                                                                                    
! Care is needed to differentiate between the sensible heat flux and the
! diffusive flux of temperature
                                                                                                                                    
gamma_t = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_t + dhdt_atm * inv_cp_air))
gamma_q = 1.0 / (1.0 - Tri_surf%dtmass * (Tri_surf%dflux_q + dedq_atm))
                                                                                                                                 
fn_t = gamma_t * (Tri_surf%delta_t + Tri_surf%dtmass * flux_t * inv_cp_air)
fn_q = gamma_q * (Tri_surf%delta_q + Tri_surf%dtmass * flux_q)
                                                                                                                                 
en_t = gamma_t * Tri_surf%dtmass * dhdt_surf * inv_cp_air
en_q = gamma_q * Tri_surf%dtmass * dedt_surf
                                                                                                                                    
!
! Note flux_sw doesn't depend on surface or lowest layer values
! Note drdt_atm is not used - should be fixed
!
alpha_t = flux_t * inv_cp_air + dhdt_atm * inv_cp_air * fn_t
alpha_q = flux_q + dedq_atm * fn_q
alpha_lw = flux_r
                                                                                                                                 
beta_t = dhdt_surf * inv_cp_air + dhdt_atm * inv_cp_air * en_t
beta_q = dedt_surf + dedq_atm * en_q
beta_lw = drdt_surf

!##########################################################################################
!                                START OCEAN QFLUX
!##########################################################################################
if (ekman_layer) then

! INITIALIZE
ocean_qflux = 0
flux_u_1d = 0
flux_m_1d = 0
t_surf_1d = 0
flux_u_replicate = 0
flux_m_replicate = 0
t_surf_replicate = 0
rad_lat_global            = 0
rad_lat_diff_global       = 0
flux_u_global             = 0
flux_m_global             = 0
t_surf_global             = 0
d_t_surf_d_lat_global     = 0
flux_h_integrand_global   = 0
flux_h_global             = 0 
d_flux_h_cos_d_lat_global = 0
div_flux_h_global         = 0 
smooth_div_flux_h_global  = 0 

call get_lat_max(lat_max)

call mpp_global_field(grid_domain, rad_lat_2d, rad_lat_global)
do j=2,(lat_max-1)
   rad_lat_diff_global(:,j)  = (rad_lat_global(:,j+1)-rad_lat_global(:,j-1))/2
enddo
coriolis_global = 2.0*omega*sin(rad_lat_global)    

! A. Fields are zonally averaged, zonally replicated, and carried globally
flux_u_1d  = sum(flux_u,1)/size(flux_u,1)
flux_m_1d  = flux_u_1d / coriolis_global(1,:) 
t_surf_1d  = sum(t_surf,1)/size(t_surf,1)

call replicate(is, ie, flux_u_replicate, flux_u_1d)
call mpp_global_field(grid_domain, flux_u_replicate, flux_u_global)

call replicate(is, ie, flux_m_replicate, flux_m_1d)
call mpp_global_field(grid_domain, flux_m_replicate, flux_m_global)

call replicate(is, ie, t_surf_replicate, t_surf_1d)
call mpp_global_field(grid_domain, t_surf_replicate, t_surf_global)

call mrdnl_gradient(t_surf_global, rad_lat_global, d_t_surf_d_lat_global)

! B. Determined latitude where surface wind changes sign in:
! a. The Northern Hemisphere
nh_lath = lat_max/2+1 
do j=(lat_max-1),(lat_max/2+1),-1
   if ((sign(1.0, flux_u_global(1,j+1))  .ne.      &
         sign(1.0, flux_u_global(1,j)))   .and.     &
             (flux_u_global(1,j)>0)         .and.     &                 
       (j .ge. ((lat_max/2+1) + lat_max/12)))  then
      nh_lath = j+1
   endif
enddo
!b. The Southern Hemisphere
sh_lath = lat_max/2
do j = 2,(lat_max/2),1
   if ((sign(1.0, flux_u_global(1,j-1)) .ne.      &
        sign(1.0, flux_u_global(1,j)))   .and.     &
            (flux_u_global(1,j)>0)         .and.     &   
      (j .le. (lat_max/2 - lat_max/12)))      then
        sh_lath = j-1
   endif
enddo

! C. Compute integrand of enthalpy flux in the tropics
do j = (sh_lath+1), (nh_lath-1)
   flux_h_integrand_global(:,j) = cp_ocean*flux_u_global(:,j)/coriolis_global(:,j)*d_t_surf_d_lat_global(:,j)
enddo

! D. Compute enthalpy flux in the tropics 
do j = (sh_lath), (lat_max/2)
    flux_h_global(:,j) =  sum(flux_h_integrand_global(:,sh_lath:j)*rad_lat_diff_global(:,sh_lath:j),2)
enddo
do j = (lat_max/2+1), (nh_lath)
    flux_h_global(:,j) = - sum(flux_h_integrand_global(:,j:nh_lath)*rad_lat_diff_global(:,j:nh_lath),2)
enddo

! E. Compute enthalpy flux divergence
call mrdnl_gradient(flux_h_global*cos(rad_lat_global), rad_lat_global, d_flux_h_cos_d_lat_global)
do j=(sh_lath), (nh_lath)
   div_flux_h_global(:,j) = 1/(radius*cos(rad_lat_global(:,j)))*d_flux_h_cos_d_lat_global(:,j) 
enddo

! F. Gaussian smooting function
rad_lat_std = 7*PI/180
if (gauss_smooth) then
   call gauss_smooth_field(div_flux_h_global, rad_lat_global, rad_lat_std, smooth_div_flux_h_global)
else
   smooth_div_flux_h_global = div_flux_h_global
endif
 
! G. Set flux_oceanq 
do j= js, je
   ocean_qflux(:,j) = smooth_div_flux_h_global(:,j)
enddo

endif
!write(*,*) 'smooth_div_flux_h', smooth_div_flux_h_global(1,:)
!write(*,*) 'sh_lath', sh_lath, 'nh_lath', nh_lath
!##########################################################################################
!                                 END OCEAN QFLUX
!##########################################################################################


!
! Implement mixed layer surface boundary condition
!
corrected_flux = - net_surf_sw_down - surf_lw_down + alpha_t * CP_AIR + alpha_lw + ocean_qflux
t_surf_dependence = beta_t * CP_AIR + beta_lw


if (evaporation) then
  corrected_flux = corrected_flux + alpha_q * HLV
  t_surf_dependence = t_surf_dependence + beta_q * HLV
endif

!
! Now update the mixed layer surface temperature using an implicit step
!
!eff_heat_capacity = depth * RHO_CP + t_surf_dependence * dt
eff_heat_capacity = depth_map * RHO_CP + t_surf_dependence * dt

if (any(eff_heat_capacity .eq. 0.0))  then 
  write(*,*) 'mixed_layer: error', eff_heat_capacity
  call error_mesg('mixed_layer', 'Avoiding division by zero',fatal)
end if

delta_t_surf = - corrected_flux  * dt / eff_heat_capacity

t_surf = t_surf + delta_t_surf
                                                                                                                                    
!
! Finally calculate the increments for the lowest atmospheric layer
!
Tri_surf%delta_t = fn_t + en_t * delta_t_surf
Tri_surf%delta_q = fn_q + en_q * delta_t_surf


!
! Note:
! When using an implicit step there is not a clearly defined flux for a given timestep
!
if(id_t_surf > 0) used = send_data(id_t_surf, t_surf, Time)
if(id_flux_t > 0) used = send_data(id_flux_t, flux_t, Time)
if(id_flux_lhe > 0) used = send_data(id_flux_lhe, HLV * flux_q, Time)
if(id_flux_oceanq > 0)   used = send_data(id_flux_oceanq, ocean_qflux, Time)

end subroutine mixed_layer

!=================================================================================================================================

subroutine mixed_layer_end(t_surf, bucket_depth)

real, intent(inout), dimension(:,:) :: t_surf
real, intent(inout), dimension(:,:,:) :: bucket_depth
integer:: unit

if (ekman_layer) then

deallocate(flux_u_replicate)
deallocate(flux_m_replicate)
deallocate(t_surf_replicate)

deallocate(rad_lat_global)
deallocate(rad_lat_diff_global)
deallocate(coriolis_global)
deallocate(flux_u_global)
deallocate(flux_m_global)
deallocate(t_surf_global)

deallocate(d_t_surf_d_lat_global)

deallocate(flux_h_integrand_global)
deallocate(flux_h_global)          
deallocate(d_flux_h_cos_d_lat_global)
deallocate(div_flux_h_global)        
deallocate(smooth_div_flux_h_global) 

endif


if(.not.module_is_initialized) return

! write a restart file for the surface temperature
call nullify_domain()
call write_data(trim('RESTART/mixed_layer.res'), 't_surf',   t_surf, grid_domain)
call write_data(trim('RESTART/mixed_layer.res'), 'bucket_depth', bucket_depth, grid_domain)

module_is_initialized = .false.

end subroutine mixed_layer_end

!=================================================================================================================================

end module mixed_layer_mod
