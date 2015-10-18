module vars

! allocatable arrays and namelist variables
  
  implicit none

! --------------------------- coordinate variables ---------------------------

  real, allocatable, dimension(:) ::                                          &
       sigma,              &   ! vertical coordinate axis (full levels)
       deg_lat,            &   ! latitude in degrees [local_num_lat]
       sin_lat,            &   ! sin(latitude)
       pot_temp_coord,        &   ! theta coodinate axis (degrees)
       d_pot_temp_coord,       &   ! spacing between theta levels (degrees)
       pk,                 &   ! pressure hybrid coefficients
       rhum_bin,           &   ! relative humidity pdf bin axis
       bk                      ! sigma hybrid coefficients

  real, allocatable, dimension(:,:) ::                                         &
       cos_lat_xy,         &   ! xy cosine(latitude) array
       sin_lat_xy,         &   ! xy   sine(latitude) array
       cos_lat_yz              ! yz cosine(latitude) array

!------------------ gridpoint variables on pressure surfaces -------------------

  real, allocatable, dimension(:,:,:) ::                                      &
       u_grid,             &   ! zonal wind 
       shum_grid,          &   ! specific humidity
       dt_shum_cond_grid,    &   ! large scale cond specific humidity tendency
       dt_shum_conv_grid,    &   ! convection specific humidity tendency
       dt_shum_diff_grid,    &  
       dt_shum_cond_counter, &   ! counter for large scale cond specific humidity tendency
       dt_shum_conv_counter, &   ! counter for convection specific humidity tendency
       dt_temp_cond_grid,    &  
       dt_temp_conv_grid,    &  
       dt_temp_diff_grid,    &  
       dt_temp_rad_grid,     &
       dt_temp_sw_grid,      &
       dt_u_diff_grid,       &  ! zonal wind tendency from the turbulent mixing;      ljj addition
       dt_v_diff_grid,       &  ! meridional wind tendency from the turbulent mixing; ljj addition
       dt_u_diss_grid,       &  ! zonal wind tendency from the spectral damping;      ljj addition
       dt_v_diss_grid,       &  ! meridional wind tendency from the spectral damping; ljj addition
       dt_temp_diss_grid,    &  ! temperature tendency from the spectral damping;     ljj addition
       v_grid,             &   ! meridional wind
       vor_grid,           &   ! vorticity
       div_grid,           &   ! divergence
       temp_grid,          &   ! temperature
       virtual_temp_grid,  &   ! virtual temperature
       vcos_grid,          &   ! cos(lat) * v
       w_grid,             &   ! vertical velocity (full levels)
       dpdt_grid,          &   ! pressure tendency (full levels)
       p_half,             &   ! pressure on half levels
       p_full,             &   ! pressure on full levels
       sfctn_grid,         &   ! streamfunction
       hor_sfctn_grid,     &   ! horizontal streamfunction
       pot_temp_grid,              &   ! potential temperature
       pot_temp_e_grid,            &   ! equivalent potential temperature
       pot_temp_e_sat_grid,        &   ! saturated equivalent potential temperature
       geopot_grid,        &   ! geopotential
       mont_grid,          &   ! Montgomery streamfunction
       dx_mont_grid,       &   ! x derivative of Montgomery streamfunction
       sat_mr_grid,        &   ! saturation mixing ratio
       sat_shum_grid,      &   ! saturation specific humidity
       rhum_grid               ! relative humidity

  real, allocatable, dimension(:,:)  ::                                       &
       ts_grid,            &   ! surface temperature
       ps_grid,            &   ! surface pressure
       geopot_sfc,         &   ! surface geopotential
       bt_sfctn_grid,      &     ! barotropic streamfunction
       conv_zon_spec_avg,  &   ! conversion as func. of wavenumber (vert. avg.)
       ucos_zon_spec_barotr,  &   ! u*cos(lat) as func. of wavenumber (vert. avg.)
       vcos_zon_spec_barotr,  &   ! v*cos(lat) as func. of wavenumber (vert. avg.)
       shum_zon_spec_barotr,  &   ! shum as func. of wavenumber (vert. avg.)
       precip_cond_grid,   &   ! large scale cond precipitation
       precip_conv_grid,   &   ! convection precipitation
       sfc_flux_lh_grid,   &
       sfc_flux_sh_grid,   &
       sfc_flux_sw_grid,   & 
       sfc_flux_lwd_grid,  &
       sfc_flux_lwu_grid,  &
       sfc_flux_ocean_grid,&
       toa_flux_sw_grid,   & 
       toa_flux_lwu_grid,  &
       drag_coeff_mo_grid, &
       drag_coeff_lh_grid, &
       drag_coeff_sh_grid, &
       bucket_depth_grid,  &     ! bucket water depth on continents
       bucket_depth_conv_grid, & ! bucket water depth on continents tendency (convection tendency)  
       bucket_depth_cond_grid, & ! bucket water depth on continents tendency (condensation tendency)  
       bucket_depth_LH_grid,   & ! bucket water depth on continents tendency (latent heat)  
       bucket_diffusion_grid,  & ! bucket water depth on continents tendency (horizontal diffusion)   
       precip_daily,       &   ! daily precipitation 
       precip_cond_daily,  &   ! daily large scale cond precipitation 
       precip_conv_daily,  &   ! daily convection precipitation 
       one_array               ! array of ones

  real, allocatable, dimension(:,:)  ::                                       &
       mont_sfc,     &   ! surface Montgomery streamfunction
       dx_mont_sfc,  &   ! derivative of surf. Mont. sf. WRT longitude
       dy_mont_sfc,  &   ! derivative of surf. Mont. sf. WRT longitude
       dx_pot_temp_sfc, &   ! derivative of surf. pot. temp. WRT longitude
       dy_pot_temp_sfc,  &           ! derivative of surf. pot. temp. WRT latitude
       bt_non_lin_spec,       &   ! non-linear interaction spectrum (barotropic)
       bt_non_lin_eddy_spec       ! exclude interaction with m=0 modes (barotropic)

  logical, allocatable, dimension(:,:)  ::                                       &
       precip_mask             ! mask for determining wet days

  real, allocatable, dimension(:, :, :) ::                                    &
       pot_temp_spec_avg,         &   ! potential temperature spectrum (all levels)
       pot_temp_e_spec_avg,       &   ! equivalent potential temperature spectrum (all levels)
       shum_spec_avg,             &   ! specific humidity spectrum (all levels)
       sat_shum_spec_avg,         &   ! sat specific humidity spectrum (all levels)
       rhum_spec_avg,             &   ! relative humidity spectrum (all levels)
       eape2eke_conv_spec_avg,    &   ! conv pe->ke spectrum (all levels)
       pot_temp_spec,         &       ! potential temperature spectrum (all levels)
       pot_temp_e_spec,       &       ! equivalent potential temperature spectrum (all levels)
       shum_spec,             &       ! specific humidity spectrum (all levels)
       sat_shum_spec,         &       ! sat specific humidity spectrum (all levels)
       rhum_spec,             &       ! relative humidity spectrum (all levels)
       eape2eke_conv_spec,    &       ! conv pe->ke spectrum (all levels)
       coriolis_spec,      &   ! spectrum related to beta effect (all levels)
       temp_spec,          &   ! Temperature variance spectrum
       non_lin_spec,       &   ! non-linear interaction spectrum
       non_lin_eddy_spec,  &     ! exclude interaction with m=0 modes
       t_non_lin_spec,       &   ! non-linear interaction spectrum (temperature)
       t_non_lin_eddy_spec       ! exclude interaction with m=0 modes (temperature)

  real, allocatable, dimension(:, :, :) ::                                    &
       rhum_pdf                ! pdf of relative humidity (lat, level, bin)

  real, allocatable, dimension(:, :) ::                                       &
       bc_sfctn_spec,      &   ! baroclinic energy spectrum
       bt_sfctn_spec           ! barotropic energy spectrum

  
!------------------- gridpoint variables on theta surfaces ---------------------

  real, allocatable, dimension(:, :, :) ::                                    &
       p_isent_grid_nw,        &   ! pressure
       dens_isent_grid_nw,     &   ! density
       vor_isent_grid,      &   ! isentropic potential vorticity (PV)
       vor_var_isent_grid,     &
       mrdnl_vor_flux_isent_grid,         &
       abs_vor_isent_grid_nw,     &   ! isentropic potential vorticity (PV)
       abs_vor_var_isent_grid, &
       mrdnl_abs_vor_flux_isent_grid,     &
       pot_vor_isent_grid,     &   ! isentropic potential vorticity (PV)
       v_pot_square_vor_isent_grid, & ! added fridoo
       pot_square_vor_isent_grid, & 
       pot_vor_var_isent_grid, &
       mrdnl_pot_vor_flux_isent_grid,     &
       u_isent_grid,           &   ! zonal wind
       u_var_isent_grid,       &
       v_isent_grid,           &   ! meridional wind
       vcos_isent_grid,        &   ! meridional wind      
       vcos_isent_grid_nw,     &   ! meridional wind      
       v_var_isent_grid,       &
       v_geostr_isent_grid,    &
       shum_isent_grid,        &   ! specific humidity
       mrdnl_shum_flux_isent_grid,       &   ! specific humidity
       rhum_isent_grid,        &   ! relative humidity
       sat_shum_isent_grid,    &   ! saturation specific humidity
       mrdnl_sat_shum_flux_isent_grid,   &   ! saturation specific humidity * meridional velocity
       sfctn_isent_grid,       &   ! streamfunction
       mont_isent_grid,        &   ! montgomery streamfunction
       mont_var_isent_grid,    &   ! montgomery streamfunction
       dx_mont_isent_grid          ! x derivative of montgomery streamfunction


  integer, allocatable, dimension(:, :) ::                                    &
       k_sfc,              &   ! index of surface potential temperature
       k_asfc                  ! index of lowest isentropic layer above surface

! --------- zonally averaged gridpoint variables on pressure surfaces ---------

  real, allocatable, dimension(:, :) ::                                       &
       p_half_zon,               &   ! pressure on half levels
       p_full_avg,               &   ! pressure on full levels
       specific_vol_avg,         &   ! specific volume
       dens_avg,                 &   ! specific volume
       u_avg,                    &   ! zonal wind
       u_var_avg,                &   ! zonal wind variance
       mrdnl_u_flux_avg,         &   ! u*cos(latitude)*v
       vrtcl_u_flux_avg,         &   ! w*u ; w =  \dot\sigma
       mrdnl_eddy_u_flux_avg,    &   ! u*cos(latitude)*v
       vrtcl_eddy_u_flux_avg,    &   ! w*u ; w =  \dot\sigma
       v_avg,                    &   ! meridional wind
       vcos_avg,                 &   ! cosine(latitude)*meridional wind
       v_var_avg,                &   ! meridional wind variance
       sfctn_avg,                &   ! streamfunction
       TEM_res_circ,             &   ! TEM residual circulation
       MTEM_res_circ,            &   ! Modified TEM residual circulation
       vort_avg,                 &   ! meridional wind
       w_avg,                    &   ! vertical wind (w =  \dot\sigma)
       w_var_avg,                &   ! vertical wind variance
       dpdt_avg,                 &   ! dp/dt, the vertical velocity in the pressure coordinate
       conv_pe_ke_avg,           &   ! conversion rate of potential energy to kinetic energy
       conv_eddy_pe_ke_avg,      &   ! eddy conversion rate of potential energy to kinetic energy
       temp_avg,                 &   ! temperature
       temp_var_avg,             &   ! temperature variance
       mrdnl_temp_flux_avg,      &   ! zonal mean [ psfc * v * cos(lat) * temp ]
       vrtcl_temp_flux_avg,      &   ! vertical temperature flux
       mrdnl_eddy_temp_flux_avg, &   ! meridional eddy temperature flux
       vrtcl_eddy_temp_flux_avg, &   ! vertical eddy temperature flux
       virtual_temp_avg,         &   ! virtual temperature
       pot_temp_avg,                 &   ! potential temperature
       pot_temp_var_avg,             &   ! potential temperature variance
       mrdnl_pot_temp_flux_avg,      &   ! zonal mean [ psfc * v * cos(lat) * pot_temp ]
       vrtcl_pot_temp_flux_avg,      &   ! vertical heat flux
       mrdnl_eddy_pot_temp_flux_avg, &  ! meridional eddy heat flux
       vrtcl_eddy_pot_temp_flux_avg, &  ! vertical eddy heat flux
       d_dy_pot_temp_avg,            &   ! derivative of pot. temp. WRT latitude
       d_dsig_pot_temp_avg,          &   ! derivative of pot. temp. WRT sigma
       d_dp_pot_temp_avg,            &   ! derivative of pot. temp. WRT pressure
       buoyancy_freq_avg,            &    ! Brunt-Vaeisaelae frequency
       z_avg,                    &   ! geopotential height
       z_var_avg,                &   ! geopotential height variance
       mrdnl_z_flux_avg,         &   ! zonal mean [ psfc * v * cos(lat) * z ]
       vrtcl_z_flux_avg,         &   ! vertical z flux 
       mrdnl_eddy_z_flux_avg,    &   ! meridional eddy z flux
       vrtcl_eddy_z_flux_avg,    &   ! vertical eddy z flux
       shum_avg,                 &   ! specific humidity
       shum_var_avg,             &   ! specific humidity variance
       mrdnl_shum_flux_avg,      &   ! specific humidity * meridional velocity
       vrtcl_shum_flux_avg,      &   ! specific humidity * vertical velocity (\dot\sigma)
       mrdnl_eddy_shum_flux_avg, &   ! specific humidity * meridional velocity
       vrtcl_eddy_shum_flux_avg, &   ! specific humidity * vertical velocity (\dot\sigma)
       sat_shum_avg,             &   ! saturation specific humidity
       sat_shum_var_avg,         &   ! saturation specific humidity variance
       rhum_avg,                 &   ! relative humidity
       pot_temp_e_avg,              &   ! equivalent potential temperature
       pot_temp_e_var_avg,          &   ! equivalent potential temperature variance
       mrdnl_pot_temp_e_flux_avg,   &   ! zonal mean [ psfc * v * cos(lat) * theta_e ]
       vrtcl_pot_temp_e_flux_avg,   &   ! vertical equiv. pot. temperature flux
       mrdnl_eddy_pot_temp_e_flux_avg,  &   ! meridional eddy equiv. pot. temp flux
       vrtcl_eddy_pot_temp_e_flux_avg,  &   ! vertical equiv. pot. temp flux
       d_dy_pot_temp_e_avg,         &   ! derivative of equiv. pot. temp. WRT latitude
       d_dsig_pot_temp_e_avg,       &   ! derivative of equiv. pot. temp. WRT sigma
       d_dp_pot_temp_e_avg,         &   ! derivative of equiv. pot. temp. WRT pressure
       pot_temp_e_sat_avg,          &   ! saturated equivalent potential temperature
       dt_shum_cond_avg,         &   ! large scale cond specific humidity tendency
       dt_shum_conv_avg,         &   ! convection specific humidity tendency
       dt_shum_diff_avg,         &  
       shum_cond_prob_avg,       &   ! large scale cond specific humidity tendency frequency
       shum_conv_prob_avg,       &  ! convection specific humidity tendency frequency
       dt_temp_cond_avg,         &
       dt_temp_conv_avg,         &
       dt_temp_diff_avg,         &
       dt_temp_rad_avg,          &
       dt_temp_sw_avg,           &
       dt_u_diff_avg,            &  ! zonal wind tendency from turbulent diffusion
       dt_v_diff_avg,            &  ! meridional wind tendency from turbulent diffusion
       dt_u_diss_avg,            &  ! zonal wind tendency from numerical dissipation (spectral damping)
       dt_v_diss_avg,            &  ! meridional wind tendency from numerical dissipation (spectral damping)
       dt_temp_diss_avg,         &  ! temperature tendency from numerical dissipation (spectral damping)
       ke_diss_zonal_avg,        &  ! zonal kinetic energy dissipation from spectral damping
       ke_diss_mrdnl_avg,        &  ! meridional kinetic energy dissipation from spectral damping
       ke_diff_zonal_avg,        &  ! zonal kinetic energy dissipation from turbulent mixing
       ke_diff_mrdnl_avg            ! meridional kinetic energy dissipation from turbulent mixing 

  real, allocatable, dimension(:) ::                                          &
       ts_avg,                &   ! surface temperature
       ps_avg,                &   ! surface pressure
       u_barotr_avg,            &   ! barotropic zonal velocity
       u_barotr_var_avg,        &   ! variance of barotropic zonal velocit
       v_barotr_avg,            &   ! barotropic meridional velocity
       v_barotr_var_avg,        &   ! variance of barotropic meridional velocity
       precip_cond_avg,       &   ! large scale cond precipiation
       precip_conv_avg,       &   ! convection precipiation
       sfc_flux_lh_avg,       &
       sfc_flux_sh_avg,       &
       sfc_flux_sw_avg,       & 
       sfc_flux_lwd_avg,      &
       sfc_flux_lwu_avg,      &
       sfc_flux_ocean_avg,    &
       toa_flux_sw_avg,       & 
       toa_flux_lwu_avg,      &
       bucket_depth_avg,      & ! bucket water depth on continents
       bucket_depth_conv_avg, & ! bucket water depth on continents (convection tendency)  
       bucket_depth_cond_avg, & ! bucket water depth on continents (condensation tendency)  
       bucket_depth_LH_avg,   & ! bucket water depth on continents (condensation tendency)  
       bucket_diffusion_avg,  & ! bucket water depth on continents (horizontal diffusion)
       drag_coeff_mo_avg,     &
       drag_coeff_lh_avg,     &
       drag_coeff_sh_avg,     &
       precip_daily_above_threshold_avg,      &  
       days_total,                            &
       days_above_threshold,                  &
       precip_daily_above_threshold_prob_avg  
! -----------------------------------------------------------------------------
! ---------- zonally averaged gridpoint variables on theta surfaces  ----------

  real, allocatable, dimension(:) ::                                          &
      mont_dx_pot_temp_sfc_avg,                                                  &
      mrdnl_geostr_eddy_pot_temp_flux_sfc_avg                                    

     
  real, allocatable, dimension(:, :) ::                                       &
       sfctn_isent_avg,        &   ! streamfunction
       pot_vor_isent_avg,          &   ! \avg{density * IPV}
       pot_vor_var_isent_avg,          &   ! \avg{density * IPV}
       pot_enstrophy_isent_avg,            &
       v_pot_enstrophy_isent_avg,         &
       mrdnl_eddy_pot_enstrophy_flux_isent_avg,            &
       mrdnl_pot_vor_flux_isent_avg,          &   ! \avg{density * IPV}
       mrdnl_eddy_pot_vor_flux_isent_avg,          &   ! \avg{density * IPV}
       p_isent_avg_nw,            &   ! pressure
       p_var_isent_avg_nw,        &   ! variance of pressure
       dens_isent_avg_nw,      &   ! density
       exner_isent_avg_nw,        &   ! exner function
       p_exner_isent_avg_nw,       &   ! pressure * (exner function)
       sfc_pot_temp_pdf,          &   ! frequency distribution of surface temperatures
       v_isent_avg,            &   ! meridional wind
       vcos_isent_avg,         &   ! meridional wind
       v_geostr_isent_avg,     & 
       u_isent_avg,            &   ! zonal wind
       shum_isent_avg,         &   ! specific humidity
       rhum_isent_avg,         &   ! relative humidity
       sat_shum_isent_avg,     &   ! saturation specific humidity
       mrdnl_shum_flux_isent_avg,        &   ! specific humidity * meridional velocity
       mrdnl_eddy_shum_flux_isent_avg,        &   ! specific humidity * meridional velocity
       mrdnl_sat_shum_flux_isent_avg,    &   ! saturation specific humidity * meridional velocity
       mrdnl_eddy_sat_shum_flux_isent_avg,    &   ! saturation specific humidity * meridional velocity
       v_var_isent_avg,           &   ! (meridional wind)^2
       u_var_isent_avg,           &   ! (zonal wind)^2
       mont_isent_avg,         &   ! montgomery streamfunction
       mont_var_isent_avg,     &   ! variance of montgomery streamfunction
       dx_mont_isent_avg            ! x derivative of montgomery streamfunction

! ---------------------------- namelist variables -----------------------------

  integer, parameter ::                                                       &
       max_len = 300           ! maximum length of namelist string variables
                    
  integer ::                                                                  &
       DataIn,             &   ! input momentum fields (u_and_v or vor_and_div)
       data_source,        &   ! source of data (FMS=1, NCEP=2, CCM3=3)
       num_fourier,        &   ! number of Fourier waves
       num_segments = 1,   &   ! number of segments (e.g., pentads)        
       num_bin,            &   ! number of bins for relative humidity pdf
       MaxIsentrLev            ! number of isentropic levels

  logical ::                                                                  &
       smooth_surface,     &   ! include topography
       is_gaussian,        &   ! .true. for data on Gaussian grid
       isentrope,          &
       moisture,           &   ! flag to analyze moisture
       virtual,            &   ! flag for virtual temp effect for z and svol
       bucket,             &   ! flag to analyse hydrology
       moist_isentropes        ! flag to analyze data on moist isentropes
  
  real ::                                                                     &
       TimeIn = 0.,        &   ! start time of analysis       
       PotTempMin,         &   ! minimum isentropic level
       PotTempMax,         &   ! maximum isentropic level
       delta_t,            &   ! timestep (seconds) for diabatic diagnostics
       roughness,          &   ! roughness length for momentum
       precip_daily_threshold        ! threshold for wet days

  character (len = max_len) ::                                                &
       InputFileName,             &   ! input file name
       SGInputFileName,           &   ! surface geopotential input file
       UVarName,                  &   ! zonal wind variable name
       VVarName,                  &   ! meridional wind variable name
       VORVarName,                &   ! vorticity variable name
       DIVVarName,                &   ! divergence variable name
       TEMPVarName,               &   ! temperature variable name
       TSVarName,                 &   ! surface temperature variable name
       PSVarName,                 &   ! surface pressure variable name
       SGVarName,                 &   ! surface geopotential variable name
       ShumVarName,               &   ! specific humidity variable name
       CondVarName,               &   ! large scale condensation tendency variable name
       ConvVarName,               &   ! convection tendency variable name
       DiffVarName,               &   ! diffusion tendency variable name
       PrecipCondVarName,         &   ! large scale condensation precipitation variable name
       PrecipConvVarName,         &   ! convection precipitation variable name
       SfcFluxLHVarName,          &
       SfcFluxSHVarName,          &
       SfcFluxSWVarName,          &
       SfcFluxLWDVarName,         &
       SfcFluxLWUVarName,         &
       SfcFluxQFVarName,          &
       ToaFluxSWVarName,          &
       ToaFluxLWUVarName,         &
       DiabCondVarName,           &
       DiabConvVarName,           &
       DiabDiffVarName,           &
       DiabRadVarName,            &
       DiabSWVarName,             &
       UMomDiffVarName,           &   ! ljj addition
       VMomDiffVarName,           &   ! ljj addition
       UMomDissVarName,           &   ! ljj addition
       VMomDissVarName,           &   ! ljj addition
       DiabDissVarName,           &   ! ljj addition
       BucketDepthVarName,        & 
       BucketDepthConvVarName,    & 
       BucketDepthCondVarName,    & 
       BucketDepthLHVarName,    & 
       BucketDiffusionVarName,    & 
       DragMOVarName,             &
       DragLHVarName,             &
       DragSHVarName,             &
       OutputFileName                 ! output file name

  namelist /main_list/                                                        &
                              DataIn,                     data_source,        &
                         num_fourier,                    MaxIsentrLev,        &
                      smooth_surface,                          TimeIn,        &
                          PotTempMin,                      PotTempMax,        &
                            UVarName,                        VVarName,        &
                          VorVarName,                      DivVarName,        &
                         TempVarName,                       TSVarName,        &
                           PSVarName,                                         &
                           SGVarName,                                         &
                         ShumVarName,                 SgInputFileName,        &
                         CondVarName,                     ConvVarName,        &
                         DiffVarName,                                         &
                   PrecipCondVarName,               PrecipConvVarName,        &
                    SfcFluxLHVarName,                SfcFluxSHVarName,        &
                    SfcFluxSWVarName,               SfcFluxLWDVarName,        &
                   SfcFluxLWUVarName,                SfcFluxQFVarName,        &
                    ToaFluxSWVarName,               ToaFluxLWUVarName,        &
                     DiabCondVarName,                 DiabConvVarName,        &
                     DiabDiffVarName,                                         &
                      DiabRadVarName,                   DiabSWVarName,        & 
                     UMomDiffVarName,                 VMomDiffVarName,        & ! ljj addition
                     UMomDissVarName,                 VMomDissVarName,        & ! ljj addition
                     DiabDissVarName,                                         & ! ljj addition
                      BucketDepthVarName,                                     & 
                      BucketDepthConvVarName,                                 & 
                      BucketDepthCondVarName,                                 & 
                      BucketDepthLHVarName,                                   & 
                      BucketDiffusionVarName,                                 & 
                       DragMOVarName,                   DragLHVarName,        &
                       DragSHVarName,                                         &
                         is_gaussian,                        moisture,        &
                    moist_isentropes,                         num_bin,        &
                           isentrope,          precip_daily_threshold,        &
                             virtual,                          bucket,        &
                        num_segments


  namelist /filename_list/                                                    &
                       InputFileName,                  OutputFileName

  contains

!------------------------------------------------------------------------------

    subroutine allocate_and_initialize(num_lat, num_lon, num_lev)

      integer, intent(in) :: num_lat, num_lon, num_lev
      integer             :: i                         ! loop over rhum pdf bins

! ------------------------------ axis variables -------------------------------

      allocate(      deg_lat(num_lat))
      allocate(      sin_lat(num_lat))
      allocate(        sigma(num_lev))
      allocate(     rhum_bin(num_bin))
      allocate(           pk(num_lev+1))
      allocate(           bk(num_lev+1))
      allocate(   sin_lat_xy(num_lon, num_lat))
      allocate(   cos_lat_xy(num_lon, num_lat))
      allocate(   cos_lat_yz(num_lat, num_lev))

      allocate(  pot_temp_coord(MaxIsentrLev))
      allocate( d_pot_temp_coord(MaxIsentrLev))


! ----------------- gridpoint variables on pressure surfaces ------------------


      allocate(            ts_grid(num_lon, num_lat))
      allocate(            ps_grid(num_lon, num_lat))
      allocate(             p_half(num_lon, num_lat, num_lev+1))
      allocate(             p_full(num_lon, num_lat, num_lev))

      allocate(     hor_sfctn_grid(num_lon, num_lat, num_lev))

      allocate(             u_grid(num_lon, num_lat, num_lev))
      allocate(             v_grid(num_lon, num_lat, num_lev))


      allocate(          vcos_grid(num_lon, num_lat, num_lev))
      allocate(         sfctn_grid(num_lon, num_lat, num_lev+1))
      allocate(           vor_grid(num_lon, num_lat, num_lev))
      allocate(           div_grid(num_lon, num_lat, num_lev))
      allocate(             w_grid(num_lon, num_lat, num_lev))
      allocate(          dpdt_grid(num_lon, num_lat, num_lev))
      allocate(          temp_grid(num_lon, num_lat, num_lev))
      allocate(  virtual_temp_grid(num_lon, num_lat, num_lev))
      allocate(     pot_temp_grid(num_lon, num_lat, num_lev))
      allocate(        geopot_grid(num_lon, num_lat, num_lev))

      allocate(    pot_temp_e_grid(num_lon, num_lat, num_lev))
      allocate(pot_temp_e_sat_grid(num_lon, num_lat, num_lev))
      allocate(          shum_grid(num_lon, num_lat, num_lev))
      allocate(    dt_shum_cond_grid(num_lon, num_lat, num_lev))
      allocate(    dt_shum_conv_grid(num_lon, num_lat, num_lev))
      allocate(    dt_shum_diff_grid(num_lon, num_lat, num_lev))
      allocate( dt_shum_cond_counter(num_lon, num_lat, num_lev))
      allocate( dt_shum_conv_counter(num_lon, num_lat, num_lev))

      allocate(        sat_mr_grid(num_lon, num_lat, num_lev))
      allocate(      sat_shum_grid(num_lon, num_lat, num_lev))
      allocate(          rhum_grid(num_lon, num_lat, num_lev))
      allocate(   precip_cond_grid(num_lon, num_lat))
      allocate(   precip_conv_grid(num_lon, num_lat))
      allocate(  precip_cond_daily(num_lon, num_lat))
      allocate(  precip_conv_daily(num_lon, num_lat))     

      allocate(        precip_mask(num_lon, num_lat))
      allocate(          one_array(num_lon, num_lat))

      allocate( bucket_depth_grid(num_lon, num_lat)) 
      allocate( bucket_depth_conv_grid(num_lon, num_lat)) 
      allocate( bucket_depth_cond_grid(num_lon, num_lat)) 
      allocate( bucket_depth_LH_grid(num_lon, num_lat)) 
      allocate( bucket_diffusion_grid(num_lon, num_lat)) 

      allocate( drag_coeff_mo_grid(num_lon, num_lat)) 
      allocate( drag_coeff_lh_grid(num_lon, num_lat)) 
      allocate( drag_coeff_sh_grid(num_lon, num_lat))

      allocate(   sfc_flux_lh_grid(num_lon, num_lat))
      allocate(   sfc_flux_sh_grid(num_lon, num_lat))
      allocate(   sfc_flux_sw_grid(num_lon, num_lat))
      allocate(  sfc_flux_lwd_grid(num_lon, num_lat))
      allocate(  sfc_flux_lwu_grid(num_lon, num_lat))
      allocate(sfc_flux_ocean_grid(num_lon, num_lat))
      allocate(   toa_flux_sw_grid(num_lon, num_lat))
      allocate(  toa_flux_lwu_grid(num_lon, num_lat))

      allocate(    dt_temp_cond_grid(num_lon, num_lat, num_lev))
      allocate(    dt_temp_conv_grid(num_lon, num_lat, num_lev))
      allocate(    dt_temp_diff_grid(num_lon, num_lat, num_lev))
      allocate(     dt_temp_rad_grid(num_lon, num_lat, num_lev))
      allocate(      dt_temp_sw_grid(num_lon, num_lat, num_lev))

      allocate(       dt_u_diff_grid(num_lon, num_lat, num_lev))  ! ljj addition
      allocate(       dt_v_diff_grid(num_lon, num_lat, num_lev))  ! ljj addition
      allocate(       dt_u_diss_grid(num_lon, num_lat, num_lev))  ! ljj addition
      allocate(       dt_v_diss_grid(num_lon, num_lat, num_lev))  ! ljj addition
      allocate(    dt_temp_diss_grid(num_lon, num_lat, num_lev))  ! ljj addition

      allocate(      bt_sfctn_grid(num_lon, num_lat))
      allocate(  conv_zon_spec_avg(num_lat, 0:num_fourier))


      allocate(  ucos_zon_spec_barotr(num_lat, 0:num_fourier))
      allocate(  vcos_zon_spec_barotr(num_lat, 0:num_fourier))
      allocate(  shum_zon_spec_barotr(num_lat, 0:num_fourier))

! --------- zonally-averaged gridpoint variables on pressure surfaces ---------
      allocate(                ts_avg(num_lat))
      allocate(                ps_avg(num_lat))
      allocate(             p_half_zon(num_lat, num_lev+1))
      allocate(             p_full_avg(num_lat, num_lev  ))
      allocate(       specific_vol_avg(num_lat, num_lev  ))
      allocate(               dens_avg(num_lat, num_lev  ))

      allocate(                  u_avg(num_lat, num_lev  ))
      allocate(              u_var_avg(num_lat, num_lev  ))


      allocate(       mrdnl_u_flux_avg(num_lat, num_lev  ))
      allocate(       vrtcl_u_flux_avg(num_lat, num_lev  ))
      allocate(  mrdnl_eddy_u_flux_avg(num_lat, num_lev  ))
      allocate(  vrtcl_eddy_u_flux_avg(num_lat, num_lev  ))
      allocate(             u_barotr_avg(num_lat))
      allocate(         u_barotr_var_avg(num_lat))

      allocate(                  v_avg(num_lat, num_lev  ))
      allocate(               vcos_avg(num_lat, num_lev  ))
      allocate(              v_var_avg(num_lat, num_lev  ))
      allocate(             v_barotr_avg(num_lat))
      allocate(         v_barotr_var_avg(num_lat))


      allocate(              sfctn_avg(num_lat, num_lev  ))
      allocate(           TEM_res_circ(num_lat, num_lev  ))
      allocate(          MTEM_res_circ(num_lat, num_lev  ))

      allocate(               vort_avg(num_lat, num_lev  ))

      allocate(                  w_avg(num_lat, num_lev  ))
      allocate(              w_var_avg(num_lat, num_lev  ))
      allocate(               dpdt_avg(num_lat, num_lev  ))
      allocate(           conv_pe_ke_avg(num_lat, num_lev  ))
      allocate(    conv_eddy_pe_ke_avg(num_lat, num_lev  ))

      allocate(               temp_avg(num_lat, num_lev  ))
      allocate(           temp_var_avg(num_lat, num_lev  ))
      allocate(     mrdnl_temp_flux_avg(num_lat, num_lev  ))
      allocate(     vrtcl_temp_flux_avg(num_lat, num_lev  ))
      allocate(mrdnl_eddy_temp_flux_avg(num_lat, num_lev  ))
      allocate(vrtcl_eddy_temp_flux_avg(num_lat, num_lev  ))
      allocate(       virtual_temp_avg(num_lat, num_lev  ))

      allocate(              pot_temp_avg(num_lat, num_lev  ))
      allocate(          pot_temp_var_avg(num_lat, num_lev  ))
      allocate(   mrdnl_pot_temp_flux_avg(num_lat, num_lev  ))
      allocate(   vrtcl_pot_temp_flux_avg(num_lat, num_lev  ))
      allocate(  mrdnl_eddy_pot_temp_flux_avg(num_lat, num_lev  ))
      allocate(  vrtcl_eddy_pot_temp_flux_avg(num_lat, num_lev  ))
      allocate(         d_dy_pot_temp_avg(num_lat, num_lev  ))
      allocate(       d_dsig_pot_temp_avg(num_lat, num_lev  ))
      allocate(      d_dp_pot_temp_avg(num_lat, num_lev  ))
      allocate(      buoyancy_freq_avg(num_lat, num_lev  ))

      allocate(                  z_avg(num_lat, num_lev  ))
      allocate(              z_var_avg(num_lat, num_lev  ))
      allocate(       mrdnl_z_flux_avg(num_lat, num_lev  ))
      allocate(       vrtcl_z_flux_avg(num_lat, num_lev  ))
      allocate(      mrdnl_eddy_z_flux_avg(num_lat, num_lev  ))
      allocate(      vrtcl_eddy_z_flux_avg(num_lat, num_lev  ))


      allocate(            pot_temp_e_avg(num_lat, num_lev  ))
      allocate(        pot_temp_e_var_avg(num_lat, num_lev  ))
      allocate( mrdnl_pot_temp_e_flux_avg(num_lat, num_lev  ))
      allocate( vrtcl_pot_temp_e_flux_avg(num_lat, num_lev  ))
      allocate(mrdnl_eddy_pot_temp_e_flux_avg(num_lat, num_lev  ))
      allocate(vrtcl_eddy_pot_temp_e_flux_avg(num_lat, num_lev  ))
      allocate(       d_dy_pot_temp_e_avg(num_lat, num_lev  ))
      allocate(     d_dsig_pot_temp_e_avg(num_lat, num_lev  ))
      allocate(       d_dp_pot_temp_e_avg(num_lat, num_lev  ))
      allocate(        pot_temp_e_sat_avg(num_lat, num_lev  ))
      allocate(         dt_shum_diff_avg(num_lat, num_lev  ))
      allocate(               shum_avg(num_lat, num_lev  ))
      allocate(           shum_var_avg(num_lat, num_lev  ))
      allocate(    mrdnl_shum_flux_avg(num_lat, num_lev  ))
      allocate(    vrtcl_shum_flux_avg(num_lat, num_lev  ))
      allocate( mrdnl_eddy_shum_flux_avg(num_lat, num_lev  ))
      allocate( vrtcl_eddy_shum_flux_avg(num_lat, num_lev  ))
      allocate(         dt_shum_cond_avg(num_lat, num_lev  ))
      allocate(         dt_shum_conv_avg(num_lat, num_lev  ))

      allocate(       shum_cond_prob_avg(num_lat, num_lev  ))
      allocate(       shum_conv_prob_avg(num_lat, num_lev  ))
      allocate(       sat_shum_var_avg(num_lat, num_lev  ))
      allocate(           sat_shum_avg(num_lat, num_lev  ))
      allocate(               rhum_avg(num_lat, num_lev  ))
      allocate(       precip_cond_avg(num_lat))
      allocate(       precip_conv_avg(num_lat))

      allocate(      precip_daily_above_threshold_avg(num_lat))
      allocate(                           days_total(num_lat))
      allocate(                 days_above_threshold(num_lat))
      allocate(precip_daily_above_threshold_prob_avg(num_lat))

      allocate(       sfc_flux_lh_avg(num_lat))
      allocate(       sfc_flux_sh_avg(num_lat))
      allocate(       sfc_flux_sw_avg(num_lat))
      allocate(      sfc_flux_lwd_avg(num_lat))
      allocate(      sfc_flux_lwu_avg(num_lat))
      allocate(    sfc_flux_ocean_avg(num_lat))
      allocate(       toa_flux_sw_avg(num_lat))
      allocate(      toa_flux_lwu_avg(num_lat))

      allocate(     pot_temp_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(   pot_temp_e_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(     temp_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(eape2eke_conv_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(      shum_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(  sat_shum_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(      rhum_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(              rhum_pdf(num_lat, num_lev, num_bin))
      allocate(    bc_sfctn_spec( num_fourier+1, num_fourier+1))
      allocate(    bt_sfctn_spec( num_fourier+1, num_fourier+1))

      allocate(    coriolis_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(     non_lin_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(non_lin_eddy_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(     t_non_lin_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(t_non_lin_eddy_spec( num_fourier+1, num_fourier+1, num_lev))
      allocate(     bt_non_lin_spec( num_fourier+1, num_fourier+1))
      allocate(bt_non_lin_eddy_spec( num_fourier+1, num_fourier+1))

      allocate(     pot_temp_spec_avg( num_fourier+1, num_fourier+1, num_lev))
      allocate(   pot_temp_e_spec_avg( num_fourier+1, num_fourier+1, num_lev))
      allocate(eape2eke_conv_spec_avg( num_fourier+1, num_fourier+1, num_lev))
      allocate(      shum_spec_avg( num_fourier+1, num_fourier+1, num_lev))
      allocate(  sat_shum_spec_avg( num_fourier+1, num_fourier+1, num_lev))
      allocate(      rhum_spec_avg( num_fourier+1, num_fourier+1, num_lev))

! -----------------------------------------------------------------------------
! ------------------- gridpoint variables on theta surfaces -------------------

      allocate(     mont_sfc(num_lon, num_lat))
      allocate(  dx_mont_sfc(num_lon, num_lat))
      allocate(  dy_mont_sfc(num_lon, num_lat))

      allocate(      dx_pot_temp_sfc(num_lon, num_lat))
      allocate(      dy_pot_temp_sfc(num_lon, num_lat))

      allocate(        geopot_sfc(num_lon, num_lat))

      allocate( mont_dx_pot_temp_sfc_avg(num_lat))
      allocate( mrdnl_geostr_eddy_pot_temp_flux_sfc_avg(num_lat))

      allocate(             mont_grid(num_lon, num_lat, num_lev))
      allocate(       mont_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(   mont_var_isent_grid(num_lon, num_lat, MaxIsentrLev))

      allocate(     dx_mont_grid(num_lon, num_lat, num_lev))
      allocate(     dx_mont_isent_grid(num_lon, num_lat, MaxIsentrLev))

      allocate( p_isent_grid_nw(num_lon, num_lat, MaxIsentrLev))
      allocate(  dens_isent_grid_nw(num_lon, num_lat, MaxIsentrLev))


      allocate(  u_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(  u_var_isent_grid(num_lon, num_lat, MaxIsentrLev))

      allocate(  v_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(  v_var_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(vcos_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(vcos_isent_grid_nw(num_lon, num_lat, MaxIsentrLev))
      allocate(v_geostr_isent_grid(num_lon, num_lat, MaxIsentrLev))

      allocate(      sfctn_isent_grid(num_lon, num_lat, MaxIsentrLev))

      allocate(           vor_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(       vor_var_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(mrdnl_vor_flux_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(    abs_vor_isent_grid_nw(num_lon, num_lat, MaxIsentrLev))
      allocate(       abs_vor_var_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(mrdnl_abs_vor_flux_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(           pot_vor_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(  v_pot_square_vor_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(    pot_square_vor_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(       pot_vor_var_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(mrdnl_pot_vor_flux_isent_grid(num_lon, num_lat,  MaxIsentrLev))

      allocate(       shum_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate( mrdnl_shum_flux_isent_grid(num_lon, num_lat, MaxIsentrLev))

      allocate(       rhum_isent_grid(num_lon, num_lat, MaxIsentrLev))

      allocate(           sat_shum_isent_grid(num_lon, num_lat, MaxIsentrLev))
      allocate(mrdnl_sat_shum_flux_isent_grid(num_lon, num_lat, MaxIsentrLev))

      allocate(             k_sfc(num_lon, num_lat              ))
      allocate(            k_asfc(num_lon, num_lat              ))

      allocate( bucket_depth_avg(num_lat)) 
      allocate( bucket_depth_conv_avg(num_lat)) 
      allocate( bucket_depth_cond_avg(num_lat)) 
      allocate( bucket_depth_LH_avg(num_lat)) 
      allocate( bucket_diffusion_avg(num_lat)) 


      allocate(     drag_coeff_mo_avg(num_lat))
      allocate(     drag_coeff_lh_avg(num_lat))
      allocate(     drag_coeff_sh_avg(num_lat))

      allocate(    dt_temp_cond_avg(num_lat, num_lev  ))
      allocate(    dt_temp_conv_avg(num_lat, num_lev  ))
      allocate(    dt_temp_diff_avg(num_lat, num_lev  ))
      allocate(     dt_temp_rad_avg(num_lat, num_lev  ))
      allocate(      dt_temp_sw_avg(num_lat, num_lev  ))

      allocate(       dt_u_diff_avg(num_lat, num_lev  ))   ! ljj addition
      allocate(       dt_v_diff_avg(num_lat, num_lev  ))   ! ljj addition
      allocate(       dt_u_diss_avg(num_lat, num_lev  ))   ! ljj addition
      allocate(       dt_v_diss_avg(num_lat, num_lev  ))   ! ljj addition
      allocate(    dt_temp_diss_avg(num_lat, num_lev  ))   ! ljj addition

      allocate(      ke_diss_zonal_avg(num_lat, num_lev))  ! ljj addition
      allocate(      ke_diss_mrdnl_avg(num_lat, num_lev))  ! ljj addition
      allocate(      ke_diff_zonal_avg(num_lat, num_lev))  ! ljj addition
      allocate(      ke_diff_mrdnl_avg(num_lat, num_lev))  ! ljj addition

! ---------- zonally averaged gridpoint variables on theta surfaces -----------

      allocate(           sfctn_isent_avg(num_lat, MaxIsentrLev))

      allocate(     pot_vor_isent_avg(num_lat, MaxIsentrLev))
      allocate(     pot_vor_var_isent_avg(num_lat, MaxIsentrLev))
      allocate( mrdnl_pot_vor_flux_isent_avg(num_lat, MaxIsentrLev))
      allocate( mrdnl_eddy_pot_vor_flux_isent_avg(num_lat, MaxIsentrLev))

      allocate(               p_isent_avg_nw(num_lat, MaxIsentrLev))
      allocate(           p_var_isent_avg_nw(num_lat, MaxIsentrLev))

      allocate(            dens_isent_avg_nw(num_lat, MaxIsentrLev))
      allocate(           exner_isent_avg_nw(num_lat, MaxIsentrLev))
      allocate(          p_exner_isent_avg_nw(num_lat, MaxIsentrLev))

      allocate(               u_isent_avg(num_lat, MaxIsentrLev))
      allocate(              u_var_isent_avg(num_lat, MaxIsentrLev))

      allocate(               v_isent_avg(num_lat, MaxIsentrLev))
      allocate(            vcos_isent_avg(num_lat, MaxIsentrLev))
      allocate(              v_var_isent_avg(num_lat, MaxIsentrLev))
      allocate(           v_geostr_isent_avg(num_lat, MaxIsentrLev))

      allocate(         sfc_pot_temp_pdf(num_lat, MaxIsentrLev))

      allocate(            mont_isent_avg(num_lat, MaxIsentrLev))
      allocate(        mont_var_isent_avg(num_lat, MaxIsentrLev))
      allocate(         dx_mont_isent_avg(num_lat, MaxIsentrLev))

      allocate(            shum_isent_avg(num_lat, MaxIsentrLev))
      allocate(      mrdnl_shum_flux_isent_avg(num_lat, MaxIsentrLev))
      allocate(  mrdnl_eddy_shum_flux_isent_avg(num_lat, MaxIsentrLev))
      allocate(            rhum_isent_avg(num_lat, MaxIsentrLev))

      allocate(        sat_shum_isent_avg(num_lat, MaxIsentrLev))
      allocate(mrdnl_sat_shum_flux_isent_avg(num_lat, MaxIsentrLev))
      allocate(mrdnl_eddy_sat_shum_flux_isent_avg(num_lat, MaxIsentrLev))

      allocate(     pot_enstrophy_isent_avg(num_lat, MaxIsentrLev))
      allocate(     v_pot_enstrophy_isent_avg(num_lat, MaxIsentrLev))
      allocate(     mrdnl_eddy_pot_enstrophy_flux_isent_avg(num_lat, MaxIsentrLev))

! -------------------------initialize variables-------------------------------


      ts_avg                  =  0.0
      ps_avg                  =  0.0
      p_full_avg              =  0.0
      specific_vol_avg        =  0.0
      dens_avg                =  0.0

      u_avg                   =  0.0
      u_var_avg               =  0.0
      mrdnl_u_flux_avg        =  0.0
      vrtcl_u_flux_avg        =  0.0
      mrdnl_eddy_u_flux_avg   =  0.0
      vrtcl_eddy_u_flux_avg   =  0.0
      u_barotr_avg              =  0.0
      u_barotr_var_avg          =  0.0

      v_avg                   =  0.0
      vcos_avg                =  0.0
      v_var_avg               =  0.0
      v_barotr_avg              =  0.0
      v_barotr_var_avg          =  0.0



      sfctn_avg               =  0.0
      TEM_res_circ            =  0.0
      MTEM_res_circ           =  0.0

      vort_avg                =  0.0

      w_avg                   =  0.0
      w_var_avg               =  0.0
      dpdt_avg            =  0.0
      conv_pe_ke_avg            =  0.0
      conv_eddy_pe_ke_avg          =  0.0

      temp_avg                =  0.0
      temp_var_avg            =  0.0
      mrdnl_temp_flux_avg     =  0.0
      vrtcl_temp_flux_avg     =  0.0
      mrdnl_eddy_temp_flux_avg     =  0.0
      vrtcl_eddy_temp_flux_avg     =  0.0
      virtual_temp_avg        =  0.0

      pot_temp_avg               =  0.0
      pot_temp_var_avg           =  0.0
      mrdnl_pot_temp_flux_avg    =  0.0
      vrtcl_pot_temp_flux_avg    =  0.0
      mrdnl_eddy_pot_temp_flux_avg    =  0.0
      vrtcl_eddy_pot_temp_flux_avg    =  0.0
      d_dy_pot_temp_avg          =  0.0
      d_dsig_pot_temp_avg        =  0.0
      d_dp_pot_temp_avg          =  0.0
      buoyancy_freq_avg          =  0.0

      z_avg                   =  0.0
      z_var_avg               =  0.0

      mrdnl_z_flux_avg        =  0.0
      vrtcl_z_flux_avg        =  0.0
      mrdnl_eddy_z_flux_avg   =  0.0
      vrtcl_eddy_z_flux_avg   =  0.0

      pot_temp_e_avg             =  0.0
      pot_temp_e_var_avg         =  0.0
      d_dy_pot_temp_e_avg        =  0.0
      d_dsig_pot_temp_e_avg      =  0.0
      d_dp_pot_temp_e_avg        =  0.0
      pot_temp_e_sat_avg         =  0.0
      mrdnl_pot_temp_e_flux_avg  =  0.0
      vrtcl_pot_temp_e_flux_avg  =  0.0
      mrdnl_eddy_pot_temp_e_flux_avg  =  0.0
      vrtcl_eddy_pot_temp_e_flux_avg  =  0.0

      shum_avg                =  0.0
      shum_var_avg            =  0.0
      mrdnl_shum_flux_avg     =  0.0
      vrtcl_shum_flux_avg     =  0.0
      mrdnl_eddy_shum_flux_avg     =  0.0
      vrtcl_eddy_shum_flux_avg     =  0.0
      dt_shum_cond_avg        =  0.0
      dt_shum_conv_avg        =  0.0
      dt_shum_diff_avg        =  0.0
      shum_cond_prob_avg   =  0.0
      shum_conv_prob_avg   =  0.0
      sat_shum_var_avg        =  0.0
      sat_shum_avg            =  0.0
      rhum_avg                =  0.0
      precip_cond_avg         =  0.0
      precip_conv_avg         =  0.0

      bucket_depth_avg = 0.0
      bucket_depth_conv_avg = 0.0
      bucket_depth_cond_avg = 0.0
      bucket_depth_LH_avg = 0.0
      bucket_diffusion_avg  = 0.0

      drag_coeff_mo_avg       = 0.0
      drag_coeff_lh_avg       = 0.0
      drag_coeff_sh_avg       = 0.0

      dt_temp_cond_avg        = 0.0
      dt_temp_conv_avg        = 0.0
      dt_temp_diff_avg        = 0.0
      dt_temp_rad_avg         = 0.0 
      dt_temp_sw_avg          = 0.0 
      
      dt_u_diff_avg           = 0.0   ! ljj addition
      dt_v_diff_avg           = 0.0   ! ljj addition
      dt_u_diss_avg           = 0.0   ! ljj addition
      dt_v_diss_avg           = 0.0   ! ljj addition
      dt_temp_diss_avg        = 0.0   ! ljj addition

      ke_diss_zonal_avg       = 0.0   ! ljj addition
      ke_diss_mrdnl_avg       = 0.0   ! ljj addition
      ke_diff_zonal_avg       = 0.0   ! ljj addition
      ke_diff_mrdnl_avg       = 0.0   ! ljj addition

      precip_daily_above_threshold_avg      =  0.0
      days_total                            =  0.0
      days_above_threshold                  =  0.0 
      precip_daily_above_threshold_prob_avg =  0.0        
      precip_cond_daily                     =  0.0
      precip_conv_daily                     =  0.0
      one_array                             =  1.0

      sfc_flux_lh_avg         =  0.0 
      sfc_flux_sh_avg         =  0.0 
      sfc_flux_sw_avg         =  0.0
      sfc_flux_lwd_avg        =  0.0
      sfc_flux_lwu_avg        =  0.0
      sfc_flux_ocean_avg      =  0.0
      toa_flux_sw_avg         =  0.0
      toa_flux_lwu_avg        =  0.0

      ucos_zon_spec_barotr       =  0.0
      vcos_zon_spec_barotr       =  0.0
      shum_zon_spec_barotr       =  0.0
      pot_temp_spec              =  0.0
      pot_temp_e_spec            =  0.0
      eape2eke_conv_spec         =  0.0
      shum_spec                  =  0.0
      sat_shum_spec              =  0.0
      rhum_spec                  =  0.0
      pot_temp_spec_avg              =  0.0
      pot_temp_e_spec_avg            =  0.0
      eape2eke_conv_spec_avg         =  0.0
      shum_spec_avg                  =  0.0
      sat_shum_spec_avg              =  0.0
      rhum_spec_avg                  =  0.0
      rhum_pdf                   =  0.0 
      ! pdf bins for relative humidity
      ! include one negative bin
      do i=1, num_bin
         rhum_bin(i) =  (float(i) - 1.5) / (float(num_bin)-1)
      enddo

      mont_sfc     = 0.0
      dx_mont_sfc  = 0.0 
      dy_mont_sfc  = 0.0

      dx_pot_temp_sfc  = 0.0 
      dy_pot_temp_sfc  = 0.0

      geopot_sfc              =  0.0
      geopot_grid             =  0.0
      k_sfc                   =  0
      k_asfc                  =  0
      sfc_pot_temp_pdf           =  0.0      
      shum_isent_avg              =  0.0
      rhum_isent_avg              =  0.0
      sat_shum_isent_avg          =  0.0
      mrdnl_shum_flux_isent_avg             =  0.0
      mrdnl_sat_shum_flux_isent_avg         =  0.0
      mrdnl_eddy_shum_flux_isent_avg             =  0.0
      mrdnl_eddy_sat_shum_flux_isent_avg         =  0.0
      sfctn_isent_avg                 =  0.0
      p_isent_avg_nw                 =  0.0
      p_var_isent_avg_nw             =  0.0
      exner_isent_avg_nw             =  0.0
      p_exner_isent_avg_nw           =  0.0
      dens_isent_avg_nw              =  0.0
      v_isent_avg                 =  0.0
      vcos_isent_avg              =  0.0
      v_geostr_isent_avg          =  0.0
      u_isent_avg                 =  0.0
      v_var_isent_avg             =  0.0
      u_var_isent_avg             =  0.0
      p_isent_grid_nw                =  0.0
      dens_isent_grid_nw             =  0.0
      v_isent_grid                =  0.0
      vcos_isent_grid             =  0.0
      vcos_isent_grid_nw          =  0.0
      u_isent_grid                =  0.0
      vor_isent_grid              =  0.0
      vor_var_isent_grid            =  0.0
      mrdnl_vor_flux_isent_grid     =  0.0
      abs_vor_isent_grid_nw         =  0.0
      abs_vor_var_isent_grid        =  0.0
      mrdnl_abs_vor_flux_isent_grid =  0.0
      pot_vor_isent_grid            =  0.0
      pot_vor_var_isent_grid        =  0.0
      v_pot_square_vor_isent_grid   = 0.0
      pot_square_vor_isent_grid   = 0.0
      mrdnl_pot_vor_flux_isent_grid     =  0.0
      mrdnl_pot_vor_flux_isent_avg      =  0.0
      mrdnl_eddy_pot_vor_flux_isent_avg =  0.0

      pot_enstrophy_isent_avg = 0.0
      v_pot_enstrophy_isent_avg = 0.0 
      mrdnl_eddy_pot_enstrophy_flux_isent_avg = 0.0 

      sfctn_isent_grid               =  0.0
      mont_isent_grid                =  0.0
      mont_var_isent_grid            =  0.0
      dx_mont_isent_avg            =  0.0
      mont_isent_avg              =  0.0
      mont_var_isent_avg          =  0.0
      mont_dx_pot_temp_sfc_avg                = 0.0
      mrdnl_geostr_eddy_pot_temp_flux_sfc_avg = 0.0

      bc_sfctn_spec      =  0.0
      bt_sfctn_spec      =  0.0

      temp_spec = 0.0
      coriolis_spec      =  0.0
      non_lin_spec       =  0.0
      non_lin_eddy_spec  =  0.0
      t_non_lin_spec       =  0.0
      t_non_lin_eddy_spec  =  0.0
      bt_non_lin_spec       =  0.0
      bt_non_lin_eddy_spec  =  0.0
      conv_zon_spec_avg  =  0.0

     end subroutine allocate_and_initialize

   end module vars
