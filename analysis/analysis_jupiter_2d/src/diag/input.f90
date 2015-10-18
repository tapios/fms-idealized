module input_mod
  
  use constants_and_switches_mod
  use netcdf

  implicit none
  private
  public ::         open_input_files,                     read_data,          &
              read_surf_geopotential,                get_dimensions,          &
                           axes_init,        pressure_variables_FMS,          &
             pressure_variables_NCEP,       pressure_variables_CCM3,          &
            pressure_variables_ERA40,   pressure_variables_FMS_init,          &
        pressure_variables_NCEP_init,  pressure_variables_CCM3_init,          &
       pressure_variables_ERA40_init,                         check

! ---------------------------input file variables--------------------------

  integer :: UncInputID,     &   ! zonal wind input file id
       VncInputID,           &   ! meridional wind input file id
       VORncInputID,         &   ! vorticity input file id
       DIVncInputID,         &   ! divergence input file id
       TEMPncInputID,        &   ! temperature input file id
       TSncInputID,          &   ! surface temperature file id
       PSncInputID,          &   ! surface pressure input file id
       SHumncInputID,        &   ! specific humidity input file id
       CondncInputID,        &   ! condensation moisture tendency input file id
       ConvncInputID,        &   ! convection moisture tendency input file id
       DiffncInputID,        &   
       PrecipCondncInputID,  &   ! condensation moisture tendency input file id 
       PrecipConvncInputID,  &   ! convection moisture tendency input file id
       SfcFluxLHncInputID,   &   
       SfcFluxSHncInputID,   &
       SfcFluxSWncInputID,   &
       SfcFluxLWDncInputID,  &
       SfcFluxLWUncInputID,  &
       SfcFluxQFncInputID,   &           
       ToaFluxSWncInputID,   &
       ToaFluxLWUncInputID,  &
       BucketDepthncInputID,    &
       BucketDepthConvncInputID,&
       BucketDepthCondncInputID,&
       BucketDepthLHncInputID,&
       BucketDiffusionncInputID,&
       DiabCondncInputID,    &
       DiabConvncInputID,    &
       DiabDiffncInputID,    &
       DiabRadncInputID,     &
       DiabSWncInputID,      &
       DragMOncInputID,      &
       DragLHncInputID,      &
       DragSHncInputID,      &
       UVarID,               &   ! zonal wind ID
       VVarID,               &   ! meridional wind ID
       VORVarID,             &   ! vorticity ID
       DIVVarID,             &   ! divergence ID
       TEMPVarID,            &   ! temperature ID
       TSVarID,              &   ! surface temperature ID
       PSVarID,              &   ! surface pressure ID
       SHumVarID,            &   ! specific humidity ID
       CondVarID,            &   ! large scale cond specific humidity tendency ID
       ConvVarID,            &   ! convection specific humidity tendency ID
       DiffVarID,            &
       PrecipCondVarID,      &   ! large scale cond precip ID
       PrecipConvVarID,      &   ! convection precip ID
       SfcFluxLHVarID,       &   
       SfcFluxSHVarID,       &
       SfcFluxSWVarID,       &
       SfcFluxLWDVarID,      & 
       SfcFluxLWUVarID,      &
       SfcFluxQFVarID,       &           
       ToaFluxSWVarID,       &
       ToaFluxLWUVarID,      &
       DiabCondVarID,        &
       DiabConvVarID,        &
       DiabDiffVarID,        &
       DiabRadVarID,         &
       DiabSWVarID,          &
       BucketDepthVarID,     &
       BucketDepthConvVarID, &
       BucketDepthCondVarID, &
       BucketDepthLHVarID,   &
       BucketDiffusionVarID, &
       DragMOVarID,          &
       DragLHVarID,          &
       DragSHVarID,          &
       ierror,               &   ! error flag for GRIB reads
       GribUnit                  ! GRIB file identifier (for ERA40 data)

  integer ::                                                                  &
       DataIn,             &   ! U and V or Vor and Div
       Data_Source,        &   ! where the data come from
       GRIBID                  ! index of GRIB ID 

  integer, dimension(nf90_max_var_dims) ::                                    &
       dimIDs                  ! dimensions of time 

  logical ::                                                                  &
       isentrope,          &
       moisture,           &    ! process moisture or not
       bucket,             &    ! bucket hydrology
       virtual                  ! use virtual temp for z and svol or not 

  character(len=128) ::                                                       &
       GRIBFile                ! GRIB filename

  ! variables for pressure routines:

  ! CCM3 and ERA40:
  real, allocatable, dimension(:) ::                                          &
       hyam,               &   ! hybrid a cofficients on layer midpoints
       hybm                    ! hybrid b cofficients on layer midpoints

  real, allocatable, dimension(:) ::                                          &
       hyai,               &   ! hybrid a coefficients on layer interfaces    &
       hybi                    ! hybrid b coefficients on layer interfaces

  real ::                                                                     &
       p0                      ! reference pressure (CCM3)

  ! FMS:
  real, allocatable, dimension(:) ::                                          &
       pk,                 &   ! hybrid pressure coefficient
       bk                      ! hybrid sigma coefficient

  real, dimension(160) ::                                                     &
       era40_lats              ! latitudes of T159 Gaussian grid used in ERA40

contains
! #############################################################################

  subroutine open_input_files(                                                &
                        DataIn_local,               Data_Source_local,        &
                      local_moisture,                    local_bucket,        & 
                       local_virtual,                 local_isentrope,        &
                            UVarName,                        VVarName,        &
                          VORVarName,                      DIVVarName,        &
                         TEMPVarName,                       TSVarName,        &
                           PSVarName,                                         &
                          SHumVarName,                                         &
                         CondVarName,                     ConvVarName,        &
                         DiffVarName,                                         &
                   PrecipCondVarName,               PrecipConvVarName,        &
                    SfcFluxLHVarName,                SfcFluxSHVarName,        &
                    SfcFluxSWVarName,               SfcFluxLWDVarName,        &
                   SfcFluxLWUVarName,                SfcFluxQFVarName,        &          
                    ToaFluxSWVarName,               ToaFluxLWUVarName,        & 
                     DiabCondVarName,                 DiabConvVarName,        &
                     DiabDiffVarName,                  DiabRadVarName,        &
                       DiabSWVarName,                                         &
                  BucketDepthVarName,          BucketDepthConvVarName,        &
              BucketDepthCondVarName,            BucketDepthLHVarName,        &            
              BucketDiffusionVarName,                                         &
                       DragMOVarName,                   DragLHVarName,        &
                       DragSHVarName,                                         &
                       InputFileName)

    integer, intent(in) ::                                                    &
         DataIn_local,     &   ! What data are input (u and v or vor and div) 
         data_source_local     ! What source the data are from

    logical, intent(in) ::                                                    &
         local_isentrope,       &
         local_moisture,        &   ! flag to analyze moisture
         local_bucket,          &   ! flag to analyze bucket (hydrology)
         local_virtual              ! flag to include virtual temp effect for z



    character(len=*), intent(in)  ::                                          &
                            UVarName,                        VVarName,        &
                          VORVarName,                      DIVVarName,        &
                         TEMPVarName,                       TSVarName,        & 
                           PSVarName,                                         &
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
                       DragMOVarName,                   DragLHVarName,        &
                       DragSHVarName,                                         &
                  BucketDepthVarName,          BucketDepthConvVarName,        &
              BucketDepthCondVarName,            BucketDepthLHVarName,        &
              BucketDiffusionVarName,                                         &
                         SHumVarName,                   InputFileName

!                           --- executable code ---
    DataIn = dataIn_local
    Data_Source = Data_Source_local

    if(data_source .eq. era40) then
       GRIBID = 2
       GRIBFile = InputFileName
    else



!!!!XL       
       call check(nf90_open(trim(InputFileName), nf90_nowrite, TEMPncInputID ) )
       call check(nf90_open(trim(InputFileName), nf90_nowrite,   PSncInputID ) )
       call check(nf90_inq_varid(TEMPncInputID, trim(TEMPVarName),   TEMPVarID ) )
       call check(nf90_inq_varid(  PSncInputID,   trim(PSVarName), PSVarID ) )

       call check(nf90_open(InputFileName, nf90_nowrite, DiabConvncInputID ) )
       call check(nf90_inq_varid(DiabConvncInputID, trim(DiabConvVarName), DiabConvVarID) )
       call check(nf90_open(InputFileName, nf90_nowrite, DiabDiffncInputID ) )
       call check(nf90_inq_varid(DiabDiffncInputID, trim(DiabDiffVarName), DiabDiffVarID) )
       call check(nf90_open(InputFileName, nf90_nowrite, DiabRadncInputID ) )
       call check(nf90_inq_varid(DiabRadncInputID, trim(DiabRadVarName), DiabRadVarID) )
!!!!XL


       if(DataIn .eq. vor_and_div) then
          call check(nf90_open(trim(InputFileName), nf90_nowrite, VORncInputID ))
          call check(nf90_open(trim(InputFileName), nf90_nowrite, DIVncInputID ))
          call check(nf90_inq_varid(VORncInputID,  trim(VORVarName),  VORVarID ) )
          call check(nf90_inq_varid(DIVncInputID,  trim(DIVVarName),  DIVVarID ) )

       elseif(DataIn .eq. u_and_v)  then
          call check(nf90_open(InputFileName, nf90_nowrite, UncInputID ) )
          call check(nf90_open(InputFileName, nf90_nowrite, VncInputID ) )
          call check(nf90_inq_varid(UncInputID, trim(UVarName), UVarID        ) )
          call check(nf90_inq_varid(VncInputID, trim(VVarName), VVarID        ) )
       endif
       
       if(local_moisture) then
          call check(nf90_open(trim(InputFileName), nf90_nowrite, TSncInputID ) )
          call check(nf90_inq_varid(  TSncInputID,   trim(TSVarName), TSVarID ) )
          call check(nf90_open(InputFileName, nf90_nowrite, SHUMncInputID ) )
          call check(nf90_inq_varid(SHUMncInputID, trim(SHUMVarName), SHUMVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, CondncInputID ) )
          call check(nf90_inq_varid(CondncInputID, trim(CondVarName), CondVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, ConvncInputID ) )
          call check(nf90_inq_varid(ConvncInputID, trim(ConvVarName), ConvVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, DiffncInputID ) )
          call check(nf90_inq_varid(DiffncInputID, trim(DiffVarName), DiffVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, PrecipCondncInputID ) )
          call check(nf90_inq_varid(PrecipCondncInputID, trim(PrecipCondVarName), PrecipCondVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, PrecipConvncInputID ) )
          call check(nf90_inq_varid(PrecipConvncInputID, trim(PrecipConvVarName), PrecipConvVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, SfcFluxLHncInputID ) )
          call check(nf90_inq_varid(SfcFluxLHncInputID, trim(SfcFluxLHVarName), SfcFluxLHVarID) )        
          call check(nf90_open(InputFileName, nf90_nowrite, SfcFluxSHncInputID ) )
          call check(nf90_inq_varid(SfcFluxSHncInputID, trim(SfcFluxSHVarName), SfcFluxSHVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, SfcFluxSWncInputID ) )
          call check(nf90_inq_varid(SfcFluxSWncInputID, trim(SfcFluxSWVarName), SfcFluxSWVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, SfcFluxLWDncInputID ) )
          call check(nf90_inq_varid(SfcFluxLWDncInputID, trim(SfcFluxLWDVarName), SfcFluxLWDVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, SfcFluxLWUncInputID ) )
          call check(nf90_inq_varid(SfcFluxLWUncInputID, trim(SfcFluxLWUVarName), SfcFluxLWUVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, SfcFluxQFncInputID ) )
          call check(nf90_inq_varid(SfcFluxQFncInputID, trim(SfcFluxQFVarName), SfcFluxQFVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, ToaFluxSWncInputID ) )
          call check(nf90_inq_varid(ToaFluxSWncInputID, trim(ToaFluxSWVarName), ToaFluxSWVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, ToaFluxLWUncInputID ) )
          call check(nf90_inq_varid(ToaFluxLWUncInputID, trim(ToaFluxLWUVarName), ToaFluxLWUVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, DiabCondncInputID ) )
          call check(nf90_inq_varid(DiabCondncInputID, trim(DiabCondVarName), DiabCondVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, DiabSWncInputID ) )
          call check(nf90_inq_varid(DiabSWncInputID, trim(DiabSWVarName), DiabSWVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, DragMOncInputID ) )
          call check(nf90_inq_varid(DragMOncInputID, trim(DragMOVarName), DragMOVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, DragLHncInputID ) )
          call check(nf90_inq_varid(DragLHncInputID, trim(DragLHVarName), DragLHVarID) )
          call check(nf90_open(InputFileName, nf90_nowrite, DragSHncInputID ) )
          call check(nf90_inq_varid(DragSHncInputID, trim(DragSHVarName), DragSHVarID) )
       endif

       ! addition fridoo sept 2012 bucket hydrology
       if(local_bucket) then
          call check(nf90_open(trim(InputFileName), nf90_nowrite,     BucketDepthncInputID  ))
          call check(nf90_inq_varid(BucketDepthncInputID,  trim(BucketDepthVarName),  BucketDepthVarID ) )
          call check(nf90_open(trim(InputFileName), nf90_nowrite,     BucketDepthConvncInputID  ))
          call check(nf90_inq_varid(BucketDepthConvncInputID,  trim(BucketDepthConvVarName),  BucketDepthConvVarID ) )
          call check(nf90_open(trim(InputFileName), nf90_nowrite,     BucketDepthCondncInputID  ))
          call check(nf90_inq_varid(BucketDepthCondncInputID,  trim(BucketDepthCondVarName),  BucketDepthCondVarID ) )
          call check(nf90_open(trim(InputFileName), nf90_nowrite,     BucketDepthLHncInputID  ))
          call check(nf90_inq_varid(BucketDepthLHncInputID,  trim(BucketDepthLHVarName),  BucketDepthLHVarID ) )
          call check(nf90_open(trim(InputFileName), nf90_nowrite,     BucketDiffusionncInputID  ))
          call check(nf90_inq_varid(BucketDiffusionncInputID,  trim(BucketDiffusionVarName),  BucketDiffusionVarID ) )
       endif
       ! end fridoo
    endif

    isentrope   = local_isentrope
    moisture    = local_moisture
    virtual     = local_virtual 
    bucket      = local_bucket


  end subroutine open_input_files


! ##############################################################################
  subroutine read_data(       nsteps,                          u_grid,        &
                              v_grid,                        vor_grid,        &
                            div_grid,                       temp_grid,        &
                             ts_grid,                                         &
                             ps_grid,                       shum_grid,        &
                   dt_shum_cond_grid,               dt_shum_conv_grid,        &
                   dt_shum_diff_grid,                                         &
                    precip_cond_grid,                precip_conv_grid,        &
                    sfc_flux_lh_grid,                sfc_flux_sh_grid,        &
                    sfc_flux_sw_grid,               sfc_flux_lwd_grid,        &
                   sfc_flux_lwu_grid,             sfc_flux_ocean_grid,        &
                    toa_flux_sw_grid,               toa_flux_lwu_grid,        &
                   dt_temp_cond_grid,               dt_temp_conv_grid,        &
                   dt_temp_diff_grid,                dt_temp_rad_grid,        &
                     dt_temp_sw_grid,                                         &
                   bucket_depth_grid,          bucket_depth_conv_grid,        &
              bucket_depth_cond_grid,            bucket_depth_lh_grid,        &
                                                bucket_diffusion_grid,        &
                  drag_coeff_mo_grid,              drag_coeff_lh_grid,        &
                  drag_coeff_sh_grid)                                      

! ------------------------------ inout variables ------------------------------

    real, intent(inout) ::                                                    &
         nsteps               ! current time

! ------------------------------ output variables -----------------------------

    real, intent(out), dimension(:,:,:) ::                                    &
                            vor_grid,                        div_grid,        &
                           temp_grid,                          u_grid,        &
                              v_grid,                       shum_grid,        &
                   dt_shum_cond_grid,               dt_shum_conv_grid,        &
                   dt_shum_diff_grid,               dt_temp_cond_grid,        &
                   dt_temp_conv_grid,               dt_temp_diff_grid,        &
                    dt_temp_rad_grid,                 dt_temp_sw_grid


    real, intent(out), dimension(:,:)  ::                                     &
                             ts_grid,                                         &
                             ps_grid,                precip_cond_grid,        &
                    precip_conv_grid,                                         &
                    sfc_flux_lh_grid,                sfc_flux_sh_grid,        &
                    sfc_flux_sw_grid,               sfc_flux_lwd_grid,        &
                   sfc_flux_lwu_grid,             sfc_flux_ocean_grid,        &
                    toa_flux_sw_grid,               toa_flux_lwu_grid,        &
                   bucket_depth_grid,          bucket_depth_conv_grid,        &
              bucket_depth_cond_grid,            bucket_depth_lh_grid,        &
                                                bucket_diffusion_grid,        &
                  drag_coeff_mo_grid,              drag_coeff_lh_grid,        &
                  drag_coeff_sh_grid                     

! ------------------------------ local variables ------------------------------

    integer ::                                                                &
         local_num_lat,    &  ! number of latitudes in current window
         local_num_lon,    &  ! number of longitudes in current window
         num_lev,          &  ! number of vertical levels 
         k                    ! level index

    integer, dimension(5) ::                                                  &
         GRIBOff               ! GRIB offsets

! ------------------------------ executable code ------------------------------


    local_num_lat = size(temp_grid, 2)
    local_num_lon = size(temp_grid, 1)
    num_lev = size(temp_grid, 3)



#ifdef GRIB 
    if(data_source .eq. era40) then

       ! Offsets for reading from ERA40 GRIB files:
       GRIBOff(1) = 3; GRIBOff(2) = 2; GRIBOff(3) = 5

       call read_grib(       GRIBID+4,                         ps_grid )
       ps_grid = exp(ps_grid)

       do k = 1, num_lev

          call read_grib(      GRIBID,                temp_grid(:,:,k) )
          if(moisture)                                                        &
              call read_grib(GRIBID+1,                shum_grid(:,:,k) )

          GRIBID = GRIBID + GRIBOff(1)

          call read_grib(      GRIBID,                 vor_grid(:,:,k) )
          GRIBID = GRIBID + GRIBOff(2)
          
          call read_grib(      GRIBID,                 div_grid(:,:,k) )

          if(k .eq. 1)  GRIBOff(2) = 1

       enddo

       GRIBID = GRIBID + 1

    else 
#endif




     if(DataIn.eq.u_and_v) then
           
       call check(nf90_get_var(UncInputID, UVarID, u_grid,                    &
            start = (/ 1, 1, 1, int(nsteps)/),                                &
            count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))
       call check(nf90_get_var(VncInputID, VVarID, v_grid,                    &
            start = (/ 1, 1, 1, int(nsteps)/),                                &
            count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))
     endif

     if(DataIn.eq.vor_and_div) then
        call check(nf90_get_var(VORncInputID, VORVarID, vor_grid,             &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))
        call check(nf90_get_var(DIVncInputID, DIVVarID, div_grid,             &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))
     endif

     if(DataIn.eq.all) then
        call check(nf90_get_var(VORncInputID, VORVarID, vor_grid,             &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))
        call check(nf90_get_var(DIVncInputID, DIVVarID, div_grid,             &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))
        call check(nf90_get_var(UncInputID, UVarID, u_grid,                   &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))
        call check(nf90_get_var(VncInputID, VVarID, v_grid,                   &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))
     endif

     call check(nf90_get_var(TEMPncInputID, TEMPVarID, temp_grid,             &
          start = (/ 1, 1, 1, int(nsteps)/),                                  &
          count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

     call check(nf90_get_var(PSncInputID, PSVarID, ps_grid,                   &
          start = (/ 1, 1, int(nsteps)/),                                     &
          count = (/ local_num_lon, local_num_lat, 1/) ))

     call check(nf90_get_var(DiabConvncInputID, DiabConvVarID, dt_temp_conv_grid,  &
          start = (/ 1, 1, 1, int(nsteps)/),                               &
          count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

     call check(nf90_get_var(DiabDiffncInputID, DiabDiffVarID, dt_temp_diff_grid,  &
          start = (/ 1, 1, 1, int(nsteps)/),                               &
          count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

     call check(nf90_get_var(DiabRadncInputID, DiabRadVarID, dt_temp_rad_grid,  &
          start = (/ 1, 1, 1, int(nsteps)/),                               &
          count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

     if(moisture) then

        call check(nf90_get_var(TSncInputID, TSVarID, ts_grid,                   &
             start = (/ 1, 1, int(nsteps)/),                                     &

             count = (/ local_num_lon, local_num_lat, 1/) ))
        call check(nf90_get_var(SHUMncInputID, SHUMVarID, shum_grid,          &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

        call check(nf90_get_var(CondncInputID, CondVarID, dt_shum_cond_grid,    &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

        call check(nf90_get_var(ConvncInputID, ConvVarID, dt_shum_conv_grid,    &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

        call check(nf90_get_var(DiffncInputID, DiffVarID, dt_shum_diff_grid,  &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

        call check(nf90_get_var(PrecipCondncInputID, PrecipCondVarID, precip_cond_grid,    &
             start = (/ 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, 1/) ))
        call check(nf90_get_var(PrecipConvncInputID, PrecipConvVarID, precip_conv_grid,    &
             start = (/ 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(SfcFluxLHncInputID, SfcFluxLHVarID, sfc_flux_lh_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                            &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(SfcFluxSHncInputID, SfcFluxSHVarID, sfc_flux_sh_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                            &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(SfcFluxSWncInputID, SfcFluxSWVarID, sfc_flux_sw_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                            &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(SfcFluxLWDncInputID, SfcFluxLWDVarID, sfc_flux_lwd_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                               &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(SfcFluxLWUncInputID, SfcFluxLWUVarID, sfc_flux_lwu_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                               &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(SfcFluxQFncInputID, SfcFluxQFVarID, sfc_flux_ocean_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                               &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(ToaFluxSWncInputID, ToaFluxSWVarID, toa_flux_sw_grid,      &
             start = (/ 1, 1, int(nsteps)/),                                               &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(ToaFluxLWUncInputID, ToaFluxLWUVarID, toa_flux_lwu_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                               &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(DiabCondncInputID, DiabCondVarID, dt_temp_cond_grid,  &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

        call check(nf90_get_var(DiabSWncInputID, DiabSWVarID, dt_temp_sw_grid,  &
             start = (/ 1, 1, 1, int(nsteps)/),                               &
             count = (/ local_num_lon, local_num_lat, num_lev, 1/) ))

        call check(nf90_get_var(DragMOncInputID, DragMOVarID, drag_coeff_mo_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                        &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(DragLHncInputID, DragLHVarID, drag_coeff_lh_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                        &
             count = (/ local_num_lon, local_num_lat, 1/) ))

        call check(nf90_get_var(DragSHncInputID, DragSHVarID, drag_coeff_sh_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                        &
             count = (/ local_num_lon, local_num_lat, 1/) ))

     endif

     ! addition fridoo sept 2012 bucket hydrology
     if(bucket) then
        call check(nf90_get_var(BucketDepthncInputID, BucketDepthVarID, bucket_depth_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                        &
             count = (/ local_num_lon, local_num_lat, 1/) ))
        call check(nf90_get_var(BucketDepthConvncInputID, BucketDepthConvVarID, bucket_depth_conv_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                        &
             count = (/ local_num_lon, local_num_lat, 1/) ))
        call check(nf90_get_var(BucketDepthCondncInputID, BucketDepthCondVarID, bucket_depth_cond_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                        &
             count = (/ local_num_lon, local_num_lat, 1/) ))
        call check(nf90_get_var(BucketDepthLHncInputID, BucketDepthLHVarID, bucket_depth_lh_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                        &
             count = (/ local_num_lon, local_num_lat, 1/) ))
        call check(nf90_get_var(BucketDiffusionncInputID, BucketDiffusionVarID, bucket_diffusion_grid,   &
             start = (/ 1, 1, int(nsteps)/),                                        &
             count = (/ local_num_lon, local_num_lat, 1/) ))
     endif


#ifdef GRIB 
  end if
#endif

end subroutine read_data

! ##############################################################################

subroutine read_surf_geopotential(                                            &
                           SGVarName,                 SGInputFileName,        &
                   surf_geopotential)
     
! ------------------------------ input arguments ------------------------------

  character(len = *), intent(in) ::                                           &
       SGVarName,          &   ! Name of surface geopotential variable (netcdf)
       SGInputFileName         ! Name of file with surface geopotential

! ------------------------------ output arguments -----------------------------

  real, intent(out), dimension(:,:) ::                                        &
       surf_geopotential       ! 2-d surface geopotential array

! ------------------------------ local variables ------------------------------

  integer :: SGncInputID, SGVarID

! ------------------------------ executable code ------------------------------
#ifdef GRIB 
  if (Data_source .eq. era40) then

     call read_grib(               1,               surf_geopotential )

    else
#endif
       call check(nf90_open(  SGInputFileName, nf90_nowrite,   SGncInputID ) )
       call check(nf90_inq_varid(SGncInputID, trim(SGVarName), SGVarID))
       call check(nf90_get_var(SGncInputID, SGVarID, surf_geopotential ))
       
       surf_geopotential = surf_geopotential * grav
       
       call check(nf90_close( SGncInputID ))
#ifdef GRIB 
    endif
#endif


  end subroutine read_surf_geopotential

! ##############################################################################

 subroutine get_dimensions( num_lat, num_lon, num_levels, num_Times)
   
! ------------------------------ input arguments ------------------------------

   integer, intent(out) ::                                                    &
        num_lon,           &   !                                              &
        num_lat,           &   ! 
        num_levels,        &
        num_Times
   
! ------------------------------ output arguments -----------------------------

   integer ::                                                                 &
        pbgtotl                ! function that determines total GRIBs

! ------------------------------ executable code ------------------------------

#ifdef GRIB 
   if(data_source .eq. era40) then
      num_lon = 320
      num_lat = 160
      num_levels = 60
      num_times = pbgtotl(GRIBFile) / ( 9*60 + 2 ) ! nine multi-, two single-level fields
   else
#endif
   ! get dimensions of temp, vor, div, and ps
      call check(nf90_Inquire_Variable (TEMPncInputID, TEMPVarID, dimids = dimIDs))
      call check(nf90_Inquire_Dimension(TEMPncInputID, dimIDs(1), len = num_lon  ))
      call check(nf90_Inquire_Dimension(TEMPncInputID, dimIDs(2), len = num_lat  ))
      call check(nf90_Inquire_Dimension(TEMPncInputID, dimIDs(3), len = num_levels))
      call check(nf90_Inquire_Dimension(TEMPncInputID, dimIDs(4), len = num_Times))
#ifdef GRIB 
   endif
#endif
 end subroutine get_dimensions

! ##############################################################################

  subroutine axes_init(      deg_lat,                         sin_lat,        &
                          cos_lat_xy,   sin_lat_xy, cos_lat_yz )

! N.B.: assumes that variable has the same name as the attribute (generally true)

! ------------------------------ output arguments -----------------------------

    real, intent(out), dimension(:) ::                                        &
         deg_lat,          &   ! latitudes in degrees
         sin_lat               ! sine(latitude)

    real, intent(out), dimension(:,:) ::                                      &
         sin_lat_xy,       &   ! xy   sine(latitude)
         cos_lat_xy,       &   ! xy cosine(latitude)
         cos_lat_yz            ! yz cosine(latitude)


! ------------------------------ local variables ------------------------------

    character(len = nf90_max_name) :: lat_name
    integer :: i, k, latVarID
  real, dimension(160) ::                                                     &
       local_era40_lats         ! latitudes of T159 Gaussian grid used in ERA40

! ------------------------------ executable code ------------------------------


    if(data_source .eq. era40) then
       

       data local_era40_lats/                                                &
         -89.142, -88.029, -86.911, -85.791, -84.670, -83.549, -82.428,      &
         -81.307, -80.185, -79.064, -77.943, -76.821, -75.700, -74.578,      &
         -73.457, -72.336, -71.214, -70.093, -68.971, -67.850, -66.728,      &
         -65.607, -64.485, -63.364, -62.242, -61.121, -60.000, -58.878,      &
         -57.757, -56.635, -55.514, -54.392, -53.271, -52.149, -51.028,      &
         -49.906, -48.785, -47.663, -46.542, -45.420, -44.299, -43.177,      &
         -42.056, -40.934, -39.813, -38.691, -37.570, -36.448, -35.327,      &
         -34.205, -33.084, -31.962, -30.841, -29.719, -28.598, -27.476,      &
         -26.355, -25.234, -24.112, -22.991, -21.869, -20.748, -19.626,      &
         -18.505, -17.383, -16.262, -15.140, -14.019, -12.897, -11.776,      &
         -10.654, -9.533,   -8.411,  -7.290,  -6.168,  -5.047,  -3.925,      &
         -2.804,  -1.682,   -0.561,   0.561,   1.682,   2.804,   3.925,      &
          5.047,   6.168,    7.290,   8.411,   9.533,  10.654,  11.776,      &
          12.897, 14.019,   15.140,  16.262,  17.383,  18.505,  19.626,      &
          20.748, 21.869,   22.991,  24.112,  25.234,  26.355,  27.476,      &
          28.598, 29.719,   30.841,  31.962,  33.084,  34.205,  35.327,      &
          36.448, 37.570,   38.691,  39.813,  40.934,  42.056,  43.177,      &
          44.299, 45.420,   46.542,  47.663,  48.785,  49.906,  51.028,      &
          52.149, 53.271,   54.392,  55.514,  56.635,  57.757,  58.878,      &
          60.000, 61.121,   62.242,  63.364,  64.485,  65.607,  66.728,      &
          67.850, 68.971,   70.093,  71.214,  72.336,  73.457,  74.578,      &
          75.700, 76.821,   77.943,  79.064,  80.185,  81.307,  82.428,      &
          83.549, 84.670,   85.791,  86.911,  88.029,  89.142 /

       deg_lat = local_era40_lats
       era40_lats = local_era40_lats

    else
       call check(nf90_Inquire_Dimension(TEMPncInputID, dimIDs(2), name=lat_name))

       ! read in latitude and longitude
       call check(nf90_inq_varid(TEMPncInputID, trim(lat_name),    latVarID) )
       call check(nf90_get_var  (TEMPncInputID, latVarID, deg_lat))
    end if

    sin_lat = sin(deg_lat*pi/180.0)
    
    do i=1, size(cos_lat_xy, 1)
       cos_lat_xy(i,:) = cos(deg_lat*pi/180.0)
    enddo

    do i=1, size(sin_lat_xy, 1)
       sin_lat_xy(i,:) = sin(deg_lat*pi/180.0)
    enddo
    
    do k=1, size(cos_lat_yz, 2)
       cos_lat_yz(:,k) = cos(deg_lat*pi/180.0)
    enddo

  end subroutine axes_init

! ##############################################################################

  subroutine pressure_variables_FMS_init(sigma, local_pk, local_bk)

    real, intent(out), dimension(:) :: sigma, local_pk, local_bk

    integer  ::  pkVarID, bkVarID, k

    allocate(pk(size(sigma)+1))
    allocate(bk(size(sigma)+1))

    ! get pk, bk from temperature input file
    call check(nf90_inq_varid(TEMPncInputID, 'pk',     pkVarID ) )
    call check(nf90_get_var  (TEMPncInputID, pkVarID, local_pk) )
    call check(nf90_inq_varid(TEMPncInputID, 'bk',     bkVarID ) )
    call check(nf90_get_var  (TEMPncInputID, bkVarID, local_bk) )
 
    pk = local_pk
    bk = local_bk
    ! full sigma levels (for now)
    do k=1,size(sigma)
       sigma(k) = .5 * (bk(k) + bk(k+1))
    enddo

  end subroutine pressure_variables_FMS_init

! ##############################################################################

  subroutine pressure_variables_FMS(                                          &
                           surface_p,                          p_half,        &
                           ln_p_half,                          p_full,        &
                           ln_p_full )
    
    real, intent (in),  dimension(:,:)   :: surface_p
    real, intent (out), dimension(:,:,:) :: p_half, ln_p_half, p_full, ln_p_full
    
    !                           --- local variables ---
    real, dimension(size(p_half,1),size(p_half,2)) :: alpha
    integer :: k
    
    !                           --- executable code ---


    do  k=1,size(p_half,3)
       p_half(:,:,k) = pk(k) + bk(k)*surface_p(:,:)
    end do
    
    if(pk(1).eq.0.0 .and. bk(1).eq.0.0) then
       
       do k=2,size(p_half,3)
          ln_p_half(:,:,k) = log(p_half(:,:,k))
       end do
       
       do k=2,size(p_half,3)-1
          alpha  = 1.0  - p_half(:,:,k)*(ln_p_half(:,:,k+1) - ln_p_half(:,:,k))/(p_half(:,:,k+1) - p_half(:,:,k))
          ln_p_full(:,:,k) = ln_p_half(:,:,k+1) - alpha
       end do
    
       ln_p_full(:,:,1) = ln_p_half(:,:,2) - 1.0
       ln_p_half(:,:,1) = 0.0

       
    else
       
       do k=1,size(p_half,3)
          ln_p_half(:,:,k) = log(p_half(:,:,k))
       end do
       
       do k=1,size(p_half,3)-1
          alpha  = 1.0  - p_half(:,:,k)*(ln_p_half(:,:,k+1) - ln_p_half(:,:,k))/(p_half(:,:,k+1) - p_half(:,:,k))
          ln_p_full(:,:,k) = ln_p_half(:,:,k+1) - alpha
       end do
       
    end if
    p_full = exp(ln_p_full)
    
    
  end subroutine pressure_variables_FMS

! ##############################################################################

  subroutine pressure_variables_NCEP_init(sigma, local_pk, local_bk)

    real, intent (out), dimension(:)  :: sigma, local_pk, local_bk

    integer  sigmaVarID, k

    ! get sigma from temperature input file
    call check(nf90_inq_varid(TEMPncInputID, 'level',     sigmaVarID ) )
    call check(nf90_get_var  (TEMPncInputID, sigmaVarID, sigma) )

    local_bk(size(local_bk)) = 1.0
    do k=size(local_bk)-1, 2, -1
       local_bk(k) = 0.5*(sigma(k-1) + sigma(k))
    enddo
    local_bk(1) = 0.5*sigma(1)
    local_pk = 0.0
  end subroutine pressure_variables_NCEP_init

! ##############################################################################

  subroutine pressure_variables_NCEP(surface_p, sigma, p_half, ln_p_half, p_full, ln_p_full)
  
    real, intent (in ), dimension(:,:  )  :: surface_p
    real, intent (out), dimension(:,:,:)  :: p_half, ln_p_half, p_full, ln_p_full
    real, intent (in), dimension(:    )  :: sigma
    
    integer :: k

    
    do  k=1,size(sigma)
       p_full(:,:,k) = sigma(k)*surface_p(:,:)
    end do
    ln_p_full = log(p_full)

    p_half(:,:,size(sigma)+1) = surface_p
    
! there is almost certainly a better way of doing this:
    do k=size(sigma), 2, -1
       p_half(:,:,k) = 0.5*(p_full(:,:,k-1) + p_full(:,:,k))
    enddo
    
    p_half(:,:,1) = 0.5*(p_full(:,:,1))
    
    ln_p_half = log(p_half)

  end subroutine pressure_variables_NCEP

!###############################################################################

  subroutine pressure_variables_CCM3_init(sigma, local_pk, local_bk)
    
    real, intent(out), dimension(:) :: sigma, local_pk, local_bk

    integer ::       &
         hyaiVarID,  & ! hybrid A coefficient at interfaces (half levels) var ID
         hybiVarID,  & ! hybrid B coefficient at interfaces (half levels) var ID
         hyamVarID,  & ! hybrid A coefficient at midpoints (full levels) var ID
         hybmVarID,  & ! hybrid B coefficient at midpoints (full levels) var ID
         p0VarID,    & ! reference pressure variable ID
         levVarID      ! level (full levels) variable ID
         
    real, dimension(size(sigma)) :: lev        ! level array

    allocate( hyam(size(sigma)))
    allocate( hybm(size(sigma)))
    allocate( hyai(size(sigma)+1))
    allocate( hybi(size(sigma)+1))
    
    ! get coefficients from temperature input file
    call check(nf90_inq_varid(TEMPncInputID, 'hyam',     hyamVarID ) )
    call check(nf90_get_var  (TEMPncInputID, hyamVarID, hyam) )
    
    call check(nf90_inq_varid(TEMPncInputID, 'hybm',     hybmVarID ) )
    call check(nf90_get_var  (TEMPncInputID, hybmVarID, hybm) )
    
    call check(nf90_inq_varid(TEMPncInputID, 'hyai',     hyaiVarID ) )
    call check(nf90_get_var  (TEMPncInputID, hyaiVarID, hyai) )
    
    call check(nf90_inq_varid(TEMPncInputID, 'hybi',     hybiVarID ) )
    call check(nf90_get_var  (TEMPncInputID, hybiVarID, hybi) )

    call check(nf90_inq_varid(TEMPncInputID, 'P0',       p0VarID ) )
    call check(nf90_get_var  (TEMPncInputID, p0VarID, p0) )
    
    call check(nf90_inq_varid(TEMPncInputID, 'lev',      levVarID ) )
    call check(nf90_get_var  (TEMPncInputID, levVarID, lev) )
    
    local_pk = hyai
    local_bk = hybi
    
! generate sigma
    sigma = lev/1000.0
    
  end subroutine pressure_variables_CCM3_init

!###############################################################################

  subroutine pressure_variables_ERA40_init(                                   &
                               sigma,                        local_pk,        &
                            local_bk )
    
    real, intent(out), dimension(:) ::                                        &
         sigma,            &   ! sigma on full levels
         local_pk,         &   ! hybrid pressure coefficient
         local_bk              ! hybrid sigma coefficient
         
! ------------------------------ local variables ------------------------------

    integer :: k ! level index

    real, dimension(60) ::                                                    &
         era40_hyam,       &   ! 
         era40_hybm  

    real, dimension(61) ::                                                    &
         era40_hyai,       &   ! 
         era40_hybi  
    
    allocate( hyam(size(sigma))) ! hybrid A coeff at midpoints (full lev)
    allocate( hybm(size(sigma))) ! hybrid B coeff at midpoints (full lev)
    allocate( hyai(size(sigma)+1)) ! hybrid A coeff at interface (half lev)
    allocate( hybi(size(sigma)+1)) ! hybrid B coeff at interface (half lev)
    
! ------------------------------ local variables ------------------------------
    
    data era40_hyam/                                                          &
          10.00000,     28.21708,     49.98032,     78.55990,    113.95865,   &
         156.40275,    206.49757,    265.36379,    334.81736,    417.65684,   &
         518.15328,    641.97983,    795.39773,    985.47705,   1220.98141,   &
        1512.77059,   1874.28858,   2322.20062,   2877.15401,   3564.72489,   &
        4416.61251,   5443.47153,   6641.50297,   8013.73416,   9547.94455,   &
       11205.26776,  12907.86276,  14563.16574,  16089.64450,  17426.43396,   &
       18534.04796,  19391.86413,  19991.39244,  20330.53822,  20412.95420,   &
       20247.54377,  19847.90837,  19231.75199,  18420.17664,  17437.18529,   &
       16309.17935,  15064.35888,  13732.06950,  12342.20552,  10924.65395,   &
        9508.77831,   8120.84150,   6793.32107,   5544.40298,   4397.38096,   &
        3370.10331,   2476.36323,   1725.33559,   1121.01763,    661.63000,   &
         339.05225,    138.24579,     36.66772,      3.68805,      0.00000 /

    data era40_hybm/                                                          &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000357072,             &
       0.0002582745,  0.0011050095,  0.0033758375,  0.0079925049,             &
       0.0157437621,  0.0271913249,  0.0426640634,  0.0623432068,             &
       0.0863148755,  0.1145456939,  0.1468931238,  0.1831184931,             &
       0.2228975265,  0.2658305349,  0.3114553888,  0.3592572773,             &
       0.4086814167,  0.4591434825,  0.5100402350,  0.5607618683,             &
       0.6107037189,  0.6592766211,  0.7059171133,  0.7501002561,             &
       0.7913508267,  0.8292546819,  0.8634693838,  0.8937345668,             &
       0.9198849561,  0.9418612920,  0.9597198053,  0.9736451215,             &
       0.9839613161,  0.9911418983,  0.9958234116,  0.9988145314 /
    
    data era40_hyai/                                                          &
           0.00000,     20.00000,     38.42530,     63.64780,     95.63700,   &
         134.48300,    180.58400,    234.77900,    298.49600,    373.97200,   &
         464.61800,    575.65100,    713.21800,    883.66000,   1094.83000,   &
        1356.47000,   1680.64000,   2082.27000,   2579.89000,   3196.42000,   &
        3960.29000,   4906.71000,   6018.02000,   7306.63000,   8765.05000,   &
       10376.12000,  12077.40000,  13775.30000,  15379.80000,  16819.50000,   &
       18045.20000,  19027.70000,  19755.10000,  20222.20000,  20429.90000,   &
       20384.50000,  20097.40000,  19584.30000,  18864.80000,  17961.40000,   &
       16899.50000,  15706.40000,  14411.10000,  13043.20000,  11632.80000,   &
       10209.50000,   8802.36000,   7438.80000,   6144.32000,   4941.78000,   &
        3850.91000,   2887.70000,   2063.78000,   1385.91000,    855.36200,   &
         467.33300,    210.39400,     65.88920,      7.36774,      0.00000,   &
           0.00000 /

    data era40_hybi/                                                          &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000000000,  0.0000000000,  0.0000000000,  0.0000000000,             &
       0.0000758235,  0.0004613950,  0.0018151600,  0.0050811200,             &
       0.0111429000,  0.0206779000,  0.0341212000,  0.0516904000,             &
       0.0735338000,  0.0996747000,  0.1300230000,  0.1643840000,             &
       0.2024760000,  0.2439330000,  0.2883230000,  0.3351550000,             &
       0.3838920000,  0.4339630000,  0.4847720000,  0.5357100000,             &
       0.5861680000,  0.6355470000,  0.6832690000,  0.7287860000,             &
       0.7715970000,  0.8112530000,  0.8473750000,  0.8796570000,             &
       0.9078840000,  0.9319400000,  0.9518220000,  0.9676450000,             &
       0.9796630000,  0.9882700000,  0.9940190000,  0.9976300000,             &
       1.0000000000 /


    hyam = era40_hyam
    hybm = era40_hybm
    hyai = era40_hyai
    hybi = era40_hybi

    do k=1, size(sigma)
       sigma(k) = hyam(k) + 101325.*hybm(k)
    enddo

! generate sigma
    sigma = sigma/101325.

    local_pk = hyai
    local_bk = hybi
    
  end subroutine pressure_variables_ERA40_init

!###############################################################################

  subroutine pressure_variables_CCM3(surface_p, p_half, ln_p_half, p_full, ln_p_full )
  
    real, intent (in ), dimension(:,:  )  :: surface_p
    real, intent (out), dimension(:,:,:)  :: p_half, p_full
    real, intent (out), dimension(:,:,:)  :: ln_p_half, ln_p_full
    
    integer k
    
    ! pressure on full levels
    do k=1,size(p_full,3)
       p_full(:,:,k) = p0*hyam(k)+surface_p(:,:)*hybm(k)
    enddo

    ln_p_full = log(p_full)

    ! pressure on half levels
    do k=1, size(p_half,3)
       p_half(:,:,k) = p0*hyai(k)+surface_p(:,:)*hybi(k)
    enddo

    do k=2,size(p_half,3)
      ln_p_half(:,:,k) = log(p_half(:,:,k))
    end do
    ln_p_half(:,:,1) = 0.0

  end subroutine pressure_variables_CCM3

!###############################################################################
  subroutine pressure_variables_ERA40(surface_p, p_half, ln_p_half, p_full, ln_p_full )
  
    real, intent (in ), dimension(:,:  )  :: surface_p
    real, intent (out), dimension(:,:,:)  :: p_half, p_full
    real, intent (out), dimension(:,:,:)  :: ln_p_half, ln_p_full
    
    integer k
    
    ! pressure on full levels
    do k=1,size(p_full,3)
       p_full(:,:,k) = hyam(k)+surface_p(:,:)*hybm(k)
    enddo

    ln_p_full = log(p_full)

    ! pressure on half levels
    do k=1, size(p_half,3)
       p_half(:,:,k) = hyai(k)+surface_p(:,:)*hybi(k)
    enddo

    do k=2,size(p_half,3)
      ln_p_half(:,:,k) = log(p_half(:,:,k))
    end do
    ln_p_half(:,:,1) = 0.0

  end subroutine pressure_variables_ERA40

!##############################################################################
  subroutine check(status)
    
    ! checks error status after each netcdf, prints out text message each time
    !   an error code is returned. 
    
    integer, intent(in) :: status
    
    if(status /= nf90_noerr) then 
       write(*, *) trim(nf90_strerror(status))
    end if
  end subroutine check

!##############################################################################

  subroutine check_error_status( emos_routine, error_code )

    character*(*) emos_routine
    integer error_code

    if ( emos_routine .eq. 'pbopen' ) then
       if ( error_code .ne. 0 ) then
          if (      error_code .eq. -1 ) then
             write(*,'(A80)') 'pbopen: could not open the file'
             stop
          else if ( error_code .eq. -2 ) then
             write(*,'(A80)') 'pbopen: invalid file name'
             stop
          else if ( error_code .eq. -3 ) then
             write(*,'(A80)') 'pbopen: invalid open mode specified'
             stop
          else
             write(*,'(A80)') 'pbopen: unspecified error code' 
             stop
          end if
       else
            return
         end if 
      else if ( emos_routine .eq. 'pbclose' ) then
         if ( error_code .ne. 0 ) then
           if (      error_code .eq. -1 ) then
             write(*,'(A80)') 'pbclose: error closing the file'
             stop
           else
             write(*,'(A80)') 'pbclose: unspecified error code'
             stop
           end if
         else
            return
         end if
      else if ( emos_routine .eq. 'pbwrite' ) then
         if ( error_code .lt. 0 ) then
           if (      error_code .eq. -1 ) then
             write(*,'(A80)') 'pbwrite: error writing to the file'
             stop
           else
             write(*,'(A80)') 'pbwrite: unspecified error code'
             stop
           end if
         else
            return
         end if
      else if ( emos_routine .eq. 'pbgrib' ) then
         if ( error_code .ne. 0 ) then
           if ( error_code .eq. -2 ) then
             write(*,'(A80)') 'pbgrib: error in file-handling'
             write(*,'(A80)') 'pbgrib: e.g. file may contain a trun-'
             write(*,'(A80)') 'cated product                        '
             stop
           else if ( error_code .eq. -3 ) then
             write(*,'(A80)') 'pbgrib: karray is not large enough'
             stop
           else
             write(*,'(A80)') 'pbgrib: unspecified error code'
             stop
           end if
         else
            return
         end if
      else if ( emos_routine .eq. 'intf' ) then
         if ( error_code .ne. 0 ) then
             write(*,'(A80)') 'intf: call not successful'
             stop
         else
            return
         end if
      else if ( emos_routine .eq. 'intuvp' ) then
         if ( error_code .ne. 0 ) then
             write(*,'(A80)') 'intuvp: call not successful'
             stop
         else
            return
         end if
      else if ( emos_routine .eq. 'intout' ) then
         if ( error_code .ne. 0 ) then
             write(*,'(A80)') 'intout: call not successful'
             stop
         else
            return
         end if
      else if ( emos_routine .eq. 'intin' ) then
         if ( error_code .ne. 0 ) then
             write(*,'(A80)') 'intin: call not successful'
             stop
         else
            return
         end if
      else
             write(*,'(A80)') 'check_error_status: ' // emos_routine
             write(*,'(A80)') 'check_error_status: unknown emos routine'
             write(*,'(A80)') 'check_error_status handles pbopen,      '
             write(*,'(A80)') 'pbclose, pbwrite, pbgrib, intf, intuvp, '
             write(*,'(A80)') 'intout, and intin                       '
             stop
      end if

      return
    end subroutine check_error_status

!###############################################################################

#ifdef GRIB 

    subroutine get_field_data( words_ioarray, iarray, rarray )
!      ------------------------------------------------------------------
!     'Implicit none':
      implicit none

      integer words_ioarray
      integer iarray(words_ioarray)
      real rarray(words_ioarray)

!     ------------------------------------------------------------------
!     Declare variables for use in gribex:
!     ------------------------------------
!     gribex( ksec0, ksec1, ksec2, psec2, ksec3, psec3,
!    +        ksec4, psec4, klenp, kgrib, kleng,
!    +        kword, hoper, kret )

!     ksec0(1): number of 'octets' in GRIB record
!     ksec0(2): GRIB edition number
!     -------------------------------------------
      integer ksec0(2)

!     Product definition section:
!     ---------------------------
      integer ksec1(1024)

!     Grid description section:
!     -------------------------
      integer ksec2(1024)
      real psec2(1024)

!     Bitmap section:
!     ---------------
      integer ksec3(2)
      real psec3(2)

!     Binary data section:
!     --------------------
      integer ksec4(1024)

!     Field data values:
!     ------------------
     real psec4(1024*1024/4)
!     integer klenp
!     parameter( klenp = words_per_mb )
!     (will use real array 'rarray' argument for 'psec4',
!     and  integer words_ioarray argument for 'klenp')

!     Grib formatted binary data:
!     ---------------------------
!     integer kgrib(words_per_mb)
!     integer kleng
!     parameter( kleng = words_per_mb )
!     (will use integer array 'iarray' argument for 'kgrib',
!     and integer words_ioarray argument for 'kleng')

!     Return value from encoding,
!     (number of elements of kgrib occupied by coded data):
!     -----------------------------------------------------
      integer kword
!     (Used only for encoding, but have to supply in any case.)

!     Coding options:
!     ---------------
      character*1 hoper

!     Error handling:
!     ---------------
      integer kret
      kret = 0
!     Abort if an error is encountered. (An initial nonzero value
!     will not cause gribex to abort even if an error is encountered.)
!     Upon return, kret &lt; 0 implies handle warning, kret &gt; 0 implies
!     handle error.

!     ------------------------------------------------------------------
!     Begin 'executable' code:
!     ------------------------------------------------------------------
!     Decode GRIB record:

      hoper = 'D'

      call gribex( ksec0, ksec1, ksec2, psec2, ksec3, psec3, &
                 ksec4, rarray, words_ioarray, iarray, words_ioarray, &
                  kword, hoper, kret )


!!$      call grprs0(ksec0)
!!$      call grprs1(ksec0,ksec1)
!!$!
!!$      if( (ksec1(5).eq.0).or.(ksec1(5).eq.64) ) then
!!$         write (*,*) 'grdemo: no section 2 in grib message.'
!!$      else
!!$         call grprs2(ksec0,ksec2,psec2)
!!$      endif
!!$ 
!!$     if( (ksec1(5).eq.0).or.(ksec1(5).eq.128) ) then
!!$         write (*,*) 'grdemo: no section 3 in grib message.'
!!$      else
!!$         call grprs3(ksec0,ksec3,psec3)
!!$      endif
!!$     call grprs4(ksec0,ksec4,psec4)
!
!     ------------------------------------------------------------------
!     If kret upon return is &gt; 0, print error message and stop program.

      if (kret .gt. 0 ) then
        write(*,'(A80)') 'gribex: error encountered'
        write(*,'(A80)') 'stopping execution in subroutine ' // &
                        'get_field_data'
        stop
      else
        return
      end if

!     ------------------------------------------------------------------
      return
    end subroutine get_field_data

! ##############################################################################
   


    subroutine read_grib(      GRIBID,                          outarr )



! ------------------------------ input variables ------------------------------

      integer, intent(in) ::                                                  &
           GRIBID

! ------------------------------ output variables -----------------------------

      real, dimension(:, :) ::                                                &
           outarr
      
! ------------------------------ local variables ------------------------------

    integer ::                                                                &
         bytes_per_kb,     &  ! GRIB (bytes per kilobyte)
         bytes_per_mb,     &  ! GRIB (bytes per megabyte)
         bytes_per_wd,     &  ! GRIB (bytes per word)
         words_per_mb         ! GRIB (words per megabyte)

!                         --- GRIB variables ---

    parameter( bytes_per_kb = 1024 )
    parameter( bytes_per_mb = bytes_per_kb * bytes_per_kb )
    parameter( bytes_per_wd = 4)
    parameter( words_per_mb = bytes_per_mb / bytes_per_wd)

    integer, dimension( words_per_mb ) ::                                     &
         iarraysphh,       &   ! GRIB integer array with harmonic information
         iarraygrid            ! GRIB integer array with grid information

    real, dimension( words_per_mb ) ::                                        &
         rarraysphh,       &   ! GRIB real array with real harmonic data
         rarraygrid            ! GRIB real array with real gridpoint data
             
    integer ::                                                                &
         words_into_grid,  &   ! words converted into gridpoint array
         ierror,           &   ! error code for GRIB functions
         lenout,           &   ! length (bytes) of GRIB product
         i,                &   ! longitude counter
         j,                &   ! latitude counter
         m                     ! grib array index
    
    integer, dimension(5) ::                                                  &
         GRIBOff               ! GRIB offsets
    
    integer ::                                                                &
         pbgget,           &   ! GRIB reading function
         intf,             &   ! GRIB interpolation initialization
         intout                ! GRIB interpolation

    integer, dimension(80) ::                                                 &
         intv                  ! GRIB integer array

    real ::                                                                   &
         realv(80),        &   ! GRIB real array
         area(4)

    character*30, dimension(80) ::                                            &
         charv                 ! GRIB character array

!                         --- executable code ---

       data area/ 4*0.0 /

       ! set up output grid (160x320 regular Gaussian grid)
       ierror = intout(       'area',                            intv,        &
                               area,                            charv )
       call check_error_status(                                               &
                            'intout',                          ierror )

       ierror = intout(                                                       &
            'user_regular_gaussian',                               80,        &
                              realv,                            charv )
       call check_error_status(                                               &
                            'intout',                          ierror )

       ierror = intout(     'g_lats',                            intv,        &
                  -era40_lats(1:80),                            charv )
       call check_error_status(                                               &
                            'intout',                          ierror )

       ierror = intout(                                                       &
                       'truncation',                              159,        &
                              realv,                            charv )
       call check_error_status(                                               &
                            'intout',                          ierror )

       ierror = pbgget(     GRIBFile,                      iarraysphh,        &
                        words_per_mb,                          GRIBID )

       words_into_grid = words_per_mb
       
       ierror = intf(     iarraysphh,                    words_per_mb,        &
                          rarraysphh,                      iarraygrid,        &
                     words_into_grid,                      rarraygrid )
       call check_error_status(                                               &
                            'intf',                            ierror )
       
       call get_field_data(                                                   &
                        words_per_mb,                      iarraygrid,        &
                          rarraygrid )
       m = 0
       do j = 160, 1, -1
          do i = 1, 320
             m = m + 1
             outarr(i, j) = rarraygrid(m)
          enddo
       enddo

     end subroutine read_grib
#endif

end module input_mod
