module output_mod

  use vars
  use typesizes, only: bytesizesOK
  use netcdf
  use constants_and_switches_mod, only: gspval

implicit none
private

public :: create_output_file, write_variables, close_output_file 

  integer, parameter ::                                                       &
       nvars = 1002,       &  ! number of output fields (max)
       theta_switch = 1,   &  ! theta variable switch
       sigma_switch = 2,   &  ! sigma variable switch
       surface_switch = 3, &  ! surface variable switch
       zon_wave_switch = 4,&  ! zonal wavenumber variable switch
       spectra_switch = 5, &  ! spectral variable switch
       pdf_switch = 6,     &  ! pdf (lat, lev, bin)
       spectra_switch_3d = 7  ! 3d spectral (wn, wn, height)

  integer ::                                                                  &
       ncFileID,           &  ! file ID
       latDimID,           &  ! latitude file ID
       sigmaDimID,         &  ! sigma dimension ID
       thetaDimID,         &  ! theta dimension ID
       zon_wavenDimID,     &  ! zonal wavenumber dimension ID 
       legendreDimID,      &  ! legendre dimension ID
       binDimID,           &  ! pdf bin dimension ID
       timeDimID,          &  ! time dimension ID
       timeVarID,          &  ! time variable ID
       latVarID,           &  ! latitude variable ID
       sigmaVarID,         &  ! sigma variable ID
       thetaVarID,         &  ! theta variable ID
       zon_wavenVarID,     &  ! zonal wavenumber variable ID
       legendreVarID,      &  ! legendre variable ID
       binVarID,           &  ! pdf bin variable ID
       statVarID(nvars)       ! array of field variable IDs 

contains
! #############################################################################

  subroutine create_output_file(                                              &
                             num_lat,                         num_lev )

! creates output file and axis variables
    
!----------------------------- input arguments --------------------------------
    
    integer, intent(in) :: num_lat, num_lev
    
!----------------------------- local variables --------------------------------

!----------------------------- executable code --------------------------------


    

    if(.not. byteSizesOK()) then
       stop 'open_and_initialize: Compiler does not appear to support required kinds of variables'
    end if

    ! -----------   create netCDF file and define variables  -------------
    ! create netcdf file
    call check(nf90_create(path = trim(OutputFileName),                       &
         cmode = nf90_share, ncid = ncFileID))

    ! global attributes
    call check(nf90_put_att(ncFileID, nf90_global, "title",                   &
         "OfflineDiag Analyses") )
       
       ! global attributes
    call check(nf90_put_att(ncFileID, nf90_global, "title",                   &
         "OfflineDiag Analyses") )


    ! Define axes
    ! -----------
    ! latitude
    call check(NF90_DEF_DIM(ncFileID, "lat", num_lat, latDimID) )
    call check(nf90_def_var(ncFileID, "lat", nf90_double, latDimID,           &
         latVarID) )  
    call check(nf90_put_att(ncFileID, latVarID,                               &
         "long_name", "latitude") )
    call check(nf90_put_att(ncFileID, latVarID,                               &
         "units", "degree") )
       
    ! sigma level
    call check(NF90_DEF_DIM(ncFileID, "sigma", num_lev, sigmaDimID) )
    call check(nf90_def_var(ncFileID, "sigma", nf90_double, sigmaDimID,       &
         sigmaVarID))
    call check(nf90_put_att(ncFileID, sigmaVarID,                             &
         "long_name", "sigma level") )
    call check(nf90_put_att(ncFileID, sigmaVarID,                             &
         "units", "none") )
    
    ! potential temperature
    call check(NF90_DEF_DIM(ncFileID, "theta", MaxIsentrLev, thetaDimID) )
    call check(nf90_def_var(ncFileID, "theta", nf90_double, thetaDimID, &
         thetaVarID))
    call check(nf90_put_att(ncFileID, thetaVarID,                             &
         "long_name", "potential temperature level") )
    call check(nf90_put_att(ncFileID, thetaVarID,                             &
         "units", "Kelvin") )
    
    ! zonal wavenumber
    call check(nf90_def_dim(ncFileID, "zon_waven", num_fourier+1,             &
         zon_wavenDimID) )
    call check(nf90_def_var(ncFileID, "zon_waven", nf90_double,               &
         zon_wavenDimID, zon_wavenVarID) )
    call check(nf90_put_att(ncFileID, zon_wavenVarID,                         &
         "long_name", "zonal wavenumber") )
    call check(nf90_put_att(ncFileID, zon_wavenVarID,                         &
         "units", "none") )
       
    ! legendre
    call check(nf90_def_dim(ncFileID, "legendre", num_fourier+1,              &
         legendreDimID) )
    call check(nf90_def_var(ncFileID, "legendre", nf90_double,                &
         legendreDimID, legendreVarID) )
    call check(nf90_put_att(ncFileID, legendreVarID,                          &
         "long_name", "legendre wavenumber") )
    call check(nf90_put_att(ncFileID, legendreVarID,                          &
         "units", "none") )

    ! pdf bin for relative humidity
    call check(NF90_DEF_DIM(ncFileID, "rhum_bin", num_bin, binDimID) )
    call check(nf90_def_var(ncFileID, "rhum_bin", nf90_double, binDimID,      &
         binVarID))
    call check(nf90_put_att(ncFileID, binVarID,                               &
         "long_name", "relative humidity bin") )
    call check(nf90_put_att(ncFileID, binVarID,                               &
         "units", "none") )
       
    ! time
    call check(NF90_DEF_DIM(ncFileID, "times", 1, timeDimID) )
    call check(nf90_def_var(ncFileID, "times", nf90_double, timeDimID,        &
         timeVarID) )  
    call check(nf90_put_att(ncFileID, timeVarID,                              &
         "long_name", "number of instants") )
       
    call initialize_variables

      
    
     end subroutine create_output_file
    
! #############################################################################

     subroutine initialize_variables

    ! Define variables
    ! ----------------

      
    call variable_init("p",                                                   &
         "Pressure on Sigma Levels",                                          &
         "Pa",                                                                &
         sigma_switch,                                                        &
         statVarID(1))

    call variable_init("p_sfc",                                               &
         "Surface Pressure",                                                  &
         "Pa",                                                                &
         surface_switch,                                                      &
         statVarID(2))

    call variable_init("temp_sfc",                                            &
         "Surface Temperature",                                               &
         "K",                                                                 &
         surface_switch,                                                      &
         statVarID(3))

    call variable_init("specific_vol",                                        &
         "Specific Volume",                                                   &
         "m.m.m/Kg",                                                          &
         sigma_switch,                                                        &
         statVarID(4))

    call variable_init("u",                                                   &
         "Zonal Wind",                                                        & 
         "m/s",                                                               &
         sigma_switch,                                                        &
         statVarID(5))

    call variable_init("u_var",                                               &
         "Zonal Wind Variance",                                               &
         "m.m/s/s",                                                           &
         sigma_switch,                                                        &
         statVarID(6))

    call variable_init("u_barotr_var",                                        &
         "Barotropic Zonal Wind Variance",                                    &
         "m.m/s/s",                                                           &
         surface_switch,                                                      &
         statVarID(7))

    call variable_init("u_mrdnl_flux",                                        &
         "Meridional Zonal Wind Flux",                                        &
         "m.m/s/s",                                                           &
         sigma_switch,                                                        &
         statVarID(8))

    call variable_init("u_vrtcl_flux",                                        &
         "Vertical Zonal Wind Flux",                                          &
         "m/s/s",                                                             &
         sigma_switch,                                                        &
         statVarID(9))

    call variable_init("u_eddy_mrdnl_flux",                                   &
         "Meridional Eddy Zonal Wind Flux",                                   &
         "m.m/s/s",                                                           &
         sigma_switch,                                                        &
         statVarID(10))

    call variable_init("u_eddy_vrtcl_flux",                                   &
         "Vertical Eddy Zonal Wind Flux",                                     &
         "m/s/s",                                                             &
         sigma_switch,                                                        &
         statVarID(11))

    call variable_init("ucos_zon_psd_barotr",                                 &
         "Zonal Wind (cos) Barotropic Zonal Power Spectral Density",          &
         "m.m/s/s",                                                           &
          zon_wave_switch,                                                    &
          statVarID(12))

    call variable_init("v",                                                   &
         "Meridional Wind",                                                   &
         "m/s",                                                               &
         sigma_switch,                                                        &
         statVarID(13))

    call variable_init("v_var",                                               &
         "Meridional Wind Variance",                                          &
         "m.m/s/s",                                                           &
         sigma_switch,                                                        &
         statVarID(14))

    call variable_init("v_barotr_var",                                        &
         "Barotropic Meridional Wind Variance",                               & 
         "m.m/s/s",                                                           &
         surface_switch,                                                      &
         statVarID(15))

    call variable_init("vcos_zon_psd_barotr",                                 &
         "Meridional Wind (cos) Barotropic Zonal Power Spectral Density",     &
         "m.m/s/s",                                                           &
         zon_wave_switch,                                                     &
         statVarID(16))

    call variable_init("w",                                                   &
         "Vertical Wind; w = \dot\sigma",                                     &
         "1/s",                                                               &
         sigma_switch,                                                        &
         statVarID(17))

    call variable_init("w_var",                                               &
         "Vertical Wind Variance",                                            &
         "1/s/s",                                                             &
         sigma_switch,                                                        &
         statVarID(18))


    call variable_init("streamfctn",                                          &
         "Eulerian-Mean Meridional Mass Flux Streamfunction",                 &
         "kg/s",                                                              &
         sigma_switch,                                                        &
         statVarID(19))

    call variable_init("streamfctn_TEM",                                      &
         "TE-Mean Meridional Mass Flux Streamfunction",                       &
         "kg/s",                                                              &
         sigma_switch,                                                        &
         statVarID(20))

    call variable_init("streamfctn_TEM_modif",                                &
         "Modified TE-Mean Meridional Mass Flux Streamfunction",              &
         "kg/s",                                                              &
         sigma_switch,                                                        &
         statVarID(21))

    call variable_init("mom_drag_coeff",                                      &
         "Aerodynamic drag coefficient for momentum flux at surface",         &
         "m/s",                                                               &
         surface_switch,                                                      &
         statVarID(221))

    call variable_init("temp",                                                &
         "Temperature",                                                       &
         "K",                                                                 &
         sigma_switch,                                                        &
         statVarID(22))

    call variable_init("temp_var",                                            &
         "Temperature Variance",                                              &
         "K.K",                                                               &
         sigma_switch,                                                        &
         statVarID(23))

    call variable_init("temp_mrdnl_flux",                                     &
         "Meridional Temperature Flux",                                       &
         "K.m/s",                                                             &
         sigma_switch,                                                        &
         statVarID(24))

    call variable_init("temp_vrtcl_flux",                                     &
         "Vertical Temperature Flux",                                         &
         "K/s",                                                               &
         sigma_switch,                                                        &
         statVarID(25))

    call variable_init("temp_eddy_mrdnl_flux",                                &
         "Meridional Eddy Temperature Flux",                                  &
         "K.m/s",                                                             &
         sigma_switch,                                                        &
         statVarID(26))

    call variable_init("temp_eddy_vrtcl_flux",                                &
         "Vertical Eddy Temperature Flux",                                    &
         "K/s",                                                               &
         sigma_switch,                                                        &
         statVarID(27))

    call variable_init("temp_drag_coeff",                                     &
         "Aerodynamic drag coefficient for sensible heat flux at surface",    &
         "m/s",                                                               &
         surface_switch,                                                      &
         statVarID(227))

    call variable_init("temp_conv_tend",                                   &
         "Temperature Tendency due to Subgrid-Scale Convection",           &
         "K/s",                                                            &
         sigma_switch,                                                     &
         statVarID(84))

    call variable_init("temp_diff_tend",                                   &
         "Temperature Tendency due to Subgrid-Scale Diffusion",            &
         "K/s",                                                            &
         sigma_switch,                                                     &
         statVarID(85))

    call variable_init("temp_rad_tend",                                    &
         "Temperature Tendency due to Radiative Fluxes",                   &
         "K/s",                                                            &
         sigma_switch,                                                     &
         statVarID(86))

    call variable_init("z",                                                   &
         "Geopotential Height",                                               &
         "m",                                                                 &
         sigma_switch,                                                        &
         statVarID(28))

    call variable_init("z_var",                                               &
         "Geopotential Height Variance",                                      &
         "m.m",                                                               &
         sigma_switch,                                                        &
         statVarID(228))

    call variable_init("z_mrdnl_flux",                                        &
         "Meridional Geopotential Height Flux",                               &
         "m.m/s",                                                             &
         sigma_switch,                                                        &
         statVarID(29))

    call variable_init("z_vrtcl_flux",                                        &
         "Vertical Geopotential Height Flux",                                 &
         "m/s",                                                               &
         sigma_switch,                                                        &
         statVarID(30))

    call variable_init("z_eddy_mrdnl_flux",                                   &
         "Meridional Eddy Geopotential Height Flux",                          &
         "m.m/s",                                                             &
         sigma_switch,                                                        &
         statVarID(31))

    call variable_init("z_eddy_vrtcl_flux",                                   &
         "Vertical Eddy Geopotential Height Flux",                            &
         "m/s",                                                               &
         sigma_switch,                                                        &
         statVarID(32))

    call variable_init("pot_temp",                                            &
         "Potential Temperature",                                             &
         "K",                                                                 &
         sigma_switch,                                                        &
         statVarID(33))

    call variable_init("pot_temp_var",                                        &
         "Potential Temperature Variance",                                    &
         "K.K",                                                               &
         sigma_switch,                                                        &
         statVarID(34))

    call variable_init("pot_temp_mrdnl_flux",                                 &
         "Meridional Potential Temperature Flux",                             &
         "K.m/s",                                                             &
         sigma_switch,                                                        &
         statVarID(35))
    
    call variable_init("pot_temp_vrtcl_flux",                                 &
         "Vertical Potential Temperature_Flux",                               &
         "K/s",                                                               &
         sigma_switch,                                                        &
         statVarID(36))

    call variable_init("pot_temp_eddy_mrdnl_flux",                            &
         "Meridional Eddy Potential Temperature Flux",                        &
         "K.m/s",                                                             &
         sigma_switch,                                                        &
         statVarID(37))

    call variable_init("pot_temp_eddy_vrtcl_flux",                            &
         "Vertical Eddy Potential Temperature Flux",                          &
         "K/s",                                                               &
         sigma_switch,                                                        &
         statVarID(38))

    call variable_init("pot_temp_mrdnl_deriv",                                &
         "Meridional Potential Temperature Gradient",                         &
         "K/m",                                                               &
         sigma_switch,                                                        &
         statVarID(39))

    call variable_init("pot_temp_vrtcl_deriv",                                &
         "Vertical Potential Temperature Gradient",                           &
         "K/Pa",                                                              &
         sigma_switch,                                                        &
         statVarID(40))


    call variable_init("pot_temp_psd",                                        &
        "Potential Temperature Horizontal Power Spectral Density",            &
        "K.K",                                                                &
        spectra_switch_3d,                                                    &
        statVarID(41))

    call variable_init("buoyancy_freq",                                       &
        "Buoyancy Frequency",                                                 &
        "1/s",                                                                &
        sigma_switch,                                                         &
        statVarID(42))

    call check(nf90_put_att(ncFileID, statVarID(42),                          &
         "_FillValue", gspval ) )

   call variable_init("eape2eke_flux_psd",                                    &
       "EAPE to EKE Conversion Horizontal Power Spectral Density",            &
       "K.m/s",                                                               &
       spectra_switch_3d,                                                     &
       statVarID(43))

    if(moisture) then

       call variable_init("shum",                                             &
            "Specific Humidity",                                              &
            "kg/kg",                                                          &
            sigma_switch,                                                     &
            statVarID(44))
       

       call variable_init("shum_var",                                         &
            "Specific Humidity Variance",                                     &
            "kg/kg kg/kg",                                                    &
            sigma_switch,                                                     &
            statVarID(45))

       call variable_init("shum_mrdnl_flux",                                  &
            "Meridional Specific Humidity Flux",                              &
            "kg.m/kg/s",                                                      &
            sigma_switch,                                                     &
            statVarID(46))

       call variable_init("shum_vrtcl_flux",                                  &
            "Vertical Specific Humidity Flux",                                &
            "kg/kg/s",                                                        &
            sigma_switch,                                                     &
            statVarID(47))

       call variable_init("shum_eddy_mrdnl_flux",                             &
            "Meridional Specific Humidity Flux",                              &
            "kg/kg.m/s",                                                      &
            sigma_switch,                                                     &
            statVarID(48))

       call variable_init("shum_eddy_vrtcl_flux",                             &
            "Vertical Specific Humidity Flux",                                &
            "kg/kg/s",                                                        &
            sigma_switch,                                                     &
            statVarID(49))

       call variable_init("shum_cond_tend",                                   &
            "Specific Humidity Tendency due to Grid-Scale Condensation",      &
            "kg/kg/s",                                                        &
            sigma_switch,                                                     &
            statVarID(50))

       call variable_init("shum_conv_tend",                                   &
            "Specific Humidity Tendency due to Subrid-Scale Convection",      &
            "kg/kg/s",                                                        &
            sigma_switch,                                                     &
            statVarID(51))


       call variable_init("shum_cond_prob",                                   &
            "Probability of Large-Scale Condensation Event",                  &
            "none",                                                           &
            sigma_switch,                                                     &
            statVarID(52))

       call variable_init("shum_diff_tend",                                   &
            "Specific Humidity Tendency due to Subgrid-Scale Diffusion",      &
            "kg/kg/s",                                                        &
            sigma_switch,                                                     &
            statVarID(152))

       call variable_init("shum_conv_prob",                                   &
            "Probability of Subgrid-Scale Convection Event",                  &
            "none",                                                           &
            sigma_switch,                                                     &
            statVarID(53))


       call variable_init("shum_zon_psd_barotr",                              &
            "Specific Humidity Barotropic Zonal Power Spectral Density",      &
            "none",                                                           &
            zon_wave_switch,                                                  &
            statVarID(54))

       call variable_init("shum_psd",                                         &
            "Specific Humidity Zonal Power Spectral Density",                 &
            "kg/kg.kg/kg",                                                    &
            spectra_switch_3d,                                                &
            statVarID(55))


       call variable_init("shum_sat",                                         &
            "Saturated Specific Humidity",                                    &
            "kg/kg",                                                          &
            sigma_switch,                                                     &
            statVarID(56))


       call variable_init("shum_sat_var",                                     &
            "Saturation Specific Humidity Variance",                          &
            "kg/kg.kg/kg",                                                    &
            sigma_switch,                                                     &
            statVarID(57))

       call variable_init("shum_drag_coeff",                                  &
            "Aerodynamic drag coefficient for moisture flux at surface",      &
            "m/s",                                                            &
            surface_switch,                                                   &
            statVarID(257))

       call variable_init("rhum",                                             &
            "Relative Humidity",                                              &
            "none",                                                           &
            sigma_switch,                                                     &
            statVarID(58))

       call variable_init("rhum_psd",                                         & 
            "Relative Humidity Horizontal Power Spectral Density",            &
            "none",                                                           &
            spectra_switch_3d,                                                &
            statVarID(59))

       call variable_init("rhum_pdf",                                         &
            "Probability Distribution for Relative Humidity",                 &
            "none",                                                           &
            pdf_switch,                                                       &
            statVarID(60))
       
       call variable_init("pot_temp_eqv",                                     &
            "Equivalent Potential Temperature",                               &
            "K",                                                              &
            sigma_switch,                                                     &
            statVarID(61))
   
       call variable_init("pot_temp_eqv_var",                                 &
            "Equivalent Potential Temperature Variance",                      &
            "K.K",                                                            &
            sigma_switch,                                                     &
            statVarID(62))
   

       call variable_init("pot_temp_eqv_mrdnl_flux",                          &
            "Meridional Equivalent Potential Temperature Flux",               &
            "K.m/s",                                                          &
            sigma_switch,                                                     &
            statVarID(63))

       call variable_init("pot_temp_eqv_vrtcl_flux",                          &
            "Vertical equivalent Potential Temperature flux",                 &
            "K/s",                                                            &
            sigma_switch,                                                     &
            statVarID(64))

       call variable_init("pot_temp_eqv_eddy_mrdnl_flux",                     &
            "Meridional Eddy Equivalent Potential Temperature Flux",          &
            "K.m/s",                                                          &
            sigma_switch,                                                     &
            statVarID(65))

       call variable_init("pot_temp_eqv_eddy_vrtcl_flux",                     &
            "Vertical Eddy Equivalent Potential Temperature Flux",            &
            "K/s",                                                            &
            sigma_switch,                                                     &
            statVarID(66))
   
       call variable_init("pot_temp_eqv_mrdnl_deriv",                         &
            "Meridional Equivalent Potential Temperature Gradient",           &
            "K/m",                                                            &
            sigma_switch,                                                     &
            statVarID(67))
   
       call variable_init("pot_temp_eqv_vrtcl_deriv",                         &
            "Vertical Equivalent Potential Temperature Gradient",             &
            "K/Pa",                                                           &
            sigma_switch,                                                     &
            statVarID(68))

       call variable_init("pot_temp_eqv_psd",                                       &
            "Equivalent Potential Temperature Horizontal Power Spectral Density",   &
            "K.K",                                                                  &
            spectra_switch_3d,                                                      &
            statVarID(69))


       call variable_init("pot_temp_sat_eqv",                                 &
            "Saturated Equivalent Potential Temperature",                     &
            "K",                                                              &
            sigma_switch,                                                     &
            statVarID(70))

       call variable_init("precip_cond",                                      &
            "Precipitation Rate due to Large-Scale Condensation",             &
            "kg/m/m/s",                                                       &
            surface_switch,                                                   &
            statVarID(71))

       call variable_init("precip_conv",                                      &
            "Precipitation Rate due to Subgrid-Scale Convection",             &
            "kg/m/m/s",                                                       &
            surface_switch,                                                   &
            statVarID(72))

       call variable_init("precip_daily_above_threshold",                     &
            "Daily Precipitation Rate for Days Above Threshold",              &
            "kg/m/m/s",                                                       &
            surface_switch,                                                   &
            statVarID(73))

       call variable_init("precip_daily_above_threshold_prob",                &
            "Probability of Daily Precipitation Rate Being Above Threshold",  &
            "none",                                                           &
            surface_switch,                                                   &
            statVarID(74))

       call variable_init("sw_flux_sfc",                                      &
            "Surface Shortwave Radiation",                                    &
            "W/m/m",                                                          &
            surface_switch,                                                   &
            statVarID(75))

       call variable_init("lwu_flux_sfc",                                     &
            "Surface Upward Longwave Radiation",                              &
            "W/m/m",                                                          &
            surface_switch,                                                   &
            statVarID(76))

       call variable_init("lwd_flux_sfc",                                     &
            "Surface Downward Longwave Radiation",                            &
            "W/m/m",                                                          &
            surface_switch,                                                   &
            statVarID(77))

       call variable_init("lh_flux_sfc",                                      &
            "Surface Latent Heat Flux",                                       &
            "W/m/m",                                                          &
            surface_switch,                                                   &
           statVarID(78))

       call variable_init("sh_flux_sfc",                                      &
            "Surface Sensible Heat Flux",                                     &
            "W/m/m",                                                          &
            surface_switch,                                                   &
            statVarID(79))

       call variable_init("ocean_qflux",                                      &
            "Ocean Q-flux",                                                   &
             "W/m/m",                                                         &
             surface_switch,                                                  &
             statVarID(80))
 
       call variable_init("sw_flux_toa",                                      &
            "Top-of-Atmosphere Shortwave Radiation",                          &
            "W/m/m",                                                          &
            surface_switch,                                                   &
            statVarID(81))

       call variable_init("lwu_flux_toa",                                     &
            "Top-of-Atmosphere Longwave Radiation",                           &
            "W/m/m",                                                          &
            surface_switch,                                                   &
            statVarID(82))

       call variable_init("temp_cond_tend",                                   &
            "Tenperature Tendency due to Grid-Scale Condensation",            &
            "K/s",                                                            &
            sigma_switch,                                                     &
            statVarID(83))

       call variable_init("temp_sw_tend",                                     &
            "Temperature Tendency due to Shortwave Radiative Fluxes",         &
            "K/s",                                                            &
            sigma_switch,                                                     &
            statVarID(87))

    endif


    if(isentrope) then

       call variable_init("p_isent_nw",                                          &
            "Pressure",                                                          &
            "Pa",                                                                &
            theta_switch,                                                        &
            statVarID(101))

       call variable_init("p_var_isent_nw",                                      &
            "Pressure Variance",                                                 &
            "Pa.Pa",                                                             &
            theta_switch,                                                        &
            statVarID(102))

       call variable_init("exner_isent_nw",                                      &
            "Exner function",                                                    &
            "J/kg/K",                                                            &
            theta_switch,                                                        &
            statVarID(103))

       call variable_init("p_exner_isent_nw",                                    &
            "Pressure * Exner function",                                         &
            "kg.m.m.m/s/s/s/s/K",                                                &
            theta_switch,                                                        &
            statVarID(104))

       call variable_init("dens_isent_nw",                                       &
            "Isentropic Density",                                                &
            "kg/m/m/K",                                                          &
            theta_switch,                                                        &
            statVarID(105))

       call variable_init("u_isent",                                             &
            "Zonal Wind",                                                        &
            "dens.m/s",                                                          &
            theta_switch,                                                        &
            statVarID(106))

       call variable_init("u_var_isent",                                         &
            "Zonal Wind Variance",                                               &
            "dens.m.m/s/s",                                                      &
            theta_switch,                                                        &
            statVarID(107))
 
       call variable_init("v_isent",                                             &
            "Meridional Wind",                                                   &
            "dens.m/s",                                                          &
            theta_switch,                                                        &
            statVarID(108))

       call variable_init("v_var_isent",                                         &
            "Meridional Wind Variance",                                          &
            "dens.m.m/s/s",                                                      &
            theta_switch,                                                        &
            statVarID(109))

       call variable_init("v_geostr_isent",                                      &
            "Meridional Geostrophic Wind",                                       &
            "m/s",                                                               &
            theta_switch,                                                        &
            statVarID(110))

       call variable_init("streamfctn_isent",                                    &
            "Streamfunction",                                                    &
            "kg/s",                                                              &
            theta_switch,                                                        &
            statVarID(112))

       call variable_init("pot_vor_isent",                                       &
            "Potential Vorticity",                                               &
            "1/s/",                                                              &
            theta_switch,                                                        &
            statVarID(113))

       call variable_init("pot_vor_mrdnl_flux_isent",                            &
            "Meridional Potential Vorticity Flux",                               &
            "m/s/s",                                                             &
            theta_switch,                                                        &
            statVarID(114))

       call variable_init("pot_vor_eddy_mrdnl_flux_isent",                       &
            "Meridional Eddy Potential Vorticity Flux",                          &
            "m/s/s",                                                             &
            theta_switch,                                                        &
            statVarID(115))

       call variable_init("mont_isent",                                          &
            "Montgomery Function",                                               &
            "dens.J/kg",                                                         &
            theta_switch,                                                        &
            statVarID(116))

       call variable_init("mont_var_isent",                                      &
            "Montgomery Function Variance",                                      &
            "dens.J.J/kg/kg",                                                    &
            theta_switch,                                                        &
            statVarID(117))

       call variable_init("pot_temp_geostr_eddy_mrdnl_flux_sfc",                 &
            "Surface Geostrophic Eddy Potential Temperature Flux ",              &
            "K.m/s",                                                             &
            surface_switch,                                                      &
            statVarID(118))

       call variable_init("pot_temp_sfc_pdf",                                     &
            "Probability Distribution for Potential Temperature at Surface",      &
            "none",                                                               &
            theta_switch,                                                         &
            statVarID(119))

        call variable_init("pot_enstrophy_isent_avg",                                &
             "Mean of Squared Potential Vorticity",                                  &
            "ISU",                                                                   &
            theta_switch,                                                            &
            statVarID(127))

        call variable_init("mrdnl_eddy_pot_enstrophy_flux_isent_avg",                 &
            "Meridional Eddy Flux of potential vorticity",                            &
            "ISU",                                                                    &
            theta_switch,                                                             &
            statVarID(128))
 

        if(moisture) then
           call variable_init("shum_isent",                                           &
                "Specific Humidity",                                                  &
                "dens.kg/kg",                                                         &
                theta_switch,                                                         &
                statVarID(120))

           call variable_init("shum_mrdnl_flux_isent",                                &
                "Meridional Specific Humidity Flux",                                  &
                "dens.kg/kg.m/s",                                                     &
                theta_switch,                                                         &
                statVarID(121))


           call variable_init("shum_eddy_mrdnl_flux_isent",                           &
                "Meridional Eddy Specific Humidity Flux",                             &
                "dens.kg/kg.m/s",                                                     &
                theta_switch,                                                         &
                statVarID(122))

           call variable_init("rhum_isent",                                           &
                "Relative Humidity",                                                  &
                "dens",                                                               &
                theta_switch,                                                         &
                statVarID(123))

           call variable_init("shum_sat_isent",                                       &
                "Saturated Specific Humidity",                                        &
                "dens.kg/kg",                                                         &
                theta_switch,                                                         &
                statVarID(124))

           call variable_init("shum_sat_mrdnl_flux_isent",                            &
                "Meridional Saturated Specific Humidity Flux",                        &
                "dens.kg/kg.m/s",                                                     &
                theta_switch,                                                         &
                statVarID(125))


           call variable_init("shum_sat_eddy_mrdnl_flux_isent",                       &
                "Meridional Eddy Saturated Specific Humidity Flux",                   &
                "dens.kg/kg.m/s",                                                     &
                theta_switch,                                                         &
                statVarID(126))
        endif  
    
    endif

     ! added by fridooo sept 2012
    if(bucket) then
       call variable_init("bucket_depth",                                     &
         "bucket depth",    &
         "m",                                                               &
         surface_switch,                                                      &
         statVarID(301))
 
       call variable_init("bucket_conv_depth",                                     &
         "bucket depth variation induced by convection",    &
         "m/s",                                                               &
         surface_switch,                                                      &
         statVarID(302))
 
       call variable_init("bucket_cond_depth",                                     &
         "bucket depth variation induced by condensation",    &
         "m/s",                                                               &
         surface_switch,                                                      &
         statVarID(303))

       call variable_init("bucket_lh_depth",                                     &
         "bucket depth variation induced by condensation",    &
         "m/s",                                                               &
         surface_switch,                                                      &
         statVarID(304))

        call variable_init("bucket_diffusion",                                     &
         "diffused surface water depth",    &
         "m/s",                                                               &
         surface_switch,                                                      &
         statVarID(305))

     endif
     ! end fridooo

    call variable_init("spec_baroclinic",                                             &
         "spectral coefficients of the baroclinic streamfunction",                    &
         "m^4/s^2",                                                                      &
         spectra_switch,                                                              &
         statVarID(140))

    call variable_init("spec_barotropic",                                             &
         "spectral coefficients of the barotropic streamfunction",                    &
         "m^4/s^2",                                                                      &
         spectra_switch,                                                              & 
         statVarID(141))

    call variable_init("non_lin_spec",                                        &
         "spectrum of non-linear interaction",                                &
         "1/s/s/s",                                                           &
         spectra_switch_3d,                                                   &
         statVarID(142))

    call variable_init("conv_zon_spec",                                       &
         "Zonal conversion spectrum",                                         &
         "m*m*m/s/kg",                                                        &
         zon_wave_switch,                                                     &
         statVarID(143))

    call variable_init("non_lin_eddy_spec",                                   &
         "spec of non-lin interaction with m\neq0",                           &
         "1/s/s/s",                                                           &
         spectra_switch_3d,                                                   &
         statVarID(144))

    call variable_init("coriolis_spec",                                       &
         "spectrum related to beta effect",                                   &
         "1/s/s/s",                                                           &
         spectra_switch_3d,                                                   &
         statVarID(145)) 

    call variable_init("bt_non_lin_spec",                                     &
         "spectrum of non-linear interaction, barotropic component",          &
         "1/s/s/s",                                                           &
         spectra_switch,                                                      &
         statVarID(146))

    call variable_init("bt_non_lin_eddy_spec",                                &
         "spec of non-lin interaction with m\neq0, barotropic component",     &
         "1/s/s/s",                                                           &
         spectra_switch,                                                      &
         statVarID(147))

    call variable_init("temp_spec",                                           &
         "spectrum of temperature variance",                                  &
         "1/K/K",                                                             &
         spectra_switch_3d,                                                   &
         statVarID(148))

    call variable_init("temp_non_lin_spec",                                   &
         "spectrum of non-linear interaction (temperature advection)",        &
         "1/s/K/K",                                                           &
         spectra_switch_3d,                                                   &
         statVarID(149))

    call variable_init("temp_non_lin_eddy_spec",                              &
         "spec of non-lin interaction with m\neq0 (temperature advection)",   &
         "1/s/K/K",                                                           &
         spectra_switch_3d,                                                   &
         statVarID(150))



    call check(nf90_enddef(ncfileID))
 
    
  end subroutine initialize_variables
  
! #############################################################################

  subroutine variable_init(name, long_name, units, switch, VarID)
  
    character(len=*), intent(in) :: name, long_name, units
    integer, intent(in) :: switch
    integer, intent(out) :: VarID

       if(switch == theta_switch) then
          call check(nf90_def_var(ncFileID, trim(name), nf90_double,          &
               (/ latDimID, thetaDimID /), VarID ) )
       elseif(switch == sigma_switch) then
          call check(nf90_def_var(ncFileID, trim(name), nf90_double,          &
               (/ latDimID, sigmaDimID /), VarID ) )
       elseif(switch == surface_switch) then
          call check(nf90_def_var(ncFileID, trim(name), nf90_double,          &
               (/ latDimID /), VarID ) )
       elseif(switch == zon_wave_switch) then
          call check(nf90_def_var(ncFileID, trim(name), nf90_double,          &
               (/latDimID, zon_wavenDimID /), VarID ) )
       elseif(switch == spectra_switch) then
          call check(nf90_def_var(ncFileID, trim(name), nf90_double,          &
               (/zon_wavenDimID, legendreDimID /), VarID ) )
       elseif(switch == spectra_switch_3d) then
          call check(nf90_def_var(ncFileID, trim(name), nf90_double,          &
               (/zon_wavenDimID, legendreDimID, sigmaDimID /), VarID ) )
       elseif(switch == pdf_switch) then
          call check(nf90_def_var(ncFileID, trim(name), nf90_double,          &
               (/latDimID, sigmaDimID, binDimID /), VarID ) )
       endif
       
       call check(nf90_put_att(ncFileID, VarID,                               &
            "long_name", trim(long_name)))
       call check(nf90_put_att(ncFileID, VarID,                               &
            "units", trim(units)))

     end subroutine variable_init

! #############################################################################

  subroutine write_variables(num_lat,                     num_fourier,        &
                             num_lev,                       num_times )

!----------------------------- input arguments --------------------------------

    integer, intent(in) ::                                                    &
         num_lat,          &   !
         num_fourier,      &   !
         num_lev               !
    
    real, intent(in) ::                                                       &
         num_times

!----------------------------- local variables --------------------------------

    real    :: zon_waven(num_fourier+1), legendre(num_fourier+1)
    integer :: i, j, k

!----------------------------- executable code --------------------------------

!------------------          write coordinate axes         --------------------

    
    ! latitudes
    call check(nf90_put_var(ncFileID, latVarID, deg_lat,                      &
         start = (/ 1 /), count = (/ num_lat /) ) )

    ! sigma level
    call check(nf90_put_var(ncFileID, sigmaVarID, sigma ) )

    ! potential temperature
    call check(nf90_put_var(ncFileID, thetaVarID, pot_temp_coord ) )

    call check(nf90_put_var(ncFileID, binVarID, rhum_bin) )

    ! zonal wave number
    do i=0, num_fourier
       zon_waven(i+1) = float(i)
       legendre(i+1) = float(i)
    enddo

    call check(nf90_put_var(ncFileID, zon_wavenVarID, zon_waven ) )
    call check(nf90_put_var(ncFileID, legendreVarID, legendre ) )

    ! times
    call check(nf90_put_var(ncFileID, timeVarID, num_times ) )


    ! --------------             write statistics             -----------------


    call check(nf90_put_var(ncFileID, statVarID(1), p_full_avg,                              &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(2), ps_avg,                                  &
         start = (/ 1 /), count = (/ num_lat /) ))


    call check(nf90_put_var(ncFileID, statVarID(4), specific_vol_avg,                        &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(5), u_avg,                                    &
         start = (/ 1, 1 /), count = (/num_lat, num_lev/) ))
        
    call check(nf90_put_var(ncFileID, statVarID(6), u_var_avg,                               &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(7), u_barotr_var_avg,                        &
         start = (/ 1 /), count = (/ num_lat /) ))

    call check(nf90_put_var(ncFileID, statVarID(8), mrdnl_u_flux_avg,                       &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(9), vrtcl_u_flux_avg,                       &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(10), mrdnl_eddy_u_flux_avg,                   &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(11), vrtcl_eddy_u_flux_avg,                   &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(12), ucos_zon_spec_barotr,                    &
            start = (/ 1, 1 /), count = (/ num_lat, num_fourier /) ))


    call check(nf90_put_var(ncFileID, statVarID(13), v_avg,                                   &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) )) 

    call check(nf90_put_var(ncFileID, statVarID(14), v_var_avg,                               &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(15), v_barotr_var_avg,                        &
         start = (/ 1 /), count = (/ num_lat /) ))

    call check(nf90_put_var(ncFileID, statVarID(16), vcos_zon_spec_barotr,                    &
            start = (/ 1, 1 /), count = (/ num_lat, num_fourier/) ))
    
    call check(nf90_put_var(ncFileID, statVarID(17), w_avg,                                   &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(18), w_var_avg,                               &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(19), sfctn_avg,                               &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(20), TEM_res_circ,                            & 
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(21), MTEM_res_circ,                           &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))


    call check(nf90_put_var(ncFileID, statVarID(22), temp_avg,                                &
         start = (/ 1, 1 /), count = (/num_lat, num_lev/) ))

    call check(nf90_put_var(ncFileID, statVarID(23), temp_var_avg,                            &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))
    
    call check(nf90_put_var(ncFileID, statVarID(24), mrdnl_temp_flux_avg,                     &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(25), vrtcl_temp_flux_avg,                     &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(26), mrdnl_eddy_temp_flux_avg,                &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(27), vrtcl_eddy_temp_flux_avg,                &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(84), dt_temp_conv_avg,                        &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(85), dt_temp_diff_avg,                        &
          start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(86), dt_temp_rad_avg,                         &
          start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(28), z_avg,                                   &
         start = (/ 1, 1 /), count = (/num_lat, num_lev/) ))

     call check(nf90_put_var(ncFileID, statVarID(228), z_var_avg,                             &
          start = (/ 1, 1 /), count = (/num_lat, num_lev/) ))

    call check(nf90_put_var(ncFileID, statVarID(29), mrdnl_z_flux_avg,                        &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(30), vrtcl_z_flux_avg,                        &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(31), mrdnl_eddy_z_flux_avg,                   &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(32), vrtcl_eddy_z_flux_avg,                   &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(33), pot_temp_avg,                            &
         start = (/ 1, 1 /), count = (/num_lat, num_lev/) ))

    call check(nf90_put_var(ncFileID, statVarID(34), pot_temp_var_avg,                        &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(35), mrdnl_pot_temp_flux_avg,                 &
         start = (/ 1, 1 /), count = (/num_lat, num_lev/) )) 

    call check(nf90_put_var(ncFileID, statVarID(36), vrtcl_pot_temp_flux_avg,                 &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(37), mrdnl_eddy_pot_temp_flux_avg,            &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(38), vrtcl_eddy_pot_temp_flux_avg,            &
         start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(39), d_dy_pot_temp_avg,                       &
         start = (/ 1, 1 /), count = (/num_lat, num_lev/) )) 

    call check(nf90_put_var(ncFileID, statVarID(40), d_dp_pot_temp_avg,                       &
         start = (/ 1, 1 /), count = (/num_lat, num_lev/) )) 

    call check(nf90_put_var(ncFileID, statVarID(41), pot_temp_spec,                           &
            start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(42), buoyancy_freq_avg,                       &
         start = (/ 1, 1 /), count = (/num_lat, num_lev/) ))

    call check(nf90_put_var(ncFileID, statVarID(43), eape2eke_conv_spec,                      &
            start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))
   
    if(moisture) then

       call check(nf90_put_var(ncFileID, statVarID(3), ts_avg,                                 &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(44), shum_avg,                             &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))
       
       call check(nf90_put_var(ncFileID, statVarID(45), shum_var_avg,                         &
            start = (/ 1, 1 /), count = (/num_lat, num_lev/) )) 

       call check(nf90_put_var(ncFileID, statVarID(46), mrdnl_shum_flux_avg,                  &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(47), vrtcl_shum_flux_avg,                  &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(48), mrdnl_eddy_shum_flux_avg,             &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(49), vrtcl_eddy_shum_flux_avg,             &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))


       call check(nf90_put_var(ncFileID, statVarID(50), dt_shum_cond_avg,                     &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(51), dt_shum_conv_avg,                     &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

        call check(nf90_put_var(ncFileID, statVarID(152), dt_shum_diff_avg,                     &
             start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(52), shum_cond_prob_avg,                   &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(53), shum_conv_prob_avg,                   &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(54), shum_zon_spec_barotr,                 &
            start = (/ 1, 1 /), count = (/ num_lat, num_fourier /) ))

       call check(nf90_put_var(ncFileID, statVarID(55), shum_spec,                            &
            start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))


       call check(nf90_put_var(ncFileID, statVarID(56), sat_shum_avg,                         &
            start = (/ 1, 1 /), count = (/num_lat, num_lev/) )) 

       call check(nf90_put_var(ncFileID, statVarID(57), sat_shum_var_avg,                     &
            start = (/ 1, 1 /), count = (/num_lat, num_lev/) )) 

       call check(nf90_put_var(ncFileID, statVarID(58), rhum_avg,                             &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))       

       call check(nf90_put_var(ncFileID, statVarID(59), rhum_spec,                            &
            start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(60), rhum_pdf,                             &
            start = (/ 1, 1 /), count = (/num_lat, num_lev, num_bin/) )) 

       call check(nf90_put_var(ncFileID, statVarID(61), pot_temp_e_avg,                       &
            start = (/ 1, 1 /), count = (/num_lat, num_lev/) ))

       call check(nf90_put_var(ncFileID, statVarID(62), pot_temp_e_var_avg,                   &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(63), mrdnl_pot_temp_e_flux_avg,            &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(64), vrtcl_pot_temp_e_flux_avg,            &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(65), mrdnl_eddy_pot_temp_e_flux_avg,       &
            start = (/ 1, 1 /), count = (/num_lat, num_lev/) )) 
   
       call check(nf90_put_var(ncFileID, statVarID(66), vrtcl_eddy_pot_temp_e_flux_avg,       &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))
      
       call check(nf90_put_var(ncFileID, statVarID(67), d_dy_pot_temp_e_avg,                  &
            start = (/ 1, 1 /), count = (/num_lat, num_lev/) )) 

       call check(nf90_put_var(ncFileID, statVarID(68), d_dp_pot_temp_e_avg,                  &
            start = (/ 1, 1 /), count = (/num_lat, num_lev/) )) 

       call check(nf90_put_var(ncFileID, statVarID(69), pot_temp_e_spec,                      &
            start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(70), pot_temp_e_sat_avg,                   &
            start = (/ 1, 1 /), count = (/num_lat, num_lev/) ))

       call check(nf90_put_var(ncFileID, statVarID(71), precip_cond_avg,                     &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(72), precip_conv_avg,                     &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(73), precip_daily_above_threshold_avg,    &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(74), precip_daily_above_threshold_prob_avg,   &
            start = (/ 1 /), count = (/ num_lat /) ))


       call check(nf90_put_var(ncFileID, statVarID(75), sfc_flux_sw_avg,                     &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(76), sfc_flux_lwu_avg,                    &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(77), sfc_flux_lwd_avg,                    &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(78), sfc_flux_lh_avg,                     &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(79), sfc_flux_sh_avg,                     &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(80), sfc_flux_ocean_avg,                  &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(81), toa_flux_sw_avg,                     &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(82), toa_flux_lwu_avg,                    &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(83), dt_temp_cond_avg,                    &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(87), dt_temp_sw_avg,                      &
            start = (/ 1, 1 /), count = (/ num_lat, num_lev /) ))

       call check(nf90_put_var(ncFileID, statVarID(221), drag_coeff_mo_avg,                     &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(227), drag_coeff_sh_avg,                      &
            start = (/ 1 /), count = (/ num_lat /) ))

       call check(nf90_put_var(ncFileID, statVarID(257), drag_coeff_lh_avg,                   &
            start = (/ 1 /), count = (/ num_lat /) ))


    endif

    if(isentrope) then
       call check(nf90_put_var(ncFileID, statVarID(101), p_isent_avg_nw,              &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))
  
       call check(nf90_put_var(ncFileID, statVarID(102), p_var_isent_avg_nw,          &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(103), exner_isent_avg_nw,          &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(104), p_exner_isent_avg_nw,        &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(105), dens_isent_avg_nw,           &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(106), u_isent_avg,                 &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(107), u_var_isent_avg,             &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(108), v_isent_avg,                 &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))         

       call check(nf90_put_var(ncFileID, statVarID(109), v_var_isent_avg,             &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(110), v_geostr_isent_avg,          &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(112), sfctn_isent_avg,             &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(113), pot_vor_isent_avg,           &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(114), mrdnl_pot_vor_flux_isent_avg,        &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(115), mrdnl_eddy_pot_vor_flux_isent_avg,   &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(116), mont_isent_avg,              &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(117), mont_var_isent_avg,          &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(118), mrdnl_geostr_eddy_pot_temp_flux_sfc_avg,   &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(119), sfc_pot_temp_pdf,                   &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

       call check(nf90_put_var(ncFileID, statVarID(127), pot_enstrophy_isent_avg,        &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))
      
       call check(nf90_put_var(ncFileID, statVarID(128), mrdnl_eddy_pot_enstrophy_flux_isent_avg,        &
            start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) )) 

       if(moisture) then

          call check(nf90_put_var(ncFileID, statVarID(120), shum_isent_avg,                  &
               start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))
 
          call check(nf90_put_var(ncFileID, statVarID(121), mrdnl_shum_flux_isent_avg,       &
               start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

          call check(nf90_put_var(ncFileID, statVarID(122), mrdnl_eddy_shum_flux_isent_avg,  &
               start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

          call check(nf90_put_var(ncFileID, statVarID(123), rhum_isent_avg,                  &
               start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

          call check(nf90_put_var(ncFileID, statVarID(124), sat_shum_isent_avg,              &
               start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))

          call check(nf90_put_var(ncFileID, statVarID(125), mrdnl_sat_shum_flux_isent_avg,   &
               start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))       

          call check(nf90_put_var(ncFileID, statVarID(126), mrdnl_eddy_sat_shum_flux_isent_avg,   &
               start = (/ 1, 1 /), count = (/ num_lat, MaxIsentrLev /) ))       

       endif

    endif   




    call check(nf90_put_var(ncFileID, statVarID(140), bc_sfctn_spec,           &
            start = (/ 1, 1 /), count = (/ num_fourier+1, num_fourier+1 /) ))

    call check(nf90_put_var(ncFileID, statVarID(141), bt_sfctn_spec,           &
            start = (/ 1, 1 /), count = (/ num_fourier+1, num_fourier+1 /) ))

    call check(nf90_put_var(ncFileID, statVarID(142), non_lin_spec,           &
         start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(143), conv_zon_spec_avg,      &
         start = (/ 1, 1 /), count = (/ num_lat, num_fourier /) ))

    call check(nf90_put_var(ncFileID, statVarID(144), non_lin_eddy_spec,      &
         start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(145), coriolis_spec,          &
         start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(146), bt_non_lin_spec,           &
         start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(147), bt_non_lin_eddy_spec,      &
         start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(148), temp_spec,           &
         start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(149), t_non_lin_spec,           &
         start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

    call check(nf90_put_var(ncFileID, statVarID(150), t_non_lin_eddy_spec,      &
         start = (/ 1, 1, 1 /), count = (/ num_fourier+1, num_fourier+1, num_lev /) ))

     

    ! added fridoo sept 2012 (hydrology)
    if(bucket) then
       
        call check(nf90_put_var(ncFileID, statVarID(301), bucket_depth_avg,                   &
         start = (/ 1/), count = (/num_lat /) ))
        call check(nf90_put_var(ncFileID, statVarID(302), bucket_depth_conv_avg,                   &
         start = (/ 1 /), count = (/num_lat /) ))
        call check(nf90_put_var(ncFileID, statVarID(303), bucket_depth_cond_avg,                   &
         start = (/ 1 /), count = (/num_lat /) ))
        call check(nf90_put_var(ncFileID, statVarID(304), bucket_depth_LH_avg,                   &
         start = (/ 1 /), count = (/num_lat /) ))
        call check(nf90_put_var(ncFileID, statVarID(305), bucket_diffusion_avg,                   &
         start = (/ 1 /), count = (/num_lat /) ))
    endif
    ! end fridoo
 
    call check(nf90_sync(ncFileID))
  


  end subroutine write_variables

! #############################################################################

  subroutine close_output_file

    call check(nf90_close(ncFileID))

  end subroutine close_output_file

! #############################################################################

  subroutine check(status)
    
    ! checks error status after each netcdf, prints out text message each time
    !   an error code is returned. 
    
    integer, intent(in) :: status
    
    if(status /= nf90_noerr) then 
       write(*, *) trim(nf90_strerror(status))
    end if
  end subroutine check 
  
end module output_mod






 


         
         










   





