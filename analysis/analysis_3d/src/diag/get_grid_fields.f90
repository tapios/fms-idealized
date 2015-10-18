module get_grid_fields_mod

  ! all subroutines that use harmonic transforms reside here
  
  use constants_and_switches_mod

  use input_mod

  use local_utilities_mod, only :                                             &
                      mrdnl_gradient,                   tetens_sat_mr_mod,    &
            interpolate_half_to_full,            compute_geopotential,        &
                       simple_sat_mr

  implicit none

  public ::         grid_fields_init,                 get_grid_fields,        &
                       surf_gradient,             zonal_spectrum_grid

  integer ::                                                                  &
       lwork,              &   ! len. single-precision work array 
       ldwork,             &   ! len. double-precision work array 
       lshags,             &   ! len. scalar analysis work array, gaussian
       lvhsgs,             &   ! len. vector synthesis work array, gaussian
       lshsgs,             &   ! len. scalar synthesis work array, gaussian
       lvhags,             &   ! len. vector analysis work array, gaussian
       lshaes,             &   ! len. scalar analysis work array, equal
       lvhses,             &   ! len. vector synthesis work array, equal
       lshses,             &   ! len. scalar synthesis work array, equal
       lvhaes,             &   ! len. vector analysis work array, equal
       num_lat,            &   ! number of latitudes
       num_lon,            &   ! number of longitudes
       num_lev,            &   ! number of sigma levels
       num_fourier,        &   ! number of fourier coefficients
       DataIn,             &   ! input data (1=u_and_v, 2=vor_and_div)
       data_source,        &   ! source of data (FMS=1, NCEP=2, CCM3=3)
       ierror                  ! Spherepack error indicator


  real, allocatable, dimension(:) ::                                          &
       work,               &   ! single-precision work array (not saved)
       wsave,              &   ! fft work array (saved)
       wshags,             &   ! scalar analysis work array, gaussian (saved)
       wvhsgs,             &   ! vector synthesis work array, gaussian (saved)
       wshsgs,             &   ! scalar synthesis work array, gaussian (saved)
       wvhags,             &   ! vector analysis work array, gaussian (saved)
       wshaes,             &   ! scalar analysis work array, equal (saved)
       wvhses,             &   ! vector synthesis work array, equal (saved)
       wshses,             &   ! scalar synthesis work array, equal (saved)
       wvhaes,             &   ! vector analysis work array, equal (saved)
       rad_lat,            &   ! latitude array
       pk_local,           &   ! pressure coefficients
       bk_local,           &   ! sigma coefficients
       fnn,                &   ! array of 1/(n*n+1) for spectral conversions
       weight,             &   ! weights for vertical averages    
       dpk,                &   ! used in four_in_one
       dbk                     ! used in four_in_one

  real, allocatable, dimension(:,:,:) ::                             &
       vor_coef_a,       &   ! vorticity cosine spectral coeff.
       vor_coef_b,       &   ! vorticity sine spectral coeff.
       div_coef_a,       &   ! divergence cosine spectral coeff.
       div_coef_b,       &   ! divergence sine spectral coeff.
       u_coef_a,         &   ! U cosine spectral coeff.
       u_coef_b,         &   ! U sine spectral coeff.
       v_coef_a,         &   ! V cosine spectral coeff.
       v_coef_b,         &   ! V sine spectral coeff.
       pot_temp_coef_a,  &   ! single-level pot. temp. cosine spectral coeff.
       pot_temp_coef_b,  &   ! single-level pot. temp. sine spectral coeff.
       shum_coef_a,      &   ! single-level spec. hum. cosine spectral coeff.
       shum_coef_b,      &   ! single-level spec. hum. sine spectral coeff.
       sat_shum_coef_a,  &   ! single-level sat. spec. hum. cos spectral coeff.
       sat_shum_coef_b,  &   ! single-level sat. spec. hum. sine spectral coeff.
       rhum_coef_a,      &   ! single-level rel. hum. cosine spectral coeff.
       rhum_coef_b,      &   ! single-level rel. hum. sine spectral coeff.
       theta_e_coef_a,   &   ! single-level equiv. pot. temp cosine spectral coeff.
       theta_e_coef_b,   &   ! single-level equiv. pot. temp sine spectral coeff.
       w_coef_a,         &   ! single-level omega cosine spectral coeff.
       w_coef_b              ! single-level omega temp. sine spectral coeff.


  double precision, allocatable, dimension(:) ::                              &
       dwork                   ! double-precision work array (not saved)
  

  logical ::                                                                  &
       is_gaussian,        &   ! true for Gaussian grid; false for equal-spaced
       moisture,           &   ! true for moisture analysis
       moist_isentropes        ! true to interpolate to moist isentropes
  
  real, allocatable, dimension(:,:) ::                                        &
       cos_lat_xy,         &   ! xy cosine(latitude) array
       surf_geopot             ! surface geopotential

contains

! ##############################################################################

  subroutine grid_fields_init(                                                &
                        local_DataIn,               local_is_gaussian,        &
                      local_moisture,                                         &
              local_moist_isentropes,               local_data_source,        &
                   local_num_fourier,                   local_num_lat,        &
                       local_num_lon,                   local_num_lev,        &
                    local_cos_lat_xy,               local_surf_geopot,        &
                               sigma,                              pk,        &
                                  bk )

!----------------------------- input arguments --------------------------------

    integer, intent(in) ::                                                    &
         local_DataIn,     &   ! input momentum fields
         local_data_source,&   ! data source (1=FMS, 2=NCEP, 3=CCM3)
         local_num_fourier,&   ! number of fourier coefficients
         local_num_lat,    &   ! number of latitudes
         local_num_lon,    &   ! number of longitudes
         local_num_lev         ! number of sigma levels

    logical, intent(in) ::                                                    &
         local_is_gaussian,&   ! switch for Gaussian grid
         local_moisture,   &   ! switch for moisture
         local_moist_isentropes    ! interpolate to equivalent potential temp.
    
    real, dimension(:,:), intent(in) ::                                       &
         local_cos_lat_xy, &   ! xy cosine(latitude) array
         local_surf_geopot     ! surface geopotential

!----------------------------- output arguments  ------------------------------

    real, dimension(:), intent(out) ::                                        &
         sigma,            &   ! sigma levels
         pk,               &   ! pressure coefficients
         bk                    ! sigma coefficients

!----------------------------- local variables  -------------------------------

    integer l1, l2, i

!----------------------------- executable code  -------------------------------

!
  ! initialize variables global to this module
  num_lat = local_num_lat
  num_lon = local_num_lon
  num_lev = local_num_lev
  num_fourier = local_num_fourier
  DataIn = local_DataIn
  is_gaussian = local_is_gaussian
  data_source = local_data_source
  moisture = local_moisture
  moist_isentropes = local_moist_isentropes

  allocate(  cos_lat_xy(num_lon, num_lat))
  allocate( surf_geopot(num_lon, num_lat))

  cos_lat_xy = local_cos_lat_xy
  surf_geopot = local_surf_geopot

  allocate(  wsave(2*num_lon+15))

  allocate(         vor_coef_a        (num_lon,num_lat,num_lev))
  allocate(         vor_coef_b        (num_lon,num_lat,num_lev))
  allocate(         div_coef_a        (num_lon,num_lat,num_lev))
  allocate(         div_coef_b        (num_lon,num_lat,num_lev))
  allocate(         u_coef_a          (num_lon,num_lat,num_lev))
  allocate(         u_coef_b          (num_lon,num_lat,num_lev))
  allocate(         v_coef_a          (num_lon,num_lat,num_lev))
  allocate(         v_coef_b          (num_lon,num_lat,num_lev))
  allocate(         pot_temp_coef_a   (num_lon,num_lat,num_lev))
  allocate(         pot_temp_coef_b   (num_lon,num_lat,num_lev))
  allocate(         shum_coef_a       (num_lon,num_lat,num_lev))
  allocate(         sat_shum_coef_a   (num_lon,num_lat,num_lev))
  allocate(         rhum_coef_a       (num_lon,num_lat,num_lev))
  allocate(         shum_coef_b       (num_lon,num_lat,num_lev))
  allocate(         sat_shum_coef_b   (num_lon,num_lat,num_lev))
  allocate(         rhum_coef_b       (num_lon,num_lat,num_lev))
  allocate(         theta_e_coef_a    (num_lon,num_lat,num_lev))
  allocate(         theta_e_coef_b    (num_lon,num_lat,num_lev))
  allocate(         w_coef_a          (num_lon,num_lat,num_lev))
  allocate(         w_coef_b          (num_lon,num_lat,num_lev))


  ! Initialize spherepack work arrays
  if(is_gaussian) then
     l1 = min0(num_lat, (num_lon+2) / 2)
     l2 = (num_lat+1) / 2
     lshags = num_lat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*num_lat-l1)-3*l1)/2+num_lon+15
     lshsgs = num_lat*(3*(l1+l2)-2)+(l1-1)*(l2*(2*num_lat-l1)-3*l1)/2+num_lon+15
     lvhags = (num_lat+1)*(num_lat+1)*num_lat/2 + num_lon + 15
     lvhsgs = l1*l2*(num_lat+num_lat-l1+1)+num_lon+15+2*num_lat
     !     
     allocate(     wshags(lshags)  )
     allocate(     wshsgs(lshsgs)  )
     allocate(     wvhsgs(lvhsgs)  )
     allocate(     wvhags(lvhags)  )
!
  !   lwork = num_lat*((2*num_lev+1)*num_lon+4*num_lev*l1+1)
     lwork  = num_lat*(2*num_lev*num_lon+max0(6*(num_lat+1)/2,num_lon))  &
          +4*(min0(num_lat,(num_lon+1)/2)*num_lat*num_lev+num_lat) * 10
     ldwork = (3*num_lat*(num_lat+3)+2)/2 * 10 
     allocate(       work(lwork)   )
     allocate(      dwork(ldwork)  )
     !
     call shagsi(            num_lat,                         num_lon,        &
                              wshags,                          lshags,        &
                                work,                           lwork,        &
                               dwork,                          ldwork,        &
                              ierror )
     if(ierror.ne.0) write(*,*) 'after shagsi, ierror=',ierror
     
     call shsgsi(            num_lat,                         num_lon,        &
                              wshsgs,                          lshsgs,        &
                                work,                           lwork,        &
                               dwork,                          ldwork,        &
                              ierror )
     if(ierror.ne.0) write(*,*) 'after shsgsi, ierror=',ierror
     
     call vhsgsi(            num_lat,                         num_lon,        &
                              wvhsgs,                          lvhsgs,        &
                               dwork,                          ldwork,        &
                              ierror )
     if(ierror.ne.0) write(*,*) 'after vhsgsi, ierror: ',ierror
     
     call vhagsi(            num_lat,                         num_lon,        &
                              wvhags,                          lvhags,        &
                                work,                           lwork,        &
                              ierror )
     if(ierror.ne.0) write(*,*) 'after vhagsi, ierror: ',ierror

  else

     if(num_lon > 3) then
        l1 = min0(num_lat,(num_lon+2)/2)
        l2 = (num_lat+1)/2 
        lshaes = (l1*l2*(num_lat+num_lat-l1+1))/2+num_lon+15
        lvhaes =  l1*l2*(num_lat+num_lat-l1+1)+num_lon+15
        lvhses = lvhaes
        lshses = lvhaes
        allocate(      wshaes(lshaes) )
        allocate(      wvhaes(lvhaes) )
        allocate(      wshses(lshses) )
        allocate(      wvhses(lvhses) )
        
        lwork = num_lat*((2*num_lev+1)*num_lon+4*num_lev*l1+1)
        ldwork = (3*num_lat*(num_lat+3)+2)/2
        allocate(       work(lwork)   )
        allocate(      dwork(ldwork)  )
     
     
        call shaesi(         num_lat,                         num_lon,        &
                              wshaes,                          lshaes,        &
                                work,                           lwork,        &
                               dwork,                          ldwork,        &
                              ierror )
        if(ierror.ne.0) write(*,*) 'after shaesi, ierror: ',ierror
        
        call vhsesi(         num_lat,                         num_lon,        &
                              wvhses,                          lvhses,        &
                                work,                           lwork,        &
                               dwork,                          ldwork,        &
                              ierror )
        if(ierror.ne.0) write(*,*) 'after vhsesi, ierror: ',ierror
     
        call vhaesi(         num_lat,                         num_lon,        &
                              wvhaes,                          lvhaes,        &
                                work,                           lwork,        &
                               dwork,                          ldwork,        &
                              ierror )
        if(ierror.ne.0) write(*,*) 'after vhaesi, ierror: ',ierror
        
        call shsesi(         num_lat,                         num_lon,        &
                              wshses,                          lshses,        &
                                work,                           lwork,        &
                               dwork,                          ldwork,        &
                              ierror )
     endif
  endif
  
  if(num_lon > 3) then
     call hrffti(num_lon, wsave)
  endif

  ! initialize pressure variables


  if(data_source == fms) then
     call pressure_variables_FMS_init(sigma, pk, bk)
  elseif(data_source == NCEP_spectral) then
     call pressure_variables_NCEP_init(sigma, pk, bk)
  elseif(data_source == CCM3) then
     call pressure_variables_CCM3_init(sigma, pk, bk)
  elseif(data_source == ERA40) then
     call pressure_variables_ERA40_init(sigma, pk, bk)
  endif
  
  allocate( pk_local(num_lev+1))
  allocate( bk_local(num_lev+1))
  
  pk_local = pk
  bk_local = bk
  
  ! generate 1/(n*(n+1))
  allocate( fnn(num_fourier+1))
  do i=2, num_fourier+1
     fnn(i) = 1./(i*(i-1))
  enddo

  ! weights for vertical averaging
  allocate( weight(num_lev))
  do i = 1, num_lev
     weight(i) = (pk(i+1) + 101325.*bk(i+1)) - (pk(i) + 101325.*bk(i))
  enddo

  weight = weight / sum(weight)

  ! four_in_one variables
  allocate ( dpk(num_lev))
  allocate ( dbk(num_lev))

  do i=1, num_lev
     dpk(i) = pk_local(i+1) - pk_local(i)
     dbk(i) = bk_local(i+1) - bk_local(i)
  enddo

  end subroutine grid_fields_init

! ##############################################################################


  subroutine get_grid_fields(                                                 &
                              u_grid,                          v_grid,        &
                            vor_grid,                        div_grid,        &
                           temp_grid,                       shum_grid,        &
                             ts_grid,                                         &
                             ps_grid,                       vcos_grid,        &
                              w_grid,                           theta,        &
                             theta_e,                     theta_e_sat,        &
                              p_half,                          p_full,        &
                           mont_grid,                     mont_x_grid,        &
                          sfctn_grid,                           sigma,        &
                          theta_spec,                       conv_spec,        &
                           shum_spec,                   sat_shum_spec,        &
                           rhum_spec,                    theta_e_spec,        &
                   virtual_temp_grid,                                         &
                         geopot_grid,                     sat_mr_grid,        &
                       sat_shum_grid,                       rhum_grid)


! This routine does the following:
!           1. Does whatever conversions are necessary so that U, V,  
!              vorticity, and divergence are all in gridpoint space.
!           2. Calculates spectra of baroclinic and barotropic streamfunctions,
!              conversion of pe -> ke, and potential temperature
!           3. Gets pressure on full and half levels (!model-specific code!)
!           4. Calculates vertical velocity
!           5. Calculates mass flux streamfunction
!
!
!                           --- Input arguments ---

    real, dimension(:,:,:), intent(in)  ::                                    &
         temp_grid,        &   ! gridpoint temperature 
         virtual_temp_grid,&   ! gridpoint virtual temperature 
         shum_grid             ! gridpoint specific humidity


    real, dimension(:, :), intent(in)  ::                                     &
         ts_grid,          &
         ps_grid               ! surface pressure
 
!                        --- Input/Output arguments ---

    real, dimension(:, :, :), intent(inout) ::                                &
         u_grid,           &   ! gridpoint zonal wind 
         v_grid,           &   ! gridpoint meridional wind
         vor_grid,         &   ! gridpoint vorticity
         div_grid              ! gridpoint divergence
    
    real, dimension(:, :, :), intent(inout) ::                                &
         theta_spec,       &   ! potential temperature spectrum
         theta_e_spec,     &   ! equiv. potential temperature spectrum
         shum_spec,        &   ! spectrum of specific humidity
         sat_shum_spec,    &   ! spectrum of saturation specific humidity
         rhum_spec,        &   ! spectrum of relative humidity
         conv_spec             ! spectrum of conversion of pe->ke

!                           --- Output arguments ---

    real, dimension(:, :, :), intent(out)  ::                                 &
         vcos_grid,        &   ! gridpoint meridional wind *cosine(latitude)
         w_grid,           &   ! gridpoint vertical wind ps * \dot\sigma (full)
         theta,            &   ! gridpoint potential temperature
         theta_e ,         &   ! gridpoint equivalent potential temperature
         theta_e_sat,      &   ! gridpoint saturated equivalent potential temperature
         p_full,           &   ! pressure on full levels
         p_half,           &   ! pressure on half levels
         sfctn_grid,       &   ! mass streamfunction
         mont_grid,        &   ! Montgomery streamfunction
         mont_x_grid,      &   ! x-derivative of Montgomery streamfunction
         geopot_grid,      &   ! geopotential
         sat_mr_grid,      &   ! saturation mixing ratio
         sat_shum_grid,    &   ! saturation specific humidity
         rhum_grid             ! relative humidity

    real, dimension(:), intent(in) ::                                        &
         sigma                 ! array of sigma levels XL!used to be out!XL


!                          --- local variables ---


    real, dimension(num_lon, num_lat) ::                                      &
         dx_ps_grid,       &   ! derivative of surface pressure WRT longitude
         dy_ps_grid,       &   ! derivative of surface pressure WRT latitude
         dx_geopot_grid,   &   ! derivative of geopotential WRT longitude
         dy_geopot_grid        ! derivative of geopotential WRT latitude


    real, dimension(num_lon, num_lat, num_lev)   ::                           &
         ln_p_full             ! log of pressure on full model levels

    real, dimension(num_lon, num_lat, num_lev)   ::                           &
         mixing_ratio,     &   ! water vapor mixing ratio
         t_lcl,            &   ! lifted condensation level temperature (K) 
         dry_pot_temp,     &   ! dry potential temperature (K) 
         vapor_pressure        ! water vapor pressure (Pa) 

    real, dimension(num_lon, num_lat, num_lev+1) ::                           &
         ln_p_half,        &   ! log of pressure on half model levels        
         w_grid_half           ! d(sigma)/dt on half model levels

    

    integer :: k               ! level counter

!                           --- executable code ---

! ensure heap variables are intitialized to zero
  vor_coef_a       = 0.0
  vor_coef_b       = 0.0
  div_coef_a       = 0.0
  div_coef_b       = 0.0
  u_coef_a         = 0.0
  u_coef_b         = 0.0
  v_coef_a         = 0.0
  v_coef_b         = 0.0
  pot_temp_coef_a  = 0.0
  pot_temp_coef_b  = 0.0
  shum_coef_a      = 0.0
  shum_coef_b      = 0.0
  sat_shum_coef_a  = 0.0
  sat_shum_coef_b  = 0.0
  rhum_coef_a      = 0.0
  rhum_coef_b      = 0.0
  theta_e_coef_a   = 0.0
  theta_e_coef_b   = 0.0
  w_coef_a         = 0.0
  w_coef_b         = 0.0


!                --- U and V  from vorticity and divergence ---

    if(DataIn .eq. vor_and_div) then     
       
       call scalar_to_spec( vor_grid,                      vor_coef_a,        &
                          vor_coef_b )

       call scalar_to_spec( div_grid,                      div_coef_a,        &
                          div_coef_b )

       call u_v_from_vor_div(                                                 &
                          vor_coef_a,                      vor_coef_b,        &
                          div_coef_a,                      div_coef_b,        &
                              u_grid,                          v_grid )
                          
       u_grid(:,:,:) = u_grid(:,:,:)*radius
       v_grid(:,:,:) = v_grid(:,:,:)*radius

    endif

!               --- vorticity and divergence from U and V  ---

    if(DataIn .eq. u_and_v)then

       call vector_to_spec(                                                   &
                              u_grid,                          v_grid,        &
                            u_coef_a,                        u_coef_b,        &
                            v_coef_a,                        v_coef_b )       

       call vor_div_from_u_v(                                                 &
                            u_coef_a,                        u_coef_b,        &
                            v_coef_a,                        v_coef_b,        &
                            vor_grid,                        div_grid)

       vor_grid(:,:,:) = vor_grid(:,:,:)/radius
       div_grid(:,:,:) = div_grid(:,:,:)/radius
       
    endif

    do k = 1, num_lev 
       vcos_grid(:,:,k) = v_grid(:,:,k)*cos_lat_xy
    end do


    call surf_gradient(      ps_grid,                      dx_ps_grid,        &
                          dy_ps_grid)


!                   --- pressure on full and half levels ---
!
!#############################################################################
!#                       Data-specific call goes here                        #
!#                                                                           #
!#                Make a subroutine that generates pressure on               #
!#               full- and half-levels and sigma on full levels              #
!#                Current routines work for GFDL's FMS (sigma                #
!#                  coordinates on half levels) and for NCEP                 #
!#              reanalysis data (sigma coordinates on full levels)           #
!#              These routines live in the module input_data_mod             #
!#                     since they read from input files.                     #
!#############################################################################
!
    if(data_source.eq.fms) then

       call pressure_variables_FMS(                                           &
                             ps_grid,                          p_half,        &
                           ln_p_half,                          p_full,        &
                           ln_p_full )

    elseif(data_source.eq.NCEP_spectral) then
       call pressure_variables_NCEP(                                          &
                             ps_grid,                           sigma,        &
                              p_half,                       ln_p_half,        &
                              p_full,                       ln_p_full )
       
    elseif(data_source.eq.CCM3) then
       call pressure_variables_CCM3(                                          &
                             ps_grid,                          p_half,        &
                           ln_p_half,                          p_full,        &
                           ln_p_full )
       
    elseif(data_source.eq.ERA40) then
       call pressure_variables_ERA40(                                         &
                             ps_grid,                          p_half,        &
                           ln_p_half,                          p_full,        &
                           ln_p_full )

       ! other model/dataset vertical coordinate defining routines here
    endif

    ! generate vertical velocity
    call vertical_velocity(                                                   &
                            div_grid,                          u_grid,        &
                              v_grid,                         ps_grid,        &
                          dx_ps_grid,                      dy_ps_grid,        &
                         w_grid_half )
    
    ! don't need num_lev+1 (w is zero there...)
    do k=1, num_lev
       w_grid_half(:,:,k) = w_grid_half(:,:,k)/ps_grid(:,:)
    enddo


    ! generate Montgomery streamfunction and its x derivative
    call compute_geopotential(                                                &
                         surf_geopot,               virtual_temp_grid,        &
                           ln_p_half,                       ln_p_full,        &
                         geopot_grid )



    mont_grid =  geopot_grid + cp * temp_grid
           
    do k = 1, num_lev
       
       call surf_gradient(                                                    &
                geopot_grid(:, :, k),                  dx_geopot_grid,        &
                      dy_geopot_grid )

       ! note this is the x derivative of mont_grid at constant theta 
       ! or equivalently the x derivative of gz at constant pressure
       mont_x_grid(:, :, k) = dx_geopot_grid +                                &
            (rdgas * virtual_temp_grid(:, :, k)) / ps_grid * dx_ps_grid
           
    enddo

    ! interpolate vertical velocity to full model levels
    w_grid = interpolate_half_to_full(w_grid_half)
    
   ! get potential temperature on full levels
    dry_pot_temp = temp_grid * (reference_sea_level_pres/p_full)**kappa 

                                                                                                                                       
    ! calculate moist quantities
    if (moisture) then

      if (data_source .eq. ERA40) then
        sat_mr_grid = tetens_sat_mr_mod(temp_grid, p_full)
      else
        sat_mr_grid = simple_sat_mr(temp_grid, p_full)
      end if

      mixing_ratio = shum_grid / (1.0 - shum_grid)

      sat_shum_grid = sat_mr_grid / (sat_mr_grid + 1.0)

      ! use ratio of vapor pressures as definition of relative humidity
      rhum_grid     = mixing_ratio / sat_mr_grid *                   &
           (1 + rvgas/rdgas * sat_mr_grid) / (1 + rvgas/rdgas * mixing_ratio)

      vapor_pressure = mixing_ratio * p_full / (rdgas/rvgas + mixing_ratio)

      ! use Bolton's (1980) formula for the saturation temperature
      t_lcl = const1 / (const2 * log(temp_grid) - log(vapor_pressure * pa_to_mb) - const3) + const4

      ! use Bolton's (1980) formula for the pseudo equiv. pot. temperature
      ! (neglecting the small change to kappa from r)
      ! this differs from the simpler formula in Holton primarily because it accounts for
      ! the pseudoadiabatic process not being at constant temperature
      theta_e = dry_pot_temp * exp(mixing_ratio * (1.0 + const5 * mixing_ratio) *   &
                                   (const6 / t_lcl - const7))

      ! with spectral advection of moisture can get negative values for specific humidity
      ! and the calculation of t_lcl will then give NANs
      where (shum_grid <= 0.0)
           theta_e = dry_pot_temp
      end where

      ! theta_e_sat is the theta_e of a hypothetical saturated atmosphere with the 
      ! same thermal structure
      theta_e_sat = dry_pot_temp * exp(sat_mr_grid * (1.0 + const5 * sat_mr_grid) *   &
                                       (const6 / temp_grid - const7))

    end if

    if(moisture.and.moist_isentropes) then
       theta      = theta_e
    else
       theta      = dry_pot_temp
    end if

   
    ! calculate spectra of potential temperature and pe->ke conversion
    pot_temp_coef_a = 0.0; pot_temp_coef_b = 0.0
    w_coef_a = 0.0; w_coef_b = 0.0
    shum_coef_a = 0.0; shum_coef_b = 0.0
    sat_shum_coef_a = 0.0; sat_shum_coef_b = 0.0
    rhum_coef_a = 0.0; rhum_coef_b = 0.0
    theta_e_coef_a = 0.0; theta_e_coef_b = 0.0


    call scalar_to_spec(    theta,                 pot_temp_coef_a,        &
                     pot_temp_coef_b ) 

    if (moisture) then                                                         
       call scalar_to_spec(                                                   &
                           shum_grid,                     shum_coef_a,        &
                           shum_coef_b ) 
       call scalar_to_spec(                                                   &
                           sat_shum_grid,             sat_shum_coef_a,        &
                           sat_shum_coef_b ) 
       call scalar_to_spec(                                                   &
                           rhum_grid,                     rhum_coef_a,        &
                           rhum_coef_b ) 
       call scalar_to_spec(                                                   &
                           theta_e,                       theta_e_coef_a,     &
                           theta_e_coef_b ) 
    end if

    call scalar_to_spec(      w_grid,                        w_coef_a,        &
                            w_coef_b ) 

    do k = 1, num_lev

       theta_spec(:, :, k) = theta_spec(:, :, k) +                      &
            pot_temp_coef_a(:, :, k)**2 + pot_temp_coef_b(:, :, k)**2
       
       if (moisture) then                                                           
            shum_spec(:, :, k) = shum_spec(:, :, k) +                         &
            shum_coef_a(:, :, k)**2 + shum_coef_b(:, :, k)**2

            sat_shum_spec(:, :, k) = sat_shum_spec(:, :, k) +                 &
            sat_shum_coef_a(:, :, k)**2 + sat_shum_coef_b(:, :, k)**2

            rhum_spec(:, :, k) = rhum_spec(:, :, k) +                         &
            rhum_coef_a(:, :, k)**2 + rhum_coef_b(:, :, k)**2

            theta_e_spec(:, :, k) = theta_e_spec(:, :, k) +                   &
            theta_e_coef_a(:, :, k)**2 + theta_e_coef_b(:, :, k)**2
       end if

       conv_spec(:, :, k) = conv_spec(:, :, k) +                              &
            pot_temp_coef_a(:, :, k) * w_coef_a(:, :, k) +                    &
            pot_temp_coef_b(:, :, k) * w_coef_b(:, :, k) 

    end do

        
    ! compute local streamfunction (the average of which would give the 
    ! Eulerian mean streamfunction)
    sfctn_grid = 0.
    do k=num_lev, 1, -1
       sfctn_grid(:, :, k) = sfctn_grid(:, :, k+1) +                          &
             vcos_grid(:, :, k) * (p_half(:, :, k+1) - p_half(:, :, k)) 
    enddo
    sfctn_grid = 2. * pi * radius / grav * sfctn_grid 

  end subroutine get_grid_fields

! #############################################################################
  
  subroutine surf_gradient(grid, dx_grid, dy_grid)

    ! uses spherepack routines to calculate the gradient for a single level
    
    real,    intent(in),  dimension(:,:) :: grid
    real,    intent(out), dimension(:,:) :: dx_grid, dy_grid
    
!                          --- local variables ---
    real, dimension(num_lon, num_lat) :: coef_a, coef_b
    real, dimension(num_lat, num_lon) :: dx_grid_trans, dy_grid_trans

!                          --- executable code ---

    if(is_gaussian) then
       call shags(           num_lat,                         num_lon,        &
                                   0,                               1,        &
                     transpose(grid),                         num_lat,        &
                             num_lon,                          coef_a,        &
                              coef_b,                         num_lon,        &
                             num_lat,                          wshags,        &
                              lshags,                            work,        &
                               lwork,                          ierror )
       if(ierror.ne.0) write(*,*) 'after shags, ierror:',ierror
   
       call gradgs(          num_lat,                         num_lon,        &
                                   0,                               1,        &
                       dy_grid_trans,                   dx_grid_trans,        &
                             num_lat,                         num_lon,        &
                              coef_a,                          coef_b,        &
                             num_lon,                         num_lat,        &
                              wvhsgs,                          lvhsgs,        &
                                work,                           lwork,        &
                              ierror )
       if(ierror.ne.0) write(*,*) 'after gradgs, ierror:',ierror

    else
       call shaes(          num_lat,                         num_lon,        &
                                  0,                               1,        &
                    transpose(grid),                         num_lat,        &
                            num_lon,                          coef_a,        &
                             coef_b,                         num_lon,        &
                            num_lat,                          wshaes,        &
                             lshaes,                            work,        &
                              lwork,                          ierror )
       if(ierror.ne.0) write(*,*) 'after shaes, ierror:',ierror
      
      
       call grades(          num_lat,                         num_lon,        &
                                   0,                               1,        &
                       dy_grid_trans,                   dx_grid_trans,        &
                             num_lat,                         num_lon,        &
                              coef_a,                          coef_b,        &
                             num_lon,                         num_lat,        &
                              wvhses,                          lvhses,        &
                                work,                           lwork,        &
                              ierror )
      if(ierror.ne.0) write(*,*) 'after grades, ierror:',ierror
      
   endif

   dx_grid = transpose(dx_grid_trans)/radius
   dy_grid = transpose(dy_grid_trans)/radius

 end subroutine surf_gradient

! #############################################################################

 function zonal_spectrum_grid(p_half_zon, grid) result(zon_spectrum)
    ! computes zonal spectrum of grid field
    ! (mass-weighted vertical average)

   real, dimension(:, :), intent(in) ::                                       &
        p_half_zon
   real, dimension(:, :, :), intent(in) ::                                    &
        grid
   real, dimension(num_lat, 0:num_fourier) ::                                 &
        zon_spectrum
   real, dimension(num_lat, num_lon) ::                                       &
        fourier
   real, dimension(num_lat) ::                                                &
        dp_top_sfc, lat_weight
   integer :: k, m, j

   ! get mass weighted vertical average spectrum
   dp_top_sfc = p_half_zon(:, num_lev+1) - p_half_zon(:, 1)
   zon_spectrum = 0.
   do k=1,num_lev
      fourier(:,:) = transpose(grid(:,:,k))
      call hrfftf(num_lat, num_lon, fourier, num_lat, wsave, work)
      fourier = fourier / num_lon
      lat_weight = (p_half_zon(:, k+1) - p_half_zon(:, k)) 
      do j=1,num_lat
         zon_spectrum(j, 0) = zon_spectrum(j, 0) +                            &
              2.0*lat_weight(j)*(fourier(j,1)**2)
         do m=2, num_fourier
            zon_spectrum(j, m-1) = zon_spectrum(j, m-1) +                     &
                 2.0*lat_weight(j)*(fourier(j,2*m-2)**2+fourier(j,2*m-1)**2)
         enddo
      enddo
   enddo
   
 end function zonal_spectrum_grid

!##############################################################################

subroutine vertical_velocity(                                                 &
                                divg,                          u_grid,        &
                              v_grid,                          p_surf,        &
                              dx_psg,                          dy_psg,        &
                                  wg )

real, intent(in),    dimension(:,:,:) :: divg, u_grid, v_grid
real, intent(in),    dimension(:,:  ) :: p_surf, dx_psg, dy_psg
real, intent(out),   dimension(:,:,:) :: wg

!  wg is dimensioned (is:ie, js:je, num_levels+1)
!  wg(:,:,k) = downward mass flux/per unit area across the K+1/2
!  cell boundary. This is the "vertical velocity" in the hybrid coordinate system.
!  When vertical coordinate is pure sigma: wg = psg*d(sigma)/dt

real, dimension(num_lon, num_lat) :: dp, dmean, dmean_tot

integer :: k

dmean_tot = 0.

do k = 1, num_lev
   dp = dpk(k) + dbk(k)*p_surf
   dmean = divg(:,:,k)*dp + dbk(k)*(u_grid(:,:,k)*dx_psg + v_grid(:,:,k)*dy_psg)
   dmean_tot = dmean_tot + dmean
   wg(:,:,k+1) = - dmean_tot
enddo

do k = 2,num_lev
  wg(:,:,k) = wg(:,:,k) + dmean_tot*bk_local(k)
enddo

wg(:,:,1           ) = 0.0
wg(:,:,num_lev+1) = 0.0

return
end subroutine vertical_velocity

!##############################################################################

subroutine streamfunction_coef(                                               &
                                mult,                       in_coef_a,        &
                           in_coef_b,                   bc_sfctn_spec,        &
                       bt_sfctn_spec,                   bt_sfctn_grid,        &
                      hor_sfctn_grid )

    ! generates spectral coefficients for streamfunction
    
!                            --- input arguments ---                           

    real, intent(in), dimension(:) ::                                         &
         mult                  ! 

    real, intent(in), dimension(:, :, :) ::                                   &
         in_coef_a,        &   ! input spectral (cosine) coefficients
         in_coef_b             ! input spectral (sine) coefficients

!                            --- inout arguments ---                         

    real, intent(inout), dimension(:, :) ::                                   &
         bc_sfctn_spec,    &   ! baroclinic streamfunction spectrum
         bt_sfctn_spec         ! barotropic streamfunction spectrum

!                            --- output arguments ---                         

    real, intent(out), dimension(:, :) ::                                     &
         bt_sfctn_grid         ! barotropic streamfunction
    
    real, intent(out), dimension(:, :, :) ::                                  &
         hor_sfctn_grid        ! horizontal streamfunction

!                            --- local variables ---  

    integer :: n, m            !  Legendre and Fourier indices, respectively

    integer :: k               ! level counter

    real, dimension(num_fourier+1, num_fourier+1) ::                          &
         sfctn_coef_a,     &   ! streamfunction a coefficient
         sfctn_coef_b,     &   ! streamfunction b coefficient
         sfctn_coef_a_avg, &   ! weighted average of streamfunction a coeffs
         sfctn_coef_b_avg      ! weighted average of streamfunction b coeffs

    real, dimension(num_lon, num_lat) ::                                      &
         bt_sfctn_coef_a,  &   ! barotropic streamfunction a coefficients
         bt_sfctn_coef_b,  &   ! barotropic streamfunction b coefficients
         hor_sfctn_coef_a, &   ! horizontal streamfunction a coefficients
         hor_sfctn_coef_b      ! horizontal streamfunction b coefficients

    real, dimension(num_lat, num_lon) ::                                      &
         bt_sfctn_grid_trans,& ! transpose of the barotropic streamfunction
         hor_sfctn_grid_trans  ! transpose of horizontal streamfunction

!                            --- executable code --- 

    sfctn_coef_a_avg = 0.0 ; sfctn_coef_b_avg = 0.0

    do k = 1, num_lev
       sfctn_coef_a = 0.0
       sfctn_coef_b = 0.0

       ! m = 0 coefficients
       do n=2, num_fourier+1
          sfctn_coef_a(1, n) = mult(n)*in_coef_a(1, n, k)
          sfctn_coef_b(1, n) = mult(n)*in_coef_b(1, n, k)
       enddo
    
       ! m > 0 coefficients

       do m=2, num_fourier+1
          do n=m, num_fourier+1
             sfctn_coef_a(m, n) = mult(n)*in_coef_a(m, n, k)
             sfctn_coef_b(m, n) = mult(n)*in_coef_b(m, n, k)
          enddo
       enddo

       ! spectral coefficients for baroclinic component
       bc_sfctn_spec = bc_sfctn_spec +                                        &
            weight(k) * (sfctn_coef_a**2 + sfctn_coef_b**2)
       
       ! weighted vertical average of streamfunction coefficients
       sfctn_coef_a_avg = sfctn_coef_a_avg + weight(k) * sfctn_coef_a
       sfctn_coef_b_avg = sfctn_coef_b_avg + weight(k) * sfctn_coef_b

       ! spherepack needs spectral arrays dimensioned (num_lon, num_lat)
       hor_sfctn_coef_a = 0.0 ; hor_sfctn_coef_b = 0.0
       hor_sfctn_coef_a(1:num_fourier+1,1:num_fourier+1) = sfctn_coef_a
       hor_sfctn_coef_b(1:num_fourier+1,1:num_fourier+1) = sfctn_coef_b

       ! horizontal barotropic streamfunction
       if(is_gaussian) then
          call shsgs(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                hor_sfctn_grid_trans,                         num_lat,        &
                             num_lon,                hor_sfctn_coef_a,        &
                    hor_sfctn_coef_b,                         num_lon,        &
                             num_lat,                          wshsgs,        &
                              lshsgs,                            work,        &
                               lwork,                          ierror )
       else
          call shses(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                hor_sfctn_grid_trans,                         num_lat,        &
                             num_lon,                hor_sfctn_coef_a,        &
                    hor_sfctn_coef_b,                         num_lon,        &
                             num_lat,                          wshses,        &
                              lshses,                            work,        &
                               lwork,                          ierror )
       endif

       call math2geos(             0,                         num_lat,        &
                             num_lon,            hor_sfctn_grid_trans,        &
             hor_sfctn_grid(:, :, k),                            work )

    enddo

    ! spherepack needs spectral arrays dimensioned (num_lon, num_lat)
    bt_sfctn_coef_a = 0.0 ; bt_sfctn_coef_b = 0.0
    bt_sfctn_coef_a(1:num_fourier+1,1:num_fourier+1) = sfctn_coef_a_avg
    bt_sfctn_coef_b(1:num_fourier+1,1:num_fourier+1) = sfctn_coef_b_avg

    ! generate barotropic streamfunction
    if(is_gaussian) then
       call shsgs(           num_lat,                         num_lon,        &
                                   0,                               1,        &
                 bt_sfctn_grid_trans,                         num_lat,        &
                             num_lon,                 bt_sfctn_coef_a,        &
                     bt_sfctn_coef_b,                         num_lon,        &
                             num_lat,                          wshsgs,        &
                              lshsgs,                            work,        &
                               lwork,                          ierror )
    else
       call shses(           num_lat,                         num_lon,        &
                                   0,                               1,        &
                 bt_sfctn_grid_trans,                         num_lat,        &
                             num_lon,                 bt_sfctn_coef_a,        &
                     bt_sfctn_coef_b,                         num_lon,        &
                             num_lat,                          wshses,        &
                              lshses,                            work,        &
                               lwork,                          ierror )
    endif

    call math2geos(                0,                         num_lat,        &
                             num_lon,             bt_sfctn_grid_trans,        &
                       bt_sfctn_grid,                            work )

    ! update barotropic streamfunction coefficients
    bt_sfctn_spec = bt_sfctn_spec + sfctn_coef_a_avg**2 + sfctn_coef_b_avg**2


  end subroutine streamfunction_coef


!##############################################################################

  subroutine scalar_to_spec(                                                  &
                            one_grid,                      one_coef_a,        &
                          one_coef_b)

!                            --- input arguments ---                           

    real, intent(in), dimension(:, :, :) ::                                   &
         one_grid              ! grid field 

!                            --- output arguments ---                         

    real, intent(in), dimension(:, :, :) ::                                   &
         one_coef_a,       &   ! spectral coefficient 'a' for grid field
         one_coef_b            ! spectral coefficient 'b' for grid field

!                            --- local variables ---                         

    integer :: k               ! level counter

    real, dimension(num_lat, num_lon) ::                                      &
     one_grid_trans            ! transpose of grid field

!                            --- executable code --- 



    do k = 1, num_lev

   ! convert from geophysical to math coordinates (spherepack)
       call geo2maths(             0,                         num_lon,        &
                             num_lat,                 one_grid(:,:,k),        &
                      one_grid_trans,                            work )

       ! convert grid field to spectral space
       if(is_gaussian) then
          call shags(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                      one_grid_trans,                         num_lat,        &
                             num_lon,             one_coef_a(:, :, k),        &
                 one_coef_b(:, :, k),                         num_lon,        &
                             num_lat,                          wshags,        &
                              lshags,                            work,        &
                               lwork,                          ierror )
          if(ierror.ne.0) write(*,*) 'after shags, ierror:',ierror
             
       else
          call shaes(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                      one_grid_trans,                         num_lat,        &
                             num_lon,             one_coef_a(:, :, k),        &
                 one_coef_b(:, :, k),                         num_lon,        &
                             num_lat,                          wshaes,        &
                              lshaes,                            work,        &
                               lwork,                          ierror )
          if(ierror.ne.0) write(*,*) 'after shegs, ierror:',ierror

       endif
    end do

  end subroutine scalar_to_spec

!##############################################################################

  subroutine vector_to_spec(                                                  &
                            one_grid,                        two_grid,        &
                          one_coef_a,                      one_coef_b,        &
                          two_coef_a,                      two_coef_b )

!                            --- input arguments ---                           

    real, intent(in), dimension(:,:,:) ::                                     &
         one_grid,         &   ! grid field one
         two_grid              ! grid field two


!                            --- output arguments ---                         

    real, intent(in), dimension(:,:,:) ::                                     &
         one_coef_a,       &   ! spectral coefficient 'a' for field one
         one_coef_b,       &   ! spectral coefficient 'b' for field one
         two_coef_a,       &   ! spectral coefficient 'a' for field two
         two_coef_b            ! spectral coefficient 'b' for field two


!                            --- local variables ---                         

    integer :: k               ! level counter

    real, dimension(num_lat, num_lon) ::                                      &
         one_grid_trans,   &   ! transpose of grid field one
         two_grid_trans        ! transpose of grid field two

!                            --- executable code --- 


    do k = 1, num_lev

   ! convert from geophysical to math coordinates (spherepack)
       call geo2mathv(             0,                         num_lon,        &
                             num_lat,                 one_grid(:,:,k),        &
                     two_grid(:,:,k),                  two_grid_trans,        &
                      one_grid_trans,                            work )


       if(is_gaussian) then

          ! transform one and two to spectral space
          call vhags(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                      two_grid_trans,                  one_grid_trans,        &
                             num_lat,                         num_lon,        &
                 two_coef_a(:, :, k),             two_coef_b(:, :, k),        &
                 one_coef_a(:, :, k),             one_coef_b(:, :, k),        &
                             num_lon,                         num_lat,        &
                              wvhags,                          lvhags,        &
                                work,                           lwork,        &
                              ierror )      
          if(ierror.ne.0) write(*,*) 'after vhags, ierror:',ierror
       else
          call vhaes(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                      one_grid_trans,                  two_grid_trans,        &
                             num_lat,                         num_lon,        &
                 two_coef_a(:, :, k),             two_coef_b(:, :, k),        &
                 one_coef_a(:, :, k),             one_coef_b(:, :, k),        &
                             num_lon,                         num_lat,        &
                              wvhaes,                          lvhaes,        &
                                work,                           lwork,        &
                              ierror )
          if(ierror.ne.0) write(*,*) 'after vhaes, ierror:',ierror
       endif
    end do

  end subroutine vector_to_spec

!##############################################################################

  subroutine u_v_from_vor_div(                                                &
                          vor_coef_a,                      vor_coef_b,        &
                          div_coef_a,                      div_coef_b,        &
                              u_grid,                          v_grid )

!                            --- input arguments ---                           

    real, intent(in), dimension(:, :, :) ::                                   &
         vor_coef_a,       &   ! vorticity coefficient a
         vor_coef_b,       &   ! vorticity coefficient b
         div_coef_a,       &   ! divergence coefficient a
         div_coef_b            ! divergence coefficient b

!                            --- output arguments ---                         

    real, intent(out), dimension(:, :, :) ::                                  &
         u_grid,           &   ! gridpoint zonal wind
         v_grid                ! gridpoint meridional wind
  

!                            --- local variables ---                         

    integer :: k               ! level counter

    real, dimension(num_lat, num_lon) ::                                      &
         u_grid_trans,     &   ! transposed zonal wind
         v_grid_trans          ! transposed meridional wind

    real :: dum1, dum2         ! dummy variables required by spherepack 

!                            --- executable code --- 

    do k = 1, num_lev

       if(is_gaussian) then   
          call idvtgs(       num_lat,                         num_lon,        &
                                   0,                               1,        &
                        v_grid_trans,                    u_grid_trans,        &
                             num_lat,                         num_lon,        &
                 div_coef_a(:, :, k),             div_coef_b(:, :, k),        &
                 vor_coef_a(:, :, k),             vor_coef_b(:, :, k),        &
                             num_lon,                         num_lat,        &
                              wvhsgs,                          lvhsgs,        &
                                work,                           lwork,        &
                                dum1,                            dum2,        &
                              ierror )
          if(ierror.ne.0) write(*,*) 'after idvtgs, ierror:',ierror

       else

          call idvtes(       num_lat,                         num_lon,        &
                                   0,                               1,        &
                        v_grid_trans,                    u_grid_trans,        &
                             num_lat,                         num_lon,        &
                 div_coef_a(:, :, k),             div_coef_b(:, :, k),        &
                 vor_coef_a(:, :, k),             vor_coef_b(:, :, k),        &
                             num_lon,                         num_lat,        &
                              wvhses,                          lvhses,        &
                                work,                           lwork,        &
                                dum1,                            dum2,        &
                              ierror )

       endif

   ! convert from math coordinates to geophysical coordinates
       call math2geov(             0,                         num_lat,        &
                             num_lon,                    v_grid_trans,        &
                        u_grid_trans,                   u_grid(:,:,k),        &
                       v_grid(:,:,k),                            work )

    end do

  end subroutine u_v_from_vor_div

!##############################################################################

  subroutine vor_div_from_u_v(                                                &
                            u_coef_a,                        u_coef_b,        &
                            v_coef_a,                        v_coef_b,        &
                            vor_grid,                        div_grid )

!                            --- input arguments ---                           

    real, intent(in), dimension(:, :, :) ::                                   &
         u_coef_a,         &   ! zonal wind coefficient a
         u_coef_b,         &   ! zonal wind coefficient b
         v_coef_a,         &   ! meridional wind coefficient a
         v_coef_b              ! meridional wind coefficient b

!                            --- output arguments ---                         

    real, intent(out), dimension(:, :, :) ::                                  &
         vor_grid,         &   ! gridpoint vorticity
         div_grid              ! gridpoint divergence
  
!                            --- local variables ---                         

    integer :: k               ! level counter

    real, dimension(num_lat, num_lon) ::                                      &
         vor_grid_trans,   &   ! transposed vorticity
         div_grid_trans        ! transposed divergence

!                            --- executable code --- 

    do k = 1, num_lev

       if(is_gaussian) then   
          call vrtgs(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                      vor_grid_trans,                         num_lat,        &
                             num_lon,               u_coef_a(:, :, k),        &
                   u_coef_b(:, :, k),                         num_lon,        &
                             num_lat,                          wshsgs,        &
                              lshsgs,                            work,        &
                               lwork,                          ierror )
          if(ierror.ne.0) write(*,*) 'after vrtgs, ierror:',ierror


          ! calculate divergence
          call divgs(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                      div_grid_trans,                         num_lat,        &
                             num_lon,               v_coef_a(:, :, k),        &
                   v_coef_b(:, :, k),                         num_lon,        &
                             num_lat,                          wshsgs,        &
                              lshsgs,                            work,        &
                               lwork,                          ierror )
          if(ierror.ne.0) write(*,*) 'after divgs, ierror:',ierror

       else

          ! calculate vorticity
          call vrtes(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                      vor_grid_trans,                         num_lat,        &
                             num_lon,               u_coef_a(:, :, k),        &
                   u_coef_b(:, :, k),                         num_lon,        &
                             num_lat,                          wshses,        &
                              lshses,                            work,        &
                               lwork,                          ierror )
          if(ierror.ne.0) write(*,*) 'after vrtes, ierror:',ierror

          call dives(        num_lat,                         num_lon,        &
                                   0,                               1,        &
                      div_grid_trans,                         num_lat,        &
                             num_lon,               v_coef_a(:, :, k),        &
                   v_coef_b(:, :, k),                         num_lon,        &
                             num_lat,                          wshses,        &
                              lshses,                            work,        &
                               lwork,                          ierror ) 
          if(ierror.ne.0) write(*,*) 'after dives, ierror:',ierror
       endif

   ! convert from math coordinates to geophysical coordinates
       call math2geos(             0,                         num_lat,        &
                             num_lon,                  div_grid_trans,        &
                   div_grid(:, :, k),                            work )

       call math2geos(             0,                         num_lat,        &
                             num_lon,                  vor_grid_trans,        &
                   vor_grid(:, :, k),                            work )
    end do

  end subroutine vor_div_from_u_v

end module get_grid_fields_mod
