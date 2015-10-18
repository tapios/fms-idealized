module local_utilities_mod

use constants_and_switches_mod, only:                                         &
                                  pi,                          radius,        &
                               omega,                            grav,        &
                               rdgas,                          gspval,        &
                               rvgas,                             hlv

implicit none

public ::             mrdnl_gradient,                  vrtcl_gradient,        &
                           vrtcl_int,                  divide_by_psfc,        &
                    multiply_by_psfc,                       mrdnl_avg,        &
                           sigma_avg,        interpolate_half_to_full,        &
                       theta_sfc_pdf,                  pot_temp_deriv,        &
                  buoyancy_frequency,            update_stat_unstable,        &
               residual_circulations,               tetens_sat_mr_mod,        &
                compute_geopotential,                         pdf_avg,        &
                       simple_sat_mr

interface mrdnl_avg
   module procedure mrdnl_avg_xyz
   module procedure mrdnl_avg_xy
end interface

interface vrtcl_gradient
   module procedure vrtcl_gradient_xyz
   module procedure vrtcl_gradient_yz
end interface

interface mrdnl_gradient
   module procedure mrdnl_gradient_yz
end interface

interface interpolate_half_to_full
   module procedure interpolate_half_to_full_yz
   module procedure interpolate_half_to_full_xyz
end interface

interface vrtcl_int
   module procedure vrtcl_int_xyz
end interface

!                            --- global variables ---                          

integer ::                                                                    &
     i,                    &   ! do-loop longitude index 
     j,                    &   ! do-loop latitude index        
     k,                    &   ! do-loop level index 
     ni,                   &   ! number of longitudes (not initialized)
     nj,                   &   ! number of latitudes (not initialized)
     nk                        ! number of levels (not initialized)

integer, dimension(1) ::                                                      &
     ibin                      ! pdf bin index

contains

! #############################################################################

  function mrdnl_gradient_yz(field_yz, deg_lat )

    ! compute meridional gradient of zonal-mean field field_yz

!                            --- input arguments ---                           

    real, dimension(:), intent(in) ::                                         &
         deg_lat               ! one-dimensional latitude array (in degrees)

    real, dimension(:, :), intent(in) ::                                      &
         field_yz              ! 2-d array to be differentiated

!                            --- output arguments ---                         

    real, dimension(size(field_yz, 1), size(field_yz, 2)) ::                  &
         mrdnl_gradient_yz     ! field_yz differentiated WRT latitude

!                            --- local variables ---                           

    real, dimension(size(field_yz, 1)) :: merid_dist 
    real :: a0, b0, c0

!                            --- executable code ---                           

    merid_dist = deg_lat * pi / 180. * radius
    nj  = size(field_yz,1)
    do j=2, nj-1
       a0 = merid_dist(j+1) - merid_dist(j)
       b0 = merid_dist(j)   - merid_dist(j-1)
       c0 = merid_dist(j+1) - merid_dist(j-1)  ! a0 + b0
       ! meridional gradient by centered differences
       mrdnl_gradient_yz(j, :) = b0 / a0 / c0 * field_yz(j+1, :)              &
                   - a0 / b0 / c0 * field_yz(j-1, :)                          &
                   + (a0 - b0) / a0 / b0 * field_yz(j, :)
    enddo
    ! one-sided differences at poles
    mrdnl_gradient_yz(1, :) =                                                 &
         ( field_yz(2, :) - field_yz(1, :) )                                  &
         / ( merid_dist(2) - merid_dist(1) )

    mrdnl_gradient_yz(nj, :) =                                                &
         ( field_yz(nj, :) - field_yz(nj-1, :) )                              &
         / ( merid_dist(nj) - merid_dist(nj-1) )

  end function mrdnl_gradient_yz

!----------------------------------------------------------------------------- 

  function vrtcl_gradient_yz(field_yz, vrtcl_coord) 

    ! compute vertical gradient of 2D field field_yz

!                            --- input arguments ---                           

    real, dimension(:, :), intent(in) ::                                      &
         field_yz              ! 3-d array to be differentiated

    real, dimension(:), intent(in)                       ::                   &
         vrtcl_coord           ! vertical coordinate array

!                            --- output arguments ---                         

    real, dimension(size(field_yz, 1), size(field_yz, 2))::                   &
         vrtcl_gradient_yz     ! vertical gradient of field_yz

!                            --- local variables ---                         

    real          :: a0, b0, c0
    
!                            --- executable code ---                         

    nk = size(field_yz, 2)
    ! centered differences in interior of domain
    do k=2, nk - 1
       a0 = vrtcl_coord(k+1) - vrtcl_coord(k)
       b0 = vrtcl_coord(k)   - vrtcl_coord(k-1)
       c0 = vrtcl_coord(k+1) - vrtcl_coord(k-1)  ! a0 + b0
       ! vertical gradient by centered differences
       vrtcl_gradient_yz(:, k) = b0 / a0 / c0 * field_yz(:, k+1)              &
            - a0 / b0 / c0 * field_yz(:, k-1)                                 &
            + (a0 - b0) / a0 / b0 * field_yz(:, k)
    enddo
    ! one-sided differences at boundaries
    vrtcl_gradient_yz(:, 1) =                                                 &
         ( field_yz(:, 2) - field_yz(:, 1) )                                  &
         / ( vrtcl_coord(2) - vrtcl_coord(1) )
    vrtcl_gradient_yz(:, nk) =                                                &
         ( field_yz(:, nk) - field_yz(:, nk-1) )                              &
         / ( vrtcl_coord(nk) - vrtcl_coord(nk-1) )

  end function vrtcl_gradient_yz

!----------------------------------------------------------------------------- 

  function vrtcl_gradient_xyz(field_xyz, vrtcl_coord) 

    ! compute vertical gradient of 3D field field_xyz

!                            --- input arguments ---                           

    real, dimension(:, :, :), intent(in) ::                                   &
         field_xyz             ! 3-d input field

    real, dimension(:), intent(in) ::                                         &
         vrtcl_coord           ! vertical coordinate to differentiate against

!                            --- output arguments ---                          

    real, dimension(size(field_xyz, 1), &
                    size(field_xyz, 2), &
                    size(field_xyz, 3)) ::                                    &
         vrtcl_gradient_xyz    ! vertical gradient of field_xyz

!                            --- executable code ---                          

    do i = 1, size(field_xyz,1)
       vrtcl_gradient_xyz(i,:,:) = vrtcl_gradient_yz(field_xyz(i,:,:),        &
            vrtcl_coord)
    enddo

  end function vrtcl_gradient_xyz

!----------------------------------------------------------------------------- 

  function divide_by_psfc(field_yz, psfc)

    ! divide 2D field by mean surface pressure

!                            --- input arguments ---                           

    real, dimension(:, :), intent(in) ::                                      &
         field_yz              ! 2-d input field

    real, dimension(:), intent(in) ::                                         &
         psfc                  ! surface pressure

!                            --- output arguments ---                         

    real, dimension(size(field_yz, 1), size(field_yz, 2)) ::                  &
         divide_by_psfc        ! field_yz divided by surface pressure

!                            --- executable code ---                         

    do k=1, size(field_yz, 2)
       divide_by_psfc(:, k) = field_yz(:, k) / psfc
    enddo

  end function divide_by_psfc

!-----------------------------------------------------------------------------
  function multiply_by_psfc(field_yz, psfc)

    ! multiply 2D field by mean surface pressure

!                            --- input arguments ---                           

    real, dimension(:, :), intent(in) ::                                      &
         field_yz              ! 2-d input field
    
    real, dimension(:), intent(in) ::                                         &
         psfc                  ! surface pressure

!                            --- output arguments ---                         

    real, dimension(size(field_yz, 1), size(field_yz, 2)) ::                  &
         multiply_by_psfc      ! field_yz multiplied by surface pressure

!                            --- executable code ---                         
        
    do k=1, size(field_yz, 2)
       multiply_by_psfc(:, k) = field_yz(:, k) * psfc
    enddo

  end function multiply_by_psfc

!-----------------------------------------------------------------------------

  function mrdnl_avg_xyz(field_xyz) 

    ! computes meridional average of 3d field in grid representation
    
!                            --- input arguments ---                           

    real, intent(in), dimension(:,:,:) ::                                     &
         field_xyz             ! 3-d input field

!                            --- output arguments ---                          

    real, dimension(size(field_xyz, 2), size(field_xyz, 3)) ::                &
         mrdnl_avg_xyz    ! meridional average of field_xyz
        
!                            --- executable code ---                           

    mrdnl_avg_xyz = sum(field_xyz, 1) / float( size(field_xyz, 1) )
    
  end function mrdnl_avg_xyz

!-----------------------------------------------------------------------------

  function mrdnl_avg_xy(field_xy) 

    ! computes meridional average of 2d field in grid representation
    
!                            --- input arguments ---                           

    real, dimension(:, :), intent(in) ::                                      &
         field_xy              ! 2-d input field

!                            --- output arguments ---                         

    real, dimension(size(field_xy, 2)) ::                                     &
         mrdnl_avg_xy     ! meridional average of field_xy
        
!                            --- executable code ---                           

    mrdnl_avg_xy = sum(field_xy, 1) / float( size(field_xy, 1) )
    
  end function mrdnl_avg_xy

!-----------------------------------------------------------------------------

  function sigma_avg(field_xyz, psfc) 
    ! computes weighted zonal average of 3d field in grid representation
    ! (averages in sigma-coordinates must be weighted by surface pressure)

!                            --- input arguments ---                           

    real, dimension(:, :, :), intent(in) ::                                   &
         field_xyz             ! 3d input field

    real, dimension(:, :), intent(in) ::                                      &
         psfc(:, :)            ! surface pressure

!                            --- output arguments ---                         

    real, dimension(size(field_xyz, 2), size(field_xyz, 3)) ::                &
         sigma_avg             ! weighted meridional average of field_xyz

!                            --- local variables ---                           

    real, dimension(size(psfc, 1), size(psfc, 2), size(field_xyz, 3)) ::      &
         field_xyz_ps          ! field_xyz * surface pressure (interim var.)
        
!                            --- executable code ---                           

    ! multiply field by surface pressure
    do k=1, size(field_xyz, 3)
       field_xyz_ps(:, :, k) = field_xyz(:, :, k) * psfc
    enddo

    sigma_avg = sum(field_xyz_ps, 1) / float( size(field_xyz, 1) )
    
  end function sigma_avg

!-----------------------------------------------------------------------------

  function pdf_avg(field_xyz, psfc, field_bin) 
     ! gives surface pressure weighted increments to a histogram for the 3d field
     ! normalizes by dividing by num_lon and multiplying by number of bins -1
     ! (elsewhere divide by ps and time also)
     ! the bin width is 1.0/(num_bin-1) as there is one negative bin 

!                            --- input arguments ---                           

    real, dimension(:, :, :), intent(in) ::                                   &
         field_xyz             ! 3d input field

    real, dimension(:), intent(in) ::                                         &
         field_bin             ! histogram bin values 

    real, dimension(:, :), intent(in) ::                                      &
         psfc(:, :)            ! surface pressure

!                            --- output arguments ---                         

    real, dimension(size(field_xyz, 2), size(field_xyz, 3), size(field_bin)) ::                &
         pdf_avg             ! change in weighted histogram of field_xyz

!                            --- executable code ---                           

    pdf_avg = 0.0

    ! loop over the grid field
    do i=1, size(field_xyz,1)
      do j=1, size(field_xyz,2)
        do k=1, size(field_xyz,3)

        ! find the appropriate bin
        ! note minloc returns a rank one array by default
        ibin = minloc(abs(field_xyz(i,j,k) - field_bin))

        ! use a surface pressure weighting
        pdf_avg(j,k,ibin(1)) = pdf_avg(j,k,ibin(1)) + psfc(i,j) 

        end do
      end do
    end do

    ! normalize with the number of longitudes (divide) and number of bins -1 (multiply)
    pdf_avg = pdf_avg / float(size(field_xyz, 1)) * float(size(field_bin)-1)

    
  end function pdf_avg

!------------------------------------------------------------------------------
  function interpolate_half_to_full_yz(field_yz)

    ! interpolate 2d field from half levels to full levels
    ! (for now, simply takes averages of neigboring levels)

!                            --- input arguments ---                           

    real, intent(in), dimension(:,:) ::                                       &
         field_yz              ! 2-d input field

!                            --- output arguments ---                          

    real, dimension( size(field_yz, 1), size(field_yz, 2)-1 ) ::              &
         interpolate_half_to_full_yz

!                            --- executable code ---                           

    nk = size(field_yz, 2) - 1
    do k=1, nk
       interpolate_half_to_full_yz(:, k)                                      &
            = .5 * ( field_yz(:, k) + field_yz(:, k+1) )
    enddo

  end function interpolate_half_to_full_yz
  
!------------------------------------------------------------------------------
  
  function interpolate_half_to_full_xyz(field_xyz)

    ! interpolate 3d field from half levels to full levels
    ! (at this point, simply takes averages of neigboring levels)

!                            --- input arguments ---                           

    real, dimension(:, :, :), intent(in) ::                                   &
         field_xyz(:, :, :)    ! 3-d input field

!                            --- output arguments ---                         

    real, dimension( size(field_xyz, 1), &
                     size(field_xyz, 2), &
                     size(field_xyz, 3) - 1 ) ::                              &
         interpolate_half_to_full_xyz    

!                            --- executable code ---                           
    
    nk = size(field_xyz, 3) - 1
    do k=1, nk
       interpolate_half_to_full_xyz(:, :, k)                                  &
            = .5 * ( field_xyz(:, :, k) + field_xyz(:, :, k+1) )
    enddo

  end function interpolate_half_to_full_xyz

!------------------------------------------------------------------------------

  subroutine theta_sfc_pdf(k_sfc, sfc_theta_pdf)
    
!                            --- input arguments ---                           

    integer, dimension(:, :), intent(in) ::                                   &
         k_sfc                 ! index of surface potential temperature

!                         --- input/output arguments ---                      

    real, dimension(:, :), intent(inout) ::                                   &
         sfc_theta_pdf         ! frequency a given level is above/below surface
    
!                            --- executable code ---                      

    do i=1, size(k_sfc, 1)
       do j=1, size(k_sfc, 2)
          sfc_theta_pdf(j, k_sfc(i,j)) = sfc_theta_pdf(j, k_sfc(i, j)) + 1.
       enddo
    enddo

  end subroutine theta_sfc_pdf

!------------------------------------------------------------------------------

  subroutine pot_temp_deriv( deg_lat,                           sigma,        &
                            pot_temp,                            psfc,        &
                       d_dy_pot_temp,                 d_dsig_pot_temp,        &
                       d_dp_pot_temp )

!                            --- input arguments ---                           

    real, dimension(:), intent(in) ::                                         &
         deg_lat,          &   ! latitude in degrees
         sigma,            &   ! sigma
         psfc                  ! time averaged surface pressure
 
    real, dimension(:, :), intent(in) ::                                      &
         pot_temp              ! time averaged potential temperature

!                            --- output arguments ---                          

    real, dimension(:, :), intent(out) ::                                     &
         d_dy_pot_temp,    &   ! derivative of pot. temp. WRT latitude
         d_dsig_pot_temp,  &   ! derivative of pot. temp. WRT sigma
         d_dp_pot_temp         ! derivative of pot. temp. WRT pressure

!                            --- executable code ---                          
    
    d_dy_pot_temp   = mrdnl_gradient(pot_temp, deg_lat)
    d_dsig_pot_temp = vrtcl_gradient(pot_temp, sigma)
    d_dp_pot_temp   = divide_by_psfc(d_dsig_pot_temp, psfc)

  end subroutine pot_temp_deriv

!------------------------------------------------------------------------------

  function buoyancy_frequency(p_full, virtual_temp_grid, pot_temp, d_dp_pot_temp)  
    ! computes mean Brunt-Vaeisaelae frequency N

!                            --- input arguments ---                           

    real, dimension(:, :), intent(in)    ::                                   &
         p_full,           &   ! pressure on full levels
         virtual_temp_grid,&   ! gridpoint virtual temperature
         pot_temp,         &   ! gridpoint potential temperature
         d_dp_pot_temp         ! vertical derviative of pot. temp. WRT pressure

!                            --- output arguments ---                          

    real, dimension(size(p_full, 1), size(p_full, 2)) ::                      &
         buoyancy_frequency    ! buoyancy frequency

!                            --- local variables ---                           

    real, dimension(size(p_full, 1), size(p_full, 2)) ::                      &
         d_dz_pot_temp         ! vertical derivative of pot. temp. in meters

!                            --- executable code ---                           

    d_dz_pot_temp      = - grav / rdgas                                       &
         * p_full / virtual_temp_grid * d_dp_pot_temp 
    buoyancy_frequency = grav * d_dz_pot_temp / pot_temp
    
    
    ! so far, buoyancy = N**2; take square root where N**2 > 0 and
    ! set N = 0 where N**2 < 0
    where (buoyancy_frequency .gt. 0)
       buoyancy_frequency = sqrt(buoyancy_frequency)
    elsewhere
       buoyancy_frequency = gspval
    endwhere

  end function buoyancy_frequency

!-------------------------------------------------------------------------- 

  subroutine update_stat_unstable(                                            &
                               sigma,                        pot_temp,        &
                  freq_stat_unstable )
    
!                            --- input arguments ---                           

    real, dimension(:, :, :), intent(in) ::                                   &
         pot_temp              ! potential temperature 

    real, dimension(:), intent(in) ::                                         &
         sigma                 ! vertical coordinate 

!                         --- input/output arguments ---                      

    real, dimension(:, :), intent(inout) ::                                   &
         freq_stat_unstable    ! frequency a gridpoint is statically unstable

!                            --- local variables ---                           

    real, dimension(size(pot_temp, 1), size(pot_temp, 2), size(sigma)) ::     &
         d_dsig_pot_temp,  &   ! vertical gradient of potential temperature 
         isunstable            ! flag for statically unstable grid points
    
!                            --- executable code ---                           

    d_dsig_pot_temp = vrtcl_gradient_xyz(pot_temp, sigma)
    where ( d_dsig_pot_temp .gt. 0 )
       isunstable = 1.
    elsewhere
       isunstable = 0.
    endwhere
    
    freq_stat_unstable  = freq_stat_unstable  + mrdnl_avg(isunstable)
    
  end subroutine update_stat_unstable

!-------------------------------------------------------------------------- 

  subroutine residual_circulations(                                           &
                          sfctn_grid,                   mrdnl_eddy_hf,        &
                       vrtcl_eddy_hf,                   d_dy_pot_temp,        &
                     d_dsig_pot_temp,                     ps_grid_avg,        &
                           res_circ1,                       res_circ2 )

!                            --- input arguments ---                           

    real, dimension(:, :), intent(in) ::                                      &
         sfctn_grid,       &   ! streamfunction
         mrdnl_eddy_hf,    &   ! meridional eddy heatflux
         vrtcl_eddy_hf,    &   ! vertical eddy heatflux
         d_dy_pot_temp,    &   ! meridional derivative of pot. temp.
         d_dsig_pot_temp       ! vertical (sigma) derivative of pot. temp.

    real, dimension(:), intent(in) ::                                         &
         ps_grid_avg           ! average (time) surface pressure

!                           --- output arguments ---                           

    real, dimension(:, :), intent(out) ::                                     &
         res_circ1,        &   ! convectional TEM
         res_circ2             ! modified TEM (see below)

!                            --- local variables ---                           

    real, dimension(size(res_circ1, 1), size(res_circ1, 2)) ::                &
         res_eddy1,                                                           &
         res_eddy2
    
!                            --- executable code ---                           

    ! conventional transformed Eulerian mean
    res_eddy1 = 2. * pi * radius / grav *                                     &
         multiply_by_psfc(mrdnl_eddy_hf, ps_grid_avg) /                       &
         d_dsig_pot_temp

    res_circ1 = sfctn_grid + res_eddy1

    ! modified residual circulation: add modified eddy term whenever 
    ! this term is smaller than the conventional TEM eddy term; add the 
    ! TEM eddy term when this gives the smaller absolute value of 
    ! the streamfunction
    res_eddy2 = - 2. * pi * radius / grav *                                   &
         multiply_by_psfc(vrtcl_eddy_hf, ps_grid_avg) /                       &
         d_dy_pot_temp
    
    where ( abs(res_eddy2) .lt. abs(res_eddy1) )
       res_circ2 = sfctn_grid + res_eddy2        
    elsewhere
       res_circ2 = sfctn_grid + res_eddy1
    endwhere
 
  end subroutine residual_circulations

  !--------------------------------------------------------------------------

  function simple_sat_mr(temp_grid, p_full) result(sat_mixing_ratio)
    
! Uses a simple exponential form for saturation vapor pressure

!                            --- input arguments ---                           

    real, dimension(:,:,:), intent(in) ::                                     &
         temp_grid,        &   ! gridpoint temperature
         p_full

!                            --- output arguments ---                          

    real, dimension(size(p_full,1),size(p_full,2),size(p_full,3)) ::          &
         sat_mixing_ratio

!                            --- local variables ---                           

    real :: gc_ratio = 0.621972 ! ratio of gas constants for dry air and water vapor

    real ::                &    ! constants used in formula:
         T0 = 273.16,      &    ! K
         e0 = 610.78            ! Pa

    
    real, dimension(size(p_full,1),size(p_full,2),size(p_full,3)) ::          &
         sat_vapor_pressure

!                            --- executable code ---                           


    sat_vapor_pressure = e0 * exp( -hlv/rvgas*(1.0/temp_grid-1.0/T0))


    sat_mixing_ratio = (gc_ratio * sat_vapor_pressure) /                         &
         (p_full - sat_vapor_pressure)


  end function simple_sat_mr

  !--------------------------------------------------------------------------

  function tetens_sat_mr_mod(temp_grid, p_full) result(sat_mixing_ratio)

  ! Uses the modified Tetens formula described by Simmons et al. (1999: QJRMS, 125,
  ! 353--386), which uses saturation over ice for temperatures less than 250 K and
  ! a quadratic interpolation between saturation over ice and over liquid water for
  ! temperatures between 250 K and 273 K.

    
!                            --- input arguments ---                           

    real, dimension(:,:,:), intent(in) ::                                     &
         temp_grid,        &   ! gridpoint temperature
         p_full

!                            --- output arguments ---                          

    real, dimension(size(p_full,1),size(p_full,2),size(p_full,3)) ::          &
         sat_mixing_ratio

!                            --- local variables ---                           

    ! gas constant ratio and coefficients in Tetens saturation 
    ! vapor pressure es = es0 * exp(a3 * (T-T0)/(T-a4))
                                                                                                                                      
    real, parameter :: gc_ratio = 0.621972   ! ratio of gas constants for dry air and water vapor
                                                                                                                                  
    real, parameter :: es0 = 611.21          ! saturation vapor pressure at T0 (Pa)
    real, parameter :: T0 = 273.16           ! (K)
    real, parameter :: Ti = 250.16           ! (K)
                                                                                                                                  
    real, parameter :: a3l = 17.502          ! liquid water (Buck 1981)
    real, parameter :: a4l = 32.19           ! (K)
                                                                                                                                  
    real, parameter :: a3i = 22.587          ! ice (Alduchov and Eskridge 1996)
    real, parameter :: a4i = -0.7            ! (K)
                                                                                                                                  

    
    real, dimension(size(p_full,1),size(p_full,2),size(p_full,3)) ::          &
         sat_vapor_pressure, esl, esi, mixed_weight

!                            --- executable code ---                           


    ! saturation vapor pressure over liquid and ice
    esl = es0 * exp(a3l * (temp_grid - T0) / (temp_grid - a4l))
    esi = es0 * exp(a3i * (temp_grid - T0) / (temp_grid - a4i))


    where (temp_grid >= T0)          ! liquid
         sat_vapor_pressure = esl
    elsewhere (temp_grid <= Ti)      ! ice
         sat_vapor_pressure = esi
    elsewhere                        ! mixed
         mixed_weight = ( (temp_grid - Ti)/(T0 - Ti) )**2
         sat_vapor_pressure = (1.0-mixed_weight) * esi + mixed_weight * esl
    end where

    sat_mixing_ratio = (gc_ratio * sat_vapor_pressure) /                      &
         (p_full - sat_vapor_pressure)


  end function tetens_sat_mr_mod

!##############################################################################

  function vrtcl_int_xyz(p_half, field_xyz) 

    ! computes mass-weighted vertical average of argument 'grid'

!                            --- input arguments ---                           

    real, dimension(:, :, :), intent(in) ::                                   &
         p_half,           &  ! pressure on half levels
         field_xyz

!                            --- output arguments ---                          

    real, dimension(size(field_xyz,1), size(field_xyz,2)) ::                  &
         vrtcl_int_xyz        ! weighted average of 'field_xyz'

!                            --- local variables ---                           

    real, dimension(size(field_xyz,1),size(field_xyz,2)) ::                   &
         weight

!                            --- executable code ---                           

    vrtcl_int_xyz = 0.0
    do k = 1, size(field_xyz,3)
       weight = (p_half(:, :, k+1) - p_half(:, :, k))
       do j = 1, size(field_xyz,2)
          do i = 1, size(field_xyz,1)
             vrtcl_int_xyz(i,j) = vrtcl_int_xyz(i,j) +                        &
                  weight(i,j) * field_xyz(i,j,k)
          enddo
       enddo
    enddo

  end function vrtcl_int_xyz
!##############################################################################

  subroutine compute_geopotential(                                            &
                         surf_geopot,               virtual_temp_grid,        &
                           ln_p_half,                       ln_p_full,        &
                         geopot_full )
                           
!                            --- input arguments ---                           

    real, intent(in), dimension(:,:) ::                                       &
         surf_geopot           ! surface geopotential

    real, intent(in), dimension(:,:,:) ::                                     &
         virtual_temp_grid,&   ! virtual temperature
         ln_p_half,        &   ! log of half-level pressures
         ln_p_full             ! log of full-level pressures

!                            --- output arguments ---                          

    real, intent(out), dimension(:,:,:) ::                                    &
         geopot_full

!                            --- local variables ---                           
    real, dimension(size(ln_p_half,1), size(ln_p_half,2), size(ln_p_half,3))::&
         geopot_half           ! geopotential on half levels

    integer ::                                                                &
         num_levels,       &   ! number of levels
         k                     ! vertical counter

!                            --- executable code ---                           

    num_levels = size(virtual_temp_grid,3)

    geopot_half(:,:,num_levels+1) = surf_geopot
    geopot_half(:,:,1) = 0.0
    
    do k = num_levels, 2, -1
       geopot_half(:,:,k)=geopot_half(:,:,k+1)+rdgas*virtual_temp_grid(:,:,k)*&
            (ln_p_half(:,:,k+1) - ln_p_half(:,:,k))
    enddo

    do k = 1, num_levels
       geopot_full(:,:,k)=geopot_half(:,:,k+1)+rdgas*virtual_temp_grid(:,:,k)*&
            (ln_p_half(:,:,k+1) - ln_p_full(:,:,k))
    end do


  end subroutine compute_geopotential
       
end module local_utilities_mod
