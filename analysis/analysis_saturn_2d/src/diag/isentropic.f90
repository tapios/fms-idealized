module isentropic_mod

use constants_and_switches_mod
use local_utilities_mod, only :  vrtcl_gradient

implicit none
private
public :: isentropic_init, isentropic_variables, sigma_to_isentropes, compute_isentropic_density

! global to this module:
real ::                                                                       &
     PotTempMin,           &   ! minimum isentropic level
     PotTempMax,           &   ! maximum isentropic level
     PotTempTrp                ! tropopause potential temperautre
     
integer ::                                                                    &
     MaxIsentrLev              ! number of isentropic levels

real, allocatable, dimension(:) ::                                            &
     theta_coord,          &   ! theta levels
     dtheta_coord,         &   ! spacing between theta levels
     scaled_theta_coord,   &   ! (theta_coord / PotTempTrp)**(-1. / kappa_dry)
     deriv_scale_theta_coord   ! 1. / grav / kappa * scaled_theta_coord / theta_coord

contains

! ##############################################################################

  subroutine isentropic_init(                                                 &
                    local_PotTempMin,                local_PotTempMax,        &
                   local_theta_coord,              local_dtheta_coord ) 

    ! set up vertical axis for isentropic coordinate system
    ! nota bene: theta_coord(1) is uppermost level with highest potential
    ! temperature, theta_coord(MaxIsentrLev) is lowermost level
    ! scaled_theta_coord is (theta_coord / PotTempTrp)**(-1. / kappa_dry), 
    ! a rescaled potential temperature coordinate.
    ! coordinate levels are equally spaced in scaled_theta_coord

!                            --- input arguments ---                          

    real,   intent(in) ::                                                     &
         local_PotTempMin, &   ! minimum isentropic level
         local_PotTempMax      ! maximum isentropic level

!                            --- output arguments ---                           
    real,    dimension(:), intent(out) ::                                     &
         local_theta_coord,&  ! theta levels (in degrees)
         local_dtheta_coord   ! spacing between theta levels (in degrees)


!                            --- local variables ---

    real    ::  scaled_theta_inc,  scaled_theta_coord_min, scaled_theta_coord_max, theta_inc
    integer :: k

!                            --- executable code ---

    PotTempMin = local_PotTempMin
    PotTempMax = local_PotTempMax
    PotTempTrp = 400.0
    MaxIsentrLev = size(local_theta_coord)

    allocate(      scaled_theta_coord(MaxIsentrLev))
    allocate( deriv_scale_theta_coord(MaxIsentrLev))
    allocate(             theta_coord(MaxIsentrLev))
    allocate(            dtheta_coord(MaxIsentrLev))


    if(.false.) then
       theta_inc = (PotTempMax - PotTempMin) / float(MaxIsentrLev - 1)
       do k = 1, MaxIsentrLev
          theta_coord(k) = PotTempMax - float(k-1) * theta_inc
       enddo
       
       do k = 1, MaxIsentrLev
          scaled_theta_coord(k) = (theta_coord(k) / PotTempTrp)**(-1. / kappa)
       enddo
       
       scaled_theta_coord_max = scaled_theta_coord(1)              ! not used
       scaled_theta_coord_min = scaled_theta_coord(MaxIsentrLev)   ! not used
    else
       scaled_theta_coord_min = (PotTempMax / PotTempTrp)**(-1. / kappa)
       scaled_theta_coord_max = (PotTempMin / PotTempTrp)**(-1. / kappa)
       
       scaled_theta_inc       = (scaled_theta_coord_max - scaled_theta_coord_min)&
            / float(MaxIsentrLev - 1)
       
       ! define coordinate levels, equally spaced in the scaled coordinate
       do k = MaxIsentrLev, 1, -1
       scaled_theta_coord(k) = scaled_theta_coord_min                         &
            + float(k-1) * scaled_theta_inc
    enddo
 endif
 


    ! coordinate levels as potential temperatures 
    local_theta_coord            = PotTempTrp * scaled_theta_coord**( -kappa)
    
    ! differenced potential temperature coordinate
    do k = MaxIsentrLev, 2, -1
       local_dtheta_coord(k) = local_theta_coord(k-1) - local_theta_coord(k)
    enddo
    local_dtheta_coord(1)    = local_dtheta_coord(2)

    ! in isentropic_variables, derivatives with respect to the
    ! scaled_theta_coord are computed. to get a derivative with
    ! respect to potential temperature, these derivatives must be
    ! multiplied by deriv_scale_theta_coord:  
    deriv_scale_theta_coord = 1. / grav /                                     &
         kappa * scaled_theta_coord / local_theta_coord

    theta_coord = local_theta_coord
    dtheta_coord = local_dtheta_coord

  end subroutine isentropic_init

! ##############################################################################

  subroutine isentropic_variables(                                            &
                             sin_lat,                   pot_temp_grid,        &
                              p_full,                            psfc,        &
                            vor_grid,                       vcos_grid,        &
                               k_sfc,                          k_asfc,        & 
                     p_isent_grid_nw,           abs_vor_isent_grid_nw,        &
                  vcos_isent_grid_nw) 

    ! interpolates variables to isentropic coordinates  
    !
    ! for sorting and interpolation, routines from the slatec library are used 

!                            --- input arguments ---                           

    real, dimension(:), intent(in) ::                                         &
         sin_lat               ! sine of latitude

    real, dimension(:,:,:), intent(in) ::                                     &
         pot_temp_grid,         &   ! grid potential temperature
         p_full,           &   ! grid pressure on full levels
         vcos_grid,        &   ! meridional wind
         vor_grid              ! grid relative vorticity

    real, dimension(:, :), intent(in) ::                                      &
         psfc                  ! grid surface pressure

!                            --- output arguments ---                         

    integer, dimension(:, :), intent(out) ::                                  &
         k_sfc,            &    ! index of surface potential temperature bin
         k_asfc                 ! index of lowest isentrope above surface

    real, dimension(:,:,:), intent(out) ::                                    &
         p_isent_grid_nw,         &   ! pressure on isentropic levels
         vcos_isent_grid_nw,      &   ! meridional wind
         abs_vor_isent_grid_nw        ! absolute vorticity

!                            --- local variables ---                           

    real, dimension(size(pot_temp_grid,1),size(pot_temp_grid,2),size(pot_temp_grid,3))  ::   &
         scaled_pot_temp,                                                     &
         abs_vor_grid

    real, dimension(size(pot_temp_grid,3)) ::                                      &
         p_column,         &
         v_column,         &
         vor_column,       &
         deriv_column,     &
         sorted_scaled_pot_temp

    integer, dimension(size(pot_temp_grid,3))   :: irank
    real, dimension(MaxIsentrLev)          :: d_p_isentr
    real, dimension(MaxIsentrLev)          :: temp_column

    integer :: ifail = 0
    integer :: i, j, k, ksfc, kasfc, num_lev
    logical :: skip
!                            --- executable code ---                           
    skip = .false.

    num_lev = size(pot_temp_grid,3)

    ! add vorticity of solid body rotation to get absolute vorticity
    do i = 1,size(pot_temp_grid,1)
       do k = 1, size(pot_temp_grid,3)
          abs_vor_grid(i, :, k)= vor_grid(i, :, k)                            &
               + 2. * omega * sin_lat
       enddo
    enddo

    ! interpolate fields columnwise to isentropic coordinates 
    ! use a scaled potential temperature in interpolation algorithms
    scaled_pot_temp = (pot_temp_grid / PotTempTrp)**(-1./kappa)

    do i = 1, size(pot_temp_grid,1)
       do j = 1, size(pot_temp_grid,2)
          ! sort scaled_pot_temp in ascending order (this changes the ordering 
          ! of indices only in regions that are statically unstable)
          sorted_scaled_pot_temp = scaled_pot_temp(i, j, :) ! vertical columns
! dpsort sorts array in ascending order and returns ranks in array irank
          call dpsort(sorted_scaled_pot_temp, num_lev, irank, 2, ifail)
          ! correspondingly rearrange pressure, meridional velocity, and
          ! vorticity in the column
          p_column      = p_full(i, j, :)
          v_column      = vcos_grid(i,j, :)
          vor_column    = abs_vor_grid(i, j, :)

! sorts based on rank vector irank
          call dpperm(p_column,      num_lev, irank, ifail)
          call dpperm(v_column,      num_lev, irank, ifail)
          call dpperm(vor_column,    num_lev, irank, ifail)

          ! fill fields in "subsurface" isentropic layers
          k          = MaxIsentrLev
          ksfc       = MaxIsentrLev
          do while(k .ge. 1 .and.                                             &
               scaled_theta_coord(k) .gt. sorted_scaled_pot_temp(num_lev))
             ! at this location, theta_coord(k) is still below the surface
             ! => fill fields with their "subsurface" values
             p_isent_grid_nw(i, j, k)            = psfc(i, j) 
             vcos_isent_grid_nw(i, j, k)         = 0.
             abs_vor_isent_grid_nw(i, j, k)      = 2. * omega * sin_lat(j)
             ksfc  = k         ! index of surface potential temperature bin
             kasfc = ksfc - 1  ! index of lowest isentropic layer above surface
             k     = k - 1
          enddo

          if (ksfc .ge. MaxIsentrLev) then
             write(*,*) sorted_scaled_pot_temp
             stop 'isentropic_circulation--Need lower isentropic coordinate level!'
          endif

          if ( abs(                                                           &
               scaled_theta_coord(ksfc-1)-sorted_scaled_pot_temp(num_lev) )   &
               .lt.                                                           &
               abs(                                                           &
               scaled_theta_coord(ksfc) - sorted_scaled_pot_temp(num_lev)) )  &
               ksfc = ksfc - 1 
          k_sfc(i, j)  = ksfc
          k_asfc(i, j) = kasfc


          ! do interpolation above the surface       
          ! interpolate using monotonicity-preserving piecewise 
          ! cubic Hermite polynomials; start with pressure and also get 
          ! vertical derivative of pressure from interpolating function

! dpchim interpolant will be monotonic on intervals where the data
! is monotonic.  At points where monotonicity switches, the derivative
! is set to zero, and there is therefore an extremum.  Uses piecewise
! cubic Hermite polynomials.

          call dpchim(num_lev, sorted_scaled_pot_temp, p_column,              &
               deriv_column, 1, ifail)         
          call dpchfe(num_lev, sorted_scaled_pot_temp, p_column,              &
               deriv_column, 1, skip, kasfc,                                  &
               scaled_theta_coord(1:kasfc), temp_column(1:kasfc),             &
               ifail)
          p_isent_grid_nw(i,j,1:kasfc) = temp_column(1:kasfc)

! vorticity
          call dpchim(num_lev, sorted_scaled_pot_temp, vor_column,            &
               deriv_column, 1, ifail)          
          call dpchfe(num_lev, sorted_scaled_pot_temp, vor_column,            &
               deriv_column, 1, skip, kasfc,                                  &
               scaled_theta_coord(1:kasfc), temp_column(1:kasfc),             &
               ifail)
          abs_vor_isent_grid_nw(i, j, 1:kasfc) = temp_column(1:kasfc)

! meridional wind
          call dpchim(num_lev, sorted_scaled_pot_temp, v_column,              &
               deriv_column, 1, ifail)          
          call dpchfe(num_lev, sorted_scaled_pot_temp, v_column,              &
               deriv_column, 1, skip, kasfc,                                  &
               scaled_theta_coord(1:kasfc), temp_column(1:kasfc),             &
               ifail)
          vcos_isent_grid_nw(i, j, 1:kasfc) = temp_column(1:kasfc)

       enddo
    enddo

  end subroutine isentropic_variables

! ##############################################################################
  
subroutine sigma_to_isentropes(                                               &
                            theta,                          p_half,           &
                     p_isent_grid,                           k_sfc,           &
                       quant_grid,                    quant_isentr,           &
                        integral )
    
    ! Interpolates gridpoint quantity 'quant_grid' to isentropic
    ! coordinates, weighting it by the isentropic density.  The
    ! interpolation is linear in a rescaled potential temperature
    ! coordinate that is proportional to 
    !         (pot_temp / 400.0)**(-1./kappa_dry).

    !                            --- input arguments ---                           

    real, dimension(:, :, :), intent(in) ::                                   &
         theta              ! grid potential temperature

    real, dimension(:, :, :), intent(in) ::                                   &
         p_half,           &   ! pressure on half levels                   
         quant_grid            ! meridional velocity in pressure coordinates   

    integer, dimension(:, :), intent(in)  ::                                  &
         k_sfc                 ! index of surface isentrope

    real, dimension(:, :, :), intent(in) ::                                   &
        p_isent_grid              ! pressure on isentropic surfaces

    logical, optional, intent(in) ::                                          &
         integral

    !                            --- output arguments ---                         

    real, dimension(:, :, :), intent(out)  ::                                 &
         quant_isentr          ! variable in isentropic coordinates

    !                            --- local variables ---  
                         
    real, dimension(size(p_isent_grid,1),size(p_isent_grid,2),size(p_isent_grid,3)) ::    &
         int_cumm_quant,   &   ! interpolated vertical integral of quantity
         scaled_quant          ! scaled "

    real, dimension(size(p_half,1), size(p_half,2), size(p_half,3)) ::        &
         cumm_quant            ! vertical integral of vshum wrt pressure

    real :: a, dp
    integer :: i, j, k, l, num_lev

    real, dimension(size(theta,1),size(theta,2),size(theta,3))  ::   &
         scaled_pot_temp                                                      

    !                            --- executable code --- 

    num_lev = size(theta,3)


    ! interpolate fields columnwise to isentropic coordinates
    ! use a scaled potential temperature in interpolation algorithms
    scaled_pot_temp = (theta / 400.0)**(-1./kappa) 

    quant_isentr = 0.0

    if (present(integral) .and. integral) then

       cumm_quant = quant_grid

    else

       ! Calculate the density-weighted vertical integral of quant_grid
       cumm_quant = 0.0
       do k = num_lev, 1, -1
          cumm_quant(:, :, k) = cumm_quant(:, :, k+1) + (quant_grid(:, :, k) *   &
               (p_half(:, :, k) -  p_half(:, :, k+1))) / grav 
       enddo

    end if

    do i = 1, size(theta,1)
       do j = 1, size(theta,2)
          
          ! Interpolate cumm_quant onto isentropic surfaces
          ! get int_cumm_quant in all isentropic levels above lowest full level
          int_cumm_quant(i, j, :) = 0.0
          l = num_lev  ! index for pressure on half levels
          do k = k_sfc(i,j), 1, -1  ! index for isentropic levels above sfc
             ! find levels that sandwich p_isentr(i, j, k)
             do while (p_half(i, j, l) .ge. p_isent_grid(i, j, k))
                l = l - 1
             enddo
             ! do interpolation             
             dp = p_half(i, j, l+1) - p_half(i, j, l)
             a  = (p_half(i, j, l+1) -  p_isent_grid(i, j, k)) / dp
             int_cumm_quant(i, j, k) = a * cumm_quant(i, j, l)                &
                  + (1. - a) * cumm_quant(i, j, l+1)
         enddo
       enddo
     enddo

     if (present(integral) .and. integral) then

        quant_isentr = int_cumm_quant

     else

        ! differentiate integral WRT scaled potential temperature
        scaled_quant(:, :, :) =  vrtcl_gradient(int_cumm_quant,                  &
             scaled_theta_coord)
        
        ! descale quantity so it's now on isentropes
        do k = MaxIsentrLev, 1, -1
           quant_isentr(:, :, k) = scaled_quant(:, :, k) *                       &
                (scaled_theta_coord(k)) / (theta_coord(k)*kappa)
        enddo
        
     end if

   end subroutine sigma_to_isentropes

! ##############################################################################
                                                                                                                  
subroutine compute_isentropic_density(p_isent_grid, dens_isent_grid)

!
! Finds the isentropic density implied by the linear interpolation routine 
! sigma_to_isentropes. 
!
                                                                                                                  
!                            --- input arguments ---
                                                                                                                  
  real, dimension(:, :, :), intent(in) ::                                       &
       p_isent_grid              ! pressure on isentropic surfaces
                                                                                                                  
!                            --- output arguments ---
                                                                                                                  
  real, dimension(:, :, :), intent(out) ::                                      &
       dens_isent_grid          ! isentropic density
                                                                                                                  
!                            --- local variables ---
                                                                                                                  
  real, dimension(size(p_isent_grid,1),size(p_isent_grid,2), size(p_isent_grid,3)) ::          &
       scaled_dens2         ! scaled density
                                                                                                                  
  integer :: k
                                                                                                                  
!                            --- executable code ---
                                                                                                                  
     dens_isent_grid = 0.0
     scaled_dens2 =  vrtcl_gradient(p_isent_grid, scaled_theta_coord)
                                                                                                                  
     ! descale density so it's now on isentropes
     do k = MaxIsentrLev, 1, -1
        dens_isent_grid(:, :, k) = scaled_dens2(:, :, k) *                          &
             (scaled_theta_coord(k)) / (theta_coord(k)*grav*kappa)
     enddo
                                                                                                                  
end subroutine compute_isentropic_density
                                                                                                                  

! ##############################################################################

end module isentropic_mod
