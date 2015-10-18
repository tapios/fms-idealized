module flux_exchange_mod

  use mpp_mod,         only: mpp_npes, mpp_pe, mpp_error, mpp_set_current_pelist
  use mpp_domains_mod, only: mpp_get_compute_domain, mpp_get_compute_domains, &
                             mpp_global_sum, mpp_redistribute, operator(.EQ.)

!model_boundary_data_type contains all model fields at the boundary.
!model1_model2_boundary_type contains fields that model2 gets
!from model1, may also include fluxes. These are declared by
!flux_exchange_mod and have private components. All model fields in
!model_boundary_data_type may not be exchanged.
!will support 3 types of flux_exchange:
!REGRID: physically distinct grids, via xgrid
!REDIST: same grid, transfer in index space only
!DIRECT: same grid, same decomp, direct copy
  use atmos_model_mod, only: atmos_data_type, land_ice_atmos_boundary_type
  use ocean_model_mod, only: ocean_data_type, ice_ocean_boundary_type
  use ice_model_mod,   only: ice_data_type, land_ice_boundary_type, &
       ocean_ice_boundary_type, atmos_ice_boundary_type
  use    land_model_mod, only:  land_data_type, atmos_land_boundary_type

  use  surface_flux_mod, only: surface_flux
  use monin_obukhov_mod, only: mo_profile     

  use xgrid_mod, only: xmap_type, setup_xmap, set_frac_area, &
       put_to_xgrid, get_from_xgrid, &
       xgrid_count, some, conservation_check, xgrid_init

  use diag_integral_mod, only:     diag_integral_field_init, &
       sum_diag_integral_field

  use     utilities_mod, only: file_exist, open_file, check_nml_error,  &
       error_mesg, FATAL, get_my_pe, close_file

  use  diag_manager_mod, only: register_diag_field,  &
       register_static_field, send_data

  use  time_manager_mod, only: time_type

  use sat_vapor_pres_mod, only: escomp

  use      constants_mod, only: rdgas, rvgas, cp

  use           soil_mod, only: send_averaged_data

  implicit none
  include 'netcdf.inc'
private

  public :: flux_exchange_init,   &
     sfc_boundary_layer,   &
     generate_sfc_xgrid,   &
     flux_down_from_atmos, &
     flux_up_to_atmos,     &
     flux_land_to_ice,     &
     flux_ice_to_ocean,    &
     flux_ocean_to_ice

!-----------------------------------------------------------------------
  character(len=128) :: version = '$Id: flux_exchange.f90,v 1.8 2002/07/16 22:47:22 fms Exp $'
  character(len=128) :: tag = '$Name: havana $'
!-----------------------------------------------------------------------
!---- exchange grid maps -----

type(xmap_type) :: xmap_sfc, xmap_runoff

integer         :: n_xgrid_sfc,  n_xgrid_runoff

!-----------------------------------------------------------------------
!-------- namelist (for diagnostics) ------

character(len=4), parameter :: mod_name = 'flux'

  integer :: id_drag_moist,  id_drag_heat,  id_drag_mom,     &
     id_rough_moist, id_rough_heat, id_rough_mom,    &
     id_land_mask,   id_ice_mask,     &
     id_u_star, id_b_star, id_u_flux, id_v_flux, id_t_surf,   &
     id_t_flux, id_q_flux, id_r_flux,                         &
     id_t_atm,  id_u_atm,  id_v_atm,  id_wind,                &
     id_t_ref,  id_rh_ref, id_u_ref,  id_v_ref,               &
     id_del_h,  id_del_m,  id_del_q,                 &
     ! + slm, Mar 25, 2002 -- add diagnostics for q_surf
     id_t_ca,   id_q_surf, id_q_atm, &
     ! - slm, Mar 25, 2002
     ! + slm, Jun 02, 2002 -- add diagnostics for reference-level values
     ! (T, RH, U, V) over land part of the cells
     id_t_ref_land, id_rh_ref_land, id_u_ref_land, id_v_ref_land
     ! - slm, Jun 02, 2002

logical :: first_static = .true.
logical :: do_init = .true.

real, parameter :: bound_tol = 1e-7

real, parameter :: d622 = rdgas/rvgas
real, parameter :: d378 = 1.0-d622

!-----------------------------------------------------------------------

  real ::  z_ref_heat =  2.,  &
           z_ref_mom  = 10.

namelist /flux_exchange_nml/ z_ref_heat, z_ref_mom

! ---- allocatable module storage --------------------------------------------
real, allocatable, dimension(:) :: &
     ! NOTE: T canopy is only differet from t_surf over vegetated land
     ex_t_surf,    &   ! surface temperature for radiation calc, degK
     ex_t_ca,      &   ! near-surface (canopy) air temperature, degK
     ex_q_surf,    &   ! near-surface (canopy) air specific humidity, kg/kg
     ex_p_surf,    &   ! surface pressure

     ex_flux_t,    &   ! sens heat flux
     ex_flux_q,    &   ! water vapor flux
     ex_flux_lw,   &   ! longwave radiation flux

     ex_dhdt_surf, &   ! d(sens.heat.flux)/d(T canopy)
     ex_dedt_surf, &   ! d(water.vap.flux)/d(T canopy)
     ex_dedq_surf, &   ! d(water.vap.flux)/d(q canopy)
     ex_drdt_surf, &   ! d(LW flux)/d(T surf)
     ex_dhdt_atm,  &   ! d(sens.heat.flux)/d(T atm)
     ex_dedq_atm,  &   ! d(water.vap.flux)/d(q atm)
     ex_albedo_fix,&
     ex_old_albedo,&   ! old value of albedo for downward flux calculations
     ex_drag_q         ! q drag.coeff.


logical, allocatable, dimension(:) :: &
     ex_avail,     &   ! true where data on exchange grid are available
     ex_land           ! true if exchange grid cell is over land
real, allocatable, dimension(:) :: &
     ex_e_t_n,      &
     ex_f_t_delt_n, &
     ex_e_q_n,      &
     ex_f_q_delt_n
           

integer :: ni_atm, nj_atm ! to do atmos diagnostic from flux_ocean_to_ice
real, dimension(3) :: ccc ! for conservation checks
!Balaji, sets boundary_type%xtype
!  REGRID: grids are physically different, pass via exchange grid
!  REDIST: same physical grid, different decomposition, must move data around
!  DIRECT: same physical grid, same domain decomposition, can directly copy data
integer, parameter :: REGRID=1, REDIST=2, DIRECT=3

contains

!#######################################################################

subroutine flux_exchange_init ( Time, Atm, Land, Ice, Ocean, &
       atmos_land_boundary, atmos_ice_boundary, land_ice_atmos_boundary, &
       land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary )

  type(time_type),                   intent(in)  :: Time
  type(atmos_data_type),             intent(in)  :: Atm
  type(land_data_type),              intent(in)  :: Land
  type(ice_data_type),               intent(in)  :: Ice
  type(ocean_data_type),             intent(in)  :: Ocean
  type(atmos_land_boundary_type),    intent(out) :: atmos_land_boundary
  type(atmos_ice_boundary_type),     intent(out) :: atmos_ice_boundary
  type(land_ice_atmos_boundary_type),intent(out) :: land_ice_atmos_boundary
  type(land_ice_boundary_type),      intent(out) :: land_ice_boundary
  type(ice_ocean_boundary_type),     intent(out) :: ice_ocean_boundary
  type(ocean_ice_boundary_type),     intent(out) :: ocean_ice_boundary
  
  integer :: unit, ierr, io, iref, i
  integer :: rcode, ncid, varid, dims(4), start(4), nread(4)
  integer :: nlon, nlat
  integer :: ndims, n_ocean_areas, n_land_areas
!    real, dimension(size(Atm%glon_bnd)) :: atmlonb
!    real, dimension(size(Atm%glat_bnd)) :: atmlatb
    real, dimension(:), allocatable :: atmlonb, atmlatb
  integer :: is, ie, js, je, kd

!-----------------------------------------------------------------------
  !------ read namelist ------

  if ( file_exist('input.nml')) then
     unit = open_file ('input.nml', action='read')
     ierr=1; 
     do while (ierr /= 0)
        read  (unit, nml=flux_exchange_nml, iostat=io, end=10)
        ierr = check_nml_error(io,'flux_exchange_nml')
     enddo
10   call close_file (unit)
  endif

  !--------- write version number and namelist ------------------

  unit = open_file ('logfile.out', action='append')
  if ( get_my_pe() == 0 ) then
     write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
     write (unit, nml=flux_exchange_nml)
  endif
  call close_file (unit)

!--------- read gridspec file ------------------
!only atmos pelists needs to do it here, ocean model will do it elsewhere

    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        rcode = nf_open('INPUT/grid_spec.nc',0,ncid)
        if (rcode/=0) call error_mesg ('flux_exchange_mod', &
             'cannot open INPUT/grid_spec.nc', FATAL)

!
! check atmosphere and grid_spec.nc have same atmosphere lat/lon boundaries
!
        rcode = nf_inq_varid(ncid, 'AREA_ATM', varid)
        if (rcode/=0) call error_mesg ('flux_exchange_mod', &
             'cannot find AREA_ATM on INPUT/grid_spec.nc', &
             FATAL)
        rcode = nf_inq_vardimid(ncid, varid, dims)
        rcode = nf_inq_dimlen(ncid, dims(1), nlon)
        rcode = nf_inq_dimlen(ncid, dims(2), nlat)
        if (nlon+1/=size(Atm%glon_bnd).or.nlat+1/=size(Atm%glat_bnd)) then
            if (mpp_pe()==0) then
                print *, 'grid_spec.nc has', nlon, 'longitudes,', nlat, 'latitudes; ', &
                     'atmosphere has', size(atmlonb)-1, 'longitudes,', &
                     size(atmlatb)-1, 'latitudes (see xba.dat and yba.dat)'
            end if
            call error_mesg ('flux_exchange_mod',  &
                 'grid_spec.nc incompatible with atmosphere resolution', FATAL)
        end if
        allocate( atmlonb(size(Atm%glon_bnd)) )
        allocate( atmlatb(size(Atm%glat_bnd)) )
        rcode = nf_inq_varid(ncid, 'xba', varid)
        start = 1; nread = 1; nread(1) = nlon+1;
        rcode = nf_get_vara_double(ncid, varid, start, nread, atmlonb)
        rcode = nf_inq_varid(ncid, 'yba', varid)
        start = 1; nread = 1; nread(1) = nlat+1;
        rcode = nf_get_vara_double(ncid, varid, start, nread, atmlatb)
        if (maxval(abs(atmlonb-Atm%glon_bnd*45/atan(1.0)))>bound_tol) then
            if (mpp_pe() == 0) then
                print *, 'GRID_SPEC/ATMOS LONGITUDE INCONSISTENCY'
                do i=1,size(atmlonb)
                   print *,atmlonb(i),Atm%glon_bnd(i)*45/atan(1.0)
                end do
            end if
            call error_mesg ('flux_exchange_mod', &
                 'grid_spec.nc incompatible with atmosphere longitudes (see xba.dat and yba.dat)'&
                 ,FATAL)
        end if
        if (maxval(abs(atmlatb-Atm%glat_bnd*45/atan(1.0)))>bound_tol) then
            if (mpp_pe() == 0) then
                print *, 'GRID_SPEC/ATMOS LATITUDE INCONSISTENCY'
                do i=1,size(atmlatb)
                   print *,atmlatb(i),Atm%glat_bnd(i)*45/atan(1.0)
                end do
            end if
            call error_mesg ('flux_exchange_mod', &
                 'grid_spec.nc incompatible with atmosphere latitudes (see xba.dat and yba.dat)'&
                 , FATAL)
        end if

        call xgrid_init

        call setup_xmap(xmap_sfc, (/ 'ATM', 'OCN', 'LND' /),   &
             (/ Atm%Domain, Ice%Domain, Land%Domain /),        &
             "INPUT/grid_spec.nc"                              )
        call generate_sfc_xgrid( Land, Ice )

        call setup_xmap(xmap_runoff, (/ 'LND', 'OCN' /),       &
             (/ Land%Domain, Ice%Domain /),                    &
             "INPUT/grid_spec.nc"             )
        n_xgrid_runoff = max(xgrid_count(xmap_runoff),1)

!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!----- initialize quantities for global integral package -----

!! call diag_integral_field_init ('prec', 'f6.3')
        call diag_integral_field_init ('evap', 'f6.3')

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----
!----- all fields will be output on the atmospheric grid -----

        call diag_field_init ( Time, Atm%axes(1:2), Land%axes )
        ni_atm = size(Atm%lon_bnd)-1 ! to dimension "diag_atm"
        nj_atm = size(Atm%lat_bnd)-1 ! in flux_ocean_to_ice

!Balaji
!allocate atmos_land_boundary
! slm Mar 20 2002 -- added dedq, drag_q, and p_surf 
        call mpp_get_compute_domain( Land%domain, is, ie, js, je )
        kd = size(Land%mask,3)
        allocate( atmos_land_boundary%t_flux(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%q_flux(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%lw_flux(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%sw_flux(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%lprec(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%fprec(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%dhdt(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%dedt(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%dedq(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%drdt(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%drag_q(is:ie,js:je,kd) )
        allocate( atmos_land_boundary%p_surf(is:ie,js:je,kd) )

!allocate atmos_ice_boundary
        call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
        kd = size(Ice%ice_mask,3)
        allocate( atmos_ice_boundary%u_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%v_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%t_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%q_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%lw_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%sw_flux(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%lprec(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%fprec(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%dhdt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%dedt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%drdt(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%coszen(is:ie,js:je,kd) )
        allocate( atmos_ice_boundary%p(is:ie,js:je,kd) )
!allocate land_ice_boundary
        allocate( land_ice_boundary%runoff(is:ie,js:je) )
        allocate( land_ice_boundary%calving(is:ie,js:je) )

!allocate land_ice_atmos_boundary
        call mpp_get_compute_domain( Atm%domain, is, ie, js, je )
        allocate( land_ice_atmos_boundary%t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%albedo(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%land_frac(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_t(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dt_q(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%u_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%v_flux(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%dtaudv(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%u_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%b_star(is:ie,js:je) )
        allocate( land_ice_atmos_boundary%rough_mom(is:ie,js:je) )
    end if
    call mpp_set_current_pelist()
!ocean_ice_boundary and ice_ocean_boundary must be done on all PES
!domain boundaries will assure no space is allocated on non-relevant PEs.
    call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
!allocate ocean_ice_boundary
    allocate( ocean_ice_boundary%u(is:ie,js:je) )
    allocate( ocean_ice_boundary%v(is:ie,js:je) )
    allocate( ocean_ice_boundary%t(is:ie,js:je) )
    allocate( ocean_ice_boundary%s(is:ie,js:je) )
!frazil and sea_level are optional, if not present they should be nullified
    allocate( ocean_ice_boundary%frazil(is:ie,js:je) )
    allocate( ocean_ice_boundary%sea_level(is:ie,js:je) )
!allocate ice_ocean_boundary
    call mpp_get_compute_domain( Ocean%domain, is, ie, js, je )
!ML ocean only requires t, q, lw, sw, fprec, calving
!AMIP ocean needs no input fields
!choice of fields will eventually be done at runtime
!via field_manager
    allocate ( ice_ocean_boundary%u_flux    (is:ie,js:je) )
    allocate ( ice_ocean_boundary%v_flux    (is:ie,js:je) )
    allocate ( ice_ocean_boundary%t_flux    (is:ie,js:je) )
    allocate ( ice_ocean_boundary%q_flux    (is:ie,js:je) )
    allocate ( ice_ocean_boundary%salt_flux (is:ie,js:je) )
    allocate ( ice_ocean_boundary%lw_flux   (is:ie,js:je) )
    allocate ( ice_ocean_boundary%sw_flux   (is:ie,js:je) )
    allocate ( ice_ocean_boundary%lprec     (is:ie,js:je) )
    allocate ( ice_ocean_boundary%fprec     (is:ie,js:je) )
    allocate ( ice_ocean_boundary%runoff    (is:ie,js:je) )
    allocate ( ice_ocean_boundary%calving   (is:ie,js:je) )
    allocate ( ice_ocean_boundary%p         (is:ie,js:je) )

    ocean_ice_boundary%xtype = REDIST
    if( Ocean%domain.EQ.Ice%domain )ocean_ice_boundary%xtype = DIRECT
    ice_ocean_boundary%xtype = ocean_ice_boundary%xtype
!---- done ----
    do_init = .false.

!-----------------------------------------------------------------------

  end subroutine flux_exchange_init

!#######################################################################

subroutine sfc_boundary_layer ( dt, Time, Atm, Land, Ice, Boundary )

  real,                  intent(in)  :: dt
  type(time_type),       intent(in)  :: Time
  type(atmos_data_type), intent(in)  :: Atm
  type(land_data_type),  intent(in)  :: Land
  type(ice_data_type),   intent(in)  :: Ice
! real, dimension(:,:),   intent(out) :: t_surf_atm, albedo_atm,    &
!                                        rough_mom_atm,             &
!                                        land_frac_atm, dtaudv_atm, &
!                                        flux_u_atm, flux_v_atm,    &
!                                        u_star_atm, b_star_atm
  type(land_ice_atmos_boundary_type), intent(out) :: Boundary

  ! ---- local vars ----------------------------------------------------------
  real, dimension(n_xgrid_sfc) :: &
       ex_albedo,     &
       ex_land_frac,  &
       ex_t_atm,      & 
       ex_q_atm,      &
       ex_z_atm,      &
       ex_p_atm,      &
       ex_u_atm, ex_v_atm,    &
       ex_gust,       &
       ex_t_surf4,    &
       ex_u_surf, ex_v_surf,  &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, &
       ex_u_star,     &
       ex_b_star, ex_q_star,  &
       ex_cd_q, ex_cd_t, ex_cd_m, &
       ex_flux_u, ex_flux_v, ex_dtaudv_atm, &
       ex_ref,        &
       ex_t_ref,      &
       ex_qs_ref,     &
       ex_del_m,      &
       ex_del_h,      &
       ex_del_q,      &
       ex_wind

  real, dimension(size(Boundary%t,1),size(Boundary%t,2)) :: diag_atm
  real, dimension(size(Land%t_ca, 1),size(Land%t_ca,2), size(Land%t_ca,3)) :: diag_land
  real    :: zrefm, zrefh
  logical :: used


  ! [1] check that the module was initialized
  if (do_init) call error_mesg ('flux_exchange_mod',  &
       'must call flux_exchange_init first', FATAL)
  
  ! [2] allocate storage for variables that are also used in flux_up_to_atmos
  allocate ( &
       ex_t_surf   (n_xgrid_sfc),  &
       ex_p_surf   (n_xgrid_sfc),  &
       ex_t_ca     (n_xgrid_sfc),  &
       ex_q_surf   (n_xgrid_sfc),  &
       ex_dhdt_surf(n_xgrid_sfc),  &
       ex_dedt_surf(n_xgrid_sfc),  &
       ex_dedq_surf(n_xgrid_sfc),  &
       ex_drdt_surf(n_xgrid_sfc),  &
       ex_dhdt_atm (n_xgrid_sfc),  &
       ex_dedq_atm (n_xgrid_sfc),  &
       ex_flux_t   (n_xgrid_sfc),  &
       ex_flux_q   (n_xgrid_sfc),  &
       ex_flux_lw  (n_xgrid_sfc),  &
       ex_drag_q   (n_xgrid_sfc),  &
       ex_avail    (n_xgrid_sfc),  &
       ex_f_t_delt_n(n_xgrid_sfc), &
       ex_f_q_delt_n(n_xgrid_sfc), &
       ex_e_t_n    (n_xgrid_sfc),  &
       ex_e_q_n    (n_xgrid_sfc),  &
       ex_land     (n_xgrid_sfc)   )

  ! [3] initialize some values on exchange grid: this is actually a safeguard
  ! against using undefined values
  ex_t_surf   = 200.

! + slm Mar 28 2002 :: commented safeguards since they are never necessary
!  ex_q_surf   = 1e-3
!  ex_drag_q   = 1e-4
! - slm Mar 28 2002
  ex_u_surf   =   0.
  ex_v_surf   =   0.

  !---- do not use if relax time /= 0 ----
  ex_cd_t = 0.0
  ex_cd_m = 0.0
  ex_cd_q = 0.0
!-----------------------------------------------------------------------
!---- put atmosphere quantities onto exchange grid ----

  ! [4] put all the qantities we need onto exchange grid
  ! [4.1] put atmosphere quantities onto exchange grid
  call put_to_xgrid (Atm%t_bot , 'ATM', ex_t_atm , xmap_sfc)
  call put_to_xgrid (Atm%q_bot , 'ATM', ex_q_atm , xmap_sfc)
  call put_to_xgrid (Atm%z_bot , 'ATM', ex_z_atm , xmap_sfc)
  call put_to_xgrid (Atm%p_bot , 'ATM', ex_p_atm , xmap_sfc)
  call put_to_xgrid (Atm%u_bot , 'ATM', ex_u_atm , xmap_sfc)
  call put_to_xgrid (Atm%v_bot , 'ATM', ex_v_atm , xmap_sfc)
  call put_to_xgrid (Atm%p_surf, 'ATM', ex_p_surf, xmap_sfc)
  call put_to_xgrid (Atm%gust,   'ATM', ex_gust,   xmap_sfc)

  ! slm, Mar 20 2002: changed order in whith the data transferred from ice and land 
  ! grids, to fill t_ca first with t_surf over ocean and then with t_ca from 
  ! land, where it is different from t_surf. It is mostly to simplify 
  ! diagnostic, since surface_flux calculations distinguish between land and 
  ! not-land anyway.

  ! [4.2] put ice quantities onto exchange grid
  ! (assume that ocean quantites are stored in no ice partition)
  ! (note: ex_avail is true at ice and ocean points)
  call put_to_xgrid (Ice%t_surf,      'OCN', ex_t_surf,      xmap_sfc)
  call put_to_xgrid (Ice%rough_mom,   'OCN', ex_rough_mom,   xmap_sfc)
  call put_to_xgrid (Ice%rough_heat,  'OCN', ex_rough_heat,  xmap_sfc)
  call put_to_xgrid (Ice%rough_moist, 'OCN', ex_rough_moist, xmap_sfc)
  call put_to_xgrid (Ice%albedo,      'OCN', ex_albedo,      xmap_sfc)
  call put_to_xgrid (Ice%u_surf,      'OCN', ex_u_surf,      xmap_sfc)
  call put_to_xgrid (Ice%v_surf,      'OCN', ex_v_surf,      xmap_sfc)
  ex_t_ca = ex_t_surf ! slm, Mar 20 2002 to define values over the ocean

  ! [4.3] put land quantities onto exchange grid ----
  ex_land = some(xmap_sfc, 'LND')

  call put_to_xgrid (Land%t_surf,     'LND', ex_t_surf,      xmap_sfc)
  call put_to_xgrid (Land%t_ca,       'LND', ex_t_ca,        xmap_sfc)
  call put_to_xgrid (Land%q_ca,       'LND', ex_q_surf,      xmap_sfc)
  call put_to_xgrid (Land%rough_mom,  'LND', ex_rough_mom,   xmap_sfc)
  call put_to_xgrid (Land%rough_heat, 'LND', ex_rough_heat,  xmap_sfc)
  call put_to_xgrid (Land%rough_heat, 'LND', ex_rough_moist, xmap_sfc)
  call put_to_xgrid (Land%albedo,     'LND', ex_albedo,      xmap_sfc)

  ex_land_frac = 0.0
  call put_logical_to_real (Land%mask,    'LND', ex_land_frac, xmap_sfc)

  ! [5] compute explicit fluxes and tendencies at all available points ---
  ex_avail = some(xmap_sfc)
  call surface_flux (&
       ex_t_atm, ex_q_atm,  ex_u_atm, ex_v_atm,  ex_p_atm,  ex_z_atm,  &
       ex_p_surf,ex_t_surf, ex_t_ca,  ex_q_surf,                       &
       ex_u_surf, ex_v_surf,                                           &
       ex_rough_mom, ex_rough_heat, ex_rough_moist, ex_gust,           &
       ex_flux_t, ex_flux_q, ex_flux_lw, ex_flux_u, ex_flux_v,         &
       ex_cd_m,   ex_cd_t, ex_cd_q,                                    &
       ex_wind,   ex_u_star, ex_b_star, ex_q_star,                     &
       ex_dhdt_surf, ex_dedt_surf, ex_dedq_surf,  ex_drdt_surf,        &
       ex_dhdt_atm,  ex_dedq_atm,   ex_dtaudv_atm,                     &
! + slm Mar 28 2002 -- it is not really necessary here
!       ex_drag_q,                                                      &
! - slm Mar 28 2002
       dt,                                                             &
       ex_land,    ex_avail                                            )

! + slm Mar 28 2002 -- additional calculation to avoid passing extra
!   argument out of surface_flux (data duplication)
  where (ex_avail) ex_drag_q = ex_wind*ex_cd_q
! - slm Mar 28 2002

  ! [6] get mean quantities on atmosphere grid
  ! [6.1] compute t surf for radiation
  ex_t_surf4 = ex_t_surf ** 4

  ! [6.2] put relevant quantities onto atmospheric boundary
  call get_from_xgrid (Boundary%t,         'ATM', ex_t_surf4  ,  xmap_sfc)
  call get_from_xgrid (Boundary%albedo,    'ATM', ex_albedo   ,  xmap_sfc)
  call get_from_xgrid (Boundary%rough_mom, 'ATM', ex_rough_mom,  xmap_sfc)
  call get_from_xgrid (Boundary%land_frac, 'ATM', ex_land_frac,  xmap_sfc)

  call get_from_xgrid (Boundary%u_flux,    'ATM', ex_flux_u,     xmap_sfc)
  call get_from_xgrid (Boundary%v_flux,    'ATM', ex_flux_v,     xmap_sfc)
  call get_from_xgrid (Boundary%dtaudv,    'ATM', ex_dtaudv_atm, xmap_sfc)
  call get_from_xgrid (Boundary%u_star,    'ATM', ex_u_star    , xmap_sfc)
  call get_from_xgrid (Boundary%b_star,    'ATM', ex_b_star    , xmap_sfc)

  Boundary%t = Boundary%t ** 0.25

  ! [6.3] save atmos albedo fix and old albedo (for downward SW flux calculations)
  ! on exchange grid
  ! allocate ( ex_old_albedo(n_xgrid_sfc)  )
  ! ex_old_albedo = ex_albedo
  
  allocate ( ex_albedo_fix(n_xgrid_sfc) )
  call put_to_xgrid (Boundary%albedo, 'ATM',  ex_albedo_fix, xmap_sfc)
  ex_albedo_fix = (1.0-ex_albedo) / (1.0-ex_albedo_fix)

  !=======================================================================
  ! [7] diagnostics section

  !------- save static fields first time only ------
  if (first_static) then

     !------- land fraction ------
     if ( id_land_mask > 0 ) then
        used = send_data ( id_land_mask, Boundary%land_frac, Time )
     endif

     first_static = .false.
  endif

  !------- drag coeff moisture -----------
  if ( id_wind > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_wind, xmap_sfc)
     used = send_data ( id_wind, diag_atm, Time )
  endif
  !------- drag coeff moisture -----------
  if ( id_drag_moist > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_cd_q, xmap_sfc)
     used = send_data ( id_drag_moist, diag_atm, Time )
  endif

  !------- drag coeff heat -----------
  if ( id_drag_heat > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_cd_t, xmap_sfc)
     used = send_data ( id_drag_heat, diag_atm, Time )
  endif
  
  !------- drag coeff momemtum -----------
  if ( id_drag_mom > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_cd_m, xmap_sfc)
     used = send_data ( id_drag_mom, diag_atm, Time )
  endif
  
  !------- roughness moisture -----------
  if ( id_rough_moist > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_rough_moist, xmap_sfc)
     used = send_data ( id_rough_moist, diag_atm, Time )
  endif
  
  !------- roughness heat -----------
  if ( id_rough_heat > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_rough_heat, xmap_sfc)
     used = send_data ( id_rough_heat, diag_atm, Time )
  endif
  
  !------- roughness momemtum -----------
  if ( id_rough_mom > 0 ) then
     used = send_data ( id_rough_mom, Boundary%rough_mom, Time )
  endif
  
  !------- friction velocity -----------
  if ( id_u_star > 0 ) then
     used = send_data ( id_u_star, Boundary%u_star, Time )
  endif
  
  !------- bouyancy -----------
  if ( id_b_star > 0 ) then
     used = send_data ( id_b_star, Boundary%b_star, Time )
  endif

  !-----------------------------------------------------------------------
  !------ diagnostics for fields at bottom atmospheric level ------
  
  if ( id_t_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_atm, xmap_sfc)
     used = send_data ( id_t_atm, diag_atm, Time )
  endif
  
  if ( id_u_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_u_atm, xmap_sfc)
     used = send_data ( id_u_atm, diag_atm, Time )
  endif
  
  if ( id_v_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_v_atm, xmap_sfc)
     used = send_data ( id_v_atm, diag_atm, Time )
  endif
  
  ! + slm, Mar 25, 2002
  if ( id_q_atm > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_q_atm, xmap_sfc)
     used = send_data ( id_q_atm, diag_atm, Time )
  endif
  ! - slm, Mar 25, 2002

  !-----------------------------------------------------------------------
  !--------- diagnostics for fields at reference level ---------
  
  if ( id_t_ref > 0 .or. id_rh_ref > 0 .or. &
       id_u_ref > 0 .or. id_v_ref  > 0 .or. &
       id_t_ref_land > 0 .or. id_rh_ref_land > 0 .or. &
       id_u_ref_land > 0 .or. id_v_ref_land  > 0 ) then
     
     zrefm = z_ref_mom
     zrefh = z_ref_heat
     !      ---- optimize calculation ----
     if ( id_t_ref <= 0 ) zrefh = zrefm
     
     call mo_profile ( zrefm, zrefh, ex_z_atm,   ex_rough_mom, &
          ex_rough_heat, ex_rough_moist,          &
          ex_u_star, ex_b_star, ex_q_star,        &
          ex_del_m, ex_del_h, ex_del_q, ex_avail  )

     !    ------- reference relative humidity -----------
     if ( id_rh_ref > 0 .or. id_rh_ref_land > 0 ) then
        ex_ref   = ex_q_surf + (ex_q_atm-ex_q_surf) * ex_del_q
        where (ex_avail) &
           ex_t_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        call escomp (ex_t_ref, ex_qs_ref)
        ex_qs_ref = d622*ex_qs_ref/(ex_p_surf-d378*ex_qs_ref)
        ex_ref    = MIN(100.,100.*ex_ref/ex_qs_ref)

        if ( id_rh_ref_land > 0 ) then
           call get_from_xgrid (diag_land,'LND', ex_ref, xmap_sfc)
           used = send_averaged_data ( id_rh_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if(id_rh_ref > 0) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_rh_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference temp -----------
     if ( id_t_ref > 0 .or. id_t_ref_land > 0 ) then
        ex_ref = ex_t_ca + (ex_t_atm-ex_t_ca) * ex_del_h
        if (id_t_ref_land > 0) then
           call get_from_xgrid (diag_land, 'LND', ex_ref, xmap_sfc)
           used = send_averaged_data ( id_t_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_t_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_t_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference u comp -----------
     if ( id_u_ref > 0 .or. id_u_ref_land > 0) then
        ex_ref = ex_u_surf + (ex_u_atm-ex_u_surf) * ex_del_m
        if ( id_u_ref_land > 0 ) then
           call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
           used = send_averaged_data ( id_u_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_u_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_u_ref, diag_atm, Time )
        endif
     endif

     !    ------- reference v comp -----------
     if ( id_v_ref > 0 .or. id_v_ref_land > 0 ) then
        ex_ref = ex_v_surf + (ex_v_atm-ex_v_surf) * ex_del_m
        if ( id_v_ref_land > 0 ) then
           call get_from_xgrid ( diag_land, 'LND', ex_ref, xmap_sfc )
           used = send_averaged_data ( id_v_ref_land, diag_land, &
                Land%tile_size, Time, mask = Land%mask )
        endif
        if ( id_v_ref > 0 ) then
           call get_from_xgrid (diag_atm, 'ATM', ex_ref, xmap_sfc)
           used = send_data ( id_v_ref, diag_atm, Time )
        endif
     endif

     !    ------- interp factor for heat ------
     if ( id_del_h > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_del_h, xmap_sfc)
        used = send_data ( id_del_h, diag_atm, Time )
     endif

     !    ------- interp factor for momentum ------
     if ( id_del_m > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_del_m, xmap_sfc)
        used = send_data ( id_del_m, diag_atm, Time )
     endif

     !    ------- interp factor for moisture ------
     if ( id_del_q > 0 ) then
        call get_from_xgrid (diag_atm, 'ATM', ex_del_q, xmap_sfc)
        used = send_data ( id_del_q, diag_atm, Time )
     endif

  endif

!=======================================================================

end subroutine sfc_boundary_layer

!#######################################################################

subroutine flux_down_from_atmos (Time, Atm, Land, Ice, &
     Atmos_boundary, Land_boundary, Ice_boundary )

  type(time_type),       intent(in) :: Time
  type(atmos_data_type), intent(in) :: Atm
  type(land_data_type),  intent(in) :: Land
  type(ice_data_type),   intent(in) :: Ice
! real, dimension(:,:),   intent(in)  :: flux_u_atm, flux_v_atm
  type(land_ice_atmos_boundary_type),intent(in) :: Atmos_boundary
! real, dimension(:,:,:), intent(out) ::                               &
!                                    flux_t_land, flux_q_land,         &
!                                    flux_lw_land, flux_sw_land,       &
!                                    dhdt_land, dedt_land, drdt_land,  &
!                                    lprec_land, fprec_land,           &
!                                    flux_t_ice, flux_q_ice,           &
!                                    flux_lw_ice, flux_sw_ice,         &
!                                    dhdt_ice, dedt_ice, drdt_ice,     &
!                                    lprec_ice , fprec_ice,            &
!                                    flux_u_ice, flux_v_ice, coszen_ice
  type(atmos_land_boundary_type),    intent(out):: Land_boundary
  type(atmos_ice_boundary_type),     intent(out):: Ice_boundary

  real, dimension(n_xgrid_sfc) :: ex_flux_sw, ex_flux_lwd, &
       ex_lprec, ex_fprec,      &
       ex_flux_u, ex_flux_v,    &
       ex_coszen, ex_ice_frac

  real, dimension(n_xgrid_sfc) :: ex_gamma  , ex_dtmass,  &
       ex_delta_t, ex_delta_q, &
       ex_dflux_t, ex_dflux_q

  real :: ice_frac (size(Ice_boundary%dhdt,1),size(Ice_boundary%dhdt,2),size(Ice_boundary%dhdt,3))
  real :: diag_atm (size(Atmos_boundary%u_flux,1),size(Atmos_boundary%u_flux,2))

  real    :: cp_inv
  logical :: used

!-----------------------------------------------------------------------
!---- put atmosphere quantities onto exchange grid ----

  call put_to_xgrid (Atm%flux_sw, 'ATM', ex_flux_sw, xmap_sfc)
  call put_to_xgrid (Atm%flux_lw, 'ATM', ex_flux_lwd, xmap_sfc)

  !  ccc = conservation_check(Atm%lprec, 'ATM', xmap_sfc)
  !  if (mpp_pe()==0) print *,'LPREC', ccc

  call put_to_xgrid (Atm%lprec,   'ATM', ex_lprec, xmap_sfc)
  call put_to_xgrid (Atm%fprec,   'ATM', ex_fprec, xmap_sfc)

  call put_to_xgrid (Atm%coszen,  'ATM', ex_coszen, xmap_sfc)

  call put_to_xgrid (Atmos_boundary%u_flux, 'ATM', ex_flux_u, xmap_sfc)
  call put_to_xgrid (Atmos_boundary%v_flux, 'ATM', ex_flux_v, xmap_sfc)

!-----------------------------------------------------------------------
!---- adjust sw flux for albedo variations on exch grid ----

  ex_flux_sw = ex_flux_sw * ex_albedo_fix
  deallocate ( ex_albedo_fix )

!----- compute net longwave flux (down-up) -----
  ! (note: lw up already in ex_flux_lw)

  ex_flux_lw = ex_flux_lwd - ex_flux_lw

!-----------------------------------------------------------------------
!----- adjust fluxes for implicit dependence on atmosphere ----


  call put_to_xgrid (Atm%Surf_Diff%dtmass , 'ATM', ex_dtmass , xmap_sfc)
  call put_to_xgrid (Atm%Surf_Diff%delta_t, 'ATM', ex_delta_t, xmap_sfc)
  call put_to_xgrid (Atm%Surf_Diff%delta_q, 'ATM', ex_delta_q, xmap_sfc)
  call put_to_xgrid (Atm%Surf_Diff%dflux_t, 'ATM', ex_dflux_t, xmap_sfc)
  call put_to_xgrid (Atm%Surf_Diff%dflux_q, 'ATM', ex_dflux_q, xmap_sfc)

  cp_inv = 1.0/cp

  where(ex_avail)

     ! temperature

     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_t + ex_dhdt_atm*cp_inv))
     ex_e_t_n      =  ex_dtmass*ex_dhdt_surf*cp_inv*ex_gamma
     ex_f_t_delt_n = (ex_delta_t + ex_dtmass * ex_flux_t*cp_inv) * ex_gamma    
     
     ex_flux_t     =  ex_flux_t        + ex_dhdt_atm * ex_f_t_delt_n 
     ex_dhdt_surf  =  ex_dhdt_surf     + ex_dhdt_atm * ex_e_t_n   

     ! moisture
     ex_gamma      =  1./ (1.0 - ex_dtmass*(ex_dflux_q + ex_dedq_atm))

! here it looks like two derivatives with different units are added together,
! but in fact they are not: ex_dedt_surf and ex_dedq_surf defined in complimentary
! regions of exchange grid, so that if one of them is not zero the other is, and
! vice versa.
     ex_e_q_n      =  ex_dtmass*(ex_dedt_surf+ex_dedq_surf) * ex_gamma

     ex_f_q_delt_n = (ex_delta_q  + ex_dtmass * ex_flux_q) * ex_gamma    
     
     ex_flux_q     =  ex_flux_q    + ex_dedq_atm * ex_f_q_delt_n 
     ex_dedt_surf  =  ex_dedt_surf + ex_dedq_atm * ex_e_q_n
     ex_dedq_surf  =  ex_dedq_surf + ex_dedq_atm * ex_e_q_n
  endwhere

!-----------------------------------------------------------------------
!---- output fields on the land grid -------

  call get_from_xgrid (Land_boundary%t_flux,  'LND', ex_flux_t,    xmap_sfc)
  call get_from_xgrid (Land_boundary%q_flux,  'LND', ex_flux_q,    xmap_sfc)
  call get_from_xgrid (Land_boundary%sw_flux, 'LND', ex_flux_sw,   xmap_sfc)
  call get_from_xgrid (Land_boundary%lw_flux, 'LND', ex_flux_lw,   xmap_sfc)
  call get_from_xgrid (Land_boundary%dhdt,    'LND', ex_dhdt_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%dedt,    'LND', ex_dedt_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%dedq,    'LND', ex_dedq_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%drdt,    'LND', ex_drdt_surf, xmap_sfc)
  call get_from_xgrid (Land_boundary%lprec,   'LND', ex_lprec,     xmap_sfc)
  call get_from_xgrid (Land_boundary%fprec,   'LND', ex_fprec,     xmap_sfc)
  call get_from_xgrid (Land_boundary%drag_q,  'LND', ex_drag_q,    xmap_sfc)
  call get_from_xgrid (Land_boundary%p_surf,  'LND', ex_p_surf,    xmap_sfc)

!-----------------------------------------------------------------------
!---- output fields on the ice grid -------

  call get_from_xgrid (Ice_boundary%t_flux,   'OCN', ex_flux_t,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%q_flux,   'OCN', ex_flux_q,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%sw_flux,  'OCN', ex_flux_sw,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%lw_flux,  'OCN', ex_flux_lw,   xmap_sfc)
  call get_from_xgrid (Ice_boundary%dhdt,     'OCN', ex_dhdt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%dedt,     'OCN', ex_dedt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%drdt,     'OCN', ex_drdt_surf, xmap_sfc)
  call get_from_xgrid (Ice_boundary%lprec,    'OCN', ex_lprec,     xmap_sfc)
  call get_from_xgrid (Ice_boundary%fprec,    'OCN', ex_fprec,     xmap_sfc)
  call get_from_xgrid (Ice_boundary%u_flux,   'OCN', ex_flux_u,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%v_flux,   'OCN', ex_flux_v,    xmap_sfc)
  call get_from_xgrid (Ice_boundary%coszen,   'OCN', ex_coszen,    xmap_sfc)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------

  !------- zonal wind stress -----------
  if ( id_u_flux > 0 ) then
     used = send_data ( id_u_flux, Atmos_boundary%u_flux, Time )
  endif

  !------- meridional wind stress -----------
  if ( id_v_flux > 0 ) then
     used = send_data ( id_v_flux, Atmos_boundary%v_flux, Time )
  endif

!=======================================================================

  end subroutine flux_down_from_atmos

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! flux_land_to_ice - translate runoff from land to ice grids                   !
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
subroutine flux_land_to_ice(Land, Ice, Boundary )
  type(land_data_type),          intent(in) :: Land
  type(ice_data_type),           intent(in) :: Ice
!real, dimension(:,:), intent(out) :: runoff_ice, calving_ice
  type(land_ice_boundary_type),  intent(out):: Boundary
  
  real, dimension(n_xgrid_runoff) :: ex_runoff, ex_calving
  real, dimension(size(Boundary%runoff,1),size(Boundary%runoff,2),1) :: ice_buf

  ! ccc = conservation_check(Land%discharge, 'LND', xmap_runoff)
  ! if (mpp_pe()==0) print *,'RUNOFF', ccc

  call put_to_xgrid ( Land%discharge,      'LND', ex_runoff,  xmap_runoff)
  call put_to_xgrid ( Land%discharge_snow, 'LND', ex_calving, xmap_runoff)
  call get_from_xgrid (ice_buf, 'OCN', ex_runoff,  xmap_runoff)
  Boundary%runoff = ice_buf(:,:,1);
  call get_from_xgrid (ice_buf, 'OCN', ex_calving, xmap_runoff)
  Boundary%calving = ice_buf(:,:,1);

end subroutine flux_land_to_ice

!#######################################################################

subroutine flux_ice_to_ocean ( Ice, Ocean, Boundary )

  type(ice_data_type),   intent(in)  :: Ice
  type(ocean_data_type), intent(in)  :: Ocean
!  real, dimension(:,:),   intent(out) :: flux_u_ocean,  flux_v_ocean,  &
!                                         flux_t_ocean,  flux_q_ocean,  &
!                                         flux_sw_ocean, flux_lw_ocean, &
!                                         lprec_ocean,   fprec_ocean,   &
!                                         runoff_ocean,  calving_ocean, &
!                                         flux_salt_ocean, p_surf_ocean
  type(ice_ocean_boundary_type), intent(inout) :: Boundary

  select case (Boundary%xtype)
  case(DIRECT)
     !same grid and domain decomp for ocean and ice    
     if( ASSOCIATED(Boundary%u_flux   ) )Boundary%u_flux    = Ice%flux_u
     if( ASSOCIATED(Boundary%v_flux   ) )Boundary%v_flux    = Ice%flux_v
     if( ASSOCIATED(Boundary%t_flux   ) )Boundary%t_flux    = Ice%flux_t
     if( ASSOCIATED(Boundary%q_flux   ) )Boundary%q_flux    = Ice%flux_q
     if( ASSOCIATED(Boundary%salt_flux) )Boundary%salt_flux = Ice%flux_salt
     if( ASSOCIATED(Boundary%sw_flux  ) )Boundary%sw_flux   = Ice%flux_sw
     if( ASSOCIATED(Boundary%lw_flux  ) )Boundary%lw_flux   = Ice%flux_lw
     if( ASSOCIATED(Boundary%lprec    ) )Boundary%lprec     = Ice%lprec
     if( ASSOCIATED(Boundary%fprec    ) )Boundary%fprec     = Ice%fprec
     if( ASSOCIATED(Boundary%runoff   ) )Boundary%runoff    = Ice%runoff
     if( ASSOCIATED(Boundary%calving  ) )Boundary%calving   = Ice%calving
     if( ASSOCIATED(Boundary%p   ) )Boundary%p    = Ice%p_surf
  case(REDIST)
     !same grid, different domain decomp for ocean and ice    
     if (ASSOCIATED(Boundary%u_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_u, Ocean%Domain, Boundary%u_flux)
     if (ASSOCIATED(Boundary%v_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_v, Ocean%Domain, Boundary%v_flux)
     if (ASSOCIATED(Boundary%t_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_t, Ocean%Domain, Boundary%t_flux)
     if (ASSOCIATED(Boundary%q_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_q, Ocean%Domain, Boundary%q_flux)
     if (ASSOCIATED(Boundary%salt_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_salt, Ocean%Domain, Boundary%salt_flux)
     if (ASSOCIATED(Boundary%sw_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_sw, Ocean%Domain, Boundary%sw_flux)
     if (ASSOCIATED(Boundary%lw_flux)) &
          call mpp_redistribute(Ice%Domain, Ice%flux_lw, Ocean%Domain, Boundary%lw_flux)
     if (ASSOCIATED(Boundary%lprec)) &
          call mpp_redistribute(Ice%Domain, Ice%lprec, Ocean%Domain, Boundary%lprec)
     if (ASSOCIATED(Boundary%fprec)) &
          call mpp_redistribute(Ice%Domain, Ice%fprec, Ocean%Domain, Boundary%fprec)
     if (ASSOCIATED(Boundary%runoff)) &
          call mpp_redistribute(Ice%Domain, Ice%runoff, Ocean%Domain, Boundary%runoff)
     if (ASSOCIATED(Boundary%calving)) &
          call mpp_redistribute(Ice%Domain, Ice%calving, Ocean%Domain, Boundary%calving)
     if (ASSOCIATED(Boundary%p)) &
          call mpp_redistribute(Ice%Domain, Ice%p_surf, Ocean%Domain, Boundary%p)
  case DEFAULT
     call mpp_error( FATAL, 'FLUX_ICE_TO_OCEAN: Boundary%xtype must be DIRECT or REDIST.' )
  end select

!-----------------------------------------------------------------------

  end subroutine flux_ice_to_ocean

!#######################################################################

subroutine flux_ocean_to_ice ( Ocean, Ice, Boundary )

  type(ocean_data_type), intent(in)  :: Ocean
  type(ice_data_type),   intent(in)  :: Ice
!  real, dimension(:,:),   intent(out) :: t_surf_ice, u_surf_ice, v_surf_ice, &
!                                         frazil_ice, s_surf_ice, sea_lev_ice
  type(ocean_ice_boundary_type), intent(inout) :: Boundary
  real, dimension(size(Boundary%t,1),size(Boundary%t,2),size(Ice%part_size,3)) &
       :: ice_frac
  real, dimension(:), allocatable :: ex_ice_frac
  real, dimension(ni_atm, nj_atm) :: diag_atm
  logical :: used
  select case (Boundary%xtype)
  case(DIRECT)
     !same grid and domain decomp for ocean and ice    
     if( ASSOCIATED(Boundary%u) )Boundary%u = Ocean%u_surf
     if( ASSOCIATED(Boundary%v) )Boundary%v = Ocean%v_surf
     if( ASSOCIATED(Boundary%t) )Boundary%t = Ocean%t_surf
     if( ASSOCIATED(Boundary%s) )Boundary%s = Ocean%s_surf
     if( ASSOCIATED(Boundary%frazil) )Boundary%frazil = Ocean%frazil
     if( ASSOCIATED(Boundary%sea_level) )Boundary%sea_level = Ocean%sea_lev
  case(REDIST)
     !same grid, different domain decomp for ocean and ice    
     if( ASSOCIATED(Boundary%u) )call mpp_redistribute(Ocean%Domain, Ocean%u_surf, Ice%Domain, Boundary%u)
     if( ASSOCIATED(Boundary%v) )call mpp_redistribute(Ocean%Domain, Ocean%v_surf, Ice%Domain, Boundary%v)
     if( ASSOCIATED(Boundary%t) )call mpp_redistribute(Ocean%Domain, Ocean%t_surf, Ice%Domain, Boundary%t)
     if( ASSOCIATED(Boundary%s) )call mpp_redistribute(Ocean%Domain, Ocean%s_surf, Ice%Domain, Boundary%s)
     if( ASSOCIATED(Boundary%frazil) )call mpp_redistribute(Ocean%Domain, Ocean%frazil, Ice%Domain, Boundary%frazil)
     if( ASSOCIATED(Boundary%sea_level) )call mpp_redistribute(Ocean%Domain, Ocean%sea_lev, Ice%Domain, Boundary%sea_level)
  case DEFAULT
     call mpp_error( FATAL, 'FLUX_OCEAN_TO_ICE: Boundary%xtype must be DIRECT or REDIST.' )
  end select
  if ( id_ice_mask > 0 ) then
     allocate ( ex_ice_frac(n_xgrid_sfc) )
     ice_frac        = 1.
     ice_frac(:,:,1) = 0.
     ex_ice_frac     = 0.
     call put_to_xgrid (ice_frac, 'OCN', ex_ice_frac, xmap_sfc)
     call get_from_xgrid (diag_atm, 'ATM', ex_ice_frac, xmap_sfc)
     used = send_data ( id_ice_mask, diag_atm, Ice%Time )
     deallocate ( ex_ice_frac )
  endif

!-----------------------------------------------------------------------

  end subroutine flux_ocean_to_ice

subroutine generate_sfc_xgrid( Land, Ice )
! subroutine to regenerate exchange grid eliminating side 2 tiles with 0 frac area
    type(land_data_type), intent(in) :: Land
  type(ice_data_type), intent(in)           :: Ice

  call set_frac_area (Ice%part_size , 'OCN', xmap_sfc)
  call set_frac_area (Land%tile_size, 'LND', xmap_sfc)
  n_xgrid_sfc = max(xgrid_count(xmap_sfc),1)
  return
end subroutine generate_sfc_xgrid

!#######################################################################

subroutine flux_up_to_atmos ( Time, Land, Ice, Boundary )

  type(time_type),      intent(in)  :: Time
  type(land_data_type), intent(in)  :: Land
  type(ice_data_type),  intent(in)  :: Ice
! real, dimension(:,:),   intent(out) :: dt_t_atm, dt_q_atm
  type(land_ice_atmos_boundary_type), intent(out) :: Boundary
  real, dimension(n_xgrid_sfc) :: ex_t_surf_new, ex_dt_t_surf,  &
       ex_dt_t, ex_dt_q, &
       ex_delta_t_n,  ex_delta_q_n, &
       ! + slm, Mar 20 2002
       ex_t_ca_new,   ex_dt_t_ca,   &
       ex_q_surf_new, ex_dt_q_surf
       ! - slm, Mar 20 2002

  real, dimension(size(Boundary%dt_t,1),size(Boundary%dt_t,2)) :: diag_atm, &
       evap_atm
  logical :: used
  !-----------------------------------------------------------------------
  !----- compute surface temperature change -----

  ex_t_surf_new = 200.0

  call put_to_xgrid (Ice%t_surf,  'OCN', ex_t_surf_new, xmap_sfc)
  ex_t_ca_new = ex_t_surf_new  ! since it is the same thing over oceans
  call put_to_xgrid (Land%t_ca,   'LND', ex_t_ca_new,   xmap_sfc)
  call put_to_xgrid (Land%t_surf, 'LND', ex_t_surf_new, xmap_sfc)

  call escomp(ex_t_surf_new, ex_q_surf_new)
  ex_q_surf_new  = d622*ex_q_surf_new/(ex_p_surf-d378*ex_q_surf_new) 
  call put_to_xgrid (Land%q_ca, 'LND', ex_q_surf_new, xmap_sfc)

  where (ex_avail)
     ex_dt_t_ca   = ex_t_ca_new   - ex_t_ca   ! changes in near-surface T
     ex_dt_t_surf = ex_t_surf_new - ex_t_surf ! changes in radiative T
     ex_dt_q_surf = ex_q_surf_new - ex_q_surf ! changes in near-surface q
  endwhere

  !-----------------------------------------------------------------------
  !-----  adjust fluxes and atmospheric increments for 
  !-----  implicit dependence on surface temperature -----

  ex_delta_t_n = 0.0
  ex_delta_q_n = 0.0

  where(ex_avail)
     ex_flux_t     = ex_flux_t  + ex_dt_t_ca   * ex_dhdt_surf
     ex_flux_lw    = ex_flux_lw - ex_dt_t_surf * ex_drdt_surf
     ex_delta_t_n  = ex_f_t_delt_n  + ex_dt_t_ca*ex_e_t_n

     ! Note that in the expressions below ex_e_q_n used to relate changes
     ! of humidity on the lower atmos level to changes of sfc temperature or 
     ! near-sfc humidity; it may seem that the units are mixed up. However
     ! it is not the case since for each exchange grid cell ex_e_q_n is
     ! defined either in terms of sfc temperature or near-sfc humidity
     
     where (ex_land) 
        ex_delta_q_n  = ex_f_q_delt_n + ex_dt_q_surf * ex_e_q_n
        ex_flux_q     = ex_flux_q     + ex_dt_q_surf * ex_dedq_surf
     elsewhere
        ! note that in this region (over ocean) ex_dt_t_surf == ex_dt_t_ca
        ex_delta_q_n  = ex_f_q_delt_n + ex_dt_t_surf * ex_e_q_n
        ex_flux_q     = ex_flux_q     + ex_dt_t_surf * ex_dedt_surf
     endwhere
  endwhere

  !-----------------------------------------------------------------------
  !---- get mean quantites on atmospheric grid ----

  call get_from_xgrid (Boundary%dt_t, 'ATM', ex_delta_t_n, xmap_sfc)
  call get_from_xgrid (Boundary%dt_q, 'ATM', ex_delta_q_n, xmap_sfc)

  !  ---- always get evaporation for diagnostic purposes ----

  call get_from_xgrid (evap_atm,      'ATM', ex_flux_q,    xmap_sfc)

  !=======================================================================
  !-------------------- diagnostics section ------------------------------

  !------- new surface temperature -----------
  if ( id_t_surf > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_surf_new, xmap_sfc)
     used = send_data ( id_t_surf, diag_atm, Time )
  endif


  ! + slm, Mar 27 2002
  ! ------ new canopy temperature --------
  !   NOTE, that in the particular case of LM2 t_ca is identical to t_surf,
  !   but this will be changed in future version of the land madel
  if ( id_t_ca > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_t_ca_new, xmap_sfc)
     used = send_data ( id_t_ca, diag_atm, Time )
  endif

  ! ---- new surface humidity ---------
  if ( id_q_surf > 0 ) then
     call get_from_xgrid(diag_atm, 'ATM', ex_q_surf_new, xmap_sfc )
     used = send_data ( id_q_surf, diag_atm, Time )
  endif
  ! - slm, Mar 27 2002

  !------- sensible heat flux -----------
  if ( id_t_flux > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_flux_t, xmap_sfc)
     used = send_data ( id_t_flux, diag_atm, Time )
  endif

!------- latent heat flux (see below) -----------

!  if ( id_q_flux > 0 ) then
!!!!  latent  = hlv
!!!!  where (ice_or_snow) latent = hlf + hlv
!     call get_from_xgrid (diag_atm, 'ATM', ex_flux_q, xmap_sfc)
!     used = send_data ( id_q_flux, diag_atm, Time )
!  endif

  !------- net longwave flux -----------
  if ( id_r_flux > 0 ) then
     call get_from_xgrid (diag_atm, 'ATM', ex_flux_lw, xmap_sfc)
     used = send_data ( id_r_flux, diag_atm, Time )
  endif

  !------- evaporation rate -----------
  
  if ( id_q_flux > 0 ) then
     used = send_data ( id_q_flux, evap_atm, Time )
  endif

  !-----------------------------------------------------------------------
  !---- accumulate global integral of evaporation (mm/day) -----

  call sum_diag_integral_field ('evap', evap_atm*86400.)

  !=======================================================================
  !---- deallocate module storage ----
  deallocate ( &
       ex_t_surf   ,  &
       ex_p_surf   ,  &
       ex_t_ca     ,  &
       ex_q_surf   ,  &
       ex_dhdt_surf,  &
       ex_dedt_surf,  &
       ex_dedq_surf,  &
       ex_drdt_surf,  &
       ex_dhdt_atm ,  &
       ex_dedq_atm ,  &
       ex_flux_t   ,  &
       ex_flux_q   ,  &
       ex_flux_lw  ,  &
       ex_drag_q   ,  &
       ex_avail    ,  &
       ex_f_t_delt_n, &
       ex_f_q_delt_n, &
       ex_e_t_n    ,  &
       ex_e_q_n    ,  &
       ex_land        )


!-----------------------------------------------------------------------

end subroutine flux_up_to_atmos

!#######################################################################

subroutine put_logical_to_real (mask, id, ex_mask, xmap)

  logical         , intent(in)    :: mask(:,:,:)
  character(len=3), intent(in)    :: id
  real            , intent(inout) :: ex_mask(:)
  type(xmap_type), intent(inout) :: xmap

  !-----------------------------------------------------------------------
  !    puts land or ice model masks (with partitions) onto the
  !    exchange grid as a real array (1.=true, 0.=false)
  !-----------------------------------------------------------------------

  real, dimension(size(mask,1),size(mask,2),size(mask,3)) :: rmask
  
  where (mask)
     rmask = 1.0
  elsewhere
     rmask = 0.0
  endwhere

  call put_to_xgrid(rmask, id, ex_mask, xmap)

end subroutine put_logical_to_real

!#######################################################################

subroutine diag_field_init ( Time, atmos_axes, land_axes )

  type(time_type), intent(in) :: Time
  integer,         intent(in) :: atmos_axes(2)
  integer,         intent(in) :: land_axes(2)

  integer :: iref
  character(len=6) :: label_zm, label_zh
  real, dimension(2) :: trange = (/  100., 400. /), &
       vrange = (/ -400., 400. /), &
       frange = (/ -0.01, 1.01 /)
!-----------------------------------------------------------------------
!  initializes diagnostic fields that may be output from this module
!  (the id numbers may be referenced anywhere in this module)
!-----------------------------------------------------------------------

  !------ labels for diagnostics -------
  !  (z_ref_mom, z_ref_heat are namelist variables)

  iref = int(z_ref_mom+0.5)
  if ( real(iref) == z_ref_mom ) then
     write (label_zm,105) iref
     if (iref < 10) write (label_zm,100) iref
  else
     write (label_zm,110) z_ref_mom
  endif

  iref = int(z_ref_heat+0.5)
  if ( real(iref) == z_ref_heat ) then
     write (label_zh,105) iref
     if (iref < 10) write (label_zh,100) iref
  else
     write (label_zh,110) z_ref_heat
  endif

100 format (i1,' m',3x)
105 format (i2,' m',2x)
110 format (f4.1,' m')

  !--------- initialize static diagnostic fields --------------------

  id_land_mask = &
       register_static_field ( mod_name, 'land_mask', atmos_axes,  &
       'fractional amount of land', 'none', &
       range=frange )
  
  !--------- initialize diagnostic fields --------------------

  id_ice_mask = &
       register_diag_field ( mod_name, 'ice_mask', atmos_axes, Time, &
       'fractional amount of sea ice', 'none',  &
       range=frange )
  
  id_wind = &
       register_diag_field ( mod_name, 'wind', atmos_axes, Time, &
       'wind speed for flux calculations', 'm/s', &
       range=(/0.,vrange(2)/) )
  
  id_drag_moist = &
       register_diag_field ( mod_name, 'drag_moist', atmos_axes, Time, &
       'drag coeff for moisture',    'none'     )
  
  id_drag_heat  = &
       register_diag_field ( mod_name, 'drag_heat', atmos_axes, Time, &
       'drag coeff for heat',    'none'     )
  
  id_drag_mom   = &
       register_diag_field ( mod_name, 'drag_mom',  atmos_axes, Time, &
       'drag coeff for momentum',     'none'     )
  
  id_rough_moist = &
       register_diag_field ( mod_name, 'rough_moist', atmos_axes, Time, &
       'surface roughness for moisture',  'm'  )

  id_rough_heat = &
       register_diag_field ( mod_name, 'rough_heat', atmos_axes, Time, &
       'surface roughness for heat',  'm'  )

  id_rough_mom  = &
       register_diag_field ( mod_name, 'rough_mom',  atmos_axes, Time, &
       'surface roughness for momentum',  'm'  )

  id_u_star     = &
       register_diag_field ( mod_name, 'u_star',     atmos_axes, Time, &
       'friction velocity',   'm/s'   )

  id_b_star     = &
       register_diag_field ( mod_name, 'b_star',     atmos_axes, Time, &
       'buoyancy scale',      'm/s2'   )

  id_u_flux     = &
       register_diag_field ( mod_name, 'tau_x',      atmos_axes, Time, &
       'zonal wind stress',     'pa'   )

  id_v_flux     = &
       register_diag_field ( mod_name, 'tau_y',      atmos_axes, Time, &
       'meridional wind stress',     'pa'   )

  id_t_surf     = &
       register_diag_field ( mod_name, 't_surf',     atmos_axes, Time, &
       'surface temperature',    'deg_k', &
       range=trange    )

  ! + slm, Mar 25, 2002 -- add diagnositcs for t_ca, q_ca, and q_atm
  id_t_ca       = &
       register_diag_field ( mod_name, 't_ca',     atmos_axes, Time, &
       'canopy air temperature',    'deg_k', &
       range=trange    )

  id_q_atm     = &
       register_diag_field ( mod_name, 'q_atm',     atmos_axes, Time, &
       'specific humidity at btm level',    'kg/kg', &
       range=trange    )

  id_q_surf     = &
       register_diag_field ( mod_name, 'q_surf',     atmos_axes, Time, &
       'surface specific humidity',    'kg/kg', &
       range=trange    )
  ! - slm, Mar 25, 2002

  id_t_flux     = &
       register_diag_field ( mod_name, 'shflx',      atmos_axes, Time, &
       'sensible heat flux',     'w/m2'    )

  id_q_flux     = &
       register_diag_field ( mod_name, 'evap',       atmos_axes, Time, &
       'evaporation rate',        'kg/m2/s'  )

  id_r_flux     = &
       register_diag_field ( mod_name, 'lwflx',      atmos_axes, Time, &
       'net (down-up) longwave flux',   'w/m2'    )

  id_t_atm      = &
       register_diag_field ( mod_name, 't_atm',      atmos_axes, Time, &
       'temperature at btm level',    'deg_k', &
       range=trange     )

  id_u_atm      = &
       register_diag_field ( mod_name, 'u_atm',      atmos_axes, Time, &
       'u wind component at btm level',  'm/s', &
       range=vrange    )

  id_v_atm      = &
       register_diag_field ( mod_name, 'v_atm',      atmos_axes, Time, &
       'v wind component at btm level',  'm/s', &
       range=vrange    )

  id_t_ref      = &
       register_diag_field ( mod_name, 't_ref',      atmos_axes, Time, &
       'temperature at '//label_zh, 'deg_k' , &
       range=trange      )

  id_rh_ref     = &
       register_diag_field ( mod_name, 'rh_ref',     atmos_axes, Time,   &
       'relative humidity at '//label_zh, 'percent' )

  id_u_ref      = &
       register_diag_field ( mod_name, 'u_ref',      atmos_axes, Time, &
       'zonal wind component at '//label_zm,  'm/s', &
       range=vrange )

  id_v_ref      = &
       register_diag_field ( mod_name, 'v_ref',      atmos_axes, Time,     &
       'meridional wind component at '//label_zm, 'm/s', &
       range=vrange )

  id_del_h      = &
       register_diag_field ( mod_name, 'del_h',      atmos_axes, Time,  &
       'ref height interp factor for heat', 'none' )
  id_del_m      = &
       register_diag_field ( mod_name, 'del_m',      atmos_axes, Time,     &
       'ref height interp factor for momentum','none' )
  id_del_q      = &
       register_diag_field ( mod_name, 'del_q',      atmos_axes, Time,     &
       'ref height interp factor for moisture','none' )

  ! + slm Jun 02, 2002 -- diagnostics of reference values over the land
  id_t_ref_land = &
       register_diag_field ( mod_name, 't_ref_land', Land_axes, Time, &
       'temperature at '//label_zh//' over land', 'deg_k' , &
       range=trange, missing_value =  -100.0)
  id_rh_ref_land= &
       register_diag_field ( mod_name, 'rh_ref_land', Land_axes, Time,   &
       'relative humidity at '//label_zh//' over land', 'percent',       &
       missing_value=-999.0)
  id_u_ref_land = &
       register_diag_field ( mod_name, 'u_ref_land',  Land_axes, Time, &
       'zonal wind component at '//label_zm//' over land',  'm/s', &
       range=vrange, missing_value=-999.0 )
  id_v_ref_land = &
       register_diag_field ( mod_name, 'v_ref_land',  Land_axes, Time,     &
       'meridional wind component at '//label_zm//' over land', 'm/s', &
       range=vrange, missing_value = -999.0 )
  ! - slm Jun 02, 2002
!-----------------------------------------------------------------------

  end subroutine diag_field_init

function in_box(i, j,is, ie, js, je)
  integer :: i, j, is, ie, js, je
  logical :: in_box

  in_box = (i>=is) .and. (i<=ie) .and. (j>=js) .and. (j<=je)
end function in_box

end module flux_exchange_mod

