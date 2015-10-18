!
!  coupler_main couples component models and controls the time integration
!
program coupler_main
! utilities used by coupler
!
  use time_manager_mod, only: time_type, set_calendar_type, set_time,  &
                              set_date, get_date, days_in_month, month_name,  &
                              operator(+), operator(-), operator (<), &
                              operator (>), operator ( /= ), operator ( / ), &
                              operator (*), thirty_day_months, julian, &
                              NOLEAP, no_calendar

  use  utilities_mod, only: open_file, file_exist, check_nml_error,  &
                            error_mesg, FATAL, WARNING,              &
                            print_version_number, get_my_pe,         &
                            utilities_init, utilities_end,           &
                            close_file, check_system_clock

  use  diag_manager_mod, only: diag_manager_init, diag_manager_end, get_base_date
!
! model interfaces used to couple the component models:
!               atmosphere, land, ice, and ocean
!
  use  atmos_model_mod, only: atmos_model_init, atmos_model_end, &
                              update_atmos_model_down,           &
                              update_atmos_model_up,             &
                              atmos_data_type, &
                              land_ice_atmos_boundary_type
  use   land_model_mod, only: land_model_init, land_model_end, &
                              land_data_type, atmos_land_boundary_type, &
                              update_land_model_fast, update_land_model_slow

  use    ice_model_mod, only: ice_model_init, ice_model_end,  &
                              update_ice_model_slow_up,          &
                              update_ice_model_fast,          &
                              update_ice_model_slow_dn,          &
                              ice_data_type, land_ice_boundary_type, &
                              ocean_ice_boundary_type, atmos_ice_boundary_type

  use  ocean_model_mod, only: update_ocean_model, ocean_model_init,  &
                              ocean_model_end, ocean_data_type, ice_ocean_boundary_type
!
! flux_ calls translate information between model grids - see flux_exchange.f90
!
  use flux_exchange_mod, only: flux_exchange_init,   &
                               sfc_boundary_layer,   &
                               generate_sfc_xgrid,   &
                               flux_down_from_atmos, &
                               flux_up_to_atmos,     &
                               flux_land_to_ice,     &
                               flux_ice_to_ocean,    &
                               flux_ocean_to_ice
  use mpp_mod, only: mpp_clock_id, mpp_clock_begin, mpp_clock_end, MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
  use mpp_mod, only: mpp_pe, mpp_npes, mpp_root_pe, stderr, mpp_set_current_pelist, mpp_declare_pelist, mpp_error
  use mpp_domains_mod, only: mpp_broadcast_domain

  implicit none

!-----------------------------------------------------------------------

  character(len=128) :: version = '$Id: coupler_main.f90,v 1.6 2002/07/16 22:47:16 fms Exp $'
  character(len=128) :: tag = '$Name: havana $'

!-----------------------------------------------------------------------
!---- model defined-types ----

  type (atmos_data_type) :: Atm
  type  (land_data_type) :: Land
  type   (ice_data_type) :: Ice
  type (ocean_data_type) :: Ocean

  type(atmos_land_boundary_type)     :: Atmos_land_boundary
  type(atmos_ice_boundary_type)      :: Atmos_ice_boundary
  type(land_ice_atmos_boundary_type) :: Land_ice_atmos_boundary
  type(land_ice_boundary_type)       :: Land_ice_boundary
  type(ice_ocean_boundary_type)      :: Ice_ocean_boundary
  type(ocean_ice_boundary_type)      :: Ocean_ice_boundary

!-----------------------------------------------------------------------
! ----- coupled model time -----

  type (time_type) :: Time, Time_init, Time_end, Time_step_ocean, &
                      Time_step_atmos, Time_step_cpld
  type(time_type) :: Time_atmos, Time_ocean
  integer :: num_ocean_calls, num_atmos_calls, no, na
  integer :: num_cpld_calls, nc

! ----- coupled model initial date -----

  logical :: ocean_seg_start
  logical :: ocean_seg_end
  integer :: date_init(6)
  integer :: calendar_type = -99

!-----------------------------------------------------------------------

  integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
  character(len=17) :: calendar = '                 '
  logical :: override = .false.  ! override restart values for date
  integer :: months=0, days=0, hours=0, minutes=0, seconds=0
  integer :: dt_atmos = 0  ! fluxes passed between atmosphere & ice/land
  integer :: dt_ocean = 0  ! ocean tracer timestep
  integer :: dt_cpld  = 0  ! fluxes passed between ice & ocean
  integer,dimension (3)           :: locmax, locmin

  integer :: atmos_pe_start=0, atmos_pe_end=0, ocean_pe_start=0, ocean_pe_end=0
  namelist /coupler_nml/ current_date, calendar, override,       &
       months, days, hours, minutes, seconds,  &
       dt_cpld, dt_atmos, dt_ocean, &
       atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end
  integer :: atmclock, ocnclock, iceclock, landclock, fluxclock, icefluxclock, landicefluxclock, fluxclockdn, fluxclockup

!#######################################################################

  call utilities_init

  call coupler_init 

  atmclock     = mpp_clock_id( 'Atmosphere',      flags=MPP_CLOCK_SYNC )
  ocnclock     = mpp_clock_id( 'Ocean',           flags=MPP_CLOCK_SYNC )
  iceclock     = mpp_clock_id( 'Ice',             flags=MPP_CLOCK_SYNC )
  landclock    = mpp_clock_id( 'Land',            flags=MPP_CLOCK_SYNC )
  fluxclock    = mpp_clock_id( 'Surface BL',      flags=MPP_CLOCK_SYNC )
  icefluxclock = mpp_clock_id( 'Ice ocean flux',        flags=MPP_CLOCK_SYNC )
  landicefluxclock = mpp_clock_id( 'Ice land flux',        flags=MPP_CLOCK_SYNC )
  fluxclockdn  = mpp_clock_id( 'Atmos flux down', flags=MPP_CLOCK_SYNC )
  fluxclockup  = mpp_clock_id( 'Atmos flux up',   flags=MPP_CLOCK_SYNC )

  call check_system_clock ('END OF INITIALIZATION')
!-----------------------------------------------------------------------
!------ ocean/slow-ice integration loop -------

  do nc = 1, num_cpld_calls

     if( Atm%pe )then
         call mpp_set_current_pelist(Atm%pelist)
         call generate_sfc_xgrid( Land, Ice )
     end if
     call mpp_set_current_pelist()
     call mpp_clock_begin(icefluxclock)
     call flux_ocean_to_ice( Ocean, Ice, Ocean_ice_boundary )
     call mpp_clock_end(icefluxclock)

!     call mpp_clock_begin(icefluxclock)
!     call flux_ice_to_ocean( Ice, Ocean, Ice_ocean_boundary )
!     call mpp_clock_end(icefluxclock)

     if( Atm%pe )then
         call mpp_set_current_pelist(Atm%pelist)
         call mpp_clock_begin(iceclock)
         call update_ice_model_slow_up( Ocean_ice_boundary, Ice )
         call mpp_clock_end(iceclock)

!-----------------------------------------------------------------------
!   ------ atmos/fast-land/fast-ice integration loop -------


         do na = 1, num_atmos_calls

!            if( mpp_pe().EQ.mpp_root_pe() )write( stderr(),'(a,2i3)' )'nc,nc=', nc,na
            Time_atmos = Time_atmos + Time_step_atmos

            call mpp_clock_begin(fluxclock)
            call sfc_boundary_layer( REAL(dt_atmos), Time_atmos, Atm, Land, Ice, &
                 Land_ice_atmos_boundary )
            call mpp_clock_end(fluxclock)

!      ---- atmosphere down ----

            call mpp_clock_begin(atmclock)
            call update_atmos_model_down( Land_ice_atmos_boundary, Atm )
            call mpp_clock_end(atmclock)

            call mpp_clock_begin(fluxclockdn)
            call flux_down_from_atmos( Time_atmos, Atm, Land, Ice, &
                 Land_ice_atmos_boundary, &
                 Atmos_land_boundary, &
                 Atmos_ice_boundary )
            call mpp_clock_end(fluxclockdn)


!      --------------------------------------------------------------

!      ---- land model ----

            call mpp_clock_begin(landclock)
            call update_land_model_fast( Atmos_land_boundary, Land )
            call mpp_clock_end(landclock)

!      ---- ice model ----


            call mpp_clock_begin(iceclock)
            call update_ice_model_fast( Atmos_ice_boundary, Ice )
            call mpp_clock_end(iceclock)

!      --------------------------------------------------------------
!      ---- atmosphere up ----

            call mpp_clock_begin(fluxclockup)
            call flux_up_to_atmos( Time_atmos, Land, Ice, Land_ice_atmos_boundary )
            call mpp_clock_end(fluxclockup)

            call mpp_clock_begin(atmclock)
            call update_atmos_model_up( Land_ice_atmos_boundary, Atm )
            call mpp_clock_end(atmclock)

!--------------

         enddo

!   ------ end of atmospheric time step loop -----
         call mpp_clock_begin(landclock)
         call update_land_model_slow(Land)
         call mpp_clock_end(landclock)
!-----------------------------------------------------------------------

!
!     need flux call to put runoff and p_surf on ice grid
!
         call mpp_clock_begin(landicefluxclock)
         call flux_land_to_ice( Land, Ice, Land_ice_boundary )

         Atmos_ice_boundary%p = 0.0 ! call flux_atmos_to_ice_slow ?
         call mpp_clock_end(landicefluxclock)

!   ------ slow-ice model ------

         call mpp_clock_begin(iceclock)
         call update_ice_model_slow_dn( Atmos_ice_boundary, Land_ice_boundary, Ice )
         call mpp_clock_end(iceclock)
         Time = Time_atmos
     end if                     !Atm%pe block
     call mpp_set_current_pelist()

     call mpp_clock_begin(icefluxclock)
     call flux_ice_to_ocean( Ice, Ocean, Ice_ocean_boundary )
     call mpp_clock_end(icefluxclock)

     if( Ocean%pe )then
         call mpp_set_current_pelist(Ocean%pelist)
         do no = 1,num_ocean_calls
            Time_ocean = Time_ocean + Time_step_ocean

            ocean_seg_start = ( no .eq. 1 )               ! could eliminate these by
            ocean_seg_end   = ( no .eq. num_ocean_calls ) ! putting this loop in
                                                          ! update_ocean_model since
                                                          ! fluxes don't change here

            call mpp_clock_begin(ocnclock)
            call update_ocean_model( Ice_ocean_boundary, Ocean, &
                 ocean_seg_start, ocean_seg_end, num_ocean_calls)
            call mpp_clock_end(ocnclock)

         enddo
!   ------ end of ocean time step loop -----
!-----------------------------------------------------------------------
         Time = Time_ocean
     end if
!--------------
  enddo

!-----------------------------------------------------------------------
  call check_system_clock ('END OF TIME LOOP')

  call coupler_end

  call diag_manager_end (Time)
  call utilities_end

!-----------------------------------------------------------------------

contains

!#######################################################################

  subroutine coupler_init

!-----------------------------------------------------------------------
!   initialize all defined exchange grids and all boundary maps
!-----------------------------------------------------------------------
    integer :: unit, log_unit, ierr, io, id, jd, kd, m, i
    integer :: date(6)
    type (time_type) :: Run_length
    character(len=9) :: month
    logical :: use_namelist
    integer :: pe, npes
!-----------------------------------------------------------------------
!----- read namelist -------

    unit = open_file ('input.nml', action='read')
    ierr=1; do while (ierr /= 0)
       read  (unit, nml=coupler_nml, iostat=io, end=10)
       ierr = check_nml_error (io, 'coupler_nml')
    enddo
10  call close_file (unit)

!----- write namelist to logfile (close log_unit later) -----

    log_unit = open_file ('logfile.out', action='append')
    if ( get_my_pe() == 0 ) then
        write (log_unit,'(/,80("="),/(a))') trim(version), trim(tag)
        write (log_unit, nml=coupler_nml)
    endif

!----- read restart file -----

    if (file_exist('INPUT/coupler.res')) then
        unit = open_file ('INPUT/coupler.res',  &
             form='native', action='read')
        read  (unit) date
        read  (unit) calendar_type
        call close_file (unit)
        use_namelist = .false.
    else
        use_namelist = .true.
    endif

!----- use namelist value (either no restart or override flag on) ---

    if ( use_namelist .or. override ) then

!----- override date with namelist values ------

        if ( sum(current_date) <= 0 ) then
            call error_mesg ('program coupler',  &
                 'no namelist value for base_date or current_date', FATAL)
        else
            date      = current_date
        endif

!----- override calendar type with namelist value -----

        if (calendar(1:6) == 'julian') then
            calendar_type = julian
        else if (calendar(1:6) == 'NOLEAP') then
            calendar_type = NOLEAP
        else if (calendar(1:10) == 'thirty_day') then
            calendar_type = thirty_day_months
        else if (calendar(1:11) == 'no_calendar') then
            calendar_type = no_calendar
        else if (calendar(1:1) /= ' ') then
            call error_mesg ('program coupler',  &
                 'invalid namelist value for calendar', FATAL)
        else
            call error_mesg ('program coupler',  &
                 'no namelist value for calendar', FATAL)
        endif

    endif

    call set_calendar_type (calendar_type)

!----- write current/initial date actually used to logfile file -----

    if ( get_my_pe() == 0 ) then
        write (log_unit,16) date(1),trim(month_name(date(2))),date(3:6)
    endif
    call close_file (log_unit)

16  format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt') 

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

    call diag_manager_init

!----- always override initial/base date with diag_manager value -----

    call get_base_date ( date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6)  )

!----- use current date if no base date ------

    if ( date_init(1) == 0 ) date_init = date

!----- set initial and current time types ------

    Time_init = set_date (date_init(1), date_init(2), date_init(3), &
         date_init(4), date_init(5), date_init(6))

    Time      = set_date (date(1), date(2), date(3),  &
         date(4), date(5), date(6))

!----- compute the ending time -----

    Time_end = Time
    do m=1,months
       Time_end = Time_end + set_time(0,days_in_month(Time_end))
    end do
    Time_end   = Time_end + set_time(hours*3600+minutes*60+seconds, days)
    Run_length = Time_end - Time

!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

    unit = open_file ('time_stamp.out', action='write')

    month = month_name(date(2))
    if ( get_my_pe() == 0 ) write (unit,20) date, month(1:3)

    call get_date (Time_end, date(1), date(2), date(3),  &
         date(4), date(5), date(6))
    month = month_name(date(2))
    if ( get_my_pe() == 0 ) write (unit,20) date, month(1:3)

    call close_file (unit)

20  format (6i4,2x,a3)

!-----------------------------------------------------------------------
!----- compute the time steps ------

    Time_step_cpld  = set_time (dt_cpld ,0)
    Time_step_ocean = set_time (dt_ocean,0)
    Time_step_atmos = set_time (dt_atmos,0)

!----- determine maximum number of iterations per loop ------

    num_cpld_calls  = Run_length      / Time_step_cpld
    num_ocean_calls = Time_step_cpld  / Time_step_ocean
    num_atmos_calls = Time_step_cpld  / Time_step_atmos

!-----------------------------------------------------------------------
!------------------- some error checks ---------------------------------

!----- initial time cannot be greater than current time -------

    if ( Time_init > Time ) call error_mesg ('program coupler',  &
         'initial time is greater than current time', FATAL)

!----- make sure run length is a multiple of ocean time step ------

    if ( num_cpld_calls * num_ocean_calls * Time_step_ocean /= Run_length )  &
         call error_mesg ('program coupler',  &
         'run length must be multiple of ocean time step', FATAL)

! ---- make sure cpld time step is a multiple of atmos time step ----

    if ( num_atmos_calls * Time_step_atmos /= Time_step_cpld )  &
         call error_mesg ('program coupler',   &
         'atmos time step is not a multiple of the cpld time step', 2)

! ---- make sure cpld time step is a multiple of ocean time step ----

    if ( num_ocean_calls * Time_step_ocean /= Time_step_cpld )  &
         call error_mesg ('program coupler',   &
         'cpld time step is not a multiple of the ocean time step', 2)
!-----------------------------------------------------------------------
!------ initialize component models ------
!------ grid info now comes from grid_spec file

!pe information
    pe = mpp_pe()
    npes = mpp_npes()
    if( atmos_pe_end.EQ.0 )atmos_pe_end = npes-1
    if( ocean_pe_end.EQ.0 )ocean_pe_end = npes-1
!    if( pe.EQ.mpp_root_pe() )write( stderr(),* ) &
!         'atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end=', &
!          atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end
    Atm%pe = atmos_pe_start.LE.pe .AND. pe.LE.atmos_pe_end
    Ocean%pe = ocean_pe_start.LE.pe .AND. pe.LE.ocean_pe_end
    allocate( Atm%pelist(atmos_pe_end-atmos_pe_start+1) )
    allocate( Ocean%pelist(ocean_pe_end-ocean_pe_start+1) )
    Atm%pelist = (/(i,i=atmos_pe_start,atmos_pe_end)/)
    Ocean%pelist = (/(i,i=ocean_pe_start,ocean_pe_end)/)
    call mpp_declare_pelist( Atm%pelist )
    call mpp_declare_pelist( Ocean%pelist )
    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
!---- atmosphere ----
        call atmos_model_init( Atm, Time_init, Time, Time_step_atmos )

!---- land ----------
        call land_model_init( Land, Time_init, Time, Time_step_atmos, Time_step_cpld )

!---- ice -----------
        call ice_model_init( Ice, Time_init, Time, Time_step_atmos, Time_step_cpld )
    end if
    if( Ocean%pe )then
        call mpp_set_current_pelist(Ocean%pelist)
!---- ocean ---------
        call ocean_model_init( Ocean, Time_init, Time, Time_step_ocean )
    end if
    call mpp_set_current_pelist()
    call mpp_broadcast_domain(Ice%domain)
    call mpp_broadcast_domain(Ocean%domain)
!-----------------------------------------------------------------------
!---- initialize flux exchange module ----

    call flux_exchange_init ( Time, Atm, Land, Ice, Ocean, &
         atmos_land_boundary, atmos_ice_boundary, land_ice_atmos_boundary, &
         land_ice_boundary, ice_ocean_boundary, ocean_ice_boundary )

    Time_atmos = Time
    Time_ocean = Time
!-----------------------------------------------------------------------
!---- open and close output restart to make sure directory is there ----

    unit = open_file ('RESTART/coupler.res',  &
         form='native', action='write')
    call close_file (unit, status='delete')

!-----------------------------------------------------------------------

  end subroutine coupler_init

!#######################################################################

  subroutine coupler_end

    integer :: unit, date(6)
!-----------------------------------------------------------------------

    if( Atm%pe )then
        call mpp_set_current_pelist(Atm%pelist)
        call atmos_model_end (Atm)
        call  land_model_end (Land)
        call   ice_model_end (Ice)
    end if
    if( Ocean%pe )then
        call mpp_set_current_pelist(Ocean%pelist)
        call ocean_model_end (Ocean)
    end if
    call mpp_set_current_pelist()

!----- compute current date ------

    call get_date (Time, date(1), date(2), date(3),  &
         date(4), date(5), date(6))

!----- check time versus expected ending time ----

    if (Time /= Time_end) call error_mesg ('program coupler',  &
         'final time does not match expected ending time', WARNING)

!----- write restart file ------

    if ( get_my_pe() /= 0 ) return

    unit = open_file ('RESTART/coupler.res',  &
         form='native', action='write')
    write (unit) date
    write (unit) calendar_type
    call close_file (unit)

!-----------------------------------------------------------------------

  end subroutine coupler_end

!#######################################################################

end program coupler_main

