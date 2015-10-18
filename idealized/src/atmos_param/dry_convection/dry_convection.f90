
module dry_convection_mod
  
  use           fms_mod, only: error_mesg, FATAL, stdlog, mpp_pe,             &
                               mpp_root_pe, open_namelist_file, close_file,   &
                               check_nml_error
  use     constants_mod, only: rdgas, cp_air
  use  time_manager_mod, only: time_type
  use  diag_manager_mod, only: register_diag_field, send_data

!-----------------------------------------------------------------------------
  
  implicit none
  private

!                             ---  namelist ---
  real :: tau, gamma

  namelist /dry_convection_nml/ tau, gamma

 !                        ---  global variables ---
 
  integer :: i, j, jit, k, num_levels
  real, parameter :: cons1 = rdgas/cp_air

  integer :: id_cape, id_cin, id_lzb, id_lcl, id_tp, id_n_tp, &
       id_dp, id_amb, id_dt

  character(len=14), parameter :: mod_name='dry_convection'
  real :: missing_value = -1.e-10
  
  public :: dry_convection_init, dry_convection

  contains

!-----------------------------------------------------------------------------

    subroutine dry_convection_init(axes, Time)

!                         ---  input arguments ---

      integer, intent(in) :: axes(4)
      type(time_type), intent(in) :: Time

!                         ---  local variables ---

      integer :: unit, ierr, io

!                         ---  executable code ---


      unit = open_namelist_file()
      ierr = 1
      do while (ierr /= 0)
         read(unit, nml=dry_convection_nml, iostat=io, end=20)
         ierr = check_nml_error (io, 'dry_convection_nml')
      enddo

20    call close_file (unit)

      if(mpp_pe() == mpp_root_pe()) write (stdlog(), nml=dry_convection_nml)

      id_cape = register_diag_field ( mod_name, 'CAPE', axes(1:2), Time, &
           'CAPE', 'J', missing_value=missing_value)

      id_cin = register_diag_field ( mod_name, 'CIN', axes(1:2), Time, &
           'CIN', 'J', missing_value=missing_value)
      
      id_dp = register_diag_field ( mod_name, 'dp', axes(1:2), Time, &
           'Pressure interval', 'Pa', missing_value=missing_value)

      id_lzb = register_diag_field ( mod_name, 'LZB', axes(1:2), Time, &
           'Level of zero buoyancy', 'index', missing_value=missing_value)

      id_lcl = register_diag_field ( mod_name, 'LCL', axes(1:2), Time, &
           'Lifting condensation level', 'index', missing_value=missing_value)

      id_tp = register_diag_field ( mod_name, 'parcel_temp', axes(1:3), &
           Time, 'Relaxation temperature', 'K', missing_value=missing_value)

      id_n_tp = register_diag_field ( mod_name, 'nonadj_parcel_temp', &
           axes(1:3), Time, 'Relaxation temperature before adjustment', 'K', &
           missing_value=missing_value)

      id_amb = register_diag_field ( mod_name, 'ambient_temp', axes(1:3), &
           Time, 'Ambient temperature', 'K', missing_value=missing_value)

      id_dt = register_diag_field ( mod_name, 'dt_tg', axes(1:3), &
           Time, 'Temperature tendency', 'K/s', missing_value=missing_value)

    end subroutine dry_convection_init

!-----------------------------------------------------------------------------
  
  subroutine dry_convection(Time, tg, p_full, p_half, dt_tg)

!                         ---  input arguments ---

    type(time_type), intent(in) :: Time

    real, intent(in), dimension(:,:,:) ::                                     &
         tg,               &   ! temperature
         p_full,           &   ! pressure on full model levels
         p_half                ! pressure on half model levels

!                         ---  output arguments ---

    real, intent(out), dimension(:,:,:) ::                                    &
         dt_tg                 ! temperature tendency

!                         ---  local variables ---

    real, dimension(size(tg,1),size(tg,2)) ::                                 &
         cape,             &   ! convectively available potential energy
         cin,              &   ! convective inhibition
         dp,               &   ! pressure interval from ground to LZB
         ener_int              ! energy integral from ground to LZB

    integer, dimension(size(tg,1),size(tg,2)) ::                              &
         lcl,              &   ! lifting condensation level (index)
         lzb,              &   ! level of zero buoyancy
         btm                   ! bottom of convecting region

    real, dimension(size(tg,1),size(tg,2), size(tg,3)) ::                     &
         tp,               &   ! parcel lifting temperature
         n_tp,             &   ! parcel lifting temperature before adjustment
         dp_half               ! spacing between half pressure levels

    logical :: used
!                         ---  executable code ---

    num_levels = size(tg,3)

    ! half-level spacings:
    do k=1, num_levels
       dp_half(:,:,k) = p_half(:,:,k+1) - p_half(:,:,k)
    end do
    
!              --- calculate lower bound of convecting region ---
!                         (now always the lowest level)

    btm = num_levels

!               --- calculate various convection quantities ---
    
    call capecalc(                tg,                          p_full,        &
                              p_half,                         dp_half,        &
                                 btm,                             lzb,        &
                                 lcl,                            cape,        &
                                 cin,                              tp )

    n_tp = tp ! save uncorrected parcel profile

!                       --- conservation of energy ---

    ener_int(:,:) = 0.0
    dp(:,:) = 0.0
    
    do i=1, size(tg,1)
       do j=1, size(tg,2)
          
          do k=1, num_levels
             if(k>=lzb(i,j).and. k<=btm(i,j)) then
                ener_int(i,j) = ener_int(i,j) +                               &
                     dp_half(i,j,k) * (tg(i,j,k) - tp(i,j,k))
                dp(i,j) = dp(i,j) + dp_half(i,j,k)
             else
                tp(i,j,k) = tg(i,j,k) ! ambient in non-convecting regions
             endif
          end do
          
          ! normalize integral 
          ener_int(i,j) = ener_int(i,j) / dp(i,j)
          
          do k=btm(i,j), lzb(i,j), -1
             tp(i,j,k) = tp(i,j,k) + ener_int(i,j)
          end do
       end do
    end do
    
    ! temperature tendency
    dt_tg = (tp - tg) / tau


    if (id_cin  > 0) used = send_data ( id_cin,  cin,        Time, 1, 1)
    if (id_cape > 0) used = send_data ( id_cape, cape,       Time, 1, 1)
    if (id_lzb  > 0) used = send_data ( id_lzb,  float(lzb), Time, 1, 1)
    if (id_lcl  > 0) used = send_data ( id_lcl,  float(lcl), Time, 1, 1)
    if (id_tp   > 0) used = send_data ( id_tp,   tp,         Time, 1, 1)
    if (id_n_tp > 0) used = send_data ( id_n_tp, n_tp,       Time, 1, 1)
    if (id_dp   > 0) used = send_data ( id_dp,   dp,         Time, 1, 1)
    if (id_amb  > 0) used = send_data ( id_amb,  tg,         Time, 1, 1)
    if (id_dt   > 0) used = send_data ( id_dt,   dt_tg,      Time, 1, 1)

  end subroutine dry_convection

! ----------------------------------------------------------------------------

  subroutine capecalc(tg, p_full, p_half, dp_half, btm, lzb, lcl, cape,       &
       cin, tp)

!                         ---  input arguments ---

    integer, intent(in), dimension(:,:) ::                                    &
         btm                   ! level parcel is lifted from

    real, intent(in), dimension(:,:,:) ::                                     &
         tg,               &   ! gridpoint temperature
         p_full,           &   ! gridpoint pressure on full model levels
         p_half,           &   ! gridpoint pressure on half model levels
         dp_half               ! spacing between half pressure levels
    
!                         ---  output arguments ---

    real, intent(out), dimension(:,:) ::                                      &
         cape,             &   ! convectively available potential energy
         cin                   ! convective inhibition

    integer, dimension(:,:) ::                                                &
         lzb,              &   ! level of zero buoyancy (top of convection)
         lcl                   ! lifting condensation level (btm of convection)
 
    real, intent(out), dimension(:,:,:) ::                                    &
         tp                    ! parcel lifting temperature

!                         ---  local variables ---

    real :: zdpkpk             ! pressure spacing
    
!                         ---  executable code ---
    
    cape(:,:) = 0.0; cin(:,:) = 0.0
    lzb(:,:) = btm; lcl(:,:) = btm

    tp(:,:,:) = tg(:,:,:)

    do i = 1, size(tg,1)
       do j = 1, size(tg,2)
          
          ! lift parcel with lapse rate given by gamma
          do k = btm(i,j)-1, 1, -1
             zdpkpk = exp( cons1 * alog(p_full(i,j,k)/p_full(i,j,k+1)))
             tp(i,j,k) = tp(i,j,k+1) +                                        &
                  gamma * (tp(i,j,k+1)*zdpkpk - tp(i,j,k+1))
          end do

          ! find LCL
          do k = btm(i,j)-1, 1, -1
             
             
             if(tp(i,j,k) > tg(i,j,k)) then ! unstable parcel

                if(lzb(i,j) == btm(i,j) ) then ! not above a lower cloud
                   ! calculate CAPE
                   cape(i,j) = cape(i,j) +                                    &
                        (tp(i,j,k)-tg(i,j,k))*dp_half(i,j,k)/p_full(i,j,k)

                   ! set LCL if parcel is stable at next lower level 
                   if(tp(i,j,k+1) < tg(i,j,k+1)) lcl(i,j) = k

                   ! set LZB if parcel is stable at next higher level 
!                   if(tp(i,j,k-1) < tg(i,j,k-1)) lzb(i,j) = k
                   ! changed by cw 3/30/04 so that convection doesn't
                   ! go through top of model
                   if( (tp(i,j,k-1) < tg(i,j,k-1)) .or. (k == 1)) lzb(i,j) = k

                else
                   tp(i,j,k) = tg(i,j,k)
                end if
                
             end if
             

             if(tp(i,j,k) <= tg(i,j,k)) then ! stable parcel
                
                if((lzb(i,j) == btm(i,j))) then ! not above a lower cloud

                   if(lcl(i,j) == btm(i,j)) then
                      ! calculate CIN (only if below LCL)
                      cin(i,j) = cin(i,j) -                                   &
                           (tp(i,j,k)-tg(i,j,k))*dp_half(i,j,k)/p_full(i,j,k)
                   end if

                else ! set parcel temp to ambient
                   tp(i,j,k) = tg(i,j,k)
                end if


             end if
             
          end do

          ! if cin > cape, turn off convection by setting 
          ! relaxation temperature to ambient temperature
          if(cin(i,j) > cape(i,j)) tp(i,j,:) = tg(i,j,:)

          ! a few checks
          if( (lcl(i,j) /= btm(i,j)) .and. (lzb(i,j) == btm(i,j)))        &
               call error_mesg ('dry_convection','LCL defined, LZB not defined', FATAL)
          if( lcl(i,j) < lzb(i,j) )                                           &
               call error_mesg ('dry_convection','LCL above LZB', FATAL)
          
          if( (lcl(i,j) == btm(i,j)) .and. (lzb(i,j) == btm(i,j))) then
             cape(i,j) = 0.0
             cin(i,j) = 0.0
          end if

       end do
    end do

  end subroutine capecalc

end module dry_convection_mod
