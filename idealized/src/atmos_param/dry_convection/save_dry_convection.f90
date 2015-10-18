

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

  integer :: id_cape, id_cin, id_lzb, id_lcl, id_tp, id_dp, id_amb, id_dt

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
           'CAPE', 'J', missing_value=missing_value,&
           verbose = .true.)

      id_cin = register_diag_field ( mod_name, 'CIN', axes(1:2), Time, &
           'CIN', 'J', missing_value=missing_value,&
           verbose = .true.)
      
      id_dp = register_diag_field ( mod_name, 'dp', axes(1:2), Time, &
           'Pressure interval', 'Pa', missing_value=missing_value,&
           verbose = .true.)

      id_lzb = register_diag_field ( mod_name, 'LZB', axes(1:2), Time, &
           'Level of zero buoyancy', 'index', missing_value=missing_value,&
           verbose = .true.)

      id_lcl = register_diag_field ( mod_name, 'LCL', axes(1:2), Time, &
           'Lifting condensation level', 'index', missing_value=missing_value,&
           verbose = .true.)

      id_tp = register_diag_field ( mod_name, 'parcel_temp', axes(1:3), &
           Time, 'Relaxation temperature', 'K', missing_value=missing_value,&
           verbose = .true.)

      id_amb = register_diag_field ( mod_name, 'ambient_temp', axes(1:3), &
           Time, 'Ambient temperature', 'K', missing_value=missing_value,&
           verbose = .true.)

      id_dt = register_diag_field ( mod_name, 'dt_tg', axes(1:3), &
           Time, 'Temperature tendency', 'K/s', missing_value=missing_value,&
           verbose = .true.)

    end subroutine dry_convection_init

!-----------------------------------------------------------------------------
  
  subroutine dry_convection(Time, tg, p_full, p_half, dt_tg)

!                         ---  input arguments ---

    type(time_type), intent(in) :: Time

    real, intent(in), dimension(:,:,:) ::                                     &
         tg,               &   ! gridpoint temperature
         p_full,           &   ! gridpoint pressure on full model levels
         p_half                ! gridpoint pressure on half model levels

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
         lzb                   ! level of zero buoyancy

    real, dimension(size(tg,1),size(tg,2), size(tg,3)) ::                     &
         tp,               &   ! parcel lifting temperature
         dp_half               ! spacing between half pressure levels

    logical :: used
!                         ---  executable code ---

    num_levels = size(tg,3)

    ! half-level spacings:
    do k=1, num_levels
       dp_half(:,:,k) = p_half(:,:,k+1) - p_half(:,:,k)
    end do
    
    call capecalc(tg, p_full, p_half, dp_half, lzb, lcl, cape, cin, tp, dp)

!                       --- conservation of energy ---

    do jit=1,2 ! two iterations of energy correction
       
       ! calculate column integral of enthaply
       ! tp is set to tg everywhere there is no convection, so calculation
       ! of integral is vectorized
       ener_int(:,:) = 0.0

!!$       do k=1, num_levels
!!$          ener_int(:,:) = ener_int(:,:) +                                     &
!!$               dp_half(:,:,k) * (tg(:,:,k) - tp(:,:,k))
!!$       end do
!!$
!!$       ! normalize integral 
!!$       ener_int(:,:) = ener_int(:,:) / dp(:,:)

       dp(:,:) = 0.0
       do i=1, size(tg,1)
          do j=1, size(tg,2)

             do k=lzb(i,j), num_levels
                ener_int(i,j) = ener_int(i,j) +                               &
                     dp_half(i,j,k) * (tg(i,j,k) - tp(i,j,k))
                dp(i,j) = dp(i,j) + dp_half(i,j,k)
             end do

             ! normalize integral 
             ener_int(i,j) = ener_int(i,j) / dp(i,j)

!!$             if(i<5.and.j<5) &
!!$             write(*,*) 'int',jit,i,j,ener_int(i,j), dp(i,j),&
!!$             p_full(i,j,lcl(i,j)),p_full(i,j,lzb(i,j))

             do k=num_levels, lzb(i,j), -1
!                tp(i,j,k) = tp(i,j,k) + ener_int(i,j)
             end do
          end do
       end do

    end do ! end of iterations


    dt_tg = (tp - tg) / tau


    if (id_cin  > 0) used = send_data ( id_cin,  cin,        Time, 1, 1)
    if (id_cape > 0) used = send_data ( id_cape, cape,       Time, 1, 1)
    if (id_lzb  > 0) used = send_data ( id_lzb,  float(lzb), Time, 1, 1)
    if (id_lcl  > 0) used = send_data ( id_lcl,  float(lcl), Time, 1, 1)
    if (id_tp   > 0) used = send_data ( id_tp,   tp,         Time, 1, 1)
    if (id_dp   > 0) used = send_data ( id_dp,   dp,         Time, 1, 1)
    if (id_amb  > 0) used = send_data ( id_amb,  tg,         Time, 1, 1)
    if (id_dt   > 0) used = send_data ( id_dt,   dt_tg,      Time, 1, 1)

  end subroutine dry_convection

! ----------------------------------------------------------------------------

  subroutine capecalc(tg, p_full, p_half, dp_half, lzb, lcl, cape, cin, tp, dp)

!                         ---  input arguments ---

    real, intent(in), dimension(:,:,:) ::                                     &
         tg,               &   ! gridpoint temperature
         p_full,           &   ! gridpoint pressure on full model levels
         p_half,           &   ! gridpoint pressure on half model levels
         dp_half               ! spacing between half pressure levels
    
!                         ---  output arguments ---

    real, intent(out), dimension(:,:) ::                                      &
         cape,             &   ! convectively available potential energy
         cin,              &   ! convective inhibition
         dp                    ! pressure interval from ground to LZB

    integer, dimension(:,:) ::                                                &
         lzb,              &   ! level of zero buoyancy (top of convection)
         lcl                   ! lifting condensation level (btm of convection)
 
    real, intent(out), dimension(:,:,:) ::                                    &
         tp                    ! parcel lifting temperature

!                         ---  local variables ---

    real :: zdpkpk             ! pressure spacing
    
!                         ---  executable code ---
    
    cape(:,:) = 0.0; cin(:,:) = 0.0
    lzb(:,:) = num_levels; lcl(:,:) = num_levels

    tp(:,:,num_levels) = tg(:,:,num_levels)
    dp(:,:) = dp_half(:,:,num_levels)

    do i = 1, size(tg,1)
       do j = 1, size(tg,2)
          
          
          ! find LCL
          do k = num_levels-1, 1, -1
             
             zdpkpk = exp( cons1 * alog(p_full(i,j,k)/p_full(i,j,k+1)))
             tp(i,j,k) = tp(i,j,k+1) +                                        &
                  gamma * (tp(i,j,k+1)*zdpkpk - tp(i,j,k+1))
             
             if(tp(i,j,k) > tg(i,j,k)) then ! unstable parcel

                if(lzb(i,j) == num_levels ) then ! not above a lower cloud
                   ! calculate CAPE
                   cape(i,j) = cape(i,j) +                                    &
                        (tp(i,j,k)-tg(i,j,k))*dp_half(i,j,k)/p_full(i,j,k)
                   ! compute pressure interval
                   dp(i,j) = dp(i,j) + dp_half(i,j,k)
!!$                   if(i<5.and.j<5) write(*,*) 'one',i,j,k,dp(i,j),        &
!!$                        cape(i,j),cin(i,j)


                   ! set LCL if parcel is stable at next lower level 
                   if((tp(i,j,k+1) < tg(i,j,k+1))) lcl(i,j) = k

                else
                   tp(i,j,k) = tg(i,j,k)
                end if
                
             end if
             

             if(tp(i,j,k) <= tg(i,j,k)) then ! stable parcel
                
                if((lzb(i,j) == num_levels)) then ! not above a lower cloud
                   ! set LZB if parcel is unstable at next lower level
                   if(tp(i,j,k+1) > tg(i,j,k+1)) lzb(i,j) = k

                   if(lcl(i,j) == num_levels) then
                      ! calculate CIN (only if below LCL)
                      cin(i,j) = cin(i,j) -                                   &
                           (tp(i,j,k)-tg(i,j,k))*dp_half(i,j,k)/p_full(i,j,k)
                      !  compute pressure interval
                      dp(i,j) = dp(i,j) + dp_half(i,j,k)
!!$                      if(i<5.and.j<5) write(*,*) 'two',i,j,k,dp(i,j),&
!!$                           cape(i,j),cin(i,j)
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
          if( (lcl(i,j) /= num_levels) .and. (lzb(i,j) == num_levels))        &
               call error_mesg ('dry_convection','LCL defined, LZB not defined', FATAL)
          if( lcl(i,j) < lzb(i,j) )                                           &
               call error_mesg ('dry_convection','LCL above LZB', FATAL)
          
          if( (lcl(i,j) == num_levels) .and. (lzb(i,j) == num_levels)) then
             cape(i,j) = 0.0
             cin(i,j) = 0.0
          end if

       end do
    end do

  end subroutine capecalc

end module dry_convection_mod
