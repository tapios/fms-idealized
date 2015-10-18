module leapfrog_mod

use fms_mod, only: mpp_pe, mpp_root_pe, error_mesg, FATAL, write_version_number

!===================================================================================================
implicit none
private
!===================================================================================================

interface leapfrog
  module procedure leapfrog_2d_complex, leapfrog_3d_complex
end interface

interface leapfrog_2level_A
  module procedure leapfrog_2level_A_2d_complex, &
                   leapfrog_2level_A_3d_complex, &
                   leapfrog_2level_A_3d_real
end interface

interface leapfrog_2level_B
  module procedure leapfrog_2level_B_2d_complex, &
                   leapfrog_2level_B_3d_complex, &
                   leapfrog_2level_B_3d_real
end interface

character(len=128), parameter :: version = '$Id leapfrog.f90 $'
character(len=128), parameter :: tagname = '$Name: latest $'

public :: leapfrog, leapfrog_2level_A, leapfrog_2level_B

logical :: entry_to_logfile_done = .false.

contains

!================================================================================

subroutine leapfrog_2level_A_3d_complex (a, dt_a, previous, current, future, delta_t, robert_coeff, raw_factor)  
! Changed by ZTAN: store the (prev-2*curr) value in dt_a
complex, intent(inout), dimension(:,:,:,:) :: a
complex, intent(inout), dimension(:,:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_factor ! raw_factor added by ZTAN

complex, dimension(size(dt_a,1), size(dt_a,2), size(dt_a,3)) :: filt     ! Added by ZTAN: prev-2*curr

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

filt(:,:,:) = a(:,:,:,previous) - 2.0*a(:,:,:,current) ! Added by ZTAN

if(previous == current) then
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))* raw_factor ! raw_factor added by ZTAN
else
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))* raw_factor ! raw_factor added by ZTAN
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
endif

dt_a (:,:,:)  = filt(:,:,:) ! Added by ZTAN

return
end subroutine leapfrog_2level_A_3d_complex

!================================================================================

subroutine leapfrog_2level_B_3d_complex (a, filt, current, future, robert_coeff, raw_factor) ! filt and raw_factor added by ZTAN
! Changed by ZTAN: the (prev-2*curr) value is in filt

complex, intent(inout), dimension(:,:,:,:) :: a
complex, intent(inout), dimension(:,:,:  ) :: filt  ! filt added by ZTAN
integer, intent(in) :: current, future
real,    intent(in) :: robert_coeff, raw_factor ! raw_factor added by ZTAN

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a(:,:,:,current) = a(:,:,:,current) + robert_coeff * a(:,:,:,future)*  raw_factor ! raw_factor added by ZTAN
a(:,:,:,future ) = a(:,:,:,future ) + robert_coeff * ( filt(:,:,:) + a(:,:,:,future))* (raw_factor-1.0)  ! Added by ZTAN

return
end subroutine leapfrog_2level_B_3d_complex

!================================================================================

subroutine leapfrog_2level_A_3d_real (a, dt_a, previous, current, future, delta_t, robert_coeff, raw_factor) ! raw_factor added by ZTAN

real, intent(inout), dimension(:,:,:,:) :: a
real, intent(inout), dimension(:,:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_factor ! raw_factor added by ZTAN

real, dimension(size(dt_a,1), size(dt_a,2), size(dt_a,3)) :: filt     ! Added by ZTAN: prev-2*curr

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

filt(:,:,:) = a(:,:,:,previous) - 2.0*a(:,:,:,current) ! Added by ZTAN

if(previous == current) then
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))* raw_factor ! raw_factor added by ZTAN
else
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))* raw_factor ! raw_factor added by ZTAN 
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
endif

dt_a (:,:,:)  = filt(:,:,:) ! Added by ZTAN

return
end subroutine leapfrog_2level_A_3d_real

!================================================================================

subroutine leapfrog_2level_B_3d_real (a, filt, current, future, robert_coeff, raw_factor) ! filt and raw_factor added by ZTAN
! Changed by ZTAN: the (prev-2*curr) value is in filt

real, intent(inout), dimension(:,:,:,:) :: a
real, intent(inout), dimension(:,:,:  ) :: filt ! filt added by ZTAN
integer, intent(in) :: current, future
real,    intent(in) :: robert_coeff, raw_factor ! raw_factor added by ZTAN

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a(:,:,:,current) = a(:,:,:,current) + robert_coeff * a(:,:,:,future)* raw_factor
a(:,:,:,future ) = a(:,:,:,future ) + robert_coeff * ( filt(:,:,:) + a(:,:,:,future))* (raw_factor-1.0)  ! Added by ZTAN

return
end subroutine leapfrog_2level_B_3d_real


!================================================================================

subroutine leapfrog_2level_A_2d_complex (a, dt_a, previous, current, future, delta_t, robert_coeff, raw_factor) ! raw_factor added by ZTAN

complex, intent(inout), dimension(:,:,:) :: a
complex, intent(inout), dimension(:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_factor ! raw_factor added by ZTAN

complex, dimension(size(a,1),size(a,2),1,size(a,3)) :: a_3d
complex, dimension(size(a,1),size(a,2),1)           :: dt_a_3d

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a_3d(:,:,1,:) = a
dt_a_3d(:,:,1) = dt_a
call leapfrog_2level_A_3d_complex(a_3d, dt_a_3d, previous, current, future, delta_t, robert_coeff, raw_factor) ! raw_factor added by ZTAN
a = a_3d(:,:,1,:)
dt_a = dt_a_3d(:,:,1)

end subroutine leapfrog_2level_A_2d_complex
!================================================================================

subroutine leapfrog_2level_B_2d_complex (a, filt, current, future, robert_coeff, raw_factor) ! filt and raw_factor added by ZTAN

complex, intent(inout), dimension(:,:,:) :: a
complex, intent(inout), dimension(:,:  ) :: filt ! added ZTAN
integer, intent(in) :: current, future
real,    intent(in) :: robert_coeff,raw_factor
complex, dimension(size(a,1),size(a,2),1,size(a,3)) :: a_3d
complex, dimension(size(a,1),size(a,2),1)           :: filt_3d ! added ZTAN
  
if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a_3d(:,:,1,:) = a
filt_3d(:,:,1) = filt ! added ZTAN
call leapfrog_2level_B_3d_complex (a_3d, filt_3d, current, future, robert_coeff, raw_factor) ! filt and raw_factor added by ZTAN
a = a_3d(:,:,1,:)
filt = filt_3d(:,:,1) ! added ZTAN

end subroutine leapfrog_2level_B_2d_complex
!================================================================================

subroutine leapfrog_3d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff, raw_factor) ! raw_factor added by ZTAN

complex, intent(inout), dimension(:,:,:,:) :: a
complex, intent(in),    dimension(:,:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_factor ! raw_factor added by ZTAN

complex, dimension(size(a,1),size(a,2),size(a,3)) :: filt ! added ZTAN


if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

filt (:,:,:) = a(:,:,:,previous) - 2.0*a(:,:,:,current)

if(previous == current) then
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current) + a(:,:,:,future ))* raw_factor ! raw_factor added by ZTAN
else
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*(a(:,:,:,previous) - 2.0*a(:,:,:,current))* raw_factor ! raw_factor added by ZTAN
  a(:,:,:,future ) = a(:,:,:,previous) + delta_t*dt_a
  a(:,:,:,current) = a(:,:,:,current ) + robert_coeff*a(:,:,:,future)* raw_factor ! raw_factor added by ZTAN
endif

a(:,:,:,future ) = a(:,:,:,future ) + robert_coeff * ( filt(:,:,:) + a(:,:,:,future))* (raw_factor-1.0) ! added by ZTAN

return
end subroutine leapfrog_3d_complex

!================================================================================

subroutine leapfrog_2d_complex(a, dt_a, previous, current, future, delta_t, robert_coeff, raw_factor) ! raw_factor added by ZTAN

complex, intent(inout), dimension(:,:,:) :: a
complex, intent(in),    dimension(:,:  ) :: dt_a
integer, intent(in) :: previous, current, future
real,    intent(in) :: delta_t, robert_coeff, raw_factor ! raw_factor added by ZTAN

complex, dimension(size(a,1),size(a,2),1,size(a,3)) :: a_3d
complex, dimension(size(a,1),size(a,2),1)           :: dt_a_3d

if(.not.entry_to_logfile_done) then
  call write_version_number(version, tagname)
  entry_to_logfile_done = .true.
endif

a_3d(:,:,1,:) = a
dt_a_3d(:,:,1) = dt_a
call leapfrog_3d_complex(a_3d, dt_a_3d, previous, current, future, delta_t, robert_coeff, raw_factor) ! raw_factor added by ZTAN
a = a_3d(:,:,1,:)

return
end subroutine leapfrog_2d_complex

!================================================================================

end module leapfrog_mod
