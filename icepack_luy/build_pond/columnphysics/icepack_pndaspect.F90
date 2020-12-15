module icepack_meltpond_cesm_pndaspect

use icepack_kinds

implicit none

!====================================================

contains

!====================================================

subroutine compute_ponds_cesm_pndaspect( istep, pndaspect)

integer(kind=int_kind) :: &
istep
real(kind=dbl_kind) :: &
pndaspect , &! ratio of pond depth to pond fraction
p1=0.0001253 , &
p2=-0.006176 , &
p3=0.1155 , &
p4=-1.016 , &
p5=4.186 , &
p6=-6.989 , &
p7=4.344

character(len=*),parameter :: subname='(compute_ponds_cesm_pndaspect)'

if ( (istep < 3072) .or. (istep >5956) ) then
  pndaspect = 0.8
else 
  pndaspect = p1*x**6 + p2*x**5 + p3*x**4 + p4*x**3 + p5*x**2 + p6*x + p7
endif

end subroutine compute_ponds_cesm_pndaspect

!====================================================

end module icepack_meltpond_cesm_pndaspect

!====================================================
