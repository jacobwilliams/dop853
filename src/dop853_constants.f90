module dop853_constants

use iso_fortran_env,    only: wp => real64

implicit none

real(wp),parameter :: uround = epsilon(1.0_wp)  !! machine \( \epsilon \)

end module dop853_constants
