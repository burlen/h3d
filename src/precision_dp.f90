module PRECISION !LAURA
integer, parameter :: SIZE_REAL = 4, SIZE_DOUBLE = 8       !LAURA
integer, parameter :: CUSTOM_REAL = SIZE_DOUBLE     !LAURA

!NOTE flux calculation requires double precision
integer, parameter :: CUSTOM_FLUX = SIZE_DOUBLE     !LAURA
end module !LAURA

module PRECISION_MPI
include "mpif.h"
integer, parameter :: CUSTOM_MPI_TYPE = MPI_DOUBLE_PRECISION
end module

