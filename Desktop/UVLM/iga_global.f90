MODULE iga_global
USE REAL_PRECISION

REAL(KIND = R8), DIMENSION(:), ALLOCATABLE :: UX, UY

!!!!!!!!!! variables for bi-splines geometry !!!!!!!!!!!!!!!!!!
INTEGER :: MX, NDX, PX
INTEGER :: MY, NDY, PY

REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: GEOM

REAL(kind=8), DIMENSION(3) :: S


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


END MODULE