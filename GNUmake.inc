PLATFORM=$(shell uname)
ifeq ($(PLATFORM),OSF1)
	RATFOR = ratfor
	FC = f90
	FOPT = -O5 -tune ev6
	LAPACK =
#	LAPACK = dsygvd.a
	BLAS = -ldxml
endif
ifeq ($(PLATFORM),AIX)
        RATFOR = ratfor
        FC = f77 
        FOPT = -O3 -qstrict  
        LAPACK = dsygvd.a
        BLAS = -lblas 
endif
ifeq ($(PLATFORM),HP-UX)
        RATFOR = ratfor
        FC = f77
        FOPT = -O3  -K +Olibcalls +Odataprefetch +Ofastaccess +Oentrysched \
          +Oregionsched +Onolimit  +Onoloop_transform 
	LAPACK =
        BLAS = +Ofastaccess -Wl,+FPD \
          -Wl,-x -Wl,-N -Wl,-n -lveclib -llapack 
endif
ifeq ($(PLATFORM),IRIX64)
        RATFOR = ratfor
        FC = f77
        FOPT = -c -O3 -mips4 -r10000 
        LAPACK = dsygvd.a
        BLAS = -lcomplib.sgimath 
endif
ifeq ($(PLATFORM),Linux)
 	RATFOR = /usr/local/bin/ratfor
	FC = ifort
  FOPT = -O3 -mtune=core2 -ip -ipo -L/opt/intel/mkl/lib/intel64 
	LAPACK = -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
	BLAS = 
endif
ifeq ($(PLATFORM),Darwin)
    RATFOR = ../ratfor/ratfor-1.01/ratfor
    FC     = gfortran
    FOPT   = -O3 -funroll-loops -pipe -funroll-all-loops -finline-functions -framework Accelerate
#   FOPT   = -O0 -g -Wall -pedantic -framework Accelerate -pipe
    LAPACK = 
    BLAS   = 
endif
