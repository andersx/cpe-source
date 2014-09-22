TARGET=dftb_dev
MAXIMA=maxima
ATOMS=500
ORBITALS=2000
EXTERNALCHARGES=100000
.SUFFIXES : .o .r .f .inc .h .f90
#
include GNUmake.inc
#
.h.inc : 
	@$(RATFOR) < $*.h > $*.inc
	@printf .
.r.o :
	@$(RATFOR) < $*.r > $*.f
	@$(FC) $(FOPT) -c $*.f
	@printf .
.f90.o :
	@$(FC) $(FOPT) -c $*.f90
.f.o :
	@$(FC) $(FOPT) -c $*.f
	@printf .
RSRC= dylcao.f gettab.f skpar.f mulliken.f \
      slkode.f vnorm.f xnorm.f geometries.f repulsiv.f fermi.f \
      output.f gradient.f vibout.f
#
INC=  commonbox.inc
OBJ=  dylcao.o eglcao.o skpar.o gettab.o xnorm.o geometries.o \
      vnorm.o slkode.o ewevge.o slktrafo.o self.o input.o \
      cgrad.o broyden.o mulliken.o gradient.o linpack.o \
      short_range.o long_range.o shift.o repulsiv.o  \
      fermi.o output.o gamma.o gammamat.o atomen.o  \
      dispersion_egr.o dispersionread.o  externalchgrad.o \
      externalshift.o  reada.o vibout.o bfgs.o bfgs_timer.o \
			bfgsgeom.o dpmeps.o gamnew.o mixer.o ranmar.o \
	    cm3charges.o fixed_fermi.o invert_eta_cpe.o cpe_scf.o cpe_gradients.o 
#
all: $(TARGET)
#
eglcao.o: $(MAXIMA).inc
#
#
$(MAXIMA).h:
ifndef ATOMS
	@echo "	** Define maximum number of atoms using "
	@echo "	** make ATOMS=100 ORBITALS=400    "
	@echo "	** (replace the numbers with the desired amount...)"
	@I_am_quitting_now,_producing_an_ugly_error_message
else # max no. of atoms defined
ifndef ORBITALS
	@echo "	** Define maximum number of orbitals using "
	@echo "	** make ATOMS=100 ORBITALS=400    "
	@echo "	** (replace the numbers with the desired amount...)"
	@I_am_quitting_now,_producing_an_ugly_error_message
else # all we need is defined...
	@echo "      integer nndim" >> $(MAXIMA).h
	@echo "      parameter( nndim= $(ATOMS))" >> $(MAXIMA).h
	@echo "      integer mdim" >> $(MAXIMA).h
	@echo "      parameter( mdim= $(ORBITALS))" >> $(MAXIMA).h
	@echo "      integer maxint" >> $(MAXIMA).h
	@echo "      parameter( maxint= 250)" >> $(MAXIMA).h
	@echo "      integer maxtyp" >> $(MAXIMA).h
	@echo "      parameter( maxtyp= 10)" >> $(MAXIMA).h
	@echo "      integer maxtab" >> $(MAXIMA).h
	@echo "      parameter( maxtab= 1000)" >> $(MAXIMA).h
	@echo "      integer ldim" >> $(MAXIMA).h
	@echo "      parameter( ldim= 9)" >> $(MAXIMA).h
	@echo "      integer maxopt" >> $(MAXIMA).h
	@echo "      parameter(maxopt=3*nndim)" >> $(MAXIMA).h
	@echo "      integer maxsiz" >> $(MAXIMA).h
	@echo "      parameter(maxsiz=mdim)" >> $(MAXIMA).h
	@echo "      integer maxtp" >> $(MAXIMA).h
	@echo "      parameter(maxtp=maxtyp)" >> $(MAXIMA).h
	@echo "      integer MXLINE" >> $(MAXIMA).h
	@echo "      parameter( MXLINE = 85)" >> $(MAXIMA).h
	@echo "      integer MAXLDT" >> $(MAXIMA).h
	@echo "      parameter( MAXLDT = 25)" >> $(MAXIMA).h
	@echo "      integer MXENTR" >> $(MAXIMA).h
	@echo "      parameter( MXENTR = 20)" >> $(MAXIMA).h
	@echo "      integer REPMAX" >> $(MAXIMA).h
	@echo "      parameter( REPMAX = 200)" >> $(MAXIMA).h
	@echo "      integer EXTDIM" >> $(MAXIMA).h
	@echo "      parameter (EXTDIM = $(EXTERNALCHARGES))" >> $(MAXIMA).h

	@printf "Compiling dftb "     
endif # ndef ORBITALS
endif # ndef ATOMS

$(TARGET): $(MAXIMA).h $(INC) $(OBJ) 
	@$(FC) $(FOPT) -o $(TARGET) $(OBJ) $(BLAS) $(LAPACK)
	@echo 
	@echo "........  Make of $(TARGET) for $(PLATFORM) completed."
	@echo "........  Required datasize limit to run it: `size $(TARGET) | \
         awk '/data/{getline; a=$$4; print int((a/1024)/1024)+1}'`M"
#	@rm -f $(RSRC) $(INC) $(OBJ) $(MAXIMA).h $(MAXIMA).inc
clean:
	@rm -f $(TARGET) $(RSRC) $(INC) $(OBJ) $(MAXIMA).h $(MAXIMA).inc 
	@echo "     ** "
	@echo "     ** $(TARGET) deleted..."
	@echo "     ** "
allclean: clean
	@rm -f dsygvd.a
	@echo "     ** "
	@echo "     **                        all clear..."
	@echo "     ** "
