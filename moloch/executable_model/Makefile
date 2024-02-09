# Following parameters are defined in main Makefile (../../Makefile) :
#
#USE_MPI = YES
#USE_RAD_ECMWF = YES
#
#FC = gfortran
#FC_MPI = mpif90
#
####FCFLAGS = Ofast -g -fbacktrace -march=native -funroll-loops -ffpe-summary=none -w
#FCFLAGS = -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -ffpe-summary=none  # debuger options
#
#FCFLAGS_MPI = -Dmpi -g -O2 -march=native -funroll-loops -w -ffpe-summary=none -finit-integer=0 -finit-local-zero -finit-real=zero
#
#DIR_SOURCES_MOLOCH = $(HOME)/sources/moloch
#DIR_SOURCES_COMMON = $(HOME)/sources/common
#DIR_RAD = /home/meteo/rad_ecmwf/executable_$(FC)

FCFLAGS += -I../
FCFLAGS_MPI += -I ../

ifeq ($(USE_MPI),NO)

FC_USED = $(FC)
FCFLAGS_USED = $(FCFLAGS)

else

FC_USED = $(FC_MPI)
FCFLAGS_USED = $(FCFLAGS_MPI)

endif

LDFLAGS_USED = $(FCFLAGS_USED)

ifeq ($(USE_RAD_ECMWF),YES)
FCFLAGS_USED += -Drad_ecmwf -I$(DIR_RAD)/module -I$(DIR_RAD)/include
endif

ifeq ($(USE_RAD_ECMWF),YES)
moloch: moloch.o pochva_module.o pochva.o $(DIR_RAD)/libecmwf_rad.a
else
moloch: moloch.o pochva_module.o pochva.o radiat_param.o
endif
	$(FC_USED) -o $@ $^ $(LDFLAGS_USED)

ifeq ($(USE_RAD_ECMWF),YES)
moloch.o: $(DIR_SOURCES_MOLOCH)/moloch.F90 pochva.o ../dimensions.inc
else
moloch.o: $(DIR_SOURCES_MOLOCH)/moloch.F90 pochva.o radiat_param.o ../dimensions.inc
endif
	$(FC_USED) -c $(FCFLAGS_USED) $<

radiat_param.o: $(DIR_SOURCES_COMMON)/radiat_param.F90
	$(FC_USED) -c $(FCFLAGS_USED) $<

pochva.o: $(DIR_SOURCES_COMMON)/pochva.F90 pochva_module.o
	$(FC_USED) -c $(FCFLAGS_USED) $<

pochva_module.o: $(DIR_SOURCES_COMMON)/pochva_module.F90
	$(FC_USED) -c $(FCFLAGS_USED) $<

clean:
	rm -f moloch *.o *.mod
