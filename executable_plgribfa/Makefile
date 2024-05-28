# Following parameters are defined in main Makefile (../../Makefile):
#DIR_ECCODES_GFORTRAN = $(HOME)/eccodes/eccodes_gfortran
#DIR_SOURCES_PLGRIBFA = $(HOME)/sources/plgribfa
#DIR_PLPLOT = /usr/lib/x86_64-linux-gnu/fortran/modules/plplot

FC = gfortran
#FCFLAGS = -g -fbounds-check
FCFLAGS = -O2
LDFLAGS = $(FCFLAGS)
DIR_ECCODES = $(DIR_ECCODES_GFORTRAN)

FCFLAGS += -I $(DIR_PLPLOT)
FCFLAGS_GRIB = $(FCFLAGS) -I $(DIR_ECCODES)/include

LDFLAGS += `pkg-config --libs plplot-fortran`
###LDFLAGS += `pkg-config --libs plplotd-f95`

LDFLAGS_GRIB = $(LDFLAGS) -L $(DIR_ECCODES)/lib -leccodes_f90 -leccodes

all: plgribfa plgribfa_work plgribfa_work_my

plgribfa: plgribfa_main_input_grib2.o plgribfa_module.o plplot_map.o
	$(FC) -o $@ $^ $(LDFLAGS_GRIB)

plgribfa_work: plgribfa_main_input_work.o plgribfa_module.o plplot_map.o
	$(FC) -o $@ $^ $(LDFLAGS_GRIB)

plgribfa_work_my: plgribfa_main_input_work_my.o plgribfa_module.o plplot_map.o
	$(FC) -o $@ $^ $(LDFLAGS_GRIB)

plplot_map.o: $(DIR_SOURCES_PLGRIBFA)/plplot_map.F90 plgribfa_module.o
	$(FC) -c $(FCFLAGS) $<

plgribfa_main_input_grib2.o: $(DIR_SOURCES_PLGRIBFA)/plgribfa_main_input_grib2.F90 plgribfa_module.o
	$(FC) -c $(FCFLAGS_GRIB) $<

plgribfa_main_input_work.o: $(DIR_SOURCES_PLGRIBFA)/plgribfa_main_input_work.F90 plgribfa_module.o
	$(FC) -c $(FCFLAGS_GRIB) $<

plgribfa_main_input_work_my.o: $(DIR_SOURCES_PLGRIBFA)/plgribfa_main_input_work_my.F90 plgribfa_module.o
	$(FC) -c $(FCFLAGS_GRIB) $<

plgribfa_module.o: $(DIR_SOURCES_PLGRIBFA)/plgribfa_module.F90
	$(FC) -c $(FCFLAGS) $<

clean:
	rm -f *.o *.mod plgribfa plgribfa_work plgribfa_work_my
