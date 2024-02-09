# Following parameters are definedi in main Makefile (../../Makefile) :
###FC = ifort
###FCFLAGS = -O2 -ip -vec-report0 -W0 -warn none -Doper
###DIR_ECCODES = ../../eccodes/eccodes_ifort
###DIR_SOURCES_BOLAM = ../../sources/bolam
###DIR_SOURCES_COMMON = ../../sources/common

LDFLAGS = $(FCFLAGS) 

ifeq ($(FC), ifort)
FCFLAGS += -assume byterecl
endif

FCFLAGS += -I$(DIR_ECCODES)/include
LDFLAGS += -L $(DIR_ECCODES)/lib -leccodes_f90 -leccodes

####$(info $$DIR_ECCODES = ${DIR_ECCODES})

prebolam: prebolam.o def_soil.o common_routines.o read_grib2_data.o
	$(FC) -o $@ $^ $(LDFLAGS)

prebolam.o: $(DIR_SOURCES_BOLAM)/prebolam.F90
	$(FC) -c $(FCFLAGS) $<

read_grib2_data.o: $(DIR_SOURCES_COMMON)/read_grib2_data.F90
	$(FC) -c $(FCFLAGS) $<

def_soil.o: $(DIR_SOURCES_COMMON)/def_soil.F90
	$(FC) -c $(FCFLAGS) $<

common_routines.o: $(DIR_SOURCES_COMMON)/common_routines.F90
	$(FC) -c $(FCFLAGS) $<

clean:
	rm -f prebolam *.o *.mod
