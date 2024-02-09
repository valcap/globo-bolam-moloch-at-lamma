# Following parameters are defined in main Makefile (../../Makefile) :
###FC = ifort
###FCFLAGS = -O2 -vec-report0 -W0 -warn none -Doper
###FCFLAGS = -g -CU -CB -fpe0 -Doper
###DIR_ECCODES = ../../eccodes/eccodes_$(FC)
###DIR_SOURCES_MOLOCH = ../../sources/moloch
###DIR_SOURCES_COMMON = ../../sources/common

LDFLAGS = $(FCFLAGS) 

FCFLAGS += -I$(DIR_ECCODES)/include
LDFLAGS += -L $(DIR_ECCODES)/lib -leccodes_f90 -leccodes

postmoloch: postmoloch.o common_routines.o\
 prepar_grib2_coding_data.o write_grib2_data.o
	$(FC) -o $@ $^ $(LDFLAGS)

postmoloch.o:  $(DIR_SOURCES_MOLOCH)/postmoloch.F90 write_grib2_data.o
	$(FC) -c $(FCFLAGS) $<

prepar_grib2_coding_data.o: $(DIR_SOURCES_COMMON)/prepar_grib2_coding_data.F90 write_grib2_data.o
	$(FC) -c $(FCFLAGS) $<

write_grib2_data.o: $(DIR_SOURCES_COMMON)/write_grib2_data.F90
	$(FC) -c $(FCFLAGS) $<

common_routines.o:  $(DIR_SOURCES_COMMON)/common_routines.F90
	$(FC) -c $(FCFLAGS) $<

clean:
	rm -f postmoloch *.o *.mod
