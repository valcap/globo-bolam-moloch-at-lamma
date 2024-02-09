# Following parameters are definedi in main Makefile (../../Makefile) :
###FC = ifort
###FCFLAGS = -O2 -vec-report0 -W0 -warn none -Doper
###DIR_ECCODES = ../../eccodes_$(FC)
###DIR_SOURCES_COMMON = ../../sources/common

LDFLAGS = $(FCFLAGS) 

FCFLAGS += -I$(DIR_ECCODES)/include
LDFLAGS += -L $(DIR_ECCODES)/lib -leccodes_f90 -leccodes

read_grib2_write_metrogram: read_grib2_write_metrogram.o common_routines.o read_grib2_data.o
	$(FC) -o $@ $^ $(LDFLAGS)

read_grib2_write_metrogram.o: $(DIR_SOURCES_COMMON)/read_grib2_write_metrogram.F90
	$(FC) -c $(FCFLAGS) $<

common_routines.o: $(DIR_SOURCES_COMMON)/common_routines.F90
	$(FC) -c $(FCFLAGS) $<

read_grib2_data.o: $(DIR_SOURCES_COMMON)/read_grib2_data.F90
	$(FC) -c $(FCFLAGS) $<

clean:
	rm -f read_grib2_write_metrogram *.o *.mod
