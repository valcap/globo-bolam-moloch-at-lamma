# Following parameters are definedi in main Makefile (../../Makefile) :
###FC = ifort
###FCFLAGS = -g -CU -CB -fpe0 -vec-report0 -W0 -warn none -Doper
###FCFLAGS = -O2 -vec-report0 -W0 -warn none -Doper
###DIR_ECCODES = /home/oper/ECCODES/ECCODES_$(FC)
###DIR_SOURCES_COMMON = ../../sources/common

LDFLAGS = $(FCFLAGS) 

FCFLAGS += -I$(DIR_ECCODES)/include
LDFLAGS += -L $(DIR_ECCODES)/lib -leccodes_f90 -leccodes

convert_shf_to_grib2: convert_shf_to_grib2.o write_grib2_data.o
	$(FC) -o $@ $^ $(LDFLAGS)

convert_shf_to_grib2.o: $(DIR_SOURCES_COMMON)/convert_shf_to_grib2.F90 write_grib2_data.o
	$(FC) -c $(FCFLAGS) $<

write_grib2_data.o: $(DIR_SOURCES_COMMON)/write_grib2_data.F90
	$(FC) -c $(FCFLAGS) $<

clean:
	rm -f convert_shf_to_grib2 *.o *.mod
