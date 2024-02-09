For conversation from grib1 format file to grib2 format file
it is necessary to integrated 2 tables.

Example for eccodes_ifort package:

cat eccodes_ifort/share/eccodes/definitions/grib1/localConcepts/ecmf/paramId.def ./grib1_add_paramId.def > table
mv table eccodes_ifort/share/eccodes/definitions/grib1/localConcepts/ecmf/paramId.def

cat eccodes_ifort/share/eccodes/definitions/grib2/localConcepts/ecmf/paramId.def ./grib2_add_paramId.def > table
mv table eccodes_ifort/share/eccodes/definitions/grib2/localConcepts/ecmf/paramId.def

Example for eccodes_gfortran package:

cat eccodes_gfortran/share/eccodes/definitions/grib1/localConcepts/ecmf/paramId.def ./grib1_add_paramId.def > table
mv table eccodes_gfortran/share/eccodes/definitions/grib1/localConcepts/ecmf/paramId.def

cat eccodes_gfortran/share/eccodes/definitions/grib2/localConcepts/ecmf/paramId.def ./grib2_add_paramId.def > table
mv table eccodes_gfortran/share/eccodes/definitions/grib2/localConcepts/ecmf/paramId.def
