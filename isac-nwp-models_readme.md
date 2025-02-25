# NWP models Bolam-Moloch-Globo (CNR-ISAC)

 Instruction for the installation, compilation and implementation
    of atmospheric numerical models software developed in CNR-ISAC:
    Numerical Weather Prediction models: Bolam, Moloch, Globo
    and software for graphical presentation of forecast outputs in grib2 format
    using PLPlot free software

(25/10/2023, Oxana Drofa, o.drofa@isac.cnr.it)

---------------------------------------------------------------------
## NWP models description: https://www.isac.cnr.it/dinamica/projects/forecasts/
---------------------------------------------------------------------
##                    Package version 24.0
---------------------------------------------------------------------

Software requested:

- fortran compiler for Linux: gfortran or ifort;
- mpi compiler for fortran compiler: mpif90 or mpiifort; 
- eccodes (ECMWF free SW);
- rad_ecwmf (ECMWF radiation processes scheme), optionally;
- plplot, libraries for fortran compiler, optionally.

## Installation
---------------------------------------------------------------------

This procedure allows to create the basic auxiliary files for preprocessing and model running procedures.

Please, read and correct the install.script file.
Execute the installation:

./install.script

The following files are created:

./bolam/dimensions.inc
./bolam/run_premodel/prebolam.inp
./bolam/run_model/bolam.inp

./moloch/dimensions.inc
./moloch/run_premodel/premoloch.inp
./moloch/run_model/moloch.inp

./globo/dimensions.inc
./globo/run_premodel/preglobo.inp
./globo/run_model/bolam.inp

./bolam/run_postmodel/grib_sample.inp
./bolam/run_convert_shf_to_grib2/grib_sample.inp
./moloch/run_postmodel/grib_sample.inp
./moloch/run_convert_shf_to_grib2/grib_sample.inp
./globo/run_postmodel/grib_sample.inp
./globo/run_convert_shf_to_grib2/grib_sample.inp

./model_chain_launch.script

Please, check and correct ALL parameters in the created files.

Use ./sources/bolam/postbolam_sample.inp and ./sources/moloch/postmoloch_sample.inp
as an example for creating auxiliary files for postprocessing procedures, that must be located in:
./bolam/run_postmodel/postbolam.inp
./moloch/run_postmodel/postmoloch.inp

If you want to create postprocessing output (for some meteorological parameters)
interpolated on defined geographical points, then you must do
the following steps:

- create auxiliary file with parameters (name, latitude (deg.),
longitude (deg.) and altitude (m)) of output points
./bolam/run_meteograms/geopoint_meteogram.txt or
./moloch/run_meteograms/geopoint_meteogram.txt using the sample file
./sources/common/geopoint_meteogram_sample.txt

- create auxiliary file with parameters of procedure
./bolam/run_meteograms/meteogram_grib2.inp or
./moloch/run_meteograms/meteogram_grib2.inp using sample file
./sources/common/meteogram_grib2_sample.inp:
it is possible to interpolate with bilinear interpolation method
(INTERP_GAUSS=0) or Gaussian interpolation method (INTERP_GAUSS=1)
with semi-radius defined by GAUSS_SEMI_RADIUS (in km),
and gauss function parameter defined by GAUSS_PAR (recommended 1).

Attention: the depth of soil levels are defined by the variable SOIL_LEV0
in ./bolam/run_premodel/prebolam.inp (or by SLT in ./moloch/run_premodel/premoloch.inp,
if you want use different soil levels in Bolam and in Moloch).
The number of defined soil levels must be equal to the parameter NLEVSOIL in
./install.script that defines parameter NLEVG in
./bolam/dimension.inc and ./moloch/dimension.inc.

Use the script ./clean_istallation.script to purge ALL files that have been created 
by installation, compilation and execution procedures.

---------------------------------------------------------------------
## Compilation
---------------------------------------------------------------------

You may use Makefile_gfortran or Makefile_ifort to compile:
cp  Makefile_gfortran  Makefile
or
cp  Makefile_ifort  Makefile

Please, set the parameters of compilation in Makefile:

Choice of sowfware to compile:

USE_MPI - MPI (YES or NO)
USE_RAD_ECMWF - ECMWF radiation processes scheme (YES or NO)
BOLAM - Bolam model (YES or NO)
MOLOCH - Moloch model (YES or NO)
GLOBO - Globo model (Bolam model on global domain) (YES or NO)
GRAPHICS_PLGRIBFA - Graphics package Plgribfa (YES or NO)

OMPI_FC - fortran compiler (gfortran or ifort)
FC – fortran compiler (gfortran or ifort)
FCFLAGS – compilation options
FC_MPI - fortran compiler for MPI
FCFLAGS_MPI - MPI compilation options
DIR_LIB – complete path of directory with Eccodes and ECMWF-Radiation installations
DIR_ECCODES -  complete path of directory with Eccodes installations
DIR_RAD -  complete path of directory with ECMWF-Radiation installation
DIR_SOURCES_BOLAM -  complete path of directory with Bolam fortran codes and auxiliary file samples
DIR_SOURCES_MOLOCH -  complete path of directory with Moloch fortran codes and auxiliary file samples
DIR_SOURCES_GLOBO -  complete path of directory with Globo (only preprocessing code) fortran codes and auxiliary file samples
DIR_SOURCES_COMMON -  complete path of directory with common and auxiliary fortran codes and various auxiliary file samples

Compilation will be performed in the directories defined by parameter SUBDIRS in depending on chosen SW:
bolam/executable_premodel
bolam/executable_model
bolam/executable_postmodel
bolam/executable_convert_shf_to_grib2
bolam/executable_meteograms
moloch/executable_premodel
moloch/executable_model
moloch/executable_postmodel
moloch/executable_convert_shf_to_grib2
moloch/executable_meteograms
globo/executable_premodel
globo/executable_model
globo/executable_postmodel
globo/executable_convert_shf_to_grib2
globo/executable_meteograms
executable_plgribfa

each of these directories contains its proper Makefile.

For compilation of all executables use the command:
make

for removal of all executables (*.o and *.mod also) use the command:
make clean

---------------------------------------------------------------------
## Implementation
---------------------------------------------------------------------

This part is useful only for users without previous experience,
and includes some explanations of execution scripts.

Basic models chain script is ./model_chain_launch.script

All implementation parameters are defined in this script. The script launches
other scripts bolam_forecast_chain.script and moloch_forecast_chain.script
for the corresponding model chain execution.

Please, read and follow the instruction included in this file.
Attention, if basic flag of the model (FLAG_BOLAM or FLAG_MOLOCH) is equal zero,
then all other flags of this model are off.

Below some specifying information is present.

DIR_GEO_DATA – complete path to directory including all physicogeographical dataset. It must include the following files:

orogr_global_latlon_1km.bin
soil_fao_global_latlon_8km.bin
worldexp.dat
veget_global_latlon_1km.bin
vegtype_GLC2000_global_latlon_1km.bin
ecmwf_glob_0_25_lsm.bin
ecmwf_ifs_1979_2014_t2_clim.bin
fveg_global_latlon_1km_day005.bin
fveg_global_latlon_1km_day015.bin
fveg_global_latlon_1km_day025.bin
fveg_global_latlon_1km_day036.bin
fveg_global_latlon_1km_day046.bin
fveg_global_latlon_1km_day056.bin
fveg_global_latlon_1km_day064.bin
fveg_global_latlon_1km_day074.bin
fveg_global_latlon_1km_day084.bin
fveg_global_latlon_1km_day095.bin
fveg_global_latlon_1km_day105.bin
fveg_global_latlon_1km_day115.bin
fveg_global_latlon_1km_day125.bin
fveg_global_latlon_1km_day135.bin
fveg_global_latlon_1km_day145.bin
fveg_global_latlon_1km_day156.bin
fveg_global_latlon_1km_day166.bin
fveg_global_latlon_1km_day176.bin
fveg_global_latlon_1km_day186.bin
fveg_global_latlon_1km_day196.bin
fveg_global_latlon_1km_day206.bin
fveg_global_latlon_1km_day217.bin
fveg_global_latlon_1km_day227.bin
fveg_global_latlon_1km_day237.bin
fveg_global_latlon_1km_day248.bin
fveg_global_latlon_1km_day258.bin
fveg_global_latlon_1km_day268.bin
fveg_global_latlon_1km_day278.bin
fveg_global_latlon_1km_day288.bin
fveg_global_latlon_1km_day298.bin
fveg_global_latlon_1km_day309.bin
fveg_global_latlon_1km_day319.bin
fveg_global_latlon_1km_day329.bin
fveg_global_latlon_1km_day339.bin
fveg_global_latlon_1km_day349.bin
fveg_global_latlon_1km_day359.bin
lai_global_latlon_1km_day005.bin
lai_global_latlon_1km_day015.bin
lai_global_latlon_1km_day025.bin
lai_global_latlon_1km_day036.bin
lai_global_latlon_1km_day046.bin
lai_global_latlon_1km_day056.bin
lai_global_latlon_1km_day064.bin
lai_global_latlon_1km_day074.bin
lai_global_latlon_1km_day084.bin
lai_global_latlon_1km_day095.bin
lai_global_latlon_1km_day105.bin
lai_global_latlon_1km_day115.bin
lai_global_latlon_1km_day125.bin
lai_global_latlon_1km_day135.bin
lai_global_latlon_1km_day145.bin
lai_global_latlon_1km_day156.bin
lai_global_latlon_1km_day166.bin
lai_global_latlon_1km_day176.bin
lai_global_latlon_1km_day186.bin
lai_global_latlon_1km_day196.bin
lai_global_latlon_1km_day206.bin
lai_global_latlon_1km_day217.bin
lai_global_latlon_1km_day227.bin
lai_global_latlon_1km_day237.bin
lai_global_latlon_1km_day248.bin
lai_global_latlon_1km_day258.bin
lai_global_latlon_1km_day268.bin
lai_global_latlon_1km_day278.bin
lai_global_latlon_1km_day288.bin
lai_global_latlon_1km_day298.bin
lai_global_latlon_1km_day309.bin
lai_global_latlon_1km_day319.bin
lai_global_latlon_1km_day329.bin
lai_global_latlon_1km_day339.bin
lai_global_latlon_1km_day349.bin
lai_global_latlon_1km_day359.bin

All this dataset files are available in https://www.isac.cnr.it/dinamica/geo_dataset/

DIR_BIN - complete path to eccodes tools directory
DIR_SOURCES_COMMON - complete path to common codes directory (used for meteograms elaboration only)

---------------------------------------------------------------------

Definition of parameters for Bolam model chain implementation 

DH_START_BOLAM=0        - Only for the case of input data created by Globo model: lag (hour) of Bolam run start with respect to Globo run start
INIDATE_BOLAM      - YYYYMMDDHH string of date and time of Bolam run start
HRUN_BOLAM           -  Hours of model run duration
DHINPUT_BOLAM     - Interval (in hours) between two boundary conditions
DHOUTPUT_BOLAM_MAIN       - Interval (in hours) between two writings of model main output (mhf)
NJUMP_BOLAM_MAIN                - Jump (time number) of mhf elaboration by postprocessing
DHOUTPUT_BOLAM_ADD         - Interval (in hours) between two writings of model additional output (shf) (surface fields)
NJUMP_BOLAM_ADD                  - Jump (time number) of shf elaboration by conversion to grib2 format

HOUR_LIMIT_PRE_BOLAM        - Time limit* (in hours) of preprocessing execution, depends on input data waiting
HOUR_LIMIT_BOLAM                 - Time limit* (in hours) of model run execution, depends on input data waiting and integration model time
HOUR_LIMIT_POST_BOLAM     - Time limit* (in hours) of postporoceesing run execution, depends on input data waiting and integration model time

*Time limit – when procedure execution time will be superior to this time limit, the procedure will stop.

DIR_INPUT_DATA_BOLAM              - Complete path to directory with input data for Bolam preprocessing (for definition of initial and boundary conditions)
DIR_PREMODEL_DATA_BOLAM    -  Complete path to directory with output data of Bolam preprocessing (the same input data of Bolam run), file (model_param_constant.bin) with time constant fields is contained in this directory 
DIR_MODEL_DATA_MHF_BOLAM - Complete path to directory with mhf_atm and mhf_soil Bolam run data
DIR_MODEL_DATA_MHF_SOIL_COPY_BOLAM -  Complete path to directory with copies of some mhf_soil files that will be used in a successive forecast run
DIR_MODEL_DATA_SHF_BOLAM                          -  Complete path to directory with shf Bolam run data
DIR_POSTMODEL_DATA_GRIB2_BOLAM            -  Complete path to directory with postprocessing output data (in grib2 format)
DIR_POSTMODEL_SHF_DATA_GRIB2_BOLAM   -  Complete path to directory with shf output data converted to grib2 format
DIR_METEOGRAMS_BOLAM                                   -  Complete path to directory with output data (ascii format) of some desired meteorological parameter interpolated to prescribed points using input data in grib2 format

FILE_FLAG_BOLAM_RUN_FINISH      -   file that plays the role of flag to indicate that Bolam run is finished
FILE_FLAG_BOLAM_RUN_ERROR       -    file that plays the role of flag to indicate that Bolam run is crashed
FILE_FLAG_PREBOLAM_FINISH       -    file that plays the role of flag to indicate that Bolam preprocessing run is finished


----------------------------------
## Bolam chain launching
----------------------------------

./bolam_forecast_chain.script   launches the following scripts:

- preprocessing

./bolam/run_premodel/premodel_gfs_data_run.script or
./bolam/run_premodel/premodel_ifs_run.script or
./bolam/run_premodel/premodel_globo_data_run.script depending on defined
input data origin, these scripts are ready to use with input data files with
conventional names:
GFS_YYYYMMDDHH+000.grib2, GFS_YYYYMMDDHH+003.grib2..., or  IFS_YYYYMMDDHH+000.grib2, IFS_YYYYMMDDHH+003.grib2…,
or globolam_atm_YYYYMMDDHH_001.mhf, globolam_soil_YYYYMMDDHH_001.mhf, globolam_atm_YYYYMMDDHH_002.mhf, globolam_soil_YYYYMMDDHH_002.mhf, ..., model_param_constant_zoom.bin,
or IFS input data there is also an option of analysis (re-analysis) data
(not forecast) for defining boundary conditions, if input data interval
is more then 3 hours, then this option is active, in this case the conventional
input data files names are:  IFS_YYYYMMDDHH1+000.grib2, IFS_YYYYMMDDHH2+000.grib2, IFS_YYYYMMDDHH3+000.grib2, ...
if you use other names for input data file, you must correct the corresponding script;

- model run execution

./bolam/run_model/bolam_run.script
Attention!!! Please, control and correct launching command (mpi),
search the string “LAUNCH” in the script;
the script launches two more scripts ./bolam/run_model/mv_mhf.script and
./bolam/run_model/mv_shf.script to move all output files to defined directories,
these scripts are useful in case of separately written instants only,
script ./bolam/run_model/bolam_run.script deletes all input local files
(that are symbolic link files) after the model run finishes;

- postprocessing

./bolam/run_postmodel/postmodel_run.script is the script for the case of separately
written instants in mhf output, the script is launched together with model run,
it waits instant mhf (atm and soil), elaborates it, waits the next mhf output, etc.;
the script moves postprocessing output in grib2 format to defined directory;

Attention!!! Please, control and correct file with postprocessing options
./bolam/run_postmodel/postbolam.inp;
the script uses input file grib_sample.inp also (see above);
 
- conversion of shf model output to grib2 format

./bolam/run_convert_shf_to_grib2/convert_shf_to_grib2_run.script this script
is analogue of ./bolam/run_postmodel/postmodel_run.script,
it is useful in case of separately written instants (in shf) only, it does not require
any configuration file, it just uses grib_sample.inp (see above);

- output data for a defined set of meteorological parameters interpolated to geographical points

./bolam/run_meteograms/read_grib2_write_meteogram.script and ./bolam/run_meteograms/meteogram_rd_wr.script
the list of output geographical points is defined by
./bolam/run_meteograms/geopoint_meteogram.txt file
(if model domain does not include some geographical points, then it is not a problem),
options of interpolation and set of output meteorological parameters are defined
by ./bolam/run_meteograms/meteogram_grib2.inp, this input file is used by
the executable called by ./bolam/run_meteograms/read_grib2_write_metrogram.script,
in this script the variable PARAM_LIST must be in agreement with the defined
meteorological parameters set;
/bolam/run_meteograms/read_grib2_write_meteogram.script is launched together
with model run and postprocessing elaboration, it waits instant output in grib2 format,
elaborates it, waits for the next output, etc., it writes output in ascii files,
each file contains data of one meteorological parameter in all points for one instant;
./bolam/run_meteograms/meteogram_rd_wr.script waits when
./bolam/run_meteograms/read_grib2_write_meteogram.script has elaborated all output
instants, reads data for all instants, writes the same in ascii format files;
one output file contains data in one geographical point for a single meteorological
parameter for all forecast instants;
in this script the variable PARAM_LIST must be in agreement with the defined
meteorological parameters set;

- cleaning of input data, cleaning or archiving of output data

./bolam/file_clean.script
Please, control with a great attention all the commands in this script;
at the moment, the period of “operational data” conservation is n=2 days,
see parameters “DAY” and “HOUR”;
the script performs the following operations:
- deletes all preprocessing output data, excluding the constant (in time) model fields;
- copies mhs_soil model run ouput file, that may be used in a successive forecast
model run: soil, vegetation and snow cover parameters,
that have been simulated by the model run and written in one output file in
mhf_soil format, may be used for the initialization of these parameters in a successive
model run by the preprocessing procedure that can read the mhf_soil format with a useful
date and time of validity, so, some mhf_soil model output must be conserved
in the special directory with a conventional name
(see variable “FILE_NAME_COPY” in ./bolam/file_clean.script), actually the copy
of this file is done each 24 hour for model run launching one time a day;
- deletes all model output data in mhf_atm and mhf_soil format that have been
created n days ago, deletes work file-flag (ascii format) of atm_mhf and
soil_mhf output writing that have been created by current model run;
- deletes all model output data in shf format that have been created n days ago,
deletes work file-flag (ascii format) of shf output that have been created by current model run;
- archives in compressed form (tar+gzip=tgz) shf output data converted in
grib2 format that have been created n days ago;
- deletes all output files created by meteograms procedures that have been created n days ago.

----------------------------------
## Moloch chain launching
----------------------------------
It is analogue of Bolam chain launching

./moloch_forecast_chain.script   launches the following scripts:

- preprocessing

./moloch/run_premodel/premodel_bolam_data_run.script, this scripts is ready
to be used with input data files with conventional names:
bolam_soil_YYYYMMDDHH_001.mhf, bolam_atm_YYYYMMDDHH_002.mhf, bolam_soil_YYYYMMDDHH_002.mhf, ..., model_param_constant.bin (of Bolam model).
If you use other names for input data file, you must correct the corresponding script;

- model run execution

./moloch/run_model/bolam_run.script
Attention!!! Please, control and correct launching command (mpi), search the string “LAUNCH” in the script;
the script launches two more scripts ./moloch/run_model/mv_mhf.script and
./moloch/run_model/mv_shf.script to move all output files to defined
directories, these scripts are useful for the case of separately written instants
only, the case of writing of mhf output with a reduced resolution
(writing every two grid points), that may be opted by the parameter
“mhfr” in moloch.inp input file, requires some special operations:
1) writing of constant (in time) model fields with the same (not full)
resolution and link of this file to directory with mhf output archive,
2) complementary writing of soil parameters output (soil_mhf) with full
resolution for use in the initialization of soil (snow and vegetation) parameters
in a successive model run and link of this file(s) to directory with mhf output archive,
script ./moloch/run_model/moloch_run.script deletes all input local files
(that are symbolic link files) after the model run finishes;

- postprocessing

./moloch/run_postmodel/postmodel_run.script is the script for the case of separately
written instants in mhf output, the script is launched together with model run,
it waits the mhf for an instant (atm and soil), elaborates it, waits the next mhf output, etc.;
the script moves postprocessing output in grib2 format to a defined directory;
Attention!!! Please, control and correct file with postprocessing options
./moloch/run_postmodel/postmoloch.inp;
the script uses also input file grib_sample.inp (see above);
 
- conversion of shf model output to grib2 format

./moloch/run_convert_shf_to_grib2/convert_shf_to_grib2_run.script this script
is analogue of ./moloch/run_postmodel/postmodel_run.script,
it is useful in case of separately written instants (in shf) only, it does not require
any configuration file, it just uses grib_sample.inp (see above);

- output data for defined set of meteorological parameters interpolated to geographical points

./moloch/run_meteograms/read_grib2_write_metrogram.script and
./moloch/run_meteograms/meteogram_rd_wr.script
the list of output geographical points is define by ./moloch/run_meteograms/geopoint_meteogram.txt file
(if model domain does not include some geographical points, then it is not problem),
options of interpolation and set of output meteorological parameters are defined by
./moloch/run_meteograms/meteogram_grib2.inp, this input file is used by
the executable called by ./moloch/run_meteograms/read_grib2_write_metrogram.script,
in this script the variable PARAM_LIST must be in agreement with the defined
meteorological parameters set;
/moloch/run_meteograms/read_grib2_write_meteogram.script is launched together
with model run and postprocessing elaboration, it waits instant output in
grib2 format, elaborates it, waits the next output, etc.,
it writes output in ascii files, each file contains data of one meteorological
parameter in all points for one instant;
./moloch/run_meteograms/meteogram_rd_wr.script waits when
./moloch/run_meteograms/read_grib2_write_meteogram.script has elaborated all output
instants, reads data for all instants, writes the same in ascii format files;
one output file contains data in one geographical point for a single
meteorological parameter for all forecast instants; in this script
the variable PARAM_LIST must be in agreement with the defined meteorological parameters set;

- cleaning of input data, cleaning or archiving of output data

./moloch/file_clean.script
Please, control with a great attention all the commands in this script;
at the moment, the period of “operational data” conservation is n=2 days, see parameters “DAY” and “HOUR”;
the script performs the following operations:
- deletes all preprocessing output data, excluding the constant (in time) model fields;
- copies mhs_soil model run ouput file, that may be used in a successive forecast
model run: soil, vegetation and snow cover parameters, that have been simulated
by the model run and written in one output file in mhf_soil format,
may be used for initialization of these parameters in a successive model run
by the preprocessing procedure that can read the mhf_soil format with a useful
date and time of validity, so, some mhf_soil model output must be conserved
in the special directory with a conventional name (see variable “FILE_NAME_COPY”
in ./moloch/file_clean.script), actually the copy of this file is done each
24 hour for model run launching one time a day;
- deletes all model output data in mhf_atm and mhf_soil format that have been
created n days ago, deletes work file-flag (ascii format) of atm_mhf and soil_mhf
output writing that have been created by current model run;
- deletes all model output data in shf format that have been created n days ago
deletes work file-flag (ascii format) of shf output that have been created by current model run;
- archives in compressed form (tar+gzip=tgz) shf output data converted in grib2 format that have been created n days ago;
- deletes all output file created by meteograms procedures that have been created n days ago.
