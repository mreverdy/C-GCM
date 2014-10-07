#!/bin/sh
unset LANG
ulimit -s unlimited

house=`pwd`
home='/homedata/mreverdy/GOCCP'
#-*-makefile-*-

### Specify your compiler
IFORT="ifort"
G95="pgf90"
COMPILO=$IFORT

### This header file is automatically included in the secondary Makefiles.
### Please tune it to your own installation

### Specify where the headers and libraries of your netCDF package reside.
# Example : 
#   if the file libnetcdf.a is located at
#   /opt/netcdf-3.5.1-IFORT-64/lib/libnetcdf.a
#   then NETCDFLIB=/opt/netcdf-3.5.1-IFORT-64/lib
#
#   if the file netcdf.mod is located at
#   /opt/netcdf-3.5.1-IFORT-64/include/netcdf.mod
#   then NETCDFINC=/opt/netcdf-3.5.1-IFORT-64/include
#
# To avoid trouble, netCDF should have been compiled with the 
# same compiler you use to build CHIMERE
# In most Linux distributions, netCDF has been compiled using g77.
# This may not be compatible with the f90 compilers required for CHIMERE.
#
NETCDFLIB=/opt/netcdf/$COMPILO/lib
NETCDFINC=/opt/netcdf/$COMPILO/include

#NETCDFLIB=/homedata/gcesana/local/ifort/netcdf-4.1.1/lib
#NETCDFINC=/homedata/gcesana/local/ifort/netcdf-4.1.1/include

### Specify where the headers and libraries of your HDF package reside.
HDFLIB=/usr/lib64/hdf
HDFINC=/usr/include/hdf



### Specify the filename
NAME=$1
instant=$2
nol=$3

### Clean the execute file
#if [ $NAME.e -e ]
#then
#./rm -f $NAME.e
#fi

#ifort $1.tmp.f90 calendar.f90 -I/usr/include/hdf -L/usr/lib64/hdf -lmfhdf -ldf -ljpeg -lz -I/opt/netcdf/ifort/include/ -L/opt/netcdf/ifort/lib -lnetcdf -o $1.e
# -I${NETCDFINC} -L${NETCDFLIB} -l${netcdf} -o ${progname}.e

### Choose your execution mode { PROD | DEVEL }
### PROD is fast, DEVEL allows for more checking and error tracking
MODE="PROD"

### If you use the Fedora Core 4 GNU/Linux distribution, you may
#   experience problems with ifort and Interprocedural Optimisation. 
#   If this is the case, you should disable it.
#   Otherwise just comment out the following line to get the maximum
#   performance of CHIMERE.
#FC4_BUG = -no-ipo
FC4_BUG=-no-ipo



if [ $MODE == "DEVEL" ]
then
echo "DEVEL"
# For debug/development
F90FLAGS1="-I${HDFINC} -L${HDFLIB} -lmfhdf -ldf -ljpeg -lz -I${NETCDFINC} -L${NETCDFLIB} -lnetcdf -fpe0 -fpp -ftrapuv -g -traceback -check all $FC4_BUG"
fi

if [ $MODE == "PROD" ]
then
# for production
F90FLAGS1="-I${HDFINC} -L${HDFLIB} -lmfhdf -ldf -ljpeg -lz -I${NETCDFINC} -L${NETCDFLIB} -lnetcdf" 
fi

VERSION=`echo $NAME | cut -c7-10`
 
if [ -z $instant ]; then

sed s/Prog_version/$NAME/g $NAME.f90 > $NAME.tmp.f90 

else

if (( $nol == 1)); then
nol2="nol"
sed s/Prog_version/$NAME/g $NAME.f90 | sed s/Num_version/$VERSION$nol2/ | sed s/instant_switch/$instant/g | sed s/nol_switch/$nol/ > $NAME.tmp.f90 

else
nol2=""
sed s/Prog_version/$NAME/g $NAME.f90 | sed s/Num_version/$VERSION$nol2/ | sed s/instant_switch/$instant/g | sed s/nol_switch/$nol/ > $NAME.tmp.f90 

fi

fi

echo ${COMPILO} $NAME.tmp.f90 calendar.f90 $F90FLAGS1 -o ${NAME}.e

${COMPILO} $NAME.tmp.f90 $home/calendar.f90 $F90FLAGS1 -o ${NAME}.e

#ifort $1.tmp.f90 calendar.f90 -I/usr/include/hdf -L/usr/lib64/hdf -lmfhdf -ldf -ljpeg -lz -I/opt/netcdf/ifort/include/ -L/opt/netcdf/ifort/lib -lnetcdf -o $1.e

# Clean trash
rm -f $NAME.tmp.f90

