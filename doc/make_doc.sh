#!/bin/bash
locpath=`pwd`
echo ${locpath}
cp ../src/brams/model/rammain.f90 ./src
cp ../src/utils/include/i8.h ./src/mpif.h
cp /opt/gnu/include/netcdf.inc ./src
$1 ${locpath}/brams.md
