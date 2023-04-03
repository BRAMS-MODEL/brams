#!/bin/bash
#script to get GFS data from nomads
if [[ $1 = "-h" || $1 = "--help" ]]
then

echo getGFS.sh YYYY MM DD DIR
echo
echo where
echo
echo YYYY - Year with 4 digits
echo   MM - Month with 2 digits
echo   DD - Day with 2 digits
echo  DIR - Output directory
echo
echo example:
echo
echo getGFS.sh 2020 09 15 ./datain/GFS
echo

exit 0

else

wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f000 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f030 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f036 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f042 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f048 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f006 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f012 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f018 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f024 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f030 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f036 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f042 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f048 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f054 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f060 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f066 
wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.$1$2$3/00/atmos/gfs.t00z.pgrb2.0p25.f072 
                                                                        atmos/
mv gfs.t00z.pgrb2.0p25.f000 $4/gfs.t00z.pgrb2.0p25.f000.$1$2$3.grib2                                        
mv gfs.t00z.pgrb2.0p25.f030 $4/gfs.t00z.pgrb2.0p25.f030.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f036 $4/gfs.t00z.pgrb2.0p25.f036.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f042 $4/gfs.t00z.pgrb2.0p25.f042.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f048 $4/gfs.t00z.pgrb2.0p25.f048.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f006 $4/gfs.t00z.pgrb2.0p25.f006.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f012 $4/gfs.t00z.pgrb2.0p25.f012.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f018 $4/gfs.t00z.pgrb2.0p25.f018.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f024 $4/gfs.t00z.pgrb2.0p25.f024.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f030 $4/gfs.t00z.pgrb2.0p25.f030.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f036 $4/gfs.t00z.pgrb2.0p25.f036.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f042 $4/gfs.t00z.pgrb2.0p25.f042.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f048 $4/gfs.t00z.pgrb2.0p25.f048.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f054 $4/gfs.t00z.pgrb2.0p25.f054.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f060 $4/gfs.t00z.pgrb2.0p25.f060.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f066 $4/gfs.t00z.pgrb2.0p25.f066.$1$2$3.grib2
mv gfs.t00z.pgrb2.0p25.f072 $4/gfs.t00z.pgrb2.0p25.f072.$1$2$3.grib2

fi