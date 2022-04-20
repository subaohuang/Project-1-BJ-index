#!bin/bash

DAILY_DIR=/home/yangsong3/data-observation/ERA5-daily/Multi-levels/
MONTH_DIR=/home/yangsong3/data-observation/ERA5-monthly/pressure/

varlist=(uwind vwind temp zg q)
size=${#varlist[*]}
size=`expr $size + 1`
# echo $size

yrStrt=1979
yrLast=2019
nyear=`expr $yrLast - $yrStrt + 1`
# echo $nyear

month=(01 02 03 04 05 06 07 08 09 10 11 12)
# echo ${month[11]}

for ((ll=0;ll<$size;ll++)) #
do
	cd $DAILY_DIR/${varlist[ll]}
	for (( yy = $yrStrt; yy < $yrLast; yy++ )); do #
		for (( i = 0; i < 12; i++ )); do
			cdo monmean ${varlist[ll]}.${yy}-${month[i]}.daily.nc ${MONTH_DIR}${varlist[ll]}.${yy}-${month[i]}.monthly.nc
		done
	done
	cd $MONTH_DIR
	cdo cat ${varlist[ll]}.*.monthly.nc ${varlist[ll]}.monthly.197901-201812.nc
	rm -f ${varlist[ll]}.*.monthly.nc
done