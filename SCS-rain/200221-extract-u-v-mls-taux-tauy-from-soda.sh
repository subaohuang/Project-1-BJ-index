#!/bin/sh
#-----------------------------------------------
#This is a shell script for extracting the U,V,OMEGA,MLS,T data 
#in soda data
#
# Good Luck!
#               creat on  2020-02-21
#
#               by Shuheng Lin  
#----------------------------------------------- 
# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
PRE_DIR_ORG=/home/yangsong3/data-observation/SODA/3.3.1/
NOW_DIR_ORG=/home/yangsong3/data-observation/SODA/3.3.1/
modelname=soda3.3.1
# modelname=cmip_10_ctl

variable_total=(temp u v wt mlt taux tauy)     #often use
# variable_total=(tauy temp)


for i in `seq 0 6` ;
 do 
 variable=${variable_total[i]}
 echo $variable
# step1 : merge the cesm data into a whole data by using cdo
  ### the prefix of data is usually CESM compet name ,alarm for time select
  if  [ ! -e  ${NOW_DIR_ORG}${modelname}_${variable}_mn_1980-2015.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    cdo select,name=${variable} ${modelname}_mn_ocean_reg_*  ${NOW_DIR_ORG}${modelname}_${variable}_mn_1980-2015.nc
  fi

done;