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
PRE_DIR_ORG=/home/yangsong3/data-observation/SODA/3.4.2/
NOW_DIR_ORG=/home/yangsong3/data-observation/SODA/3.4.2/
modelname=soda3.4.2
# modelname=cmip_10_ctl

variable_total=(salt prho mlp temp u v wt mlt taux tauy)     #often use
# variable_total=(tauy temp)


for i in `seq 0 9` ;
 do 
 variable=${variable_total[i]}
# step1 : merge the cesm data into a whole data by using cdo
  ### the prefix of data is usually CESM compet name ,alarm for time select
  if  [ ! -e  ${NOW_DIR_ORG}${modelname}_${variable}_mn_1980-2015.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    cdo select,name=${variable} ${modelname}_mn_ocean_reg_*  ${NOW_DIR_ORG}${modelname}_${variable}_mn_1980-2015.nc
  fi

  # step2 : interpolate the data from hybird level to pressure level
  if [ ! -e ${NOW_DIR_ORG}${modelname}_${variable}_mn_1980-2015_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
    echo "don't exit chazhi file, procecing..."
    cd /home/ys17-19/lsh/Project/SCS-rain/
    pwd
    ncl  -nQ inpath=\"${NOW_DIR_ORG}${modelname}_${variable}_mn_1980-2015.nc\" \
         outpath=\"${NOW_DIR_ORG}${modelname}_${variable}_mn_1980-2015_chazhi.nc\" \
         var_need=\"${variable}\" \
       ./200306-soda3.4.2-chazhi.ncl
    echo "finish soda chazhi"
  fi 
  
done;


