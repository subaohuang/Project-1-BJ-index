#!/bin/sh
#-----------------------------------------------
#This is a shell script for extracting the U,V,OMEGA,MLS,T data 
#in soda data
#
# Good Luck!
#               creat on  2020-03-08
#
#               by Shuheng Lin  
#----------------------------------------------- 
# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
PRE_DIR_ORG=/home/yangsong3/data-observation/SODA/2.2.4/
NOW_DIR_ORG=/home/yangsong3/data-observation/SODA/2.2.4/
modelname=soda2.2.4
# modelname=cmip_10_ctl

variable_total=(temp u v w taux tauy)     #often use
# variable_total=(tauy temp)


for i in `seq 0 5` ;
 do 
 variable=${variable_total[i]}
  # step2 : interpolate the data from hybird level to pressure level
  if [ ! -e ${NOW_DIR_ORG}${modelname}_mn_ocean_reg_${variable}_1957-2008_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
    echo "don't exit chazhi file, procecing..."
    cd /home/ys17-19/lsh/Project/SCS-rain/
    pwd
    ncl  -nQ inpath=\"${NOW_DIR_ORG}${modelname}_mn_ocean_reg_${variable}_1957-2008.nc\" \
         outpath=\"${NOW_DIR_ORG}${modelname}_mn_ocean_reg_${variable}_1957-2008_1dgree.nc\" \
         var_need=\"${variable}\" \
       ./200308-soda2.2.4-chazhi.ncl
    echo "finish soda chazhi"
  fi 
  
done;


