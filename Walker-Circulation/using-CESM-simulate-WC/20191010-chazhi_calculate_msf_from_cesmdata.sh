#!/bin/sh
#-----------------------------------------------
#   This is a shell script for calculating the
# mass streamfunction by using some NCL script, 
# and mainly focus on CESM ouput data,this shell script call
# ncl function script : 
#2019925-CESM-data-chazhi.ncl (function:cesm_hybird2pre(infilepath,outfilepath))
#2019926-mass_streamfunciton_cesm.ncl(functuon:get_msf(path,filename,timestart,timeend,outpath,outputname))
#step1:interpolate the CESM output data hybird level to pressure level ,step2:cauculate massstreamfuction
# You should set the basic parameters as below. 
# Good Luck!
#               creat on  2019-10-10
#
#               by Shuheng Lin  
#-----------------------------------------------

# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F_1850-PDAY/
STEP=3
echo test 
#step1 : merge the cesm data into a whole data by using cdo
  ## the prefix of data is usually CESM compet name ,alarm for time select
if  [ $STEP == 3 ] ; then
  echo test
    # cd $PRE_DIR_ORG 
    # cdo select,name=U,V,OMEGA,PRECL,PSL,PS F_1850-PDAY.cam.h0.* F_1850-PDAY.cam.h0.197901-201212.nc 
    # echo "merge cesm data done"
fi

#step2 : interpolate the data from hybird level to pressure level
 
if  [ $STEP == 3 ] ; then
    cd /home/ys17-19/lsh/Project/Walker-Circulation/using-CESM-simulate-WC/ 
    ncl  inpath='"~/lsh/CESM-data/f19_f19_FAMIP_rmindian/f19_f19_FAMIP_rmindian.cam.h0.197901-200611.nc"' \
      outpath='"~/lsh/CESM-data/f19_f19_FAMIP_rmindian/f19_f19_FAMIP_rmindian.cam.h0.197901-200611_chazhi.nc"' \
       ./2019925-CESM-data-chazhi.ncl 
   echo "finish step2"
fi
# step3 calculate mass stream function

if  [ $STEP == 3 ] ; then
     ## 斜杠用来添加"，不然"会被默认成输出变量用的符号 
   ncl path=\"${PRE_DIR_ORG}\"               \
       filename='"F_1850-PDAY.cam.h0.197901-201212.nc"' \
     ./2019926-mass_streamfunciton_cesm.ncl 
   echo "finish step3"
fi




