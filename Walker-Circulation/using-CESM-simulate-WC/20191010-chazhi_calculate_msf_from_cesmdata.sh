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
PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F_1850-PDAY-rmatlantic/

STEP=3
# step1 : merge the cesm data into a whole data by using cdo
  ### the prefix of data is usually CESM compet name ,alarm for time select
if  [ ! -e  ${PRE_DIR_ORG}F_1850-PDAY-rmatlan.cam.h0.197801-201212.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    rm F_1850-PDAY-rmatlan.cam.h0.197801-201212.nc
    cdo select,name=U,V,OMEGA,PRECL,PSL,PS F_1850-PDAY-rmatlan.cam.h0.* F_1850-PDAY-rmatlan.cam.h0.197801-201212.nc 
fi

# step2 : interpolate the data from hybird level to pressure level
  if [ ! -e ${PRE_DIR_ORG}F_1850-PDAY-rmatlan.cam.h0.197801-201212_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
    echo "don't exit chazhi file, procecing..."
    cd /home/ys17-19/lsh/Project/Walker-Circulation/using-CESM-simulate-WC/ 
    ncl  inpath='"~/lsh/CESM-data/F_1850-PDAY-rmatlantic/F_1850-PDAY-rmatlan.cam.h0.197801-201212.nc"' \
         outpath='"~/lsh/CESM-data/F_1850-PDAY-rmatlantic/F_1850-PDAY-rmatlan.cam.h0.197801-201212_chazhi.nc"' \
       ./190925-CESM-data-chazhi.ncl
    echo "finish CESM chazhi"
  fi 
# step3 calculate mass stream function

  if  [ ! -e /home/lsh/data/wc-result/msf_F_1850-PDAY-rmatlantic_197912-201212.nc ] ; then
     ## 斜杠用来添加"，不然"会被默认成输出变量用的符号 
   echo "don't exit msf file, procecing..."
   ncl inpath=\"${PRE_DIR_ORG}\"               \
       filename='"F_1850-PDAY-rmatlan.cam.h0.197801-201212_chazhi.nc"' \
       outputpath='"~/lsh/data/wc-result/"' \
       outputname='"msf_F_1850-PDAY-rmatlan_197912-201212.nc"' \
     ./190926-mass_streamfunciton_cesm.ncl
   echo "finish CESM mass_stream function"
  fi

#-----------------------------------------------------------



