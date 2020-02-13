#!/bin/sh
#-----------------------------------------------
#   This is a shell script for calculating the
# mass streamfunction by using some NCL script, 
# and mainly focus on CESM ouput data,this shell script call
# ncl function script : 
#191209-CESM-data-chazhi.ncl (function:cesm_hybird2pre(infilepath,outfilepath))
#191209-mass_streamfunciton_cesm-F2000.ncl(functuon:get_msf(path,filename,timestart,timeend,outpath,outputname))
#step1:interpolate the CESM output data hybird level to pressure level ,step2:cauculate massstreamfuction
# You should set the basic parameters as below. 
# Good Luck!
#               creat on  2020-02-13
#
#               by Shuheng Lin  
#----------------------------------------------- 

# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F/F_2000_IPO/
# PRE_DIR_ORG=/home/yangsong3/data-model/MC-ztt/amip_10_noland/atm/monthly/
# PRE_DIR_ORG=/home/yangsong3/data-model/MC-ztt/cmip_10_noland/atm/monthly/
# NOW_DIR_ORG=/home/ys17-19/lsh/CESM-data/B/cmip_10_noland/
PRE_DIR_ORG=/home/yangsong3/data-model/MC-ztt/cmip_10_ctl/atm/monthly/
NOW_DIR_ORG=/home/ys17-19/lsh/CESM-data/B/cmip_10_ctl/

# modelname=cmip_10_noland
modelname=cmip_10_ctl


variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q    


# step1 : merge the cesm data into a whole data by using cdo
  ### the prefix of data is usually CESM compet name ,alarm for time select

  if  [ ! -e  ${NOW_DIR_ORG}${modelname}.cam.h0.197901-200512.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    cdo select,name=${variable} ${modelname}.cam.h0.*  ${NOW_DIR_ORG}${modelname}.cam.h0.197901-200512.nc
  fi

# step2 : interpolate the data from hybird level to pressure level
  if [ ! -e ${NOW_DIR_ORG}${modelname}.cam.h0.197901-200512_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
    echo "don't exit chazhi file, procecing..."
    cd /home/ys17-19/lsh/Project/Walker-Circulation/
    pwd
    ncl  -nQ inpath=\"${NOW_DIR_ORG}${modelname}.cam.h0.197901-200512.nc\" \
         outpath=\"${NOW_DIR_ORG}${modelname}.cam.h0.197901-200512_chazhi.nc\" \
       ./191209-CESM-data-chazhi.ncl
    echo "finish CESM chazhi"
  fi 
# step3 calculate mass stream function

  if  [ ! -e /home/ys17-19/lsh/data/wc-result/msf_${modelname}_197901-200512.nc ] ; then
     ## 斜杠用来添加"，不然"会被默认成输出变量用的符号 
   echo "don't exit msf file, procecing..."
   cd  /home/ys17-19/lsh/Project/Walker-Circulation/
   pwd
   ncl -nQ inpath=\"${NOW_DIR_ORG}\"               \
       filename=\"${modelname}.cam.h0.197901-200512_chazhi.nc\" \
       outputpath='"~/lsh/data/wc-result/"' \
       outputname=\"msf_${modelname}_197901-200512.nc\" \
     ./191209-mass_streamfunciton_cesm-F2000.ncl
   echo "finish CESM mass_stream function"
  fi
  
   echo "finish this script"
#-----------------------------------------------------------



