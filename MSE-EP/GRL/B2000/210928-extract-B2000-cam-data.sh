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
#               creat on  2019-12-19
#
#               by Shuheng Lin  
#----------------------------------------------- 

# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
PRE_DIR_ORG=/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_alt_north_year_obsclm/
STEP=3
modelname=lsh_B2000_alt_north_year_obsclm
# variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q,UBOT,VBOT,TREFHT
variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q,TREFHT,T,
# variable=hyam


#step1 : merge the cesm data into a whole data by using cdo
  ## the prefix of data is usually CESM compet name ,alarm for time select

  if  [ ! -e  ${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    rm ${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280.nc
    cdo_old select,name=${variable} ${modelname}.cam.h0.* ${modelname}.cam.h1.0251-0280.nc
  fi

# step2 : interpolate the data from hybird level to pressure level
  if [ ! -e ${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
    echo "don't exit chazhi file, procecing..."
    ncl  -nQ inpath=\"${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280.nc\" \
         outpath=\"${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280_chazhi.nc\" \
      /home/ys17-19/lsh/Project/MSE-EP/F2000/210423-CESM-data-chazhi-finalver.ncl
    echo "finish CESM chazhi"
  fi 

   echo "finish this script"
#-----------------------------------------------------------


