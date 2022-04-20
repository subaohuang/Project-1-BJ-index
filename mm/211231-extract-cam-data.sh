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
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F/F_2000_IhPO/
PRE_DIR_ORG=/home/yangsong3/CMIP6/linshh/CESM-data/F_mm/F2000_TPP_indNTA_JAwarming/
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F/F_2000_GHG/
STEP=2
modelname=F2000_TPP_indNTA_JAwarming
# modelname=F_2000_ghg
# variable=PSL,U,V,PRECL,PRECC,PS
variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q,TREFHT,T
# variable=hyam


#step1 : merge the cesm data into a whole data by using cdo
  ## the prefix of data is usually CESM compet name ,alarm for time select

if [ ${STEP} -eq 1 ]; then 
  rm  ${PRE_DIR_ORG}${modelname}.cam.h1.0001-0030.nc
  if  [ ! -e  ${PRE_DIR_ORG}${modelname}.cam.h1.0001-0030.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    rm ${PRE_DIR_ORG}${modelname}.cam.h1.0001-0030.nc
    cdo_old -L -select,name=${variable} ${modelname}.cam.h0.*  ${modelname}.cam.h1.0001-0030.nc
  fi
fi

if [ ${STEP} -eq 2 ]; then 
# step2 : interpolate the data from hybird level to pressure level
  if [ ! -e ${PRE_DIR_ORG}${modelname}.cam.h1.0001-0030_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
    echo "don't exit chazhi file, procecing..."
    pwd
    ncl  -nQ inpath=\"${PRE_DIR_ORG}${modelname}.cam.h1.0001-0030.nc\" \
         outpath=\"${PRE_DIR_ORG}${modelname}.cam.h1.0001-0030_chazhi.nc\" \
      /home/ys17-19/lsh/Project/function_source/CESM-data-chazhi-finalver.ncl
    echo "finish CESM chazhi"
  fi 
fi

   echo "finish this script"
#-----------------------------------------------------------


