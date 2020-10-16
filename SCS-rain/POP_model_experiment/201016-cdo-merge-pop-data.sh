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
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F/F_2000_IPO/
PRE_DIR_ORG=/home/yangsong3/data-observation/linshh/CESM-data/C/lsh_C_CTRL
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F/F_2000_addallocean_tropical/

STEP=3
modelname=lsh_C_CTRL
# variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q

 cd $PRE_DIR_ORG
 pwd
  if  [ ! -e  ${PRE_DIR_ORG}${modelname}.TEMP_TAUX_TAUY.h0.0001-0030.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    rm ${PRE_DIR_ORG}${modelname}.TEMP_TAUX_TAUY.h0.0001-0030.nc
    cdo select,name=TAUX,TAUY,TEMP,level=500.0  ${modelname}.pop.h.*  ${modelname}.TEMP_TAUX_TAUY.h0.0001-0030.nc
  fi
  
   echo "finish this script"
#-----------------------------------------------------------



