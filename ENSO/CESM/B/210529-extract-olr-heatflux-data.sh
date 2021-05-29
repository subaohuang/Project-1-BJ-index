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
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_spring_fixmonth/
PRE_DIR_ORG=/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F/F_2000_addtropical_winter/

STEP=3
modelname=B2000_f19g16_CP_CTRL
# modelname=F_2000_tro_winter
variable=FLUT,FSNS,FLNS,LHFLX,SHFLX
# variable=TAUX,TAUY
# variable=FLUT


# step1 : merge the cesm data into a whole data by using cdo
  ### the prefix of data is usually CESM compet name ,alarm for time select

  if  [ ! -e  ${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280_OLR_heatflux.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    rm ${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280_OLR_heatflux.nc
    cdo select,name=${variable} ${modelname}.cam.h0.{0251..0280}* ${modelname}.cam.h1.0251-0280_OLR_heatflux.nc
  fi


   echo "finish this script"
#-----------------------------------------------------------



