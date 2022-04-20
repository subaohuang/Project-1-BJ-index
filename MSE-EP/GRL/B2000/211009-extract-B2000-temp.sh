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
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_year/
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/B/B2000_alt_north_winter_fixmonth/
PRE_DIR_ORG=/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_alt_north_db_year_80_20/
# PRE_DIR_ORG=/home/yangsong3/data-model/lsh/CESM/B/B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer/
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F/F_2000_addallocean_tropical/

STEP=3
modelname=lsh_B2000_alt_north_db_year_80_20

variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q
# variable=U,V

# wantyear={0251..0280}

 cd $PRE_DIR_ORG
 pwd
  if  [ ! -e  ${PRE_DIR_ORG}${modelname}.TEMP.h1.0251-0280.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    rm ${PRE_DIR_ORG}${modelname}.TEMP.h1.0251-0280.nc
    cdo_old select,name=TEMP,level=500.0  ${modelname}.pop.h.{0251..0280}*  ${modelname}.TEMP.h1.0251-0280.nc
  fi

   echo "finish this script"
#-----------------------------------------------------------



