#!/bin/sh
#-----------------------------------------------
#   This is a shell script for calculating the
# mass streamfunction by using some NCL script, 
# and mainly focus on CESM ouput data,this shell script call
# ncl function script : 
#2019925-CESM-data-chazhi.ncl (function:cesm_hybird2pre(infilepath,outfilepath))
#2019926-mass_streamfunciton_cesm.ncl(functuon:get_msf(path,filename,startyear,endyear,outpath,outputname))
#step1:interpolate the CESM output data hybird level to pressure level ,step2:cauculate massstreamfuction
# You should set the basic parameters as below. 
# Good Luck!
#               creat on  2019-12-19
#
#               by Shuheng Lin  
#----------------------------------------------- 

# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
PRE_DIR_ORG=/home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/
PRE_DIR_OUT=/home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/

STEP=3
modelname=B2000_f19g16_CTRL

startyear=0250
endyear=0325
yearselect=${startyear}..${endyear}
echo $yearselect

 cd $PRE_DIR_ORG
 pwd
  rm ${PRE_DIR_OUT}${modelname}.sst.${startyear}-${endyear}.nc

  if  [ ! -e  ${PRE_DIR_OUT}${modelname}.sst.${startyear}-${endyear}.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG 
    cdo_old select,name=TEMP,level=500.0  ${modelname}.pop.h.{0250..0325}*  ${PRE_DIR_OUT}${modelname}.sst.${startyear}-${endyear}.nc
    # cdo_old select,name=TEMP,level=500.0  ${modelname}.pop.h.`{$select_year}`*  ${PRE_DIR_OUT}${modelname}.sst.${startyear}-${endyear}.nc
  fi


  if  [ ! -e  ${PRE_DIR_OUT}${modelname}.$Var_name.fvgrid.${startyear}-${endyear}.nc ] ; then
     ## 斜杠用来添加"，不然"会被默认成输出变量用的符号 
   echo "don't exit chazhi file, procecing..."
   cd $PRE_DIR_OUT
   pwd
   ncl -nQ infilepath=\"${PRE_DIR_OUT}${modelname}.sst.${startyear}-${endyear}.nc\" \
       outfilepath=\"${PRE_DIR_OUT}${modelname}.sst.${startyear}-${endyear}_fvgrid.nc\" \
     /home/ys17-19/lsh/Project/function_source/210527-POP2fvgrid.ncl
   echo "finish POP to fv grid interpolation"
  fi


   echo "finish this script"
#-----------------------------------------------------------



