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
PRE_DIR_ORG=/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/
# PRE_DIR_ORG=/home/yangsong3/data-model/lsh/CESM/B/B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer/
# PRE_DIR_ORG=/home/ys17-19/lsh/CESM-data/F/F_2000_addallocean_tropical/

STEP=3
modelname=B2000_f19g16_CP_CTRL

variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q
# variable=U,V

# wantyear={0251..0280}

#step1 : merge the cesm data into a whole data by using cdo
  # the prefix of data is usually CESM compet name ,alarm for time select

  # if  [ ! -e  ${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280.nc ] ; then
  #   echo "don't exit merge file, procecing..."
  #   cd $PRE_DIR_ORG
  #   rm ${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280.nc
  #   # cdo select,name=${variable} ${modelname}.cam.h0.* ${modelname}.cam.h1.0251-0280.nc
  #   cdo_old select,name=${variable} ${modelname}.cam.h0.{0251..0280}* ${modelname}.cam.h1.0251-0280.nc
  # fi

  # if  [ ! -e  ${PRE_DIR_ORG}${modelname}.cam.h0.0251-0280_UV.nc ] ; then
  #   echo "don't exit merge file, procecing..."
  #   cd $PRE_DIR_ORG
  #   rm ${PRE_DIR_ORG}${modelname}.cam.h0.0251-0280_UV.nc
  #   cdo select,name=${variable} ${modelname}.cam.h0.* ${modelname}.cam.h0.0251-0280_UV.nc
  # fi

# #step2 : interpolate the data from hybird level to pressure level
#   if [ ! -e ${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
#     echo "don't exit chazhi file, procecing..."
#     cd /home/ys17-19/lsh/Project/Walker-Circulation/using-CESM-simulate-WC/F_2000/
#     pwd
#     ncl  -nQ inpath=\"${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280.nc\" \
#          outpath=\"${PRE_DIR_ORG}${modelname}.cam.h1.0251-0280_chazhi.nc\" \
#        /home/ys17-19/lsh/Project/SCS-rain/annual/191209-CESM-data-chazhi.ncl
#     echo "finish CESM chazhi"
#   fi 

 cd $PRE_DIR_ORG
 pwd
  if  [ ! -e  ${PRE_DIR_ORG}${modelname}.TEMP.h1.0251-0280.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    rm ${PRE_DIR_ORG}${modelname}.TEMP.h1.0251-0280.nc
    cdo_old select,name=TEMP,level=500.0  ${modelname}.pop.h.{0210..0340}*  ${modelname}.TEMP.h1.0210-0340.nc
  fi


  # if  [ ! -e  ${PRE_DIR_ORG}${modelname}.TEMP_taux_tauy.h0.0251-0280.nc ] ; then
  #   echo "don't exit merge file, procecing..."
  #   cd $PRE_DIR_ORG
  #   rm ${PRE_DIR_ORG}${modelname}.TEMP_taux_tauy.h0.0251-0280.nc
  #   cdo select,name=TEMP,TAUX,TAUY,level=500.0  ${modelname}.pop.h.*  ${modelname}.TEMP_taux_tauy.h0.0251-0280.nc
  # fi


# # step3 calculate mass stream function

#   if  [ ! -e /home/ys17-19/lsh/data/wc-result/msf_${modelname}_0101-4012.nc ] ; then
#      ## 斜杠用来添加"，不然"会被默认成输出变量用的符号 
#    echo "don't exit msf file, procecing..."
#    cd /home/ys17-19/lsh/Project/Walker-Circulation/using-CESM-simulate-WC/F_2000/
#    pwd
#    ncl -nQ inpath=\"${PRE_DIR_ORG}\"               \
#        filename=\"${modelname}.cam.h0.0101-4012_chazhi.nc\" \
#        outputpath='"~/lsh/data/wc-result/"' \
#        outputname=\"msf_${modelname}_0101-4012.nc\" \
#      ./191209-mass_streamfunciton_cesm-F2000.ncl
#    echo "finish CESM mass_stream function"
#   fi


  # if  [ ! -e  ${PRE_DIR_ORG}${modelname}.mld.h0.0251-0280.nc ] ; then
  #   echo "don't exit merge file, procecing..."
  #   cd $PRE_DIR_ORG
  #   rm ${PRE_DIR_ORG}${modelname}.mld.h0.0251-0280.nc
  #   cdo select,name=HMXL  ${modelname}.pop.h.*  ${modelname}.mld.h0.0251-0280.nc
  # fi
  
   echo "finish this script"
#-----------------------------------------------------------



