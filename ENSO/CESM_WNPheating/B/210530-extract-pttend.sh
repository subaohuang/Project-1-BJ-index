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
# PRE_DIR_ORG=/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/
PRE_DIR_ORG=/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_Idl/
STEP=3
modelname=lsh_B2000_WNP_heating_Idl
# modelname=F_2000_ghg
# variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q,UBOT,VBOT,TREFHT
variable=DTCOND,PTENDT,PS
# variable=U,V,PS

# variable=hyam

#step1 : merge the cesm data into a whole data by using cdo
  ## the prefix of data is usually CESM compet name ,alarm for time select

  if  [ ! -e  ${PRE_DIR_ORG}${modelname}.PTENDT.cam.h1.0251-0280.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd $PRE_DIR_ORG
    cdo_old select,name=${variable} ${modelname}.cam.h0.* ${modelname}.PTENDT.cam.h1.0251-0280.nc
  fi

# step2 : interpolate the data from hybird level to pressure level
  if [ ! -e ${PRE_DIR_ORG}${modelname}.PTENDT.cam.h1.0251-0280_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
    echo "don't exit chazhi file, procecing..."
    cd /home/ys17-19/lsh/Project/ENSO/CESM/B/
    pwd
    ncl  -nQ inpath=\"${PRE_DIR_ORG}${modelname}.PTENDT.cam.h1.0251-0280.nc\" \
         outpath=\"${PRE_DIR_ORG}${modelname}.PTENDT.cam.h1.0251-0280_chazhi.nc\" \
     ./210423-CESM-data-chazhi-finalver.ncl
    echo "finish CESM chazhi"
  fi 
# # step3 calculate mass stream function

#   if  [ ! -e /home/ys17-19/lsh/data/wc-result/msf_${modelname}_0251-0280.nc ] ; then
#      ## 斜杠用来添加"，不然"会被默认成输出变量用的符号 
#    echo "don't exit msf file, procecing..."
#    cd /home/ys17-19/lsh/Project/Walker-Circulation/using-CESM-simulate-WC/F_2000/
#    pwd
#    ncl -nQ inpath=\"${PRE_DIR_ORG}\"               \
#        filename=\"${modelname}.cam.h0.0251-0280_chazhi.nc\" \
#        outputpath='"~/lsh/data/wc-result/"' \
#        outputname=\"msf_${modelname}_0251-0280.nc\" \
#      ./191209-mass_streamfunciton_cesm-F2000.ncl
#    echo "finish CESM mass_stream function"
#   fi


# # # step4 calculate Local WC and HC
#    if [ ! -e ${PRE_DIR_ORG}${modelname}.cam.h0.0251-0280_local_wk_hc_500.nc ] ; then    ####判断差值的文件是否已经存在
#    echo "don't exit WC and HC file, procecing..."
#    cd /home/ys17-19/lsh/Project/Walker-Circulation/using-CESM-simulate-WC/F_2000/
#    pwd
#    ncl -nQ inpath=\"${PRE_DIR_ORG}\"               \
#        filename=\"${modelname}.cam.h0.0251-0280_chazhi.nc\" \
#        outputpath=\"${PRE_DIR_ORG}\"               \
#        outputname=\"${modelname}.cam.h0.0251-0280_local_wk_hc_500.nc\" \
#      ./200902-local_wk_hc_cesm-F2000.ncl
#    echo "finish CESM mass_stream function"
#   fi
#    echo "finish this script"
# #-----------------------------------------------------------


# # # step5 calculate vp and sf
#    if [ ! -e ${PRE_DIR_ORG}${modelname}.PTENDT.cam.h1.0251-0280_vp_sf.nc ] ; then    ####判断差值的文件是否已经存在
#    echo "don't exit vp and sf file, procecing..."
#    cd /home/ys17-19/lsh/Project/SCS-rain/F2000_model_experiment/
#    pwd
#    ncl -nQ infilepath=\"${PRE_DIR_ORG}${modelname}.PTENDT.cam.h1.0251-0280_chazhi.nc\" \
#        outfilepath=\"${PRE_DIR_ORG}${modelname}.PTENDT.cam.h1.0251-0280_vp_sf.nc\"   \
#      ./200909-F2000-cal-write-era-velocity-potensial-streamfuc.ncl
#    echo "finish CESM VP and SF calculate"
#   fi

# # # # step6 calculate vr and ur
#    if [ ! -e ${PRE_DIR_ORG}${modelname}.cam.h0.0251-0280_ur_vr.nc ] ; then    ####判断差值的文件是否已经存在
#    echo "don't exit rotation wind file, procecing..."
#    cd /home/ys17-19/lsh/Project/SCS-rain/F2000_model_experiment/
#    pwd
#    ncl -nQ inpath=\"${PRE_DIR_ORG}\"                           \
#         filename=\"${modelname}.PTENDT.cam.h1.0251-0280_chazhi.nc\" \
#         outputpath=\"${PRE_DIR_ORG}\"                          \
#        outputname=\"${modelname}.PTENDT.cam.h1.0251-0280_ur_vr.nc\"   \
#      ./201103-cal-F2000-ur-vr.ncl
#    echo "finish CESM rotation wind calculate"
#   fi
   echo "finish this script"
#-----------------------------------------------------------


