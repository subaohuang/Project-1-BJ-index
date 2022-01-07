#!/bin/sh
#-----------------------------------------------
#merge and calculate ensemble mean 
               ###by Shuheng Lin  
#----------------------------------------------- 

# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
modelname=year286_enso_wind_forcing_heatsen_5m
PRE_DIR_ORG=/home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/${modelname}/
STEP=1
# variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q,UBOT,VBOT,TREFHT
variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q,PTENDT
# variable=hyam
ensemble_num=10
ensemble_start=1
ensemble_end=26
# module load cdo/1.9.9rc7-icc15
if [ ${STEP} -eq 1 ]; then 
# for ((i=1;i<$ensemble_num+1;i++)) ;
for i in `seq $ensemble_start $ensemble_end`
do 
   echo ${i}
#step1 : merge the cesm data into a whole data by using cdo
  ## the prefix of data is usually CESM compet name ,alarm for time select
   # rm ${PRE_DIR_ORG}Ensemble_${i}/${modelname}.cam.h1.0286-0287.nc
  if  [ ! -e  ${PRE_DIR_ORG}Ensemble_${i}/${modelname}.cam.h1.0286-0287.nc ] ; then
    echo "don't exit merge file, procecing..."
    cd ${PRE_DIR_ORG}Ensemble_${i}
    cdo_old  select,name=${variable} ${modelname}.cam.h0.* ${modelname}.cam.h1.0286-0287.nc
  fi
done 
fi

# if [ ${STEP} -eq 2 ]; then 
# for ((i=1;i<$ensemble_num+1;i++)) ;
# do 
# # step2 : interpolate the data from hybird level to pressure level
#   if [ ! -e ${PRE_DIR_ORG}Ensemble_${i}/${modelname}.cam.h1.0286-0287_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
#     echo "don't exit chazhi file, procecing..."
#     ncl  -nQ inpath=\"${PRE_DIR_ORG}Ensemble_${i}/${modelname}.cam.h1.0286-0287.nc\" \
#          outpath=\"${PRE_DIR_ORG}Ensemble_${i}/${modelname}.cam.h1.0286-0287_chazhi.nc\" \
#     /WORK/sysu_hjkx_ys/linshh/B/ensemble_run/ncl_for_monsoon_enso_experiment/210423-CESM-data-chazhi-finalver.ncl
#     echo "finish CESM chazhi"
#   fi 
# done 
# fi 


if [ ${STEP} -eq 1 ]; then 
cd ${PRE_DIR_ORG}
rm ${modelname}.cam.h1.0286-0287_ensemblemean.nc 
# step3 : ensemble mean 
# cdo ensavg ${PRE_DIR_ORG}Ensemble_{1..20}/${modelname}.cam.h1.0286-0287.nc  ${modelname}.cam.h1.0286-0287_ensemblemean.nc 
cdo_old ensavg ${PRE_DIR_ORG}Ensemble_{1..26}/${modelname}.cam.h1.0286-0287.nc  ${modelname}.cam.h1.0286-0287_ensemblemean.nc 
#-----------------------------------------------------------
fi


if [ ${STEP} -eq 4 ]; then 
rm ${PRE_DIR_ORG}${modelname}.cam.h1.0286-0287_ensemblemean_chazhi.nc
for ((i=1;i<$ensemble_num+1;i++)) ;
do 
  if [ ! -e ${PRE_DIR_ORG}${modelname}.cam.h1.0286-0287_ensemblemean_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
    echo "don't exit chazhi file, procecing..."
    ncl  -nQ inpath=\"${PRE_DIR_ORG}${modelname}.cam.h1.0286-0287_ensemblemean.nc\" \
         outpath=\"${PRE_DIR_ORG}${modelname}.cam.h1.0286-0287_ensemblemean_chazhi.nc\" \
    /WORK/sysu_hjkx_ys/linshh/B/ensemble_run/ncl_for_monsoon_enso_experiment/210423-CESM-data-chazhi-finalver.ncl
    echo "finish CESM chazhi"
  fi 
done 
fi 

echo "finish this script"

