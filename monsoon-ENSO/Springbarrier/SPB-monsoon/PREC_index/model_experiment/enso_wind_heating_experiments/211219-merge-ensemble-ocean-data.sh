#!/bin/sh
#-----------------------------------------------
#merge and calculate ensemble mean 
               ###by Shuheng Lin  
#----------------------------------------------- 

# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
modelname=year286_enso_wind_forcing_dbcoolsen_5m
# modelname=year286_enso_wind_forcing_heatsen_5m
PRE_DIR_ORG=/WORK/sysu_hjkx_ys/linshh/B/ensemble_run/${modelname}/ensemble/
STEP=1
# variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q,UBOT,VBOT,TREFHT
variable=TEMP
# variable=hyam
ensemble_num=20

ensemble_start=21
ensemble_end=30

# # module load cdo/1.9.9rc7-icc15
# if [ ${STEP} -eq 1 ]; then 
# # for ((i=1;i<$ensemble_num+1;i++)) ;
# for i in `seq $ensemble_start $ensemble_end`
# do 
#   echo ${i}
# #step1 : merge the cesm data into a whole data by using cdo
#   ## the prefix of data is usually CESM compet name ,alarm for time select
#   rm ${PRE_DIR_ORG}Ensemble_${i}/${modelname}.sst.0286-0287.nc
#   if  [ ! -e  ${PRE_DIR_ORG}Ensemble_${i}/${modelname}.sst.0286-0287.nc ] ; then
#     echo "don't exit merge file, procecing..."
#     cd ${PRE_DIR_ORG}Ensemble_${i}
#     cdo  -select,name=${variable},level=500.0 ${modelname}.pop.h.0* ${modelname}.sst.0286-0287.nc
#   fi
# done 
# fi


# if [ ${STEP} -eq 2 ]; then 
cd ${PRE_DIR_ORG}
rm ${modelname}.sst.0286-0287_ensemblemean.nc 
# step3 : ensemble mean 
cdo ensavg ${PRE_DIR_ORG}Ensemble_{1..10}/${modelname}.sst.0286-0287.nc  ${modelname}.sst.0286-0287_ensemblemean_1-10.nc 
#-----------------------------------------------------------
# fi

# rm ${PRE_DIR_ORG}${modelname}.sst.0286-0287_ensemblemean_chazhi.nc
# if [ ${STEP} -eq 4 ]; then 
# for ((i=1;i<$ensemble_num+1;i++)) ;
# do 
#   if [ ! -e ${PRE_DIR_ORG}${modelname}.sst.0286-0287_ensemblemean_chazhi.nc ] ; then    ####判断差值的文件是否已经存在
#     echo "don't exit chazhi file, procecing..."
#     ncl  -nQ inpath=\"${PRE_DIR_ORG}${modelname}.sst.0286-0287_ensemblemean.nc\" \
#          outpath=\"${PRE_DIR_ORG}${modelname}.sst.0286-0287_ensemblemean_chazhi.nc\" \
#     /WORK/sysu_hjkx_ys/linshh/B/ensemble_run/ncl_for_monsoon_enso_experiment/210423-CESM-data-chazhi-finalver.ncl
#     echo "finish CESM chazhi"
#   fi 
# done 
# fi 

echo "finish this script"

