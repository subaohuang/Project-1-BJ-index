#!/bin/sh
#!/bin/bash
## written by Shuheng Lin 2021.03.29
## interpolate CMIP6 DATA to 2.5 grid
###select model level and interpolate to 5m intervals
function getArrItemIdx(){
local arr=$1
local item=$2
local index=0
local int
for i in ${arr[*]}
do
  int=`echo $i |awk '{print int($0)}'`  ##取整，因为if只能判断整数

  if [  $int -gt $item ]
    then
    # echo $index
    return $index 
    # retrun $(( $index + 1 ))
  fi
  index=$(( $index + 1 ))
done
}

# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
modelname=year286_enso_wind_forcing_ctrl_5m
# modelname=year286_enso_wind_forcing_heatsen_5m
PRE_DIR_ORG=/home/yangsong3/CMIP6/linshh/CESM-data/B/enso_monsoon_ensemble/${modelname}/
STEP=1
# variable=U,V,OMEGA,PRECL,PRECC,PSL,PS,Z3,Q,UBOT,VBOT,TREFHT
variable=TEMP
# variable=hyam
ensemble_start=1
ensemble_end=30
# year1=0257
# year2=0258

year1=0286
year2=0287

cd /home/yangsong3/CMIP6/linshh/CESM-data/B/B2000_f19g16_CTRL/
cdo -s -showlevel pop.temp_upperlevel.0250-0300.nc > lelinfo.txt  ##print level info 
cat lelinfo.txt | xargs -n1  > levinfo_new.txt ##convert to colum

declare -a lev
for line in `cat levinfo_new.txt`
do 
   lev[$c]=${line}","    ### 加逗号
   ((c++))
done 

rm lelinfo.txt
rm  levinfo_new.txt

echo ${lev[*]} > test_lev.txt ##挑取0-400的深度

# getArrItemIdx "${lev[*]}" 40000  ##调用寻找深度为400的位置的函数
# ind=`echo $?`   ## 返回 functon return的值
# ind=`expr $ind + 1`  ##多取小于400的下一层
# echo ${lev[*]:0:$ind} > test_lev.txt ##挑取0-400的深度
levselect=`cat test_lev.txt | sed s/[[:space:]]//g` ###去掉空格
# exit
#  z_t=500,1500,2500,3500,4500,5500,6500,7500,8500,9500,10500,11500,12500,13500,14500,15500,16509.84,17547.9,18629.13,19766.03,\
# 20971.14,22257.83,23640.88,25137.02,26765.42,28548.37,30511.92,32686.8,35109.35,37822.76,40878.46
z_t=$levselect
echo $z_t


cd ${PRE_DIR_ORG}

# module load cdo/1.9.9rc7-icc15
if [ ${STEP} -eq 1 ]; then 
# for ((i=1;i<$ensemble_num+1;i++)) ;
for i in `seq $ensemble_start $ensemble_end`
do 
  echo ${i}
#step1 : merge the cesm data into a whole data by using cdo
  ## the prefix of data is usually CESM compet name ,alarm for time select
  # if  [ -e ${PRE_DIR_ORG}Ensemble_fixinitial_${i}/${modelname}.fixinitial.sst.${year1}-${year2}.nc ] ; then
  # rm ${PRE_DIR_ORG}Ensemble_${i}/${modelname}.upperlevel.temp.${year1}-${year2}.nc
  # fi
  if  [ ! -e  ${PRE_DIR_ORG}Ensemble_${i}/${modelname}.upperlevel.temp.${year1}-${year2}.nc ] ; then
  	echo "ensemble '${i}'"
    echo "don't exit merge file, procecing..."
    cd ${PRE_DIR_ORG}Ensemble_${i}
    cdo  -select,name=${variable},level=$z_t ${modelname}.pop.h.0* ${modelname}.upperlevel.temp.${year1}-${year2}.nc
  fi
done 
fi

# if [ ${STEP} -eq 2 ]; then 
cd ${PRE_DIR_ORG}
# rm ${modelname}.fixinitial.sst.${year1}-${year2}_ensemblemean.nc 
# step3 : ensemble mean 
cdo ensavg ${PRE_DIR_ORG}Ensemble_{1..30}/${modelname}.upperlevel.temp.${year1}-${year2}.nc  ${modelname}.upperlevel.temp.${year1}-${year2}_ensemblemean_1-30.nc

echo "finish script"
