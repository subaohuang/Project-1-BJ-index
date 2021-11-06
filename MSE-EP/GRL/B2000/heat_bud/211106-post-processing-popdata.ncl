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
  fi
  index=$(( $index + 1 ))
done
}
# Path of the original data
# Caution: DO NOT DELETE /" IN STRING!
## select : TEMP(0:300),UVEL,VVEL,WVEL,SHF,TAUX,TAUY 
PRE_DIR_ORG=/home/yangsong3/data-model/lsh/CESM/B/B2000_alt_north_year_CTRL/
# PRE_DIR_ORG=/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_db/

STEP=3
modelname=B2000_alt_north_year_CTRL
# modelname=lsh_B2000_WNP_heating_db

timestart=0251
timeend=0280

variable=UVEL,VVEL,WVEL,SHF,TAUX,TAUY
# variable=U,V

#step1 : merge the cesm data :UVEL,VVEL,WVEL
echo "step1 : merge the cesm data :UVEL,VVEL,WVEL"
# the prefix of data is usually CESM compet name ,alarm for time select
levintp=500.0,1500.0,2500.0,3500.0,4500.0,5500.0,6500
levintf=0.0,1000.0,2000.0,3000.0,4000.0,5000.0,6000.0
Var_Name_total=("UVEL" "VVEL" "TEMP")
cd $PRE_DIR_ORG
for Var_name in  ${Var_Name_total[*]}
    do 
      if  [ ! -e  ${PRE_DIR_ORG}${modelname}.$Var_name.pop.${timestart}-${timeend}.nc ] ; then
        echo "don't exit merge file, procecing..."
        echo $Var_name
        # rm ${PRE_DIR_ORG}${modelname}.$Var_name.pop.${timestart}-${timeend}.nc 
        cdo  -intlevel,$levintf -select,name=$Var_name,level=$levintp ${modelname}.pop.h.{0251..0280}*  ${modelname}.$Var_name.pop.${timestart}-${timeend}.nc
      fi
    done

levintp=0.0,1000.0,2000.0,3000.0,4000.0,5000.0
Var_Name_total=("WVEL")
cd $PRE_DIR_ORG
for Var_name in  ${Var_Name_total[*]}
    do 
      if  [ ! -e  ${PRE_DIR_ORG}${modelname}.$Var_name.pop.${timestart}-${timeend}.nc ] ; then
        echo "don't exit merge file, procecing..."
        echo $Var_name
        # rm ${PRE_DIR_ORG}${modelname}.$Var_name.pop.${timestart}-${timeend}.nc 
        cdo_old select,name=$Var_name,level=$levintp ${modelname}.pop.h.{0251..0280}* ${modelname}.$Var_name.pop.${timestart}-${timeend}.nc
      fi
    done
echo "finish STEP1"

#step2 : merge the cesm data :SHF,TAUX,TAUY
echo "step2 : merge the cesm data :SHF,TAUX,TAUY"

Var_Name_total=("SHF" "TAUX" "TAUY")
for Var_name in  ${Var_Name_total[*]}
    do
      if  [ ! -e  ${PRE_DIR_ORG}${modelname}.$Var_name.pop.${timestart}-${timeend}.nc ] ; then
        echo "don't exit merge file, procecing..."
        cd $PRE_DIR_ORG
        # rm ${PRE_DIR_ORG}${modelname}.$Var_name.pop.${timestart}-${timeend}.nc 
        cdo_old select,name=$Var_name ${modelname}.pop.h.{0251..0280}* ${modelname}.$Var_name.pop.${timestart}-${timeend}.nc
      fi
    done
echo "finish STEP2"


#step4 : interp the cesm POP data using PopLatLon(data,"gx1v6","fv1.9x2.5","bilin","da","090206")
Var_Name_total=("UVEL" "VVEL" "WVEL" "SHF" "TAUX" "TAUY" "TEMP")
# Var_Name_total=("UVEL")
for Var_name in  ${Var_Name_total[*]}
    do
   if  [ ! -e  ${PRE_DIR_ORG}${modelname}.$Var_name.fvgrid.${timestart}-${timeend}.nc ] ; then
     ## 斜杠用来添加"，不然"会被默认成输出变量用的符号 
   echo "don't exit chazhi file, procecing..."
   cd /home/ys17-19/lsh/Project/ENSO/CESM_WNPheating/B/
   pwd
   ncl -nQ infilepath=\"${PRE_DIR_ORG}${modelname}.$Var_name.pop.${timestart}-${timeend}.nc\" \
       outfilepath=\"${PRE_DIR_ORG}${modelname}.$Var_name.fvgrid.${timestart}-${timeend}.nc\" \
     ./210527-POP2fvgrid.ncl
   echo "finish POP to fv grid interpolation"
  fi
 done
echo "finish STEP3"

echo "finish this script"
#-----------------------------------------------------------



