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
# PRE_DIR_ORG=/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/
PRE_DIR_ORG=/home/yangsong3/CMIP6/linshh/CESM-data/B/lsh_B2000_WNP_heating_double/

STEP=3
# modelname=B2000_f19g16_CP_CTRL
modelname=lsh_B2000_WNP_heating_double

variable=UVEL,VVEL,WVEL,SHF,TAUX,TAUY
# variable=U,V
# wantyear={0251..0280}

#step1 : merge the cesm data :UVEL,VVEL,WVEL
echo "step1 : merge the cesm data :UVEL,VVEL,WVEL"
# the prefix of data is usually CESM compet name ,alarm for time select
levintp=500.0,1500.0,2500.0,3500.0,4500.0,5500.0
Var_Name_total=("UVEL" "VVEL")
cd $PRE_DIR_ORG
for Var_name in  ${Var_Name_total[*]}
    do 
      if  [ ! -e  ${PRE_DIR_ORG}${modelname}.$Var_name.pop.0251-0280.nc ] ; then
        echo "don't exit merge file, procecing..."
        echo $Var_name
        # rm ${PRE_DIR_ORG}${modelname}.$Var_name.pop.0251-0280.nc 
        cdo_old select,name=$Var_name,level=$levintp ${modelname}.pop.h.{0251..0280}* ${modelname}.$Var_name.pop.0251-0280.nc
      fi
    done

levintp=0.0,1000.0,2000.0,3000.0,4000.0,5000.0
Var_Name_total=("WVEL")
cd $PRE_DIR_ORG
for Var_name in  ${Var_Name_total[*]}
    do 
      if  [ ! -e  ${PRE_DIR_ORG}${modelname}.$Var_name.pop.0251-0280.nc ] ; then
        echo "don't exit merge file, procecing..."
        echo $Var_name
        # rm ${PRE_DIR_ORG}${modelname}.$Var_name.pop.0251-0280.nc 
        cdo_old select,name=$Var_name,level=$levintp ${modelname}.pop.h.{0251..0280}* ${modelname}.$Var_name.pop.0251-0280.nc
      fi
    done
echo "finish STEP1"



#step2 : merge the cesm data :SHF,TAUX,TAUY
echo "step2 : merge the cesm data :SHF,TAUX,TAUY"

Var_Name_total=("SHF" "TAUX" "TAUY")
for Var_name in  ${Var_Name_total[*]}
    do
      if  [ ! -e  ${PRE_DIR_ORG}${modelname}.$Var_name.pop.0251-0280.nc ] ; then
        echo "don't exit merge file, procecing..."
        cd $PRE_DIR_ORG
        # rm ${PRE_DIR_ORG}${modelname}.$Var_name.pop.0251-0280.nc 
        cdo_old select,name=$Var_name ${modelname}.pop.h.{0251..0280}* ${modelname}.$Var_name.pop.0251-0280.nc
      fi
    done
echo "finish STEP2"


echo "step3 : merge the cesm data: TEMP ,select vertical range of 0-300m"
#step3 : merge the cesm data: TEMP ,select vertical range of 0-300m
    ### pick 0-300 level 
     cdo -showlevel /home/yangsong3/data-model/lsh/CESM/B/B2000_alt_north_year/temp.pop.h.0251-01.nc > lelinfo.txt  ##print level info 
        cat lelinfo.txt | xargs -n1  > levinfo_new.txt ##convert to colum

        declare -a lev
        for line in `cat levinfo_new.txt`
        do 
           lev[$c]=${line}","    ### 加逗号
           ((c++))
        done 
      rm lelinfo.txt
      rm levinfo_new.txt

      getArrItemIdx "${lev[*]}" 40000  ##调用寻找深度为400的位置的函数
      ind=`echo $?`   ## 返回 functon return的值
      ind=`expr $ind + 1`  ##多取小于400的下一层

      echo ${lev[*]:0:$ind} > test_lev.txt ##挑取0-400的深度
      levselect=`cat test_lev.txt | sed s/[[:space:]]//g` ###去掉空格
      rm test_lev.txt
      echo $levselect
      ####垂直坐标差值范围的选取
      lev1=`seq -s , 500.0 1000.0 3500.0`
      lev2=`seq -s , 5000.0 1500.0 8000.0`
      lev3=`seq -s , 10000.0 2000.0 20000.0`
      lev4=`seq -s , 22500.0 2500.0 30000.0`
      lev5=`seq -s , 35000.0 5000.0 40000.0`

      levintp=${lev1},${lev2},${lev3},${lev4},${lev5}
      # echo $levintp
      if  [ ! -e  ${PRE_DIR_ORG}${modelname}.TEMP.pop.0251-0280.nc ] ; then
        echo "don't exit merge file, procecing..."
        cd $PRE_DIR_ORG
        cdo  -intlevel,$levintp -select,name=TEMP,level=$levselect ${modelname}.pop.h.{0251..0280}* ${modelname}.TEMP.pop.0251-0280.nc
      fi
echo "finish STEP3"


#step4 : interp the cesm POP data using PopLatLon(data,"gx1v6","fv1.9x2.5","bilin","da","090206")
Var_Name_total=("UVEL" "VVEL" "WVEL" "SHF" "TAUX" "TAUY" "TEMP")
# Var_Name_total=("UVEL")
for Var_name in  ${Var_Name_total[*]}
    do
   if  [ ! -e  ${PRE_DIR_ORG}${modelname}.$Var_name.fvgrid.0251-0280.nc ] ; then
     ## 斜杠用来添加"，不然"会被默认成输出变量用的符号 
   echo "don't exit chazhi file, procecing..."
   cd /home/ys17-19/lsh/Project/ENSO/CESM/B
   pwd
   ncl -nQ infilepath=\"${PRE_DIR_ORG}${modelname}.$Var_name.pop.0251-0280.nc\" \
       outfilepath=\"${PRE_DIR_ORG}${modelname}.$Var_name.fvgrid.0251-0280.nc\" \
     ./210527-POP2fvgrid.ncl
   echo "finish POP to fv grid interpolation"
  fi
 done
echo "finish STEP4"

echo "finish this script"
#-----------------------------------------------------------



