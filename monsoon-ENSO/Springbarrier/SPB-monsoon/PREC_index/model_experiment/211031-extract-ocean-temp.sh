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

cd /home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/

Var_Name_total=("TEMP")
modelname=B2000_f19g16_CP_CTRL

echo "********start  $Model_Name  *************************"
OUT_DIR=/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/
IN_DIR=/home/yangsong3/data-model/CESM_CTRL/B2000_F19G16_CAM4_CTRL/pre/


# ncl_filedump -v z_t  B2000_f19g16_CP_CTRL.pop.h.0268-11.nc > lelinfo.txt  ##print level info 
# cat lelinfo.txt | xargs -n1  > levinfo_new.txt ##convert to colum
 

# declare -a lev
# for line in `cat levinfo_new.txt`
# do 
#    lev[$c]=${line}","    ### 加逗号
#    ((c++))
# done 
# rm lelinfo.txt
# rm  levinfo_new.txt
# getArrItemIdx "${lev[*]}" 40000  ##调用寻找深度为400的位置的函数
# ind=`echo $?`   ## 返回 functon return的值
# ind=`expr $ind + 1`  ##多取小于400的下一层
# echo ${lev[*]:0:$ind} > test_lev.txt ##挑取0-400的深度
# levselect=`cat test_lev.txt | sed s/[[:space:]]//g` ###去掉空格
# rm test_lev.txt
# echo $levselect

#  z_t=500,1500,2500,3500,4500,5500,6500,7500,8500,9500,10500,11500,12500,13500,14500,15500,16509.84,17547.9,18629.13,19766.03,\
# 20971.14,22257.83,23640.88,25137.02,26765.42,28548.37,30511.92,32686.8,35109.35,37822.76,40878.46

 z_t=165

cdo  -select,name=${Var_Name_total},level=$z_t ${modelname}.pop.h.0{250..349}* ${modelname}.pop.temp_upperlevel.0250-0349nc
        # rename .nc _2.5.nc $OUT_DIR/$file
echo "finish script"
