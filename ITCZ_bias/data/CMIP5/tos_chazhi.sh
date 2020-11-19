#!/bin/bash
## written by shuheng Lin 2020.11.19
## interpolate CMIP5 tos DATA to 2.5 grid
cd /home/yangsong3/data-model/CMIP5-wzq/ocean/mon/tos
Model_Name_total=`ls `
echo $Model_Name_total
# Model_Name_total=("CMCC-CESM" )

for Model_Name in ${Model_Name_total[*]}
do 
  echo "********start  $Model_Name  *************************"
  OUT_DIR=/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/tos_chazhi/$Model_Name/
  IN_DIR=/home/yangsong3/data-model/CMIP5-wzq/ocean/mon/tos/$Model_Name/r1i1p1/

  if [ ! -d $OUT_DIR ]; then
    mkdir $OUT_DIR
    cd $IN_DIR
    file_test=`ls $IN_DIR/*`
    if [  ${#file_test[@]} -gt 0 ]; then
    for file in `ls *`
    do
    cdo  remapbil,global_2.5 $IN_DIR/$file $OUT_DIR/$file
    rename .nc _2.5.nc $OUT_DIR/$file
    done
    fi
   echo "********finish  $Model_Name  *************************"
  fi 
done

  echo "finish script"

