#!/bin/bash

#do
  OUT_DIR=/home/ys17-19/lsh/Project/ITCZ_bias/data/CMIP5/v_chazhi/
#  if [ ! -d $OUT_DIR ]; then
#    mkdir $OUT_DIR
#  fi

#  cd $OUT_DIR
  IN_DIR=/home/yangsong3/data-model/CMIP5-hes/va/historical/
  echo $IN_DIR

  for file in `ls $IN_DIR`
  do
    cdo remapbil,global_2.5 $IN_DIR/$file $OUT_DIR/$file
  done
#done
  echo "finish"