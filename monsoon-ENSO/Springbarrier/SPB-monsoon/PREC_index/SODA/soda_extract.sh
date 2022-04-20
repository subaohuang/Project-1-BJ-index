#!/bin/sh
variable=temp
cdo -select,name=${variable},level=5.0  "http://iridl.ldeo.columbia.edu/SOURCES/.CARTON-GIESE/.SODA/.v2p2p4/temp/dods" soda.sst.nc