;;********************************************************
;;将模式的海温数据的印度洋地区的海温去趋势，并且其他海温不变
;;得到一个新的海温数据，注意在文件中的其他数据要保留不能变化
;;linshuheng 2019.9.22 
;**********************************************************
;test
begin
  ;;;----------------------------------------------------------读入数据
    fils1 = "~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_c130411.nc"
    f1 = addfile(fils1, "r")
    time = f1->time
    lon = f1->lon
    lat = f1->lat

    SST_cpl = f1->SST_cpl(:,:,:)
    SST_cpl_prediddle = f1->SST_cpl_prediddle(:,:,:)

    ice_cov  = f1->ice_cov
    ice_cov_prediddle = f1->ice_cov_prediddle

  ;;;------去掉太平洋的海温增温趋势，海盆的大小定义为50N-50S ,150E-280E---------
    lat1 = -90
    lat2 = 90
    lon1 = 0
    lon2 = 360

     SST_cpl_pacific = SST_cpl(:,{lat1:lat2},{lon1:lon2})
     SST_cpl_prediddle_pacific = SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2})

     printVarSummary(SST_cpl_pacific)

     ;;-------cesm中的海温在陆地上也有数据，所以这里不能将陆地的温度也去掉--------
       fils2 = "~/lsh/data/landsea.nc"
       f2 = addfile(fils2, "r")
       lsdata = f2->LSMASK
       lsm = landsea_mask(lsdata,SST_cpl_pacific&lat,SST_cpl_pacific&lon) 
       lsm1 = landsea_mask(lsdata,SST_cpl_prediddle_pacific&lat,SST_cpl_prediddle_pacific&lon)
       
       SST_cpl_pacific_mask = SST_cpl_pacific
       SST_cpl_pacific_mask = mask(SST_cpl_pacific, lsm.eq.1, False)

       SST_cpl_prediddle_pacific_mask = SST_cpl_prediddle_pacific
       SST_cpl_prediddle_pacific_mask = mask(SST_cpl_prediddle_pacific, lsm.eq.1, False)

       printVarSummary(SST_cpl_pacific_mask)

       SST_pacific_de = dtrend_msg_n(ispan(1,dimsizes(SST_cpl_pacific_mask&time),1),SST_cpl_pacific_mask,False,True,0)
       SST_prediddle_pacific_de = dtrend_msg_n(ispan(1,dimsizes(SST_cpl_prediddle_pacific_mask&time),1),SST_cpl_prediddle_pacific_mask,False,True,0)


       printVarSummary(SST_pacific_de)
       print(SST_pacific_de@_FillValue)
       SST_pacific_de = where(ismissing(SST_pacific_de) ,SST_cpl_pacific, SST_pacific_de)                         ;;;把陆地mask掉的海温再还原
       SST_prediddle_pacific_de = where(ismissing(SST_prediddle_pacific_de) ,SST_cpl_prediddle_pacific, SST_prediddle_pacific_de)                         ;;;把陆地mask掉的海温再还原
  ;;;----------------------------------------------------------------------
       

  ;;;-----将原来的数据中的印度洋数据用去趋势的海温代替--------------  
    
    SST_cpl(:,{lat1:lat2},{lon1:lon2}) = (/SST_pacific_de/)
    printVarSummary(SST_cpl)
    
    SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2}) = (/SST_prediddle_pacific_de/)
    printVarSummary(SST_cpl_prediddle)
      


    ;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f ~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_alltrend.nc")   ; remove any pre-existing file
    ncdf = addfile("~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_alltrend.nc" ,"c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications

    ;===================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
       ncdf->SST_cpl   =  SST_cpl                           ; 4D               
       ncdf->SST_cpl_prediddle  = SST_cpl_prediddle                        ; 3D

       ncdf->ice_cov   =  ice_cov
       ncdf->ice_cov_prediddle = ice_cov_prediddle

       ncdf->lat   =  lat
       ncdf->lon = lon 
       ncdf->time = time
 

end 

     






