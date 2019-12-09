;;********************************************************
;;将模式的海温数据的印度洋地区的海温增加，增加的幅度是通过计算sst_HadOIBl_bc_1.9x2.5_1978_2012_c130411.nc这个文件的趋势
;;得到一个新的海温数据，注意在文件中的其他数据要保留不能变化
;;用来跑F2000实验
;;lnshuheng 2019.12.07 
;**********************************************************
;test
begin
;;;---------------读入变化的海温数据----------------------------
    fils1 = "/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_1978_2012_c130411.nc"
    f1 = addfile(fils1, "r")
    time = f1->time
    date = cd_calendar(time, 1) 
    lon = f1->lon
    lat = f1->lat
    timeind = ind(date.le.201212.and.date.ge.197901)

    SST_cpl = f1->SST_cpl(timeind,:,:)
    SST_cpl_prediddle = f1->SST_cpl_prediddle(timeind,:,:)

;;;------只保留印度的海温增温趋势，其他海盆的温度趋势去掉---------
     
 
      SST_cpl_indian = SST_cpl(:,{-20:20},{30:100})
      SST_cpl_prediddle_indian = SST_cpl_prediddle(:,{-20:20},{30:100})

      printVarSummary(SST_cpl_indian)

      printVarSummary(SST_cpl_indian)

       ;;-------cesm中的海温在陆地上也有数据，所以这里不能将陆地的温度也去掉--------
       fils2 = "~/lsh/data/landsea.nc"
       f2 = addfile(fils2, "r")
       lsdata = f2->LSMASK
       lsm = landsea_mask(lsdata,SST_cpl_indian&lat,SST_cpl_indian&lon) 
       lsm1 = landsea_mask(lsdata,SST_cpl_prediddle_indian&lat,SST_cpl_prediddle_indian&lon)
       
       SST_cpl_indian_mask = SST_cpl_indian
       SST_cpl_indian_mask = mask(SST_cpl_indian, lsm.eq.1, False)

       SST_cpl_prediddle_indian_mask = SST_cpl_prediddle_indian
       SST_cpl_prediddle_indian_mask = mask(SST_cpl_prediddle_indian, lsm.eq.1, False)

       printVarSummary(SST_cpl_indian_mask)

       SST_indian_rc = regCoef_n(ispan(1,dimsizes(SST_cpl_indian_mask&time),1),SST_cpl_indian_mask,0,0)
       SST_prediddle_indian_rc = regCoef_n(ispan(1,dimsizes(SST_cpl_prediddle_indian_mask&time),1),SST_cpl_prediddle_indian_mask,0,0)

       printVarSummary(SST_indian_rc)
       print(SST_indian_rc@_FillValue)
     
       ssta_indian = SST_indian_rc * dimsizes(SST_cpl_indian_mask&time)
       ssta_indian_prediddle = SST_indian_rc * dimsizes(SST_cpl_indian_mask&time)
     
       ssta_indian = where(ismissing(ssta_indian) ,0, ssta_indian)                         ;;;把陆地mask掉的海温再还原
       ssta_indian_prediddle = where(ismissing(ssta_indian_prediddle) ,0, ssta_indian) ;;;把陆地mask掉的海温再还原

       delete(SST_cpl_prediddle)
       delete(SST_cpl)
       delete(time)
       delete(lon)
       delete(lat)
       delete(date)
   
       ;;;;
  ;;;;
  
;;;---------------读入气候态海温数据----------------------------
    fils3 = "/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031.nc"
    f3 = addfile(fils3, "r")
    lon = f3->lon
    lat = f3->lat
    time= f3->time
    date = f3->date
    datesec = f3->datesec
    
    SST_cpl = f3->SST_cpl(:,:,:)
    SST_cpl_prediddle = f3->SST_cpl_prediddle(:,:,:)
    
    ice_cov  = f3->ice_cov
    ice_cov_prediddle = f3->ice_cov_prediddle

;;;
    ssta_indian_time = conform_dims(dimsizes(SST_cpl(:,{-20:20},{30:100})), ssta_indian, (/1,2/))
    ssta_indian_prediddle_time = conform_dims(dimsizes(SST_cpl_prediddle(:,{-20:20},{30:100})), ssta_indian_prediddle, (/1,2/))

    SST_cpl(:,{-20:20},{30:100}) = (/SST_cpl(:,{-20:20},{30:100}) + ssta_indian_time/)
    printVarSummary(SST_cpl)

    SST_cpl_prediddle(:,{-20:20},{30:100}) = (/SST_cpl_prediddle(:,{-20:20},{30:100}) + ssta_indian_prediddle_time/)
    printVarSummary(SST_cpl_prediddle)


;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f /home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_addindianwarming.nc")   ; remove any pre-existing file
    ncdf = addfile("/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_addindianwarming.nc" ,"c")  ; open output netCDF file

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
       ncdf->date = date
       ncdf->datesec = datesec
 

end 

     






