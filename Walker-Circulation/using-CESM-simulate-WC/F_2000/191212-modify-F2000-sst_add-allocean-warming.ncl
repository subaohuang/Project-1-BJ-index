;;********************************************************
;;将整个热带地区海温的海温趋势，变化的幅度是通过计算sst_HadOIBl_bc_1.9x2.5_1978_2012_c130411.nc这个文件的趋势
;;得到一个新的海温数据，注意在文件中的其他数据要保留不能变化
;;用来跑F2000实验
;;lnshuheng 2019.12.12 
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

;;;------保留整个热带地区海温的异常---------
    lat1 = -30
    lat2 = 30
    lon1 = 0
    lon2 = 360
 
      SST_cpl_allocean = SST_cpl(:,{lat1:lat2},{lon1:lon2})
      SST_cpl_prediddle_allocean = SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2})

      printVarSummary(SST_cpl_allocean)

      printVarSummary(SST_cpl_allocean)

       ;;-------cesm中的海温在陆地上也有数据，所以这里不能将陆地的温度也去掉--------
       fils2 = "~/lsh/data/landsea.nc"
       f2 = addfile(fils2, "r")
       lsdata = f2->LSMASK
       lsm = landsea_mask(lsdata,SST_cpl_allocean&lat,SST_cpl_allocean&lon) 
       lsm1 = landsea_mask(lsdata,SST_cpl_prediddle_allocean&lat,SST_cpl_prediddle_allocean&lon)
       
       SST_cpl_allocean_mask = SST_cpl_allocean
       SST_cpl_allocean_mask = mask(SST_cpl_allocean, lsm.eq.1, False)

       SST_cpl_prediddle_allocean_mask = SST_cpl_prediddle_allocean
       SST_cpl_prediddle_allocean_mask = mask(SST_cpl_prediddle_allocean, lsm.eq.1, False)

       printVarSummary(SST_cpl_allocean_mask)

       SST_allocean_rc = regCoef_n(ispan(1,dimsizes(SST_cpl_allocean_mask&time),1),SST_cpl_allocean_mask,0,0)
       SST_prediddle_allocean_rc = regCoef_n(ispan(1,dimsizes(SST_cpl_prediddle_allocean_mask&time),1),SST_cpl_prediddle_allocean_mask,0,0)

       printVarSummary(SST_allocean_rc)
       print(SST_allocean_rc@_FillValue)
     
       ssta_allocean = SST_allocean_rc * dimsizes(SST_cpl_allocean_mask&time)
       ssta_allocean_prediddle = SST_allocean_rc * dimsizes(SST_cpl_allocean_mask&time)
     
       ssta_allocean = where(ismissing(ssta_allocean) ,0, ssta_allocean)                         ;;;把陆地mask掉的海温再还原
       ssta_allocean_prediddle = where(ismissing(ssta_allocean_prediddle) ,0, ssta_allocean) ;;;把陆地mask掉的海温再还原

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
    ssta_allocean_time = conform_dims(dimsizes(SST_cpl(:,{lat1:lat2},{lon1:lon2})), ssta_allocean, (/1,2/))
    ssta_allocean_prediddle_time = conform_dims(dimsizes(SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2})), ssta_allocean_prediddle, (/1,2/))

    SST_cpl(:,{lat1:lat2},{lon1:lon2}) = (/SST_cpl(:,{lat1:lat2},{lon1:lon2}) + ssta_allocean_time/)
    printVarSummary(SST_cpl)

    SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2}) = (/SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2}) + ssta_allocean_prediddle_time/)
    printVarSummary(SST_cpl_prediddle)


;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f /home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_addalloceanwarming.nc")   ; remove any pre-existing file
    ncdf = addfile("/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_addalloceanwarming.nc" ,"c")  ; open output netCDF file

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

     






