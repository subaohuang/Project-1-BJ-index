;;********************************************************
;;将模式的海温数据的大西洋地区的海温增加，增加的幅度是1摄氏度
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

    
    SST_cpl = lonFlip(f1->SST_cpl(:,:,:))
    SST_cpl_prediddle =lonFlip(f1->SST_cpl_prediddle(:,:,:))
   

  
;;;

;;;------只保留印度的海温增温趋势，其他海盆的温度趋势去掉---------
      ; lat1 = 0
      ; lat2 = 30
      lat1 = -30
      lat2 = 0
      lon1 = -70
      lon2 = 20
 
      ;; 将0-360E的数据先换为-180-180的数据  
     SST_cpl_atlantic = SST_cpl(:,{lat1:lat2},{lon1:lon2})
     SST_cpl_prediddle_atlantic = SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2})

      printVarSummary(SST_cpl_atlantic)

      printVarSummary(SST_cpl_atlantic)

       ;;-------cesm中的海温在陆地上也有数据，所以这里不能将陆地的温度也去掉--------
       fils2 = "~/lsh/data/landsea.nc"
       f2 = addfile(fils2, "r")
       lsdata = f2->LSMASK
       lsm = landsea_mask(lsdata,SST_cpl_atlantic&lat,SST_cpl_atlantic&lon) 
       lsm1 = landsea_mask(lsdata,SST_cpl_prediddle_atlantic&lat,SST_cpl_prediddle_atlantic&lon)
       
       SST_cpl_atlantic_mask = SST_cpl_atlantic
       SST_cpl_atlantic_mask = mask(SST_cpl_atlantic, lsm.eq.1, False)

       SST_cpl_prediddle_atlantic_mask = SST_cpl_prediddle_atlantic
       SST_cpl_prediddle_atlantic_mask = mask(SST_cpl_prediddle_atlantic, lsm.eq.1, False)

       printVarSummary(SST_cpl_atlantic_mask)

       SST_atlantic_rc = regCoef_n(ispan(1,dimsizes(SST_cpl_atlantic_mask&time),1),SST_cpl_atlantic_mask,0,0)
       SST_prediddle_atlantic_rc = regCoef_n(ispan(1,dimsizes(SST_cpl_prediddle_atlantic_mask&time),1),SST_cpl_prediddle_atlantic_mask,0,0)

       printVarSummary(SST_atlantic_rc)
       print(SST_atlantic_rc@_FillValue)
     
       ssta_atlantic = SST_atlantic_rc * dimsizes(SST_cpl_atlantic_mask&time)
       ssta_atlantic_prediddle = SST_atlantic_rc * dimsizes(SST_cpl_atlantic_mask&time)
     
       ssta_atlantic = where(ismissing(ssta_atlantic) ,0, 1.0)                         ;;;把陆地mask掉的海温再还原
       ssta_atlantic_prediddle = where(ismissing(ssta_atlantic_prediddle) ,0, 1.0) ;;;把陆地mask掉的海温再还原

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
    
    SST_cpl = lonFlip(f3->SST_cpl(:,:,:))
    SST_cpl_prediddle =lonFlip(f3->SST_cpl_prediddle(:,:,:))

    ice_cov  = f3->ice_cov
    ice_cov_prediddle = f3->ice_cov_prediddle

;;;
    ssta_atlantic_time = conform_dims(dimsizes(SST_cpl(:,{lat1:lat2},{lon1:lon2})), ssta_atlantic, (/1,2/))
    ssta_atlantic_prediddle_time = conform_dims(dimsizes(SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2})), ssta_atlantic_prediddle, (/1,2/))

    SST_cpl(:,{lat1:lat2},{lon1:lon2}) = (/SST_cpl(:,{lat1:lat2},{lon1:lon2}) + ssta_atlantic_time/)
    printVarSummary(SST_cpl)

    SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2}) = (/SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2}) + ssta_atlantic_prediddle_time/)
    printVarSummary(SST_cpl_prediddle)


;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f /home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_addatlanticwarming_1degree_south.nc")   ; remove any pre-existing file
    ncdf = addfile("/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_addatlanticwarming_1degree_south.nc" ,"c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications

    ;===================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
       ncdf->SST_cpl   =  lonFlip(SST_cpl)                           ; 4D               
       ncdf->SST_cpl_prediddle  = lonFlip(SST_cpl_prediddle)                        ; 3D

       ncdf->ice_cov   =  ice_cov
       ncdf->ice_cov_prediddle = ice_cov_prediddle

       ncdf->lat   =  lat
       ncdf->lon = lon 
       ncdf->time = time
       ncdf->date = date
       ncdf->datesec = datesec

end 

     






