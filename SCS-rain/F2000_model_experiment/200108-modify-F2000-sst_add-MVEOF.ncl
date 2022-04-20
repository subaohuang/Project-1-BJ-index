;;********************************************************
;;将MVEOF得到的第一模态和第二模态的海温异常添加到气候态的海温异常中
;;F2000实验，修改文件sst_HadOIBl_bc_1.9x2.5_1978_2012_c130411.nc
;;得到一个新的海温数据，注意在文件中的其他数据要保留不能变化
;;用来跑F2000实验
;;lnshuheng 2020.01.08
;;

;**********************************************************
;test
begin
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

;;;-------------读入MV-EOF的海温异常--------------------------
    fils4 = "/home/ys17-19/lsh/Project/SCS-rain/ev_sst_ghg.nc"
    f4= addfile(fils4, "r")
    
    ev_sst = f4->ev_sst
    ev_sst = where(ev_sst.ge.2.0, 2.0, ev_sst)
    ev_sst = where(ev_sst.le.-2.0, -2.0, ev_sst)
    ev_sst = where(ismissing(ev_sst), 0, ev_sst)

    printVarSummary(ev_sst)
   

    opt = True
    opt@fout = "ev_sst"
    write_matrix (ev_sst, "144f7.2", opt)


    ;;;-----添加MV-EOF的异常海温的区域---------
    lat1 = -60
    lat2 = 60
    lon1 = 0
    lon2 = 360
;;;
 
    
    ev_sst_time = conform_dims(dimsizes(SST_cpl(:,{lat1:lat2},{lon1:lon2})), ev_sst, (/1,2/))
    ev_sst_time = conform_dims(dimsizes(SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2})), ev_sst, (/1,2/))

    SST_cpl(:,{lat1:lat2},{lon1:lon2}) = (/SST_cpl(:,{lat1:lat2},{lon1:lon2}) + ev_sst_time/)
    printVarSummary(SST_cpl)

    SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2}) = (/SST_cpl_prediddle(:,{lat1:lat2},{lon1:lon2}) + ev_sst_time/)
    printVarSummary(SST_cpl_prediddle)
;;;



;;-----------------------写入nc文件---------------------------
    system("rm /home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_mvEOF-ghg.nc")   ; remove any pre-existing file
    ncdf = addfile("/home/ys17-19/lsh/CESM-data/modsst/sst_HadOIBl_bc_1.9x2.5_clim_c061031_mvEOF-ghg.nc" ,"c")  ; open output netCDF file

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

     






