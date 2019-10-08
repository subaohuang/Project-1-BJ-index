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

  ;;;------去掉印度的海温增温趋势，海盆的大小定义为20N-20S ,30E-100E---------

     SST_cpl_indian = SST_cpl(:,{-20:20},{30:100})
     SST_cpl_prediddle_indian = SST_cpl_prediddle(:,{-20:20},{30:100})

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

       SST_indian_de = dtrend_msg_n(ispan(1,dimsizes(SST_cpl_indian_mask&time),1),SST_cpl_indian_mask,False,True,0)
       SST_prediddle_indian_de = dtrend_msg_n(ispan(1,dimsizes(SST_cpl_prediddle_indian_mask&time),1),SST_cpl_prediddle_indian_mask,False,True,0)


       printVarSummary(SST_indian_de)
       print(SST_indian_de@_FillValue)
       SST_indian_de = where(ismissing(SST_indian_de) ,SST_cpl_indian, SST_indian_de)                         ;;;把陆地mask掉的海温再还原
       SST_prediddle_indian_de = where(ismissing(SST_prediddle_indian_de) ,SST_cpl_prediddle_indian, SST_prediddle_indian_de)                         ;;;把陆地mask掉的海温再还原
  ;;;----------------------------------------------------------------------
       

  ;;;-----将原来的数据中的印度洋数据用去趋势的海温代替--------------  
    
    SST_cpl(:,{-20:20},{30:100}) = (/SST_indian_de/)
    printVarSummary(SST_cpl)
    
    SST_cpl_prediddle(:,{-20:20},{30:100}) = (/SST_prediddle_indian_de/)
    printVarSummary(SST_cpl_prediddle)
      


    ;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f ~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_indiantrend.nc")   ; remove any pre-existing file
    ncdf = addfile("~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_indiantrend.nc" ,"c")  ; open output netCDF file

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

     






