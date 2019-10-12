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
    
    SST_cpl = lonFlip(f1->SST_cpl(:,:,:))
 
    SST_cpl_prediddle =lonFlip(f1->SST_cpl_prediddle(:,:,:))

    ice_cov  =lonFlip(f1->ice_cov)
    ice_cov_prediddle =lonFlip(f1->ice_cov_prediddle)
    

    
  ;;;------去掉大西洋的海温增温趋势，海盆的大小定义为30N-30S ,80w(280E)-20E---------

     ;; 将0-360E的数据先换为-180-180的数据  
     SST_cpl_atlantic = SST_cpl(:,{-30:60},{-80:20})
     SST_cpl_prediddle_atlantic = SST_cpl_prediddle(:,{-30:60},{-80:20})

      printVarSummary(SST_cpl_atlantic)

     ;;-------cesm中的海温在陆地上也有数据，所以这里不能将陆地的温度也去掉--------
       fils2 = "~/lsh/data/landsea.nc"
       f2 = addfile(fils2, "r")
      
       lsdata = f2->LSMASK

       lsm = landsea_mask(lsdata,SST_cpl_atlantic&lat,SST_cpl_atlantic&lon) 
       lsm1 = landsea_mask(lsdata,SST_cpl_prediddle_atlantic&lat,SST_cpl_prediddle_atlantic&lon)
       printVarSummary(lsm)
       printVarSummary(SST_cpl_atlantic)
       
       SST_cpl_atlantic_mask = SST_cpl_atlantic
       SST_cpl_atlantic_mask = mask(SST_cpl_atlantic, lsm.eq.1, False)
   

       SST_cpl_prediddle_atlantic_mask = SST_cpl_prediddle_atlantic
       SST_cpl_prediddle_atlantic_mask = mask(SST_cpl_prediddle_atlantic, lsm.eq.1, False)

       printVarSummary(SST_cpl_atlantic_mask)


       SST_atlantic_de = dtrend_msg_n(ispan(1,dimsizes(SST_cpl_atlantic_mask&time),1),SST_cpl_atlantic_mask,False,True,0)
       SST_prediddle_atlantic_de = dtrend_msg_n(ispan(1,dimsizes(SST_cpl_prediddle_atlantic_mask&time),1),SST_cpl_prediddle_atlantic_mask,False,True,0)


       printVarSummary(SST_atlantic_de)
       print(SST_atlantic_de@_FillValue)
       SST_atlantic_de = where(ismissing(SST_atlantic_de) ,SST_cpl_atlantic, SST_atlantic_de)                         ;;;把陆地mask掉的海温再还原
       SST_prediddle_atlantic_de = where(ismissing(SST_prediddle_atlantic_de) ,SST_cpl_prediddle_atlantic, SST_prediddle_atlantic_de) ;;;把陆地mask掉的海温再还原
       
      

  ;;;----------------------------------------------------------------------
       

    ;;;-----将原来的数据中的大西洋数据用去趋势的海温代替--------------  
    
    SST_cpl(:,{-30:60},{-80:20}) = (/SST_atlantic_de/)
    printVarSummary(SST_cpl)
    
    SST_cpl_prediddle(:,{-30:60},{-80:20}) = (/SST_prediddle_atlantic_de/)
    printVarSummary(SST_cpl_prediddle)
      


    ;;-----------------------写入nc文件---------------------------
    system("/bin/rm -f ~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_atlantictrend.nc")   ; remove any pre-existing file
    ncdf = addfile("~/lsh/sst_HadOIBl_bc_1.9x2.5_1978_2012_remove_atlantictrend.nc" ,"c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications

    ;===================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
       ncdf->SST_cpl   =  lonFlip(SST_cpl)                                         
       ncdf->SST_cpl_prediddle  = lonFlip(SST_cpl_prediddle)

       ncdf->ice_cov   =  lonFlip(ice_cov)
       ncdf->ice_cov_prediddle =lonFlip(ice_cov_prediddle)

       ncdf->lat   =  lat
       ncdf->lon = lon 
       ncdf->time = time
 

end 

     






