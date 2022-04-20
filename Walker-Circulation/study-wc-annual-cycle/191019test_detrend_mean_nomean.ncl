;;********************************************************
;;test detrend mean and no mean 
;
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
    SST_cpl_prediddle =f1->SST_cpl_prediddle(:,:,:)

    ice_cov  =f1->ice_cov
    ice_cov_prediddle =f1->ice_cov_prediddle
  ;;;

  ;;;------去掉大西洋的海温增温趋势，海盆的大小定义为30N-30S ,80w(280E)-20E---------
     ;; 将0-360E的数据先换为-180-180的数据  
    SST_cpl_atlantic = SST_cpl(:,{-30:30},{-80:20})
    SST_cpl_prediddle_atlantic = SST_cpl_prediddle(:,{-30:30},{-80:20})

    printVarSummary(SST_cpl_atlantic)
    
    SST_atlantic_de = dtrend_msg_n(ispan(1,dimsizes(SST_cpl_atlantic&time),1),SST_cpl_atlantic,False,True,0)
    SST_atlantic_de_mean = dtrend_msg_n(ispan(1,dimsizes(SST_cpl_atlantic&time),1),SST_cpl_atlantic,True,True,0)

    SST_atlantic_mean = dim_avg_n_Wrap(SST_cpl_atlantic, 0)

    SST_atlantic_mean_con = conform_dims(dimsizes(SST_atlantic_de), SST_atlantic_mean, (/1,2/))
    printVarSummary(SST_atlantic_mean_con)
    exit()

    SST_atlantic_de_mean = (/SST_atlantic_de_mean + SST_atlantic_mean_con/)
    
    diff = SST_atlantic_de - SST_atlantic_de_mean

    print(diff)
    
    
 

end 

     






