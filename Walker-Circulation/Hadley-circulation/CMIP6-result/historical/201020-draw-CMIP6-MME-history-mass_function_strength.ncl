function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
end


begin
;读入数据

; model_name = (/"ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","BCC-CSM2-MR", "BCC-ESM1",\
;               "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
;     "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
;     "GFDL-CM4", "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
;         "HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
;         "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM"/)
model_name = (/"ACCESS-CM2","ACCESS-ESM1-5", "BCC-ESM1",\
              "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
   , "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
  "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
        ,"HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
        "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM","SAM0-UNICON","UKESM1-0-LL"/)

; model_name = (/"ACCESS-ESM1-5" , "CanESM5"/)

dim_model = dimsizes(model_name)

plot = new(1,graphic)
iplot= plot
plot1 = plot


a=6371000
g=9.8
pi=get_pi("double")
dfi=(pi*10)/180

diri_output="/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/"




;;;;
   
    yearstart = 190001
    yearend = 201412
    fils1 :=  (/"/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/mme/msf_mme_history_r1i1p1f1_185001-201412.nc"/)

    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    msf = f1->msf(timeind,:,:,:)

    printVarSummary(msf)

  
    ; msf_mean = dim_avg_n_Wrap(msf(:,:,:,{240:270}), (/3/))
    msf_mean = dim_avg_n_Wrap(msf(:,:,:,{150:270}), (/3/))

    msf_mean_annual = month_to_annual(msf_mean, 1)


    ;;;;;区域平均
    
    msf_mean_annual_vermean = dim_avg_n_Wrap(msf_mean_annual(:,:,{-10:10}), (/1,2/))
   
    msf_mean_annual_vermean = (/-1*msf_mean_annual_vermean/)

    msf_mean_annual_vermean = (/dim_standardize_n(msf_mean_annual_vermean, 1, 0)/)

    msf_mean_annual_9 =  runave_n(msf_mean_annual_vermean, 11, 0, 0)
;;;;

  
    fils2 :=  (/"/home/ys17-19/lsh/Project/Walker-Circulation/Hadley-circulation/CMIP6-result/historical/Vep_timeseries_1900-2014.nc"/)

    f2 := addfile(fils2, "r")

    time:=f2->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    vep_all = f2->vep_all(:,timeind)

    vep_all_annual = month_to_annual(vep_all(time|:,model_name|:), 0)

    printVarSummary(vep_all_annual)

    vep_annual = dim_avg_n_Wrap(vep_all_annual, 1)

    vep_annual = (/dim_standardize_n(vep_annual, 1, 0)/)
  
;;;;;;

    fils3 :=  (/"/home/ys17-19/lsh/Project/Walker-Circulation/Hadley-circulation/CMIP6-result/historical/sstg_timeseries_1900-2014.nc"/)

    f3 := addfile(fils3, "r")

    time:=f3->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)

    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    sstg_all = f3->sstg_all(:,timeind)

    sstg_all_annual = month_to_annual(sstg_all(time|:,model_name|:), 0)
   
    printVarSummary(sstg_all_annual)


    print("num="+num(ismissing(sstg_all_annual)))
    
    sstg_annual = dim_avg_n(sstg_all_annual, 1)
    ; sstg_annual = sstg_all_annual(:,6)


    sstg_annual = (/dim_standardize_n(sstg_annual, 1, 0)/)


;;;


;;;计算相关系数
   corr1 = escorc_n(vep_annual, msf_mean_annual_vermean, 0, 0)
   corr2 = escorc_n(vep_annual, sstg_annual, 0, 0)
   
   print("corr1="+corr1)
   print("corr2="+corr2)

;;;


;;;计算趋势
   
    x = ispan(1900,2014,1)

    dim = dimsizes(msf_mean_annual_vermean)

    rc_msf = regCoef_n(x,msf_mean_annual_vermean,0,0)
 
    y_msf = rc_msf*x + rc_msf@yintercept

    rc_vep = regCoef_n(x,vep_annual,0,0)
 
    y_vep = rc_vep*x + rc_vep@yintercept

    rc_sstg = regCoef_n(x,sstg_annual,0,0)

    y_sstg = rc_sstg*x + rc_sstg@yintercept
;;;;


;;;计算解释趋势

   
    rc_msf_vep = regCoef_n(msf_mean_annual_vermean,vep_annual,0,0)

    rc_realte_msf = rc_msf_vep*rc_msf

    percent = rc_realte_msf/rc_vep

    ; print("解释百分比= "+percent)
 ;     print( "HCtrend="+sprintf("%4.3f",rc_msf*10))
 ; print( "VEPtrend="+sprintf("%4.3f",rc_vep*10))
 ;  print( "解释trend="+sprintf("%4.3f",rc_realte_msf*10))
 ;  exit()



    ; rc_sstg_vep = regCoef_n(sstg_annual,vep_annual,0,0)

    ; rc_realte_sstg = rc_sstg_vep*rc_sstg

    ; percent = rc_realte_sstg/rc_vep

    print("解释百分比= "+percent)

exit()
    ; print( "SSTGtrend="+sprintf("%4.3f",rc_sstg*10))
    ; print( "VEPtrend="+sprintf("%4.3f",rc_vep*10))
    ; print( "解释trend="+sprintf("%4.3f",rc_realte_sstg*10))




;;;;
    wks = gsn_open_wks("pdf", "150-90w的MME的CMIP6模式东太平洋哈德来环流强度的时间序列")
    res  = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@vpHeightF = 0.3
    res@vpWidthF = 0.7
    ; res@tiYAxisString = "Intensity"
    res@tiMainString  = "Time series of HC intensity,CEW index and Meridional SST gradients"   ; title'
    res@tiMainFontHeightF = 0.01
    res@trXMinF = 1900
    res@trXMaxF = 2014



    res@xyLineColors      = (/"black","blue","red","black","blue","red"/)  ; line colors
    res@xyLineThicknesses = (/2.0,2.0,2.0,1.0,1.0,1.0/)        ; line thicknesses
    res@xyDashPatterns    = (/0.0,0.0,0.0,1.0,1.0,1.0/)        ; line patterns
    
    ; resL = res
    ; resL@xyLineColors      := (/"red","red"/)  ; line colors
    ; resL@xyLineThicknesses := (/2.0,1.0/)        ; line thicknesses
    ; resL@xyDashPatterns    := (/0.0,1.0/)        ; line patterns
    
    
    plot = gsn_csm_xy(wks, x, (/msf_mean_annual_vermean,vep_annual,sstg_annual, y_msf,y_vep,y_sstg/), res)
    ;;;  

    

    ;;添加legend 
    gres = True
    gres@YPosPercent = 95.    ; expressed as %, 0->100, sets position of top border of legend 
                ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
    gres@XPosPercent = 8    ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)

    lineres = True
    ; lineres@lgLineColors = (/"black"，/) ; line colors
    lineres@lgLineColors       = (/"black","Blue","red"/) ; colors chosen
    lineres@lgLineThicknesses = 4.0                       ; line thicknesses
    lineres@LineLengthPercent = 9.                         ; expressed as %, 0->100, length of line

    textres = True

    textres@lgLabels = (/"HC intensity","CEW index","Meridional SST gradients"/)
    textres@lgLabelFontHeights = (/0.01,0.01,0.01/)
    ; textres@lgLabelOpacity = (/0/) 

    mid = simple_legend(wks,plot,gres,lineres,textres)   


    ;;添加相关系数

    resxy = True

    resxy@txFontHeightF = 0.01

    txid = gsn_add_text(wks,plot,"Corr(CEW,HC)=" +sprintf("%3.2f",corr1), 2000.0, 2.4 ,resxy)
    txi2 = gsn_add_text(wks,plot,"Corr(CEW,SSTG)=" +sprintf("%3.2f",corr2), 2000.0, 2.0 ,resxy)
    
    draw(plot)
    frame(wks)
   

end 



