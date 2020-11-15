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
wks = gsn_open_wks("pdf", "MME的CMIP6模式-全球平均的哈德来环流的趋势")



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
    ; msf_mean = dim_avg_n_Wrap(msf(:,:,:,{180:270}), (/3/))
    msf_mean = dim_avg_n_Wrap(msf(:,:,:,:), (/3/))



    msf_mean_annual = month_to_annual(msf_mean, 1)
    msf_mean_clm = dim_avg_n_Wrap(msf_mean_annual, 0)


    dim = dimsizes(msf_mean_annual)

    ;;;计算趋势
    ; rc = msf_mean_annual(1,:,:)
    rc := regCoef_n(ispan(1,dim(0),1),msf_mean_annual(:,:,{-60:60}),0,0)

    copy_VarCoords(msf_mean_annual(1,:,{-60:60}), rc)

    rc  = (/rc*dim(0)/)


    prob := rc_test(rc)

    copy_VarCoords(msf_mean_annual(1,:,{-60:60}), prob)
    prob&level@units = "hPa"
    


;;;;;draw 
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLinesOn=True
 
    res@cnFillOn =True
   
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    ; res@cnMinLevelValF = -1.0
    ; res@cnMaxLevelValF = 1.0
    ; res@cnLevelSpacingF = 0.25

    res@cnMinLevelValF = -0.5
    res@cnMaxLevelValF = 0.5
    res@cnLevelSpacingF = 0.1


    res@lbLabelBarOn =True
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    


    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0
  
    
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
 
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(-60,60,10)

  
    res@tmYRMode="Automatic"
     
    res@vpHeightF=0.3
    res@vpWidthF=0.7

    ;;;;
    res@gsnLeftString = "Linear trend of Regional Hadley Circulation"
    res@gsnLeftStringFontHeightF = 0.015
    plot = gsn_csm_pres_hgt(wks,rc,res)


    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

   
    iplot=gsn_csm_contour(wks,prob,res_pro)


    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.05      ; make dots larger

    iplot = gsn_contour_shade(iplot, 0,0.1, opt)
    ;;;  
    overlay(plot, iplot)


    ;;;画气候态的等值线
    res@gsnRightString=""
    res@gsnLeftString=""
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = ispan(-8,8,1)
    res@cnLineLabelsOn=True
    res@cnLineLabelBackgroundColor="White"
    res@cnLineThicknessF = 2.
    res@cnFillOn = False
    res@gsnContourNegLineDashPattern=1
    res@gsnContourZeroLineThicknessF=1.0

    plot1 = gsn_csm_pres_hgt(wks, msf_mean_clm, res)
    ;;;  
    overlay(plot, plot1)

    draw(plot)
    frame(wks)
end 



