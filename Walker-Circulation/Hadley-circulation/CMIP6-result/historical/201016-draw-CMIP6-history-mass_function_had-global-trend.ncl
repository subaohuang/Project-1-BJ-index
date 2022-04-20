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

plot = new(dim_model,graphic)
iplot= plot
plot1 = plot


a=6371000
g=9.8
pi=get_pi("double")
dfi=(pi*10)/180

diri_output="/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/"

wks = gsn_open_wks("pdf", "30个CMIP模式东太平洋平均的哈德来环流的趋势")


do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
    yearstart = 190001
    yearend = 201412
    fils1 :=  (/"/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/1850-1949/msf_"+model+"_history_r1i1p1f1_185001-194912.nc",\
               "/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/1950-2014/msf_"+model+"_history_r1i1p1f1_195001-201412.nc"/)

    f1 := addfiles(fils1, "r")
    
    time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    msf := f1[:]->msf(timeind,:,:,:)
    msf = (msf/10^10)
    printVarSummary(msf)
   

    msf_mean = dim_avg_n_Wrap(msf(:,:,:,{240:270}), (/3/))
    ; msf_mean = dim_avg_n_Wrap(msf(:,:,:,{150:270}), (/3/))
    ; msf_mean = dim_avg_n_Wrap(msf(:,:,:,:), (/3/))


    msf_mean_annual = month_to_annual(msf_mean, 1)
    msf_mean_clm = dim_avg_n_Wrap(msf_mean_annual, 0)


    dim = dimsizes(msf_mean_annual)

   ;;;计算趋势
    ; rc = msf_mean_annual(1,:,:)
    rc := regCoef_n(ispan(1,dim(0),1),msf_mean_annual(:,:,{-30:30}),0,0)
    
     copy_VarCoords(msf_mean_annual(1,:,{-30:30}), rc)

     rc  = (/rc*dim(0)/)


    prob := rc_test(rc)

    copy_VarCoords(msf_mean_annual(1,:,{-30:30}), prob)
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
    res@cnMinLevelValF = -1.2
    res@cnMaxLevelValF = 1.2
    res@cnLevelSpacingF = 0.2
    ;  res@cnMinLevelValF = -0.8
    ; res@cnMaxLevelValF = 0.8
    ; res@cnLevelSpacingF = 0.1

    res@lbLabelBarOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    res@cnLineThicknessF = 2.0

    res@tmXBMajorThicknessF=1.5
    res@tmYLMajorThicknessF=1.5
    res@tmBorderThicknessF=1.5
   
    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02

    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
 
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(-30,30,5)

  
    res@tmYRMode="Automatic"
     
    res@vpHeightF=0.3
    res@vpWidthF=0.7

    ;;;;
    res@gsnLeftString = model
    res@gsnLeftStringFontHeightF = 0.02
    plot(k) = gsn_csm_pres_hgt(wks,rc,res)


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

   
    iplot(k)=gsn_csm_contour(wks,prob,res_pro)


    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.05      ; make dots larger

    iplot(k) = gsn_contour_shade(iplot(k), 0,0.1, opt)
    ;;;  
    overlay(plot(k), iplot(k))


    ;;;画气候态的等值线
    res@gsnRightString=""
    res@gsnLeftString=""
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = ispan(-8,8,2)
    res@cnLineLabelsOn=True
    res@cnLineLabelBackgroundColor="Transparent"
    res@cnLineThicknessF = 2.
    res@cnFillOn = False
    res@gsnContourNegLineDashPattern=1
    res@gsnContourZeroLineThicknessF=2.2

    plot1(k) = gsn_csm_pres_hgt(wks, msf_mean_clm, res)
  
 
    ;;;  
    overlay(plot(k), plot1(k))
 print("********************************finish calculate"+" "+model+"************************************")
 end do 


pres = True
pres@gsnMaximize = True
pres@gsnPanelLabelBar = True
pres@pmLabelBarHeightF = 0.05
pres@pmLabelBarWidthF  = 0.6
pres@pmLabelBarOrthogonalPosF = -0.02
; pres@lbOrientation = "Vertical"
; pres@lbLabelFontHeightF = 0.01
gsn_panel(wks, plot, (/6,5/), pres)
end 



