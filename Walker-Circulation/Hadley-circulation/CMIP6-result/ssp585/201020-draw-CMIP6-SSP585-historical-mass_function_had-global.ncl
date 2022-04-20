
function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
end 

begin
;读入数据

; model_name = (/"ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","BCC-CSM2-MR", "BCC-ESM1",\
;               "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
;     "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
;     "GFDL-CM4", "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
;         "HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
;         "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM"/)
; model_name = (/"ACCESS-CM2","ACCESS-ESM1-5", "BCC-ESM1",\
;               "CAMS-CSM1-0","CanESM5","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-ESM2-1"\
;    , "E3SM-1-0", "EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L",  "FGOALS-g3", \
;   "GFDL-ESM4",    "GISS-E2-1-G",  "GISS-E2-1-G-CC",   "GISS-E2-1-H"\
;         ,"HadGEM3-GC31-LL","IPSL-CM6A-LR",   "MCM-UA-1-0","MIROC6", \
;         "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3","NorCPM1",  "NorESM2-LM","SAM0-UNICON","UKESM1-0-LL"/)

model_name = (/ "CESM2","CESM2-WACCM","CNRM-ESM2-1"\
   , "EC-Earth3","EC-Earth3-Veg", "GFDL-ESM4","HadGEM3-GC31-LL","IPSL-CM6A-LR","MIROC6", \
        "MIROC-ES2L","MPI-ESM1-2-HR","MRI-ESM2-0","NESM3",  "NorESM2-LM","UKESM1-0-LL"/)
; model_name = (/"ACCESS-ESM1-5" , "CanESM5"/)

dim_model = dimsizes(model_name)

plot = new(dim_model,graphic)
iplot= plot
plot1 = plot


a=6371000
g=9.8
pi=get_pi("double")
dfi=(pi*10)/180

wks = gsn_open_wks("pdf", "CMIP6模式洋哈德来环流SSP - Historical")


do k = 0,dim_model-1
;;;;
    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
    yearstart1 = 190001
    yearend1 = 195012
    fils1 :=  (/"/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/1850-1949/msf_"+model+"_history_r1i1p1f1_185001-194912.nc",\
               "/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/1950-2014/msf_"+model+"_history_r1i1p1f1_195001-201412.nc"/)
    
   
    fils2 :=  (/"/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/ssp585/msf_"+model+"_ssp585_r1i1p1f1_201501-210012.nc"/)
    


    f1 := addfiles(fils1, "r")
    
    time:=f1[:]->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend1.and.date.ge.yearstart1)
    delete(date)     


    msf_his := f1[:]->msf(timeind,:,:,:)
    msf_his = (msf_his/10^10)
    printVarSummary(msf_his)


    f2 := addfile(fils2, "r")

    time := f2->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)

    yearstart2 = 204501
    yearend2 = 210012

    timeind:=ind(date.le.yearend2.and.date.ge.yearstart2)
    delete(date)     

    msf_ssp585 := f2->msf(timeind,:,:,:)
    msf_ssp585 = (msf_ssp585/10^10)
    printVarSummary(msf_ssp585)



    msf_his_mean = dim_avg_n_Wrap(msf_his(:,:,:,{240:270}), (/3/))
    msf_his_mean_annual = month_to_annual(msf_his_mean, 1)
    msf_his_mean_annual_clm = dim_avg_n_Wrap(msf_his_mean_annual,0)
    ; msf_mean_clm = dim_avg_n_Wrap(msf_his_mean_annual, 0)


    msf_ssp585_mean = dim_avg_n_Wrap(msf_ssp585(:,:,:,{240:270}), (/3/))
    msf_ssp585_mean_annual = month_to_annual(msf_ssp585_mean, 1)
    msf_ssp585_mean_annual_clm = dim_avg_n_Wrap(msf_ssp585_mean_annual,0)


   ; ;;;计算趋势
   ;  ; rc = msf_mean_annual(1,:,:)
   ;  rc := regCoef_n(ispan(1,dim(0),1),msf_mean_annual(:,:,{-60:60}),0,0)
    
   ;   copy_VarCoords(msf_mean_annual(1,:,{-60:60}), rc)

   ;   rc  = (/rc*dim(0)/)


   ;  prob := rc_test(rc)

   ;  copy_VarCoords(msf_mean_annual(1,:,{-60:60}), prob)
   ;  prob&level@units = "hPa"
  
    msf_diff = msf_his_mean_annual_clm(:,{-30:30})
    msf_diff = (/msf_ssp585_mean_annual_clm(:,{-30:30}) - msf_his_mean_annual_clm(:,{-30:30})/)

 prob = t_test(msf_his_mean_annual,msf_ssp585_mean_annual)

;;;;;draw 
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLinesOn=True
 
    res@cnFillOn =True
   
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -2.0
    res@cnMaxLevelValF = 2.0
    res@cnLevelSpacingF = 0.5

    res@lbLabelBarOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False

    res@tmXBMajorThicknessF=2.0
    res@tmYLMajorThicknessF=2.0
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
     
    res@vpHeightF=0.4
    res@vpWidthF=0.6

    ;;;;
    res@gsnLeftString = model
    res@gsnLeftStringFontHeightF = 0.02
    plot(k) = gsn_csm_pres_hgt(wks,msf_diff,res)


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
    res@gsnContourZeroLineThicknessF=2.0

    plot1(k) = gsn_csm_pres_hgt(wks, msf_his_mean_annual_clm, res)
  
 
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



