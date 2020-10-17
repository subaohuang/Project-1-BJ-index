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
; model_name = (/ "GISS-E2-1-H",  "MIROC6", "NorCPM1", "UKESM1-0-LL",\
; "ACCESS-ESM1-5",  "CanESM5" , "E3SM-1-0"    ,"HadGEM3-GC31-LL" , "MIROC-ES2L",  "NorESM2-LM" ,\
;  "CESM2",  "EC-Earth3"  , "GFDL-ESM4"  ,"MPI-ESM1-2-HR" ,\
;  "CESM2-WACCM" , "EC-Earth3-Veg",  "GISS-E2-1-G", "IPSL-CM6A-LR" ,"MRI-ESM2-0","SAM0-UNICON",\
; "BCC-ESM1","CNRM-CM6-1","FGOALS-f3-L","GISS-E2-1-G-CC","MCM-UA-1-0","NESM3"/)

model_name = (/"ACCESS-ESM1-5" , "CanESM5"/)

dim_model = dimsizes(model_name)

plot = new(dim_model,graphic)


a=6371000
g=9.8
pi=get_pi("double")
dfi=(pi*10)/180

diri_output="/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/"

    wks = gsn_open_wks("pdf", "东太平洋哈德来环流的趋势")
do k = 0,dim_model-1
;;;;

    model = model_name(k)

    print("*******************start calculate"+" "+model+"**********************************************")
    yearstart = 195001
    yearend = 201412
    fils1 :=  "/home/yangsong3/data-observation/linshh/CMIP6-derived/msf/msf_"+model+"_history_r1i1p1f1_195001-201412.nc"
    f1 := addfile(fils1, "r")
    
    time:=f1->time    ;;;这里用了20cr的time，因为jra的time是一个个来的
    date:=cd_calendar(time, 1)
    delete(time)
    timeind:=ind(date.le.yearend.and.date.ge.yearstart)
    delete(date)     
    msf := f1->msf(timeind,:,:,:)
    msf = (msf/10^10)
    printVarSummary(msf)
   

    msf_mean = dim_avg_n_Wrap(msf(:,:,:,{240:270}), (/3/))


    msf_mean_annual = month_to_annual(msf_mean, 1)
    

  dim = dimsizes(msf_mean_annual)

   ;;;计算趋势
    ; rc = msf_mean_annual(1,:,:)
    rc := regCoef_n(ispan(1,dim(0),1),msf_mean_annual,0,0)
    
     copy_VarCoords(msf_mean_annual(1,:,:), rc)

     rc  = (/rc*dim(0)/)


    prob = rc_test(rc)

;;;;;draw 
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlRe"

    res@cnLinesOn=True
 
    res@cnFillOn =True
   
    res@cnLinesOn=True
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -1.2
    res@cnMaxLevelValF = 1.2
    res@cnLevelSpacingF = 0.2

    res@lbLabelBarOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
  
    
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
 
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(-90,90,10)

  
    res@tmYRMode="Automatic"
     
    res@vpHeightF=0.3
    res@vpWidthF=0.7

    ;;;;

    
    res@gsnLeftString = "Spring"
    plot(k) = gsn_csm_pres_hgt(wks,rc,res)
 print("********************************finish calculate"+" "+model+"************************************")
 end do 
pres = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot, (/2,1/), pres)
end 



