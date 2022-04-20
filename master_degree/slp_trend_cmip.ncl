function decade_trend(data)
    local data_a, data_a_new, rc
    begin
    time=ispan(1,dimsizes(data&time),1)
    rc = regCoef_n(time,data,0,0)
    rc!0="lat"
    rc!1="lon"
    rc&lat=data&lat
    rc&lon=data&lon
    rc&lat@units = "degrees_north"
    rc&lon@units = "degrees_east"
    ; print(rc)
return(rc)
end
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

;************************************************
; 画20cr,ncep,era的msf的变化趋势。
;************************************************
begin
    diri="/mnt/f/share/"
    model="GFDL-CM3"

    fils1=systemfunc("ls /mnt/j/CMIP5/history/psl/"+model+"/r1i1p1/psl_Amon_"+model+"_historical_r1i1p1_"+"*"+".nc")
        
    ; fils2=systemfunc("ls /mnt/j/CMIP5/history_GHG/va/"+model+"/r1i1p1/va_Amon_"+model+"_historicalGHG_r1i1p1_"+"*"+".nc")
    print(fils1)
   
    yearstart=195801
    yearend=200512
    wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/"+model+"slp_history_的趋势_"+yearstart+"-"+yearend)
    f1=addfiles(fils1, "join")
    
    date=cd_calendar(f1[:]->time, 1)
    timeind=ind(date.ge.yearstart.and.date.le.yearend)
    
    slp_20cr=f1[:]->psl(timeind,{20:-20},:)
    slp_20cr=(/slp_20cr/100.0/)
    
    slp_20cr&time=date(timeind)
  delete(f1)


;;;计算趋势 

    rc1=decade_trend(slp_20cr(:,:,:))
    rc1=(/rc1*120/)
;;;显著性检验
    prob1=rc_test(rc1)
;;;


;;;画图


    gsn_define_colormap(wks, "GMT_polar")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=False
    res@gsnLeftString="20CR"
    res@gsnRightString=yearstart+"-"+yearend
    
    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnFillMode="RasterFill"
    res@cnRasterSmoothingOn=True
    
 
    res@mpDataBaseVersion = "mediumRes"
    res@mpFillAreaSpecifiers=(/"Land","water"/)      ;确定需要填充的具体位置
    res@mpSpecifiedFillColors=(/"transparent","transparent"/) ;填充位置相应的具体颜色
    res@mpOutlineOn = True; turn off outline 



    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =   -0.3
    res@cnMaxLevelValF       =   0.3
    res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = min(rc1&lat)
    res@mpMaxLatF             = max(rc1&lat)
    res@mpMinLonF             = min(rc1&lon)
    res@mpMaxLonF             = max(rc1&lon)
    res@mpCenterLonF          =avg(rc1&lon)


    res@tmYLMode="Explicit"
    res@tmYLValues=ispan(-20,20,10)
    res@tmYLLabels=(/"20~S~o~N~S","10~S~o~N~S","0~S~o~N","10~S~o~N~N","20~S~o~N~N"/)
    res@tmYLLabelFontHeightF=0.01

    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(0,360,30)
    res@tmXBLabelFontHeightF=0.01
    res@tmXBLabelFontThicknessF=0.02
    res@tmXBMinorOn=True
    res@tmXBLabels=(/"0","30~S~o~N~E","60~S~o~N~E","90~S~o~N~E","120~S~o~N~E","150~S~o~N~E",\
                "180","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W","60~S~o~N~W","30~S~o~N~W","0"/)

    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.8
    res@vpHeightF = 0.2

   
 
    
    plot=gsn_csm_contour_map(wks, rc1, res)
 


    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(1,graphic)

    iplot=gsn_csm_contour(wks,prob1,res_pro)



    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17         ; stipple
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.2         ; add extra density
    opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger


    iplot= gsn_contour_shade(iplot, 0,0.05, opt)


    overlay(plot, iplot)
    
    draw(plot)
    frame(wks)
    end