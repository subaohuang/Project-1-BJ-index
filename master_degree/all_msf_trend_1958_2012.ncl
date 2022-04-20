function decade_trend(data)
    local data_a, data_a_new, rc
    begin

    ; data_new=doubletofloat(dim_standardize_n_Wrap(data(::-1,:,:),1,(/0,1/)))       ;标准化
    ; print(data_new)
    ; copy_VarMeta(data(::-1,:,:), data_new)
    ;-----------------------------------------------去掉年循环
    ; data_ano=rmMonAnnCycLLT(data_a_new)
    ; printVarSummary(data_ano)
    ;-----------------------------------------------
    rc = regCoef_n(data&time,data,0,2)
    rc!0="lev"
    rc!1="lon"
    rc&lev=data&level
    rc&lon=data&lon
    rc&lev@units = "hpa"
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
; 画20cr,ncep,era40及jRA-55的msf的变化趋势。
;************************************************
undef("draw_msf_trend")
function draw_msf_trend(fils1,wks,LeftString)
begin
    diri="/mnt/f/share/"
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.195801.and.date.le.201212)
    
    msf_time_20cr=f1->msf_time(:,:,timeind)
    ud_time_20cr=f1->ud_time(:,:,timeind)
    w_time_20cr=f1->w_time(:,:,timeind)
    msf_time_20cr=(/msf_time_20cr/10^9/)

    delete(msf_time_20cr&time)
    msf_time_20cr&time=date(timeind)
;;;
   msf_clm_20cr=dim_avg_n_Wrap(msf_time_20cr, 2)
  
;;;计算趋势 

    ; rc1=decade_trend(msf_time_20cr(:,:,{195801:197812}))

    ; rc1=decade_trend(msf_time_20cr(:,:,{197901:201212}))

    rc1=decade_trend(msf_time_20cr(:,:,{195801:201212}))
    prob1=rc_test(rc1)
    rc1=(/rc1*120/)
    ; rc2=(/rc2*120/)
    ; rc3=(/rc3*120/)

;;;

;;;显著性检验
    prob1=rc_test(rc1)
    ; prob2=rc_test(rc2)
    ; prob3=rc_test(rc3)
;;;

;;;画图
    ; gsn_define_colormap(wks,"GMT_polar")
    cmap=read_colormap_file("MPL_PuOr")
    
    cmap=(/cmap(::-1,0:3)/)
    
    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString=LeftString
    res2@gsnRightString="1958-2012"
   
    res2@cnLinesOn=False

    res2@cnLevelSelectionMode ="ManualLevels"
    res2@cnMinLevelValF = -.08
    res2@cnMaxLevelValF = 0.08
    res2@cnLevelSpacingF = 0.01
    res2@cnFillPalette=cmap

    res2@cnLineLabelFontHeightF=0.01
    res2@cnFillOn = True 
    
    res2@cnLineLabelsOn=False
    res2@cnInfoLabelOn=False
    res2@lbLabelBarOn=False

    res2@tiYAxisString=""
    res2@tmXTOn=False
    res2@tmYROn=False
    res2@tmXBMode="Explicit"
  
    res2@tmXBValues=(/0,90,180,270,357.5/)
    res2@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res2@tmXBMinorOn=True
    res2@tmXBMinorValues=ispan(0,360,10)
    res2@tmYRMode="Automatic"
     
    res2@vpHeightF=0.3
    res2@vpWidthF=0.7
   
    plot=new(1,graphic)
    
    plot=gsn_csm_pres_hgt(wks,rc1,res2)
  
;;;gsn_contour_shade
        res_pro=True
        res_pro@gsnDraw=False
        res_pro@gsnFrame=False
        res_pro@gsnRightString=""
        res_pro@gsnLeftString=""
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
    iplot=gsn_csm_pres_hgt(wks,prob1,res_pro)

        opt                      = True
        opt@gsnShadeFillType     = "pattern"
        opt@gsnShadeMid     = 17         ; stipple
        ; opt@gsnShadeHigh   = "green"
        opt@gsnShadeFillScaleF   = 0.2         ; add extra density
        opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

        iplot = gsn_contour_shade(iplot(0), 0,0.05, opt)

 ;;;
        res2@cnLevelSelectionMode = "ExplicitLevels"
        res2@cnLinesOn=True
        res2@cnLevels = ispan(-6,6,1)
        res2@cnLineLabelsOn=True
        res2@cnLineLabelBackgroundColor="Transparent"
        res2@cnLineThicknessF = 1.
        res2@cnFillOn = False
        res2@gsnContourNegLineDashPattern=1
        res2@gsnContourZeroLineThicknessF=2.0
        res2@gsnLeftString=""
        

        plot1 = gsn_csm_pres_hgt(wks, msf_clm_20cr, res2)
        
        overlay(plot,plot1)
        overlay(plot,iplot)
return(plot)
end

;***********************************************************
;主程序
;***********************************************************         
begin

fils1="ncep/20CR/msf_time_20cr_global_6--6.nc"
fils2="ncep/ncep1/msf_time_ncep1_global.nc"
fils3="era/era40/msf_time_era40_in_global.nc"
fils4="JRA-55/msf_time_jra55_global.nc"

wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/4个再分析资料计算的msf的趋势")

plot=new(4,graphic)
plot(0)=draw_msf_trend(fils1,wks,"20CR")
plot(1)=draw_msf_trend(fils2,wks,"NCPE1")
plot(2)=draw_msf_trend(fils3,wks,"ERA")
plot(3)=draw_msf_trend(fils4,wks,"JRA-55")

;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 8.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=1
   ;---Attach the polyline
    ploylevel=1000
    ploylon1=fspan(9, 39, 25)
    ploylon2=fspan(100, 149, 10)
    ploylon3=fspan(280, 315, 25)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
 
    dum1=new(4,graphic)
    dum2=dum1
    dum3=dum1

    do i=0,3
    ployres@gsLineDashPattern=0
    dum1(i)= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
    end do 

    pres=True
    pres@gsnPanelLabelBar=True
    pres@pmLabelBarWidthF=0.8
    gsn_panel(wks, plot, (/2,2/), pres)
end



