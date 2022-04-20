function decade_trend(data)
    local data_a, data_a_new, rc
    begin
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
;;****************************************************************
;;main program
;;*************************************************************

begin
 diri="/mnt/f/share/CMIP5/msf_history/"
 model="IPSL-CM5A-LR"

 do  i=0,dimsizes(xi_model)-1
 fils1=systemfunc("ls "+diri+"msf_time_"+model+"_history_r1i1p1_197901-200512.nc")
 fils2=systemfunc("ls "+diri+"msf_time_"+model+"_history_r1i1p1_197901-200512.nc")


    fils=new(2,string)
    fils(0)=fils1
    fils(1)=fils2
    f=addfile(fils, "r")
    date=cd_calendar(f[:]->time, 1)
    print(date)

    yearstart=197901
    yearend=201212

    timeind=ind(date.le.yearend.and.date.ge.yearstart)

    msf_time=f1->msf_time(:,:,timeind)
    msf_time=(/msf_time/10^9/)
    msf_time&level@units="hpa"
    ; msf_time&time=ispan(0, dimsizes(timeind)-1, 1)
    delete(msf_time&time)
    msf_time&time=date(timeind)
    ; msf_clm=dim_avg_n_Wrap(msf_time, 2)

    ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    msf_son=month_to_season(msf_time(time|:,level|:,lon|:), "SON")
    msf_djf=month_to_season(msf_time(time|:,level|:,lon|:), "DJF")

    msf_son&time=ispan(1979,2005,1)
    msf_djf&time=ispan(1979,2005,1)  ;;;因为这里第一个是1958 1月和2月的平均，第二个是195812,195901,195902的平均
    ;;注意点
    msf_djf(0:dimsizes(msf_djf&time)-2,:,:)=(/msf_djf(1:,:,:)/)
    msf_djf(dimsizes(msf_djf&time)-1,:,:)=(/msf_time(:,:,{200512})/)
    ;;;

    msf_clm=dim_avg_n((/msf_son(:,:,:)+msf_djf(:,:,:)/)*0.5,0)
    copy_VarCoords(msf_son(1,:,:), msf_clm)
;;;




 ;;season
 ; msf_clm=dim_avg_n_Wrap(month_to_season(msf_time(time|:,level|:,lon|:),"SON"), 0)
  ;;
 rc1=decade_trend(msf_time(:,:,:))
 prob1=rc_test(rc1)
 rc1=(/rc1*120/)
 
 
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/"+model+"GHG_msf_"+yearstart+"-"+yearend+"")
;********************************************************
;;draw
;********************************************************

    cmap=read_colormap_file("MPL_PuOr")
    cmap=(/cmap(::-1,0:3)/)

    res2 = True
    res2@gsnDraw=False
    res2@gsnFrame=False
    res2@gsnLeftString=model
    res2@gsnRightString=""+yearstart+"-"+yearend
   
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
    res2@lbLabelBarOn=True

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
        res2@gsnContourZeroLineThicknessF=4.0
        res2@gsnLeftString=""
    
        plot1 = gsn_csm_pres_hgt(wks, msf_clm, res2)

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
 
    dum1=new(1,graphic)
    dum2=dum1
    dum3=dum1

     i=0
    ployres@gsLineDashPattern=0
    dum1(i)= gsn_add_polyline(wks,plot(i),ploylon1,ploylevel1,  ployres)
    ployres@gsLineDashPattern=16
    dum2(i)= gsn_add_polyline(wks,plot(i),ploylon2,ploylevel2,  ployres)
    ployres@gsLineDashPattern=0
    dum3(i)= gsn_add_polyline(wks,plot(i),ploylon3,ploylevel3,  ployres)
     

        overlay(plot,plot1)
        overlay(plot,iplot)
        draw(plot)
        frame(wks)
  end 