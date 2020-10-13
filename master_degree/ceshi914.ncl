;***************************************************************
  ;变化趋势平面图
;***************************************************************
function decade_trend(data,units,long_name)
    local data_a, data_a_new, rc, data_ano
    begin

    if (data@long_name.eq."Sea surface temperature") then
    data_a=data@scale_factor*data+data@add_offset-274.15
    else 
    data_a=data@scale_factor*data+data@add_offset
    end if
    
    ; data_a2=dim_standardize_n_Wrap(data_a,1,0)        ;标准化
    ; data_a1=doubletofloat(data_a2)

    ; data_a=data
    copy_VarMeta(data, data_a)
    ; data_a1&time=ispan(1,456,1)
    data_a&time=ispan(1,dimsizes(data&time),1)
    data_a@units=units
    data_a@long_name=long_name
    data_a_new = data_a(latitude|:,longitude|:,time|:)
    printVarSummary(data_a_new)
    ;-----------------------------------------------去掉年循环
    data_ano=rmMonAnnCycLLT(data_a_new)
    printVarSummary(data_ano)
    ;-----------------------------------------------
    rc = regCoef(data_ano&time,data_ano)
    rc!0="lat"
    rc!1="lon"
    rc&lat=data_ano&latitude
    rc&lon=data_ano&longitude
    rc&lat@units = "degrees_north"
    rc&lon@units = "degrees_east"
    
    tval = onedtond(rc@tval , dimsizes(rc))
    df   = onedtond(rc@nptxy, dimsizes(rc)) - 2
    b = tval    ; b must be same size as tval (and df)
    b = 0.5
    prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   
    copy_VarCoords(rc,prob )
    printVarSummary(prob)
    rc=where(prob.lt.0.05, rc, rc@_FillValue)

return(rc)
end
;*************************************************************************
 ;draw plot
;************************************************************************

function  draw_plot(rc,wks,cnMinLevelValF,cnMaxLevelValF,cnLevelSpacingF,vpYF)
    local res, plot
    begin

    res = True
    res@gsnAddCyclic=False
    res@gsnDraw=False
    res@gsnFrame=False

    gsn_define_colormap(wks,"NCV_blue_red")
    res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
    res@cnFillOn = True 
    res@cnLinesOn = False
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       =  cnMinLevelValF
    res@cnMaxLevelValF       =   cnMaxLevelValF
    res@cnLevelSpacingF      =  cnLevelSpacingF
    res@cnFillDrawOrder="PreDraw"
    res@lbOrientation = "Vertical"
    res@lbBoxEndCapStyle="TriangleBothEnds"
    ; res@lbLabelStrings=ispan(-4,4,1)
    ; res@lbTitleFontHeightF=0.015              ;给colorbar添加坐标单位，，u10和sst需要
    ; res@lbTitleString="*10~S~-3~N"
    ; res@lbTitleExtentF=0.01

    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(100,270,40)
    res@tmYLTickSpacingF=10
    res@tmYLTickStartF=-30
    res@tmYLTickEndF=30
    res@tmYROn=False
    res@tmXTOn=False
    res@tmYLLabelFontHeightF=0.01
    res@tmXBLabelFontHeightF=0.01
    res@tmXBLabelFontThicknessF=0.02
    res@tmXBMinorOn=True
    res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W"/)
    res@tmYLLabels=(/"30~S~o~N~S","20~S~o~N~S","10~S~o~N~S","0~S~o~N","10~S~o~N~N","20~S~o~N~N","30~S~o~N~N"/)

    res@vpHeightF=0.2249
    res@vpWidthF=0.2249
    res@vpXF=0.04
    res@vpYF=vpYF

    ; res@vpKeepAspect =False

    ; res@mpDataSetName = "Earth..4"
    res@mpDataBaseVersion = "mediumRes"
    res@mpFillAreaSpecifiers=(/"Land","water"/)      ;确定需要填充的具体位置
    res@mpSpecifiedFillColors=(/"white","transparent"/) ;填充位置相应的具体颜色

    res@mpOutlineOn = True; turn off outline 
    ; res@mpGeophysicalLineThicknessF= 4.          ; double the thickness of geophysical boundaries 
    res@mpNationalLineThicknessF= 2.               ; double the thickness of national boundaries 
    res@mpOutlineSpecifiers = (/"Land"/)
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = min(rc&lat)
    res@mpMaxLatF             = max(rc&lat)
    res@mpMinLonF             = min(rc&lon)
    res@mpMaxLonF             = max(rc&lon)
    res@mpCenterLonF          =avg(rc&lon)
   
    plot = gsn_csm_contour_map(wks, rc, res)

return(plot)
end 

function area_avg(data)
 local data_a,data_a1,data_avg
 begin 
 if (data@long_name.eq."Sea surface temperature") then
 data_a=data@scale_factor*data+data@add_offset-274.15
 else 
 data_a=data@scale_factor*data+data@add_offset
 end if

 data_a1=doubletofloat(data_a)
 copy_VarMeta(data,data_a1 )
 data_a_new = data_a1(latitude|:,longitude|:,time|:)

 data_a_new&time=ispan(1,432,1)
 
 

  data_mask=where(ismissing(data_a_new(ind(data_a_new&latitude.eq.10.0):ind(data_a_new&latitude.eq.5.0),ind(data_a_new&longitude.eq.260.0):ind(data_a_new&longitude.eq.290.0),:)).ne.True,data_a_new@missing_value,(data_a_new(ind(data_a_new&latitude.eq.10.0):ind(data_a_new&latitude.eq.5.0),ind(data_a_new&longitude.eq.260.0):ind(data_a_new&longitude.eq.290.0),:)))
  copy_VarMeta((data_a_new(ind(data_a_new&latitude.eq.10.0):ind(data_a_new&latitude.eq.5.0),ind(data_a_new&longitude.eq.260.0):ind(data_a_new&longitude.eq.290.0),:)), data_mask)
  printVarSummary(data_mask)
  data_a_new(ind(data_a_new&latitude.eq.10.0):ind(data_a_new&latitude.eq.5.0),ind(data_a_new&longitude.eq.260.0):ind(data_a_new&longitude.eq.290.0),:)=data_mask

 data_avg_area1=dim_avg_n_Wrap(data_a_new(:,0:ind(data_a_new&longitude.eq.195.0),:), (/0,1/))
 data_avg_area2=dim_avg_n_Wrap(data_a_new(:,(ind(data_a_new&longitude.eq.195.0)+1):380,:), (/0,1/))

 data_avg=new((/2,dimsizes(data_a_new&time)/),double)

 data_avg(0,:)=data_avg_area1
 data_avg(1,:)=data_avg_area2

 data_avg_new=doubletofloat(data_avg)
 data_avg_new!0 = "area"
 data_avg_new!1 = "time"
 data_avg_new&area= (/"dong","xi"/)
 data_avg_new&time=data_a_new&time
 ;平滑处理
  data_avg_new_a=runave_n_Wrap(data_avg_new, 5, 0, 1)


return(data_avg_new_a)
end 

function gain_data(data)
 begin
 sst_avg_area=area_avg(data)

 sst_rc1=regCoef(sst_avg_area&time,sst_avg_area(0,:))

 sst_rc2=regCoef(sst_avg_area&time,sst_avg_area(1,:))

 data_fuc=new((/2,3,dimsizes(sst_avg_area&time)/),float)


 data_fuc(0,0,:) = sst_avg_area(0,:)
 data_fuc(0,1,:) = sst_rc1*sst_avg_area&time+sst_rc1@yintercept
 data_fuc(0,2,:) =avg(sst_avg_area(0,:))

 data_fuc(1,0,:) = sst_avg_area(1,:)
 data_fuc(1,1,:) = sst_rc2*sst_avg_area&time+sst_rc2@yintercept
 data_fuc(1,2,:) =avg(sst_avg_area(1,:))
return(data_fuc)
end 
function draw_panel(wks,time,data1,data2,title,vpYF)
 begin
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@tmXBMode="Explicit"
    resxy@trXMaxF=432
    resxy@tmXBValues=ispan(13,432,60)      
    resxy@tmXBLabels=ispan(1980,2015,5)

    resxy@xyLineThicknesses = (/ 3.0,2.0,2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"black","red","black"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines","lines","Lines"/)
    resxy@xyExplicitLegendLabels = (/"Anual change", "5-point-smoothing", "Linear trend"/)   

    resxy@xyDashPatterns    =(/0,0,2/)


    resxy@tmXBLabelFontAspectF=2.0
    resxy@tmYLLabelFontAspectF=2.0
    resxy@tmXTOn=False
    resxy@tmYROn=False
    resxy@tiYAxisString = title

    
    resxy@vpWidthF=0.1
    resxy@vpHeightF=0.1
    resxy@vpKeepAspect=False
    ; resxy@vpXF=0.54
    ; resxy@vpYF=vpYF
  
   
    ; resxy@tmLabelAutoStride = True                    ; nice stride on labels
    contour=new(2,graphic)
   
    contour(0)=gsn_csm_xy(wks,time,data1,resxy)
    contour(1)=gsn_csm_xy(wks, time,data2,resxy)
    ; contour@res=resxy

    ; contour(0)=gsn_csm_xy(wks,time,data1,resxy)
    ; contour(1)=gsn_csm_xy(wks, time,data2,resxy)

return contour
end 

;******************************************************************************************************
;******************************************************************************************************
begin 


    diri="/mnt/f/share/era/data/"
    fils1=systemfunc("ls "+diri+"era_u_sst_slp.nc")
    f1=addfile(fils1, "r")
    fils2=systemfunc("ls "+diri+"era_rq_q.nc")
    f2=addfile(fils2, "r")
    fils3=systemfunc("ls "+diri+"era_div_200.nc")
    f3=addfile(fils3, "r")

    sst=f1->sst(0:431,:,:)
    msl=f1->msl(0:431,:,:)
    u10=f1->u10(0:431,:,:)  
    rq_r=f2->r(0:431,:,:)
      
    sst_time=f1->sst(0:431,40:80,:)
    msl_time=f1->msl(0:431,40:80,:)
    u10_time=f1->u10(0:431,40:80,:)

    sst_time&time=ispan(1,432,1)
    msl_time&time=ispan(1,432,1)
    u10_time&time=ispan(1,432,1)
    u10_time_new=-1*u10_time
    copy_VarMeta(u10_time, u10_time_new)
    printVarSummary(u10_time_new)

    q=f2->q(0:431,:,:)
    d=f3->d(0:431,:,:)
    ;------------------------------------------------去掉年循环

    sst_rc=decade_trend(sst,"℃","sst")
    msl_rc=decade_trend(msl,"Pa","mslp")
    u10_rc=decade_trend(u10,"ms**-1","u")
    ; q_rc=decade_trend(q,"%","q")
    ; rq_rc=decade_trend(rq,"%","rq")
    ; d_rc=decade_trend(d,"s**-1","d")
  
    wks=gsn_open_wks("png", "/mnt/f/share/era/era_rc_total")
   
    u10_rc_new=-u10_rc
    copy_VarMeta(u10_rc, u10_rc_new)

    u10title="*10~S~-3~N"
    plot=new(6,graphic) 
    plot(0)=draw_plot(u10_rc_new,wks,-0.005,0.005,0.001,0.9433)
    plot(2)=draw_plot(msl_rc,wks,-0.4,0.4,0.1,0.67)
    plot(4)=draw_plot(sst_rc,wks,-0.004,0.004,0.001,0.41)
    

    data_u10= gain_data(u10_time_new)
    data_msl= gain_data(msl_time)
    data_sst= gain_data(sst_time)

    pres1=True
    pres1@gsnMaximize=False

    pres2=True
    pres2@gsnMaximize=False

    contour_u10=draw_panel(wks,sst_time&time,data_u10(0,:,:),data_u10(1,:,:),"u10/m*s-1" ,0.706391)
    
    plot1= gsn_attach_plots(contour_u10(0),contour_u10(1),pres1,pres2)
    
    contour_msl=draw_panel(wks,sst_time&time,data_u10(0,:,:),data_u10(1,:,:),"slp/hpa" ,0.565919)
    
    plot2= gsn_attach_plots(contour_msl(0),contour_msl(1),pres1,pres2)
  
     contour_sst=draw_panel(wks,sst_time&time,data_sst(0,:,:),data_sst(1,:,:),"sst/k",0.425446 )
    
    plot3= gsn_attach_plots(contour_sst(0),contour_sst(1),pres1,pres2)
    
    draw(contour_u10(0))
    plot(1)=contour_u10(0)
    plot(3)=contour_msl(0)
    plot(5)=contour_sst(0)

     delete(pres1)
     delete(pres2)
frame(wks)
   


end