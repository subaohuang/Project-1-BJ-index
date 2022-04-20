;***************************************************************
  ;变化趋势平面图
;***************************************************************

function decade_trend(data,units,long_name)
    local data_a, data_a_new, rc
    begin

    if (data@long_name.eq."Sea surface temperature") then
    data_a=data@scale_factor*data+data@add_offset-274.15
    else 
    data_a=data@scale_factor*data+data@add_offset
    end if

                                           ;二维线性插值
    copy_VarMeta(data, data_a)
    data_a&time=ispan(1,456,1)
    data_a@units=units
    data_a@long_name=long_name
    printVarSummary(data_a)
    lat=fspan(-30, 30, 25)
    LAT=fspan(-30.0,30.0,121)
    LON=fspan(100.0,290.0,381)
    data_inter=linint2_Wrap(data_a&lon,lat, data_a(:,::-1,:), False, LON, LAT, 0)
   
    data_a_new = data_inter(LAT|:,LON|:,time|:)
    printVarSummary(data_a_new) 
    rc = regCoef(data_a_new&time,data_a_new)
    rc!0="lat"
    rc!1="lon"
    rc&lat=data_a_new&LAT
    rc&lon=data_a_new&LON
    rc&lat@units = "degrees_north"
    rc&lon@units = "degrees_east"
return(rc)
end
;*************************************************************************
 ;draw plot
;************************************************************************

function  draw_plot(rc,wks,cnMinLevelValF,cnMaxLevelValF,cnLevelSpacingF)
    local res, plot
    begin
    res = True
    res@gsnAddCyclic=False
    res@gsnDraw=False
    res@gsnFrame=False

    ; res@cnFillPalette  = "temp_diff_18lev"
    ; res@cnFillPalette="BlueDarkRed18"
    res@cnFillPalette="NCV_blu_red"
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


    res@pmLabelBarWidthF=0.1     ;colorbar的宽度


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


    res@mpDataSetName = "Earth..4"
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
;******************************************************************************************************
;******************************************************************************************************
begin 


    diri="/media/sf_share/ncep/"
    fils1=systemfunc("ls "+diri+"rhum.mon.mean.nc")
    fils2=systemfunc("ls "+diri+"uwnd.mon.mean.nc")

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
    rq=f1->rhum(0:455,9,24:48,40:116)
    printVarSummary(rq)
    rq_rc=decade_trend(rq,"%","rq")

    wks=gsn_open_wks("png", "/media/sf_share/ncep/necp_rq_rc_month_inter")
    ; u10_rc_new=-u10_rc
    ; copy_VarMeta(u10_rc, u10_rc_new)

    u10title="*10~S~-3~N"




    plot=draw_plot(rq_rc,wks,-0.02,0.02,0.005)
    draw(plot)
    frame(wks)

end