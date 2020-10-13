;***************************************************************
  ;变化趋势平面图
;***************************************************************

function decade_trend(data,units,long_name)
    local data_a, data_a_new, rc
    begin

    if (data@long_name.eq."Monthly mean v wind") then
    data_a=data
    else 
    data_a=data@scale_factor*data+data@add_offset
    end if

    copy_VarMeta(data, data_a)
    data_a&time=ispan(1,dimsizes(data&time),1)
    data_a@units=units
    data_a@long_name=long_name
    data_a_new=data_a(lat|:,lon|:,time|:)
    printVarSummary(data_a_new)
 ;-----------------------------------------------去掉年循环
    data_ano=rmMonAnnCycLLT(data_a_new)
    
;-----------------------------------------------
  
;********************************************************************
   ; ; 插值 
   ;  lat=fspan(-30, 30, 25)
   ;  LAT=fspan(-30.0,30.0,121)
   ;  LON=fspan(100.0,290.0,381)
   ;  ; data_inter=linint2_Wrap(data_a&lon,lat, data_a(:,::-1,:), False, LON, LAT, 0)    
   ;  data_inter=linint2_Wrap(data_a&lon,lat, data_a(:,:,:), False, LON, LAT, 0)        ;注意速度势这里的纬度是从小到大
   ;  data_a_new = data_inter(LAT|:,LON|:,time|:)
;********************************************************************************


    rc = regCoef(data_ano&time,data_ano)
    rc!0="lat"
    rc!1="lon" 
    rc&lat=data_ano&lat
    rc&lon=data_ano&lon
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
    ; rc&lat=data_a_new&LAT
    ; rc&lon=data_a_new&LON

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
    ; res@cnFillPalette="NCV_blu_red"
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

    txid_tr = new(1 , graphic)
    amid_tr = new(1, graphic)
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.02
    txres@txFontAspectF         =1.5
    
    amres_tr                  = True
  ; amres_tr@amJust      =TopRight
    amres_tr@amParallelPosF   =  -0.48    ; This is the right edge of the plot.
    amres_tr@amOrthogonalPosF =  -0.38    ; This is the top edge of the plot.

    letter = (/"i"/)

;---Create text strings
    tr_label = letter 

    txres@txBackgroundFillColor = "white"
    txres@txPerimOn = True
    txid_tr= gsn_create_text(wks, tr_label, txres)

;---Attach text strings to plot
    amid_tr= gsn_add_annotation(plot, txid_tr, amres_tr)

return(plot)
end 
;******************************************************************************************************
;******************************************************************************************************
begin 


    diri="/mnt/f/share/ncep/ncep1/"
    fils1=systemfunc("ls "+diri+"rhum.mon.mean.nc")
     fils2=systemfunc("ls "+diri+"uwnd.mon.mean.nc")
     fils3=systemfunc("ls "+diri+"vwnd.mon.mean.nc")
    fils4=systemfunc("ls "+diri+"olr.mon.mean.nc")
    f1=addfile(fils1, "r")
    ; rq=f1->rhum(371:827,9,24:48,40:116)
    f2=addfile(fils2, "r")
    f3=addfile(fils3, "r")
    f4=addfile(fils4, "r")
    timeolr=f4->time
    date=cd_calendar(timeolr, 1)
    ti1=ind(date.eq.197901)
    ti2=ind(date.eq.201301)
    olr=f4->olr(ti1:ti2-1,24:48,40:116)
   printVarSummary(olr)
   ;------------------------------------
    ; uwnd=f2->uwnd(371:827,9,:,:)
    ; vwnd=f3->vwnd(371:827,9,:,:)
    ; uwnd_new=uwnd(:,::-1,:)
    ; copy_VarMeta(uwnd, uwnd_new)
    ; uwnd_new&lat=fspan(-90,90,73)

    ; vwnd_new=vwnd(:,::-1,:)
    ; copy_VarMeta(vwnd, vwnd_new)
    ; vwnd_new&lat=fspan(-90,90,73)

    ; result=uv2sfvpF(uwnd_new, vwnd_new)
    ; printVarSummary(result)
    ; vp=result(1,:,:,:)
    ; copy_VarMeta(vwnd_new, vp)
    ; vp_new=vp(:,24:48,40:116)
    ;---------------------------------------
    ; rq_rc=decade_trend(rq,"%","rq")
    olr_rc=decade_trend(olr,"%","olr")
    ; vp_rc=decade_trend(vp_new,"%","vp")
    wks=gsn_open_wks("png", "/mnt/f/share/era/olr")
    ; u10_rc_new=-u10_rc
    ; copy_VarMeta(u10_rc, u10_rc_new)

    u10title="*10~S~-3~N"
    
    plot=draw_plot(olr_rc,wks,-0.02,0.02,0.005)
    draw(plot)
    frame(wks)

end