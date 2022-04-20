;***********************************************************
;分区域，东西区域进行平均后进行5点平滑
;*************************************************************


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

 data_a_new&time=ispan(1,456,1)
 
 

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
  data_avg_new_a=runave_n_Wrap(data_avg_new, 11, 0, 1)


return(data_avg_new_a)
end 

;---------------------------------------------------
;获取拟合曲线，平均值直线和年际变化曲线
;------------------------------------------------------
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


;***********************************************************
;画几个图，为了能在gsn_attch上使用
;***********************************************************
function draw_panel(wks,time,data1,data2,title)
 begin
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@tmXBMode="Explicit"
    resxy@trXMaxF=456
    resxy@tmXBValues=ispan(13,456,60)
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

    ; resxy@tmLabelAutoStride = True                    ; nice stride on labels
    contour=new(2,graphic)
   
    contour(0)=gsn_csm_xy(wks,time,data1,resxy)
    contour(1)=gsn_csm_xy(wks, time,data2,resxy)


    ; contour(0)=gsn_csm_xy(wks,time,data1,resxy)
    ; contour(1)=gsn_csm_xy(wks, time,data2,resxy)

return contour
end 
;*************************************************************************
;主程序
;*************************************************************************
begin 

    diri="/media/sf_share/era/"     ;读取数据
    fils1=systemfunc("ls "+diri+"era_u_sst_slp.nc")
    f1=addfile(fils1, "r")
    fils2=systemfunc("ls "+diri+"era_w_rq.nc")
    f2=addfile(fils2, "r")
 


    sst=f1->sst(0:455,40:80,:)
    msl=f1->msl(0:455,40:80,:)
    u10=f1->u10(0:455,40:80,:)
    rq=f2->r(0:455,0,40:80,:)

    sst&time=ispan(1,456,1)

    data= gain_data(rq)
    wks2=gsn_open_wks("png", "/media/sf_share/era/era_rq_series")
    pres1=True
    pres2=True
    pres1@gsnMaximize=True

    contour=draw_panel(wks2,sst&time,data(0,:,:),data(1,:,:),"rq/%" )
    plot = gsn_attach_plots(contour(0),contour(1),pres1,pres2)
    draw(contour(0))
    frame(wks2)

end
