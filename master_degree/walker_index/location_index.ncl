function gain_data_an(data1)
 begin
 
 intesion_an=runave_n_Wrap(dim_rmvmean(month_to_annual(data1(0:455), 1)),3,0,0)
 intesion_an!0="time"
 intesion_an&time=ispan(1979,2016,1)


 rc1=regCoef(intesion_an&time,intesion_an)

 data_fuc=new((/3,dimsizes(intesion_an&time)/),double)


 data_fuc(0,:) = intesion_an(:)
 data_fuc(1,:) = rc1*intesion_an&time+rc1@yintercept
 data_fuc(2,:) =avg(intesion_an)

return(data_fuc)
end 

function gain_data_de(data)
 begin
 data1=dtrend_n(data, True, 0)
 intesion_an=dim_rmvmean(month_to_annual(data1(0:455), 1))
 intesion_an!0="time"
 intesion_an&time=ispan(1979,2016,1)


 ; rc1=regCoef(intesion_an&time,intesion_an)

 data_fuc=new((/2,dimsizes(intesion_an&time)/),double)


 data_fuc(0,:) = intesion_an(:)
 ; data_fuc(1,:) = rc1*intesion_an&time+rc1@yintercept
 data_fuc(1,:) =avg(intesion_an)

return(data_fuc)
end 


begin
 diri="/mnt/f/share/era/data/"
 diri_pi="/mnt/f/share/era/"
 fils0=systemfunc("ls "+diri+"msf_era_level_location1979-2017.nc")
 f0=addfile(fils0, "r")
 fils1=systemfunc("ls "+diri+"msf_era_gravity_center_location1979-2017.nc")
 f1=addfile(fils1, "r")
 x=f1->x
 locationxi_a=f0->locationxi_a
 locationdong_a=f0->locationdong_a
 printVarSummary(locationxi_a)
 data=gain_data_de(locationxi_a)
 data1=gain_data_de(x)
 time=ispan(1979,2016,1)

;********************************************************************************************
    resxy=True
    resxy@gsnDraw=False
    resxy@gsnFrame=False
    resxy@gsnLeftString="(a)msf_zeroline"
    ; resxy@gsnRightString=""+slp_data(3,1)+""
    resxy@gsnLeftStringFontHeightF=0.01
    resxy@gsnRightStringFontHeightF=0.01
    ; resxy@trXMaxF=dimsizes(intesion&time)
    ; resxy@tmXBValues=ispan(13,dimsizes(intesion&time),60)
    ; resxy@tmXBLabels=ispan(1980,2015,5)

    resxy@tmXBMode="Explicit"
    resxy@trXMaxF=2016
    resxy@trXMinF=1979
    resxy@tmXBValues=ispan(1979,2016,5)
    resxy@tmXBLabels=ispan(1979,2016,5)
    resxy@tmXBMinorOn=True
    resxy@tmXBMinorValues=ispan(1979,2016,1)

    resxy@tmYLMinorOn=False

    resxy@vpHeightF=0.13
    resxy@vpWidthF=0.7
    resxy@vpXF=0.15
    resxy@vpYF=0.95

    resxy@xyLineThicknesses = (/ 3.0,1.5,2.0/)          ; make second line thicker
    resxy@xyLineColors      = (/"black","red","black"/)          ; change line color
    resxy@xyMarkLineModes   = (/"lines","lines","Lines"/)
    ; resxy@xyExplicitLegendLabels = (/"Anual change",  "Linear trend"/)     ;加线的标注

    resxy@xyDashPatterns    =(/0,0,2/)
    
    
    resxy@tmXBLabelFontAspectF=1.0
    resxy@tmXBLabelFontHeightF=0.008
    resxy@tmXBLabelFont=30     ;黑体
    resxy@tmXBLabelFontThicknessF=2.0
    resxy@tmXBMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmXBMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
    ; resxy@tmYLLabelFontAspectF=2.0      ;字体的高除以长的比率

    resxy@tmXBMinorLengthF=0.003
    resxy@tmXBMinorOutwardLengthF=0.003

    resxy@tmYLMajorLengthF=0.005        ;調整刻度的长度
    resxy@tmYLMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
    resxy@tmYLLabelFont=30
    resxy@tmYLLabelFontHeightF=0.008
    resxy@tmXTOn=False
    resxy@tmYROn=False





;***************************************************************************************************************8
 wks=gsn_open_wks("png", ""+diri+"zeroline_an")

 wks1=gsn_open_wks("png", ""+diri+"gravity_an")
 wks2=gsn_open_wks("png", ""+diri_pi+"location_index_panel")

 plot=gsn_csm_xy(wks,time, data, resxy)
 resxy@gsnLeftString="(b)gravity_center"
 plot1=gsn_csm_xy(wks1,time,data1 , resxy)

 plot_panel=new(2,graphic)
 resxy@gsnLeftString="(a)msf_zeroline"
 plot_panel(0)=gsn_csm_xy(wks2,time, data, resxy)
 resxy@gsnLeftString="(b)gravity_center"
 plot_panel(1)=gsn_csm_xy(wks2,time,data1 , resxy)


 
 enino=(/1982,1983,1986,1987,1991,1992,1993,1994,1997,1998/)
 
    ilon=new(dimsizes(enino),integer)
    do i=0,dimsizes(enino)-1
    ilon(i)=ind(data&time.eq.enino(i))
  end do
    itime=data&time(ilon)
   print(itime)
    ilev=new((/2,dimsizes(enino)/),double)
    locationxi_an=data(0,:)
    x_an=data1(0,:)
    ilev(0,:)=locationxi_an(ilon)
    ilev(1,:)=x_an(ilon)
    
    
   polyres                = True
   polyres@gsMarkerIndex     = 16          ; polymarker style
   polyres@gsMarkerSizeF     = 5.         ; polymarker size
   polyres@gsMarkerColor     ="red"
   polyres@tfPolyDrawOrder="PreDraw"
   dum1=new(2,graphic)
   do i=0,1
    dum1(i)=gsn_add_polymarker(wks, plot_panel(i), itime, ilev(i,:), polyres)
   end do 

    
 pres=True
 pres@gsnMaximize=True
 gsn_panel(wks2, plot_panel, (/2,1/), pres)

end