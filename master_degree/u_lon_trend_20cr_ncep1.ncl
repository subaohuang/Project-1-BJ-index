
function decade_trend(data)
    local data_a, data_a_new, rc
    begin
    time=ispan(1,dimsizes(data&time),1)
    rc = regCoef_n(time,data,0,0)
    rc!0="lon"
    rc&lon=data&lon
    rc&lon@units = "degrees_east"
    ; print(rc)
return(rc)
end

begin
  diri="/mnt/f/share/"
  fils1=systemfunc("ls "+diri+"ncep/20CR/uwnd.sig995.mon.mean.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.195801.and.date.le.201212)

  u_20cr=f1->uwnd(timeind,{-5:5},:)
  
  ; delete(msf_time_20cr&time)
  ; msf_mean_20cr&time=date(timeind)
  u_20cr_mean=dim_avg_n_Wrap(u_20cr, 1)
	delete(fils1)
	delete(f1)
	delete(date)
	delete(timeind)
;;;



;;;
  fils1=systemfunc("ls "+diri+"ncep/ncep1/uwnd.sur.mon.mean.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.195801.and.date.le.201212)


  u_ncep1=f1->uwnd(timeind,{5:-5},:)
  u_ncep1_mean=dim_avg_n_Wrap(u_ncep1, 1)
   
   ; delete(msf_time_ncep1&time)
   ; msf_mean_ncep1&time=date(timeind)

	delete(fils1)
	delete(f1)
	delete(date)
	delete(timeind)
;;;


;;;;

    rc_20cr=decade_trend(u_20cr_mean)
    rc_ncep1=decade_trend(u_ncep1_mean)

    rc_20cr=(/rc_20cr*120/)
    rc_ncep1=(/rc_ncep1*120/)
;;;


;;;气候态
    u_20cr_clm=dim_avg_n_Wrap(u_20cr_mean, 0)
    u_ncep1_clm=dim_avg_n_Wrap(u_ncep1_mean, 0)
;;;


;;;画图
wks=gsn_open_wks("pdf", "/mnt/f/share/era/data/20CR&necp1_u随经度的变化和趋势变化")


; resources for "left" variable
  resL                        = True 
  resL@gsnDraw=False
  resL@gsnFrame=False
  resL@gsnRightString="20CR"
  resL@xyLineThicknesses      = 2.                  ; thicker line
  resL@xyLineColors      = "red"                   ; line color
  resL@gsnYRefLine=0
  resL@tiYAxisString = "u trend"  ; axis string
   
  resL@tmYLMode="Explicit"
  resL@trYMaxF=0.5
  resL@trYMinF=-0.5
  resL@tmYLValues=fspan(-0.5, 0.5, 11)
  resL@tmYLLabels=fspan(-0.5, 0.5, 11)

	resL@tmXTOn=False
	resL@tmXBMode="Explicit"
	resL@trXMaxF=360.0
	resL@tmXBValues=ispan(0,360,30)

	resL@tmXBLabelFontHeightF=0.01
	resL@tmXBLabelFontThicknessF=0.02
	resL@tmXBMinorOn=True
	resL@tmXBLabels=(/"0","30~S~o~N~E","60~S~o~N~E","90~S~o~N~E","120~S~o~N~E","150~S~o~N~E",\
	        "180","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W","60~S~o~N~W","30~S~o~N~W","0"/)

    resL@vpHeightF=0.3
    resL@vpWidthF=0.6
     
; resources for "right" variable
	resR                        = True   
	resR@gsnDraw=False
	resR@gsnFrame=False                   
	resR@xyDashPatterns         = 1                   ; dashed line for 2nd
	resR@xyLineThicknesses      = 2                   ; thicker line
	resR@xyLineColors      = "blue"                   ; line color
	resR@tiYAxisString = "Mean State"  ; axis string
  
  plot=new(2,graphic)
  plot(0)= gsn_csm_xy2(wks,u_20cr&lon,rc_20cr,u_20cr_clm,resL,resR)
  resL@gsnRightString="ncep1"
  plot(1)= gsn_csm_xy2(wks,u_ncep1&lon,rc_ncep1,u_ncep1_clm,resL,resR)
  
;;;添加地形
   ployres                  = True          ; polyline mods desired
      ployres@gsLineThicknessF = 15.0            ; line thickness   
      ployres@gsLineColor      = "black"         ; color of lines
      ployres@gsLineDashPattern=15
 ;---Attach the polyline
      ploylevel=-0.5
      ploylon1=fspan(9, 39, 25)
      ploylon2=fspan(100, 149, 25)
      ploylon3=fspan(280, 315, 25)
      ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
      ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
      ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)
      printVarSummary(ploylevel1)
  
     dum1=new(2,graphic)
     dum2=new(2,graphic) 
     dum3=new(2,graphic) 

    ployres@gsLineDashPattern=0
    do i=0,1
    dum1(i)= gsn_add_polyline(wks,plot(i),  ploylon1,ploylevel1,  ployres)
    end do
    ployres@gsLineDashPattern=16
    do i=0,1
    dum2(i)= gsn_add_polyline(wks,plot(i),  ploylon2,ploylevel2,  ployres)
    end do
    ployres@gsLineDashPattern=0
    do i=0,1
    dum3(i)= gsn_add_polyline(wks,plot(i),  ploylon3,ploylevel3,  ployres)
    end do

  pres=True
  gsn_panel(wks, plot, (/2,1/), pres)


  end


    


