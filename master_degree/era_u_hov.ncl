begin

diri="/mnt/f/share/era/data/"
fils1=systemfunc("ls "+diri+"era_pre_global.nc")
f1=addfile(fils1, "r")
date=cd_calendar(f1->time, 1)
timeind=ind(date.le.201012.and.date.ge.197901)
u_r=f1->u(timeind,{200},{10:-10},:)

u=u_r*u_r@scale_factor+u_r@add_offset
copy_VarCoords(u_r, u)
u_annual=month_to_annual(u, 1)
printVarSummary(u_annual)

u_annual_avg=dim_avg_n_Wrap(u_annual,1)

u_annual_avg=(/runave_n(u_annual_avg, 3, 0, 0)/)

wks  = gsn_open_wks ("pdf", ""+diri+"era u的时间-经度图" )           ; send graphics to PNG file
 gsn_define_colormap(wks,"BlAqGrWh2YeOrReVi22")
 ; gsn_define_colormap(wks, "BlRe")
  ; gsn_define_colormap(wks, "GMT_paired")

	res            = True                         ; plot mods desired
	res@gsnDraw    = False                        ; don't draw
	res@gsnFrame   = False                        ; don't advance frame
    res@gsnLeftString=""
    res@gsnRightString=""

    res@lbLabelBarOn=True
    res@lbOrientation="Vertical"

	res@cnFillOn=True
	res@cnFillDrawOrder="PreDraw"
	res@cnLinesOn=False
	res@cnFillMode="RasterFill"
	res@cnRasterSmoothingOn=True
    res@cnInfoLabelOn=False
    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
	res@cnMinLevelValF       = -12.0               ; min level
	res@cnMaxLevelValF       =  12.0              ; max leve
	res@cnLevelSpacingF      =1.0
	 
 ;    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
	; res@cnMinLevelValF       = -2.0               ; min level
	; res@cnMaxLevelValF       =  2.0              ; max leve
	; res@cnLevelSpacingF      =0.25

 
	res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
	res@tmYLLabelsOn=True
	res@tmYUseLeft=True

	res@tmYLValues    = ispan(0,31,2)
	res@tmYLLabels    = ispan(1979,2010,2)

	res@tmXTOn=False
	res@tmYROn=False

	res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
	res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)
    
     plot = gsn_csm_hov(wks, u_annual_avg, res)   ; create plot

 ;;;添加地形
    ployres                  = True          ; polyline mods desired
    ployres@gsLineThicknessF = 5.0            ; line thickness   
    ployres@gsLineColor      = "black"         ; color of lines
    ployres@gsLineDashPattern=15
   ;---Attach the polyline
    ploylevel=0.0
    ploylon1=fspan(9, 39, 25)
    ploylon2=fspan(100, 149, 25)
    ploylon3=fspan(280, 315, 25)
    ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
    ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
    ploylevel3=conform_dims(dimsizes(ploylon3), ploylevel, 0)  
      
	dum1= gsn_add_polyline(wks,plot,ploylon1,ploylevel1,  ployres)
	ployres@gsLineDashPattern=16
	dum2= gsn_add_polyline(wks,plot,ploylon2,ploylevel2,  ployres)
	ployres@gsLineDashPattern=0
	dum3= gsn_add_polyline(wks,plot,ploylon3,ploylevel3,  ployres)

    draw(plot)
    frame(wks)
end