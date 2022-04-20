begin
diri="/mnt/f/share/era/data/"
fils1=systemfunc("ls "+diri+"msf_time_era_global.nc")
f1=addfile(fils1, "r")
fils2=systemfunc("ls /mnt/f/share/ncep/20CR/omega.mon.mean.nc") 
f2=addfile(fils2, "r")
date=cd_calendar(f2->time, 1)
timeind=ind(date.le.201012.and.date.ge.195001)
w=f2->omega(timeind,:,{30:-30},:)
delete(w&time)
w&time=date(timeind)
w=(/w*-100/)

w_avg=dim_avg_n_Wrap(w(:,{850:200},{10:-10},:), (/1,2/))
 ; w_avg=dim_avg_n_Wrap(w(:,{500},{10:-10},:), (/1/))
w_yr=month_to_annual(w_avg, 1)


w_clm=dim_avg_n_Wrap(w_yr, 0)

; do i=0,60
;   w_yr(i,:)=(/w_yr(i,:)-w_clm/)   ;;;求距平场
; end do

; w_yr_sm=runave_n_Wrap(w_yr, 9, 1, 0)
w_yr_sm=w_yr

;************************************************
; 读入太平洋和印度洋的气候态位置
;************************************************
 paindex=asciiread("/mnt/f/share/era/data/pa_index_1950_2010.ascii", 732, "float")
 idindex=asciiread("/mnt/f/share/era/data/id_index_1950_2010.ascii", 732, "float")
 
 idindex!0="time"
 idindex&time=ispan(1,732,1)
 copy_VarMeta(idindex, paindex)

 paindex_year=month_to_annual(paindex, 1)
 idindex_year=month_to_annual(idindex, 1)

 idindex_year&year=ispan(1,61,1)
 copy_VarMeta(idindex_year, paindex_year)


;;;年变化
 rc1=regCoef(paindex_year&year,paindex_year)
 rc2=regCoef(idindex_year&year,idindex_year)

 paindex_rc =rc1*paindex_year&year+rc1@yintercept
 idindex_rc =rc2*idindex_year&year+rc2@yintercept
;;;

;;;月变化
 rc1_month=regCoef(paindex&time,paindex)
 rc2_month=regCoef(idindex&time,idindex)


 paindex_rc_month =rc1_month*paindex&time+rc1_month@yintercept
 idindex_rc_month =rc2_month*idindex&time+rc2_month@yintercept

;;;



 pa_clm=dim_avg_n(paindex, 0)
 id_clm=dim_avg_n(idindex, 0)

  paindex_year_sm=paindex_year
  idindex_year_sm=idindex_year 
 ; paindex_year_sm=runave_n_Wrap(paindex_year, 3, 1, 0)
 ; idindex_year_sm=runave_n_Wrap(idindex_year, 3, 1, 0)


;************************************************
; 画时间-经度图
;************************************************

 wks  = gsn_open_wks ("pdf", ""+diri+"20cr的w的月数据时间-经度图" )           ; send graphics to PNG file
 ; gsn_define_colormap(wks,"BlAqGrWh2YeOrReVi22")

 gsn_define_colormap(wks, "BlAqGrYeOrReVi200")
 

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
  res@cnMinLevelValF       = -7.0               ; min level
  res@cnMaxLevelValF       =  7.0              ; max leve
  res@cnLevelSpacingF      =0.5

  ;    res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
  ; res@cnMinLevelValF       = -2.0               ; min level
  ; res@cnMaxLevelValF       =  2.0              ; max leve
  ; res@cnLevelSpacingF      =0.25


  res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
  res@tmYLLabelsOn=True
  res@tmYUseLeft=True

  ; res@tmYLValues    = ispan(0,60,5)
  ; res@tmYLLabels    = ispan(1950,2010,5)

  res@tmYLValues    = ispan(0,731,60)
  res@tmYLLabels    = ispan(1950,2010,5)


  res@tmXTOn=False
  res@tmYROn=False

  res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
  res@tmXBValues=(/0,90,180,270,357.5/)
  res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
  res@tmXBMinorOn=True
  res@tmXBMinorValues=ispan(0,360,10)



  res@vpHeightF=0.5
  res@vpWidthF=0.8



  plot = gsn_csm_hov(wks, w_avg, res)   ; create plot

   ;;;添加0线
	res@cnLevelSelectionMode = "ExplicitLevels"
	res@cnLinesOn=True
	res@cnLevels = 0
	res@cnLineLabelsOn=False
	res@cnLineThicknessF = 2.
	res@cnFillOn        = False
   plot1 = gsn_csm_hov(wks, w_avg, res)   ; create plot

  

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

;;;添加气候态的位置
    lineres=True
    lineres@gsLineColor="green"
    lineres@gsLineDashPattern=16
    lineres@gsLineThicknessF=2.0
     time=ispan(0,731,5)     
     linelon1=conform_dims(dimsizes(time), pa_clm, 0)
     linelon2=conform_dims(dimsizes(time), id_clm, 0)
    
	dum4= gsn_add_polyline(wks,plot,linelon1,time,  lineres)
	dum5= gsn_add_polyline(wks,plot,linelon2,time,  lineres)

 ;;;添加位置变化
    line=True
    line@gsLineColor="black"
    line@gsLineDashPattern=0
    line@gsLineThicknessF=2.0

    dum6= gsn_add_polyline(wks,plot,paindex, time,  line)
	  dum7= gsn_add_polyline(wks,plot,idindex, time,  line)
    
    dum8= gsn_add_polyline(wks,plot,paindex_rc_month, time,  line)
  	dum9= gsn_add_polyline(wks,plot,idindex_rc_month, time,  line)

 ; overlay(plot, plot1)
  draw(plot)
  frame(wks)
end 