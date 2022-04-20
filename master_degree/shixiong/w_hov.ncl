begin
diri="/mnt/f/share/era/data/shixiong/"
fils=systemfunc("ls "+diri+"era_iod_pacific_w_u.nc")
f1=addfile(fils, "r")
time=f1->time
date=cd_calendar(time, 1)
tind=ind(date.le.201412.and.date.ge.197901)
w_r=f1->w(tind,{500},{10:-10},{100:290})
w=(w_r@scale_factor*w_r+w_r@add_offset)*(-1)
copy_VarCoords(w_r, w)
w_avg=dim_avg_n_Wrap(w, 1)
printVarSummary(w_avg)
w_avg&time=ispan(1, 432, 1)
diri3=diri+"locationxi.ascii"
diri4=diri+"locationdong.ascii"

locationxi=asciiread(diri3, 432, "float") 
locationdong=asciiread(diri4, 432, "float")
locationxi!0="time"
locationdong!0="time"
locationxi&time=ispan(1,432,1)
locationdong&time=ispan(1,432,1)
rc1=regCoef(locationxi&time,locationxi)
rc2=regCoef(locationdong&time,locationdong)
y1= rc1*locationxi&time+rc1@yintercept
y2= rc2*locationdong&time+rc2@yintercept

 wks  = gsn_open_wks ("pdf", ""+diri+"w_hov" )           ; send graphics to PNG file
 gsn_define_colormap(wks,"hotcolr_19lev")


 res            = True                         ; plot mods desired
 res@gsnDraw    = False                        ; don't draw
 res@gsnFrame   = False                        ; don't advance frame
 ; res@cnLineLabelsOn =False
 res@cnFillOn=True
 res@cnFillDrawOrder="PreDraw"
 res@cnLinesOn=False
 ; res@cnMonoLineColor  = False           ; Tells NCL not to draw contour lines in one color
 res@cnLevelSelectionMode = "ManualLevels"     ; manually select contours
 res@cnMinLevelValF       = -0.15               ; min level
 res@cnMaxLevelValF       =  0.15               ; max leve
 res@cnLevelSpacingF      =0.03
 res@lbOrientation         ="Vertical"
 res@gsnLeftString         ="(a)"
  res@tmYLMode      = "Explicit"                ; Define own tick mark labels.
  res@tmYLLabelsOn=True
  res@tmYUseLeft=True
  res@tmYLValues    = ispan(13,432,60)
  res@tmYLLabels    = ispan(1980,2010,5)
  res@tmXTOn=False
  res@tmYROn=False

  res@tmXBMode      = "Explicit"                ; Define own tick mark labels.
  res@tmXBValues    = ispan(100,290,40)
  res@tmXBLabels    = (/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W",\
  	"100~S~o~N~W"/)
  res@tmXBMinorOn=False
  res@vpHeightF=0.6
  res@vpWidthF=0.4
 
  ; res@lbBoxLinesOn=False

  plot = gsn_csm_hov(wks, w_avg, res)   ; create plot


  plres=True
  plres@gsLineThicknessF=5.0
  plres@gsLineDashPattern=0
  plres@tfPolyDrawOrder="PostDraw"
  plot_line1=gsn_add_polyline(wks, plot,  runave_n(locationxi,13,0,0), locationxi&time,plres)
  plres@gsLineColor="red"
  plot_line2=gsn_add_polyline(wks, plot, y1,locationxi&time, plres)
  plres@gsLineColor="black"
  plot_line3=gsn_add_polyline(wks, plot,  runave_n(locationdong,13,0,0),locationxi&time, plres)
  plres@gsLineColor="red"
  plot_line4=gsn_add_polyline(wks, plot,  y2,locationxi&time, plres)

   ; resxy=True
   ; resxy@gsnDraw=False
   ; resxy@gsnFrame=False
   ; resxy@tfPolyDrawOrder="PostDraw"
   ; plot_line1=gsn_csm_xy(wks, locationxi&time, runave_n(locationxi,13,0,0), resxy)

   ; overlay(plot, plot_line1)

  draw(plot)
  frame(wks)

  end