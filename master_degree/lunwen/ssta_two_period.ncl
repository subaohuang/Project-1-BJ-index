begin
  diri="/mnt/f/share/ersst/"
  fils1=systemfunc("ls "+diri+"sst.mnmean.v4.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201312.and.date.ge.197901)
  sst=f1->sst(timeind,{30:-30},:)
  ; sst=f1->sst(timeind,:,:)
  
  sst!0="time"
  sst!1="lat"
  sst!2="lon"
  delete(date)
  delete(f1)

;;;处理成异常的冬季平均场
  ; ssta=(/rmMonAnnCycTLL(sst(:,:,{30:300}))/)
  ; copy_VarCoords(sst(:,:,{30:300}), ssta)
  
  ssta=(/rmMonAnnCycTLL(sst(:,:,:))/)
  copy_VarCoords(sst(:,:,:), ssta)

  ssta_djf=month_to_season(ssta, "DJF")
  ssta_son=month_to_season(ssta, "SON")

  ssta_s=(/ssta_djf(1:,:,:)/)*1.0
    ; ssta_s=(/ssta_son(0:33,:,:)/)
    copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&time)-2,:,:), ssta_s)
   delete(ssta_s&time)
   ssta_s&time=ispan(1979, 2012, 1)  
 ;;;两个时间段的海温异常
 ssta_s1=dim_avg_n_Wrap(ssta_s({1979:1994},:,:),0)
 ssta_s2=dim_avg_n_Wrap(ssta_s({1995:},:,:),0)
 ;;;

;;;处理成冬季平均场
  sst_djf=month_to_season(sst, "DJF")
  sst_son=month_to_season(sst, "SON")

  sst_s=(/sst_djf(1:,:,:)/)*1.0
    ; sst_s=(/sst_son(0:33,:,:)/)
  copy_VarCoords(sst_djf(0:dimsizes(sst_djf&time)-2,:,:), sst_s)
   delete(sst_s&time)
   sst_s&time=ispan(1979, 2012, 1)  
 ;;;两个时间段的海温异常
 sst_s1=dim_avg_n_Wrap(sst_s({1979:1994},:,:),0)
 sst_s2=dim_avg_n_Wrap(sst_s({1995:},:,:),0)
 sst_s3=sst_s2 - sst_s1
 copy_VarCoords(sst_s1, sst_s3)

 ;;;


wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/ 2个时间段的海温异常")
wks2=gsn_open_wks("pdf", "/mnt/f/share/era/data/ 2个时间段的海温1")

plot=new(2,graphic)
plot1=new(3,graphic)


; colormap=read_colormap_file("GMT_jet")
; printVarSummary(colormap)
; cmap=new((/18,4/),typeof(colormap))
; printVarSummary(colormap)

; cmap(0:6,:)=colormap(2:92:15,:)
; ; cmap(0:6,:)=colormap(2:50:8,:)
; cmap(7:17,:)=colormap(105:255:15,:)

  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  ; res@cnFillPalette="cmocean_balance"
  res@cnFillPalette="MPL_Spectral"
  res@gsnSpreadColorStart=5
  res@gsnSpreadColorEnd=115
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLineLabelsOn=False
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  ; res@cnMinLevelValF       =  14
  ; res@cnMaxLevelValF       =  30
  ; res@cnLevelSpacingF      =  1
  res@lbLabelBarOn=False

  res@lbOrientation = "Horizontal"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(ssta&lat)
  res@mpMaxLatF             = max(ssta&lat)
  res@mpMinLonF             = min(ssta&lon)
  res@mpMaxLonF             = max(ssta&lon)
  res@mpCenterLonF          =avg(ssta&lon)
  res@mpShapeMode   ="FreeAspect"
  res@vpHeightF=0.35
  res@vpWidthF=0.65
  
  res@tmXTOn=False
  res@tmYROn=False
  
  plot(0)=gsn_csm_contour_map(wks1,ssta_s1,res)
  plot(1) =gsn_csm_contour_map(wks1,ssta_s2,res)
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =  14
  res@cnMaxLevelValF       =  30
  res@cnLevelSpacingF      = 0.5
  res@gsnLeftString        =  "1979-1994 mean"
  plot1(0)=gsn_csm_contour_map(wks2,sst_s1,res)
   res@gsnLeftString        =  "1995-2012 mean"
  plot1(1) =gsn_csm_contour_map(wks2,sst_s2,res)
  delete( res@gsnSpreadColorStart)
  delete(res@gsnSpreadColorEnd)
    res@cnFillPalette="matlab_jet"
   res@cnMinLevelValF       =  -1
  res@cnMaxLevelValF       =  1
  res@cnLevelSpacingF      =  0.1
    res@lbLabelBarOn=True
  res@vpXF=0.217047
  res@vpYF=0.32
  res@vpWidthF=0.488313
  res@vpHeightF=0.262938
  res@lbOrientation="Vertical"
  res@pmLabelBarWidthF=0.15
  res@gsnLeftString        =  "(b)-(a)"
  plot1(2) =gsn_csm_contour_map(wks2,sst_s3,res)
 
 ;---Create text strings
    txres                       = True
    txres@txPerimOn             = True
    txres@txFontHeightF         = 0.013

    amres_tr                  = True
    amres_tr@amParallelPosF   =  0.46    ; This is the right edge of the plot.
    amres_tr@amOrthogonalPosF =  0.43    ; This is the top edge of the plot.
     tr_label = (/"c"/)
     txres@txBackgroundFillColor = "White"
     txres@txPerimOn = True
     txid_tr = gsn_create_text(wks2, tr_label, txres)
 ;---Attach text strings to plot
     amid_tr = gsn_add_annotation(plot1(2), txid_tr, amres_tr)

   delete(res@gsnLeftString)
  res@cnLevelSelectionMode="ExplicitLevels"
  res@cnLevels=28.75
  res@cnFillOn=False
  res@cnLineThicknessF =3.0
  res@cnLinesOn=True
  res@cnInfoLabelOn=False
  
  plot3=gsn_csm_contour(wks2,sst_s1,res)
  plot4 =gsn_csm_contour(wks2,sst_s2,res)
  overlay(plot1(0), plot3)
  overlay(plot1(1), plot4)

  pres=True
  pres@gsnFrame=False
  pres@gsnPanelLabelBar=True
  pres@lbOrientation="Vertical"
  pres@pmLabelBarWidthF=0.15
  pres@pmLabelBarHeightF=0.5
  pres@gsnPanelDebug=True
  pres@gsnPanelBottom=0.36
  pres@pmLabelBarOrthogonalPosF=-0.01
  pres@gsnPanelFigureStrings=(/"a","b"/)


  gsn_panel(wks1, plot, (/2,1/), pres)
  gsn_panel(wks2, (/plot1(0),plot1(1)/), (/2,1/), pres)
  

  draw(plot1(2))
  frame(wks2) 
 end