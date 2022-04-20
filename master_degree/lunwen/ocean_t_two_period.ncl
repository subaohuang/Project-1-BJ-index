begin

  ;;;读入MERRA的数据
  diri="/mnt/f/share/ORAS4/"
  fils1=systemfunc("ls "+diri+"ORAS4_1979-2012.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->TIME, 1)
  timeind=ind(date.le.201212.and.date.ge.197901)
  sst=f1->THETAO(timeind,{5:400},{-5:5},:)
  printVarSummary(sst)

  sst_m=dim_avg_n_Wrap(sst, 2)


  ssta_djf=month_to_season(sst_m, "DJF")
  ssta_s=(/ssta_djf(:,:,:)/)
  copy_VarCoords(ssta_djf(0:dimsizes(ssta_djf&TIME)-1,:,:), ssta_s)
  printVarSummary(ssta_s)

  ; delete(ssta_s&TIME)

  ssta_s&TIME=ispan(1979, 2012, 1)  
  ;;;两个时间段的海温异常
  ssta_s1=dim_avg_n_Wrap(ssta_s({1979:1994},:,:),0)
  ssta_s2=dim_avg_n_Wrap(ssta_s({1995:},:,:),0)
  ;;;

  ssta_s3=ssta_s2 - ssta_s1
  copy_VarCoords(ssta_s1, ssta_s3)

 ;;;


  wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/ 2个时间段的下层海温")

  plot1=new(3,graphic)



  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  res@cnFillPalette="MPL_Spectral"
  res@gsnSpreadColorStart=5
  res@gsnSpreadColorEnd=115
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn =False
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation = "Horizontal"
  res@gsnYAxisIrregular2Linear=True
  res@trYReverse=True
  res@cnInfoLabelOn=False

  res@vpHeightF=0.35
  res@vpWidthF=0.65
  
  res@tmXTOn=False
  res@tmYROn=False
  
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =  10
  res@cnMaxLevelValF       =  30
  res@cnLevelSpacingF      = 1.0
  res@gsnLeftString        =  "1979-1994 mean"
  plot1(0)=gsn_csm_contour(wks1,ssta_s1,res)
   res@gsnLeftString        =  "1995-2012 mean"
  plot1(1) =gsn_csm_contour(wks1,ssta_s2,res)
  delete( res@gsnSpreadColorStart)
  delete(res@gsnSpreadColorEnd)
  res@cnFillPalette="matlab_jet"
  res@cnMinLevelValF       =  -1
  res@cnMaxLevelValF       =  1
  res@cnLevelSpacingF      =  0.1
  res@lbLabelBarOn=True
  res@vpXF=0.233191
  res@vpYF=0.32
  res@vpWidthF=0.491055
  res@vpHeightF=0.264414
  res@lbOrientation="Vertical"
  res@pmLabelBarWidthF=0.15
  res@gsnLeftString        =  "(b)-(a)"
  plot1(2) =gsn_csm_contour(wks1,ssta_s3,res)
 
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
     txid_tr = gsn_create_text(wks1, tr_label, txres)
 ;---Attach text strings to plot
     amid_tr = gsn_add_annotation(plot1(2), txid_tr, amres_tr)

   delete(res@gsnLeftString)
  res@cnLevelSelectionMode="ExplicitLevels"
  res@cnLevels=20.0
  res@cnFillOn=False
  res@cnLineThicknessF =3.0
  res@cnLinesOn=True
  res@cnInfoLabelOn=False
  
  plot3=gsn_csm_contour(wks1,ssta_s1,res)
  plot4 =gsn_csm_contour(wks1,ssta_s2,res)
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

  gsn_panel(wks1, (/plot1(0),plot1(1)/), (/2,1/), pres)
   draw(plot1(2))
  frame(wks1) 

  
 end