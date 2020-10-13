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
  
  ssta_djf1=month_to_season(sst_m, "DJF")
  sst_clm=dim_avg_n_Wrap(ssta_djf1, 0)

  sst_m=(/rmMonAnnCycTLL(sst_m)/)
 
  


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
  

  ;;;两个时间段的海温异常
  nino3=dim_avg_n(ssta_s(:,0,{210:270}),1)
  nino3_std=dim_standardize_n(nino3, 1, 0)
  
  ;;;挑选年份 canino  enino year
  caind=ind(nino3_std.ge.1)
  laind=ind(nino3_std.le.-1)
  midind=ind(nino3_std.le.1.and.nino3_std.ge.-1)
  
  ssta_e=dim_avg_n_Wrap(ssta_s(caind,:,:),0)
  ssta_l=dim_avg_n_Wrap(ssta_s(laind,:,:),0)
  ssta_n=dim_avg_n_Wrap(ssta_s(midind,:,:),0)



  ; wks1=gsn_open_wks("pdf", "/mnt/f/share/era/data/ ENSO年份的下层海温")
     wks_type = "pdf"
   ;   wks_type@wkWidth = 8000
   ; wks_type@wkHeight = 8000
  wks1 = gsn_open_wks(wks_type,"/mnt/f/share/era/data/ enino年份的下层海温")          ; send graphics to PNG file
  wks2 = gsn_open_wks(wks_type,"/mnt/f/share/era/data/ lanina年份的下层海温")          ; send graphics to PNG file
  wks3 = gsn_open_wks(wks_type,"/mnt/f/share/era/data/ 正常年份-气候态")          ; send graphics to PNG file
  
  plot1=new(3,graphic)



  res=True
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnAddCyclic=False
  ; res@cnFillPalette="MPL_Spectral"
    res@cnFillPalette="matlab_jet"
  ; res@gsnSpreadColorStart=5
  ; res@gsnSpreadColorEnd=115
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
  res@tmXBMode="Explicit"
  res@tmXBValues=fspan(0.5, 359.5,7 )
  res@tmXBLabels=(/"0~S~o","60~S~o~N~E","120~S~o~N~E","180~S~o"\
                ,"120~S~o~N~W","60~S~o~N~W","0~S~o"/)
  res@tmXBMinorOn=True
  res@tmXBMinorValues=fspan(0.5,359.5,13)


  res@tmXTOn=False
  res@tmYROn=False
  
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =  -5
  res@cnMaxLevelValF       =  5
  res@cnLevelSpacingF      =  0.5



  ; res@gsnLeftString        =  "Enino "

  plot1(0)=gsn_csm_contour(wks1,ssta_e,res)
   ; res@gsnLeftString        =  "Lanina"
  plot1(1) =gsn_csm_contour(wks2,ssta_l,res)
  ; delete( res@gsnSpreadColorStart)
  ; delete(res@gsnSpreadColorEnd)
  res@cnFillPalette="matlab_jet"
  res@cnMinLevelValF       =  -1
  res@cnMaxLevelValF       =  1
  res@cnLevelSpacingF      =  0.1
   ; res@cnMinLevelValF       =  10
  ; res@cnMaxLevelValF       =  1
  ; res@cnLevelSpacingF      =  30
  res@lbLabelBarOn=True
  res@vpXF=0.233191
  res@vpYF=0.32
  res@vpWidthF=0.491055
  res@vpHeightF=0.264414
  res@lbOrientation="Vertical"
  res@pmLabelBarWidthF=0.15
  ; res@gsnLeftString        =  "(b)-(a)"
  ; plot1(2) =gsn_csm_contour(wks3,ssta_n,res)
  plot1(2) =gsn_csm_contour(wks3,ssta_n,res)
 
 ; ;---Create text strings
 ;    txres                       = True
 ;    txres@txPerimOn             = True
 ;    txres@txFontHeightF         = 0.013

 ;    amres_tr                  = True
 ;    amres_tr@amParallelPosF   =  0.46    ; This is the right edge of the plot.
 ;    amres_tr@amOrthogonalPosF =  0.43    ; This is the top edge of the plot.
 ;     tr_label = (/"c"/)
 ;     txres@txBackgroundFillColor = "White"
 ;     txres@txPerimOn = True
 ;     txid_tr = gsn_create_text(wks1, tr_label, txres)
 ; ;---Attach text strings to plot
 ;     amid_tr = gsn_add_annotation(plot1(2), txid_tr, amres_tr)

   ; delete(res@gsnLeftString)

  ; pres=True
  ; ; pres@gsnFrame=False
  ; pres@gsnPanelLabelBar=True
  ; pres@lbOrientation="Vertical"
  ; pres@pmLabelBarWidthF=0.15
  ; pres@pmLabelBarHeightF=0.5
  ; pres@gsnPanelDebug=True
  ; pres@gsnPanelBottom=0.36
  ; pres@pmLabelBarOrthogonalPosF=-0.01
  ; pres@lbLabelStride =2
  ; pres@gsnPanelFigureStrings=(/"a","b"/)

  ; gsn_panel(wks1, (/plot1(0),plot1(1)/), (/2,1/), pres)
   
   draw(plot1(0))
   draw(plot1(1))
   draw(plot1(2))
   frame(wks1) 
   frame(wks2)
   frame(wks3)

  
 end