 ;;hadisst_eof分析
begin
  diri="/mnt/f/share/ersst/"
  diri_pi="/mnt/f/share/era/data/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  
  timeind=ind(date.le.201012.and.date.ge.197901)

  sst=lonFlip(f1->sst(timeind,{50:-50},:))
 

  sst=(/rmMonAnnCycTLL(sst)/)
  sst_clm=dim_avg_n_Wrap(sst, 0)

  sst_annual=month_to_annual(sst, 1)
  sst_annual&year=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  printVarSummary(sst_annual)
  
  ; sst_sm=sst_annual
  sst_sm=runave_n_Wrap(sst_annual, 5, 1, 0)
  ; sst_sm=sst_annual
  ; sst_sm=sst_annual
  sst_sm_clm=dim_avg_n_Wrap(sst_sm, 0)
  printVarSummary(sst_sm)

  neval  = 30                   ; calculate 2 EOFs out of 7 ;specifies the number of eigenvalues and eigenvectors to be returned
  ev_sst=eofunc_n_Wrap(sst_sm,neval,True,0)    ;;;得到的标准化后的场，要使得不标准化
  printVarSummary(ev_sst)
   
  

  optETS      = True
  optETS@jopt = 0
  do ne=0,neval-1
     ev_sst(ne,:,:) = ev_sst(ne,:,:)*sqrt(ev_sst@eval(ne) ) ;;;使得量级能够相同
  end do
  ; ev_timef=eofunc_ts_n_Wrap(sst_sm, ev_sst, False,0)
  ev_timef=eofunc_ts_n_Wrap(sst_sm, ev_sst, False,0)
  printVarSummary(ev_timef)
 
   ; do n=0,neval-1
   ;    ev_timef(n,:) = ev_timef(n,:) + ev_timef@ts_mean(n)  ; add time series mean
   ; end do         
   printVarSummary(ev_sst)
   
  
   ev_timef_std=new((/3,dimsizes(ev_timef&year)/), double)
   eval=ev_sst@eval

  ; tsmean=ev_timef@ts_mean
  
  ev_timef_std(0,:)=ev_timef(0,:)/eval(0)
  ev_timef_std(1,:)=ev_timef(1,:)/eval(1)
  ev_timef_std(2,:)=ev_timef(2,:)/eval(2)
   

  ; opt1=diri_pi+"hadsst_timef1_1958_2012.ascii"
  ; opt2=diri_pi+"hadsst_timef2_1958_2012.ascii"
  ; asciiwrite(opt1, sprintf ("%6.2f,", ev_timef_std(0,:)))
  ; asciiwrite(opt2, sprintf ("%6.2f,", ev_timef_std(1,:)))
  
  ; ev_sst=(/1000*ev_sst/)

  wks=gsn_open_wks("png", ""+diri_pi+"hadisst_eof")
  wks1=gsn_open_wks("pdf", ""+diri_pi+"hadisst_eof_enso")
  ; wks2=gsn_open_wks("pdf", ""+diri_pi+"hadisst_eof2_recon")

  res = True
  res@gsnAddCyclic=False
  res@gsnDraw=False
  res@gsnFrame=False

  ; res@cnFillPalette  = "temp_diff_18lev"
  res@cnFillPalette="BlueRed"
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn = False

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -0.3
  res@cnMaxLevelValF       =   0.3
  res@cnLevelSpacingF      =  0.05
  
  res@lbLabelBarOn=True

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -50.0
  res@mpMaxLatF             =  50.0
  res@mpMinLonF             = min(ev_sst&longitude)
  res@mpMaxLonF             = max(ev_sst&longitude)
  res@mpCenterLonF          =avg(ev_sst&longitude)


  plot=new(2,graphic)
  
   
 
  plot(0)=gsn_csm_contour_map(wks1,ev_sst(1,:,:), res)
  




  ; plot(0)=gsn_csm_contour_map(wks, ev_rot(0,:,:), res)
  ; plot(1)=gsn_csm_contour_map(wks, ev_rot(1,:,:), res)
  
  lineres=True 
  lineres@gsnDraw=False
  lineres@gsnFrame=False
  lineres@gsnYRefLine=0.0
  lineres@xyLineColors=(/"blue","red"/)
  lineres@xyDashPatterns=(/0,0/)
 
  lineres@tmXTOn=False
  lineres@tmYROn=False
  lineres@tmXBMode="Explicit"
  lineres@tmXBMinorValues=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  lineres@tmXBValues=ispan(1,dimsizes(sst_annual(:,1,1)),10)
  lineres@tmXBLabels=ispan(1958,2012,10)
  lineres@tmYLLabelFontHeightF=0.01
  lineres@tmXBLabelFontHeightF=0.01
  lineres@gsnYRefLine=0.0
  
  lineres@vpHeightF=0.24
  lineres@vpWidthF=0.75
     

  plot(1)=gsn_csm_xy(wks1, ispan(1,dimsizes(sst_annual(:,1,1)),1), (/ev_timef_std(0,:),ev_timef_std(1,:)/), \
                     lineres)
  
  ; draw(plot(2))
  ; frame(wks1)
  pres=True
  pres@gsnPanelLabelBar=True
  gsn_panel(wks1, plot, (/2,1/), pres)



end