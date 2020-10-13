;;MV-EOF分析
begin
;;;读入hadisst
  diri="/mnt/f/share/ersst/"
  diri_pi="/mnt/f/share/era/data/"
  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  
  timeind=ind(date.le.201212.and.date.ge.192001)
  sst=lonFlip(f1->sst(timeind,{60:-60},:))
  sst_annual=month_to_annual(sst, 1)
  sst_annual=(/dim_rmvmean_n(sst_annual, 0)/)
  sst_annual&year=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  printVarSummary(sst_annual)
  sst_sm=runave_n_Wrap(sst_annual, 3, 1, 0)
  printVarSummary(sst_sm)

; =================================================================
; create weights:  sqrt(cos(lat))   [or sqrt(gw) ]
; =================================================================
  rad    = 4.*atan(1.)/180.
  clat   = sst_sm&latitude           
  clat   = sqrt( cos(rad*clat) )                 ; gw for gaussian grid

; =================================================================
; weight all observations 
; =================================================================
  wsst_sm   = sst_sm                                 ; copy meta data
  wsst_sm   = sst_sm*conform(sst_sm, clat, 1)
  ; =================================================================
; weight all observations 
; =================================================================

 
  sst_sm_lin=wsst_sm
  sst_sm_lin!1="lat"
  sst_sm_lin!2="lon"


  printVarSummary(sst_sm_lin)
  delete(diri)  
  delete(fils1)
  delete(f1)
  delete(date)
  delete(timeind)

 


;;;
diri="/mnt/f/share/era/data/"
sst_p=sst_sm_lin


;;;
  optEOF      = True
  optEOF@jopt = 1           ;;计算协方差振，当输入的矩阵是标准化的时候，是相关系数振
  neval  =2                   ; calculate 2 EOFs out of 7 ;specifies the number of eigenvalues and eigenvectors to be returned
  ev_sst=eofunc_n_Wrap(sst_p,neval,optEOF,0)    ;;;得到的标准化后的场，要使得不标准化
 
  optETS      = True
  optETS@jopt = 1     ;;;选择1是程序会帮你算标准化的数据，再计算t，这里输入的矩阵已经是标准化的了，无所谓
  ev_timef=eofunc_ts_n_Wrap(sst_p, ev_sst,optETS,0)    
  printVarSummary(ev_timef)
  ; ev_timef_std=new((/2,dimsizes(ev_timef&year)/), double)
  ; eval=ev_sst@eval
  ; tsmean=ev_timef@ts_mean
  ; ev_timef_std(0,:)=ev_timef(0,:)/eval(0)
  ; ev_timef_std(1,:)=ev_timef(1,:)/eval(1)
   ev_timef=(/dim_standardize_n(ev_timef, 1, 1)/) 
   do ne=0,neval-1
     ev_sst(ne,:,:) = ev_sst(ne,:,:)*sqrt(ev_sst@eval(ne)) ;;;使得量级能够相同
  end do
  ev_sst&lat@units="degrees_north"
  ev_sst&lon@units="degrees_east"

   ev_timef_std=ev_timef
; =================================================================

 
  printVarSummary(ev_sst)
   print(ev_timef_std(0,:))
   
  


  opt1=diri+"sst_timef1_1920_2012.ascii"
  opt2=diri+"sst_timef2_1920_2012.ascii"
  asciiwrite(opt1, sprintf ("%6.2f,", ev_timef_std(0,:)))
  asciiwrite(opt2, sprintf ("%6.2f,", ev_timef_std(1,:)))

  
  wks=gsn_open_wks("pdf", ""+diri_pi+"hadisst_eof")
  res = True
  res@gsnAddCyclic=False
  res@gsnDraw=False
  res@gsnFrame=False
  res@cnFillPalette="BlueRed"
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn = False
 
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -0.4
  res@cnMaxLevelValF       =   0.4
  res@cnLevelSpacingF      =  0.08



  res@lbLabelBarOn=True
  res@lbOrientation   ="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = min(ev_sst&lat)
  res@mpMaxLatF             = max(ev_sst&lat)
  res@mpMinLonF             = min(ev_sst&lon)
  res@mpMaxLonF             = max(ev_sst&lon)
  res@mpCenterLonF          =avg(ev_sst&lon)
  
  res@tmXTOn=False
  res@tmYROn=False
  res@tmYLMode="Explicit"
  res@tmYLValues=ispan(-60,60,30)
  res@tmYLMinorValues=ispan(-60,60,10)
  res@tmYLMinorOn=True
  res@tmYLLabels=(/"60~S~o~N~S","30~S~o~N~S","0","30~S~o~N~N","60~S~o~N~N"/)
  res@tmYLLabelFontHeightF=0.01

  res@tmXBMode="Explicit"
  res@tmXBValues=fspan(0.5,358,13)
  res@tmXBMinorOn=False  
  res@tmXBLabelFontHeightF=0.01
  res@tmXBLabelFontThicknessF=0.02
  res@tmXBLabels=(/"0","30~S~o~N~E","60~S~o~N~E","90~S~o~N~E","120~S~o~N~E","150~S~o~N~E",\
            "180","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W","60~S~o~N~W","30~S~o~N~W","0"/)


  plot=new(3,graphic)


  plot(0)=gsn_csm_contour_map(wks, ev_sst(0,:,:), res)
  res@cnFillPalette="precip_diff_12lev"
  plot(1)=gsn_csm_contour_map(wks, ev_sst(1,:,:), res)
  
  lineres=True 
  lineres@gsnDraw=False
  lineres@gsnFrame=False
  lineres@gsnYRefLine=0.0
  lineres@xyLineColors="blue"
  lineres@xyDashPatterns=0
 
  lineres@tmXTOn=False
  lineres@tmYROn=False
  lineres@tmXBMode="Explicit"
  lineres@tmXBMinorValues=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  lineres@tmXBValues=ispan(1,dimsizes(sst_annual(:,1,1)),10)
  lineres@tmXBLabels=ispan(1920,2012,10)
  lineres@tmYLLabelFontHeightF=0.01
  lineres@gsnYRefLine=0.0
  
  lineres@vpHeightF=0.25
  lineres@vpWidthF=0.75
  

  plot(2)=gsn_csm_xy(wks, ispan(1,dimsizes(sst_annual(:,1,1)),1), ev_timef_std(1,:), \
                     lineres)
  
  pres=True
  pres@gsnPanelLabelBar=False
  pres@gsnPanelDebug=True
  gsn_panel(wks, plot, (/3,1/), pres)




end