
;;********************************************************
;;将NOAA的降水场，Had的SST以及计算的ud进行MV-EOF分解
;;
;;linshuheng 2020.01.07
;**********************************************************
   
begin
;;;读入hadisst
 
  diri="/home/ys17-19/lsh/data/sst/"


  fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  
  timeind=ind(date.le.201212.and.date.ge.195801)
  sst=lonFlip(f1->sst(timeind,{60:-60},:))
  sst=(/rmMonAnnCycTLL(sst)/)
  sst_annual=month_to_annual(sst, 1)
  sst_annual&year=ispan(1,dimsizes(sst_annual(:,1,1)),1)
  printVarSummary(sst_annual)
  sst_sm=runave_n_Wrap(sst_annual, 5, 1, 0)
  printVarSummary(sst_sm)
  
  ;;插值到2.5乘2.5网格
   yo=fspan(-59.5,58,48)
   xo=fspan(0.5, 358, 144)
  ;;
   
  sst_sm_lin=linint2_Wrap(sst_sm&longitude, sst_sm&latitude(::-1), sst_sm(:,::-1,:), \
                          True, xo, yo, 0)
 

  printVarSummary(sst_sm_lin)

  sst_sm_lin_std=dim_stddev_n_Wrap(sst_sm_lin, 0)
  sst_sm_lin_std!0="lat"
  sst_sm_lin_std!1="lon"
  printVarSummary(sst_sm_lin_std)

  sst_sm_lin=(/dim_standardize_n(sst_sm_lin, 1, 0)/)

  sst_sm_lin!1 = "lat"
  sst_sm_lin!2 = "lon"
  
    sst_sm_lin = (/SqrtCosWeight(sst_sm_lin)/)

  printVarSummary(sst_sm_lin)
  delete(diri)  
  delete(fils1)
  delete(f1)
  delete(date)
  delete(timeind)
  delete(sst)

;;;读入降水
  fils1=systemfunc("ls /home/yangsong3/data-observation/PREC/precip.mon.anom.nc") 
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201212.and.date.ge.195801)
  p=f1->precip(timeind,{60:-60},:)
  
  p=(/rmMonAnnCycTLL(p)/)
  p_annual=month_to_annual(p, 1)
  p_annual&year=ispan(1,dimsizes(p_annual(:,1,1)),1)
  printVarSummary(p_annual)
  p_sm=runave_n_Wrap(p_annual, 5, 1, 0)
  printVarSummary(p_sm)

  p_sm_lin=linint2_Wrap(p_sm&lon, p_sm&lat(::-1), p_sm(:,::-1,:), \
                        True, xo, yo, 0)
  
  p_sm_lin_std=dim_stddev_n_Wrap(p_sm_lin, 0)
  p_sm_lin_std!0="lat"
  p_sm_lin_std!1="lon"
  printVarSummary(p_sm_lin_std)

  p_sm_lin!1 = "lat"
  p_sm_lin!2 = "lon"
  
  p_sm_lin = (/SqrtCosWeight(p_sm_lin)/)
; 
  p_sm_lin=(/dim_standardize_n(p_sm_lin, 1, 0)/)
   
  

  delete(fils1)
  delete(f1)
  delete(date)
  delete(timeind)
;;;

;;;读入ud数据

  fils1=systemfunc("ls /home/ys17-19/lsh/data/MRE/ud_mre_global.nc") 
  f1=addfile(fils1, "r")
  date=cd_calendar(f1->time, 1)
  timeind=ind(date.le.201212.and.date.ge.195801)
  ud_r=f1->ud_new({850},{60:-60},:,timeind)
  ud=ud_r(time|:,lat|:,lon|:)
  
  delete(ud_r)
  ud=(/rmMonAnnCycTLL(ud(time|:,lat|:,lon|:))/)
  ud_annual=month_to_annual(ud, 1)
  ud_annual&year=ispan(1,dimsizes(ud_annual(:,1,1)),1)
  printVarSummary(ud_annual)
  ud_sm=runave_n_Wrap(ud_annual, 5, 1, 0)
  printVarSummary(ud_sm)

  ud_sm_lin=linint2_Wrap(ud_sm&lon, ud_sm&lat(::-1), ud_sm(:,::-1,:), \
                        True, xo, yo, 0)
 
  ud_sm_lin_std=dim_stddev_n_Wrap(ud_sm_lin, 0)
  ud_sm_lin_std!0="lat"
  ud_sm_lin_std!1="lon"
  printVarSummary(ud_sm_lin_std)
  ud_sm_lin=(/dim_standardize_n(ud_sm_lin, 1, 0)/)
     ud_sm_lin!1 = "lat"
  ud_sm_lin!2 = "lon"
    ud_sm_lin = (/SqrtCosWeight(ud_sm_lin)/)

;;;

 printVarSummary(p_sm_lin)

  sst_p=new((/dimsizes(p_sm_lin&year),dimsizes(p_sm_lin&lat)\
           ,dimsizes(p_sm_lin&lon),3/),typeof(p_sm_lin))

  sst_p(:,:,:,0)=sst_sm_lin
  sst_p(:,:,:,1)=p_sm_lin
  sst_p(:,:,:,2)=ud_sm_lin

  ; printVarSummary(sst_p)
  sst_p!3="lev"

;;;
  
  optEOF      = True
  optEOF@jopt = 0            ;;计算协方差振，当输入的矩阵是标准化的时候，是相关系数振
  neval  = 30                   ; calculate 2 EOFs out of 7 ;specifies the number of eigenvalues and eigenvectors to be returned
  ev_sst=eofunc_n_Wrap(sst_p,neval,optEOF,0)    ;;;得到的标准化后的场，要使得不标准化
  printVarSummary(ev_sst)
  ev_sst&lon@units="degrees_east"
  ev_sst&lat@units="degrees_north"

  do ne=0,neval-1
     ev_sst(ne,:,:,:) = ev_sst(ne,:,:,:)*sqrt(ev_sst@eval(ne)) ;;;使得量级能够相同
  end do

  
  ; ;;;将得到的eof进行旋转分析
  ;    ev_rot = eofunc_varimax_Wrap(ev_sst(0:9,:,:,:), -1 )
  ;    eofunc_varimax_reorder(ev_rot)
  optETS      = True
  optETS@jopt = 0          ;;;选择1是程序会帮你算标准化的数据，再计算t，这里输入的矩阵已经是标准化的了，无所谓
  ev_timef=eofunc_ts_n_Wrap(sst_p, ev_sst,optETS ,0)    
  printVarSummary(ev_timef)
  ev_timef_std=new((/2,dimsizes(ev_timef&year)/), double)
  eval=ev_sst@eval
  tsmean=ev_timef@ts_mean
  ev_timef_std(0,:)=ev_timef(0,:)/eval(0)
  ev_timef_std(1,:)=ev_timef(1,:)/eval(1)



  ;;;乘各自模态的标准差，以得到原来的距平场
     ev_sst(0,:,:,0)=(/ev_sst(0,:,:,0)*sst_sm_lin_std/)
     ev_sst(1,:,:,0)=(/ev_sst(1,:,:,0)*sst_sm_lin_std/)
     
     ev_sst(0,:,:,1)=(/ev_sst(0,:,:,1)*p_sm_lin_std/)
     ev_sst(1,:,:,1)=(/ev_sst(1,:,:,1)*p_sm_lin_std/)
  ;;;
 





  ; opt1=diri+"sst_p_ud_timef1_1958_2012.ascii"
  ; opt2=diri+"sst_p_ud_timef2_1958_2012.ascii"
  ; asciiwrite(opt1, sprintf ("%6.2f,", ev_timef_std(0,:)))
  ; asciiwrite(opt2, sprintf ("%6.2f,", ev_timef_std(1,:)))
  



  mode=0

  wks=gsn_open_wks("pdf", "hadisst_p_ud_mveof"+(/mode+1/))
  res = True
  res@gsnAddCyclic=False
  res@gsnDraw=False
  res@gsnFrame=False
  res@gsnLeftString="EOF"+(/mode+1/)
  res@gsnRightString=""+ev_sst@pcvar(mode)+""
  ; res@cnFillPalette  = "temp_diff_18lev"
  res@cnFillPalette="BlueRed"
  ; res@cnFillPalette="precip_diff_12lev"
  res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
  res@cnFillOn = True 
  res@cnLinesOn = False
 
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  res@cnMinLevelValF       =   -0.4
  res@cnMaxLevelValF       =   0.4
  res@cnLevelSpacingF      =  0.05

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
  ; res@cnMinLevelValF       =   -0.40
  ; res@cnMaxLevelValF       =   0.40
  ; res@cnLevelSpacingF      =  0.08

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


  plot(0)=gsn_csm_contour_map(wks, ev_sst(mode,:,:,0), res)
  res@cnFillPalette="precip_diff_12lev"

  res@cnMinLevelValF       =   -0.40
  res@cnMaxLevelValF       =   0.40
  res@cnLevelSpacingF      =  0.08

  plot(1)=gsn_csm_contour_map(wks, ev_sst(mode,:,:,1), res)
  
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
  lineres@tmXBLabels=ispan(1958,2012,10)
  lineres@tmYLLabelFontHeightF=0.01
  lineres@gsnYRefLine=0.0
  
  lineres@vpHeightF=0.25
  lineres@vpWidthF=0.75
  

  plot(2)=gsn_csm_xy(wks, ispan(1,dimsizes(sst_annual(:,1,1)),1), ev_timef_std(mode,:), \
                     lineres)
  
  pres=True
  pres@gsnPanelLabelBar=False
  pres@gsnPanelDebug=True
  gsn_panel(wks, plot, (/3,1/), pres)



end