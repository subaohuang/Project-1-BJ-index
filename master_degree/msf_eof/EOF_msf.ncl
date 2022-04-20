;******************************************************************
;EOF分解msf
;*******************************************************************


 function marker(data,time)
 begin
 mono=1
   markerx=data(0,time)
   markery=data(1,time)
   markerii=dim_pqsort_n(markerx,mono,0)

   markerxx=markerx(markerii)
   markeryy=markery(markerii)
   marker2=new((/2,dimsizes(markerxx)/) ,typeof(markerxx))
   marker2(0,:)=markerxx
   marker2(1,:)=markeryy
 return marker2
end


function band_pass(data)
begin
 ihp   = 2                            ; band_pass
  sigma = 1.0                           ; Lanczos sigma
  nWgt  = 25                           ; loose 60 months each end
  fcb   = 1./24.                        ; 2 years
  fca   = 1./84.                         ;7年
  wgts  = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )

   msf_500_band=wgt_runave_Wrap(data, wgts, 0)
  return msf_500_band
  end 
begin

diri="/mnt/f/share/era/"
fils0=systemfunc("ls "+diri+"era_nino3_4.ascii")
nino3_4=asciiread(fils0, 456, "float")
nino3_4_de=dtrend_n(nino3_4(0:407), True,0)
printVarSummary(nino3_4_de)
nino3_4_de!0="time"
nino3_4_de&time=ispan(1,408,1)

; enino=ind(nino3_4_nor.gt.1.0)
; lanila=ind(nino3_4_nor.lt.-1.0)
; Neutral=ind(nino3_4_nor.le.1.0.and.nino3_4_nor.ge.-1.0)
  ; nino3_4_nor=dim_standardize_n(nino3_4_de, 1, 0)
  ; str_enino=ind(nino3_4_nor.gt.1.0)
  ; str_lanina=ind(nino3_4_nor.lt.-1.0)
  ; weak_enino=ind(nino3_4_nor.gt.0.5.and.nino3_4_nor.lt.1.0)
  ; weak_lanina=ind(nino3_4_nor.gt.-1.0.and.nino3_4_nor.lt.-0.5)
  ; Neutral=ind(nino3_4_nor.gt.-0.5.and.nino3_4_nor.lt.0.5)

str_enino=ind(nino3_4_de.gt.1.0)
str_lanina=ind(nino3_4_de.lt.-1.0)
weak_enino=ind(nino3_4_de.gt.0.5.and.nino3_4_de.lt.1.0)
weak_lanina=ind(nino3_4_de.gt.-1.0.and.nino3_4_de.lt.-0.5)
Neutral=ind(nino3_4_de.gt.-0.5.and.nino3_4_de.lt.0.5)


diri1="/mnt/f/share/era/data/"
diri_pi="/mnt/f/share/era/"
fils1=systemfunc("ls "+diri1+"msf_time_era.nc")
f1=addfile(fils1, "r")
msf_time=f1->msf_time(:,:,0:407)     ;1979年到2012年12月份
; nino3_4_de_3d=conform_dims(dimsizes(msf_time), nino3_4_de, 2)
; printVarSummary(nino3_4_de_3d)

msf_time_de=dtrend_n(rmMonAnnCycLLT(msf_time),True,2)/10^9
; msf_time_a=dtrend_n(msf_time,True,2)/10^9
copy_VarCoords(msf_time, msf_time_de)
; copy_VarCoords(msf_time, msf_time_a)
; print("msf_time_de="+msf_time_de(:,:,1))
neval  = 30                    ; calculate 2 EOFs out of 7 ;specifies the number of eigenvalues and eigenvectors to be returned
ev_msf=eofunc_n_Wrap(msf_time_de,neval,True,2)
printVarSummary(ev_msf)
optETS      = True
optETS@jopt = 1
ev_timef=eofunc_ts_n_Wrap(msf_time_de, ev_msf, False,2)

ev_timef_std=new((/3,dimsizes(ev_timef&time)/), double)
eval=ev_msf@eval
tsmean=ev_timef@ts_mean


ev_timef_std(0,:)=ev_timef(0,:)/eval(0)
ev_timef_std(1,:)=ev_timef(1,:)/eval(1)
ev_timef_std(2,:)=ev_timef(2,:)/eval(2)

ev_sum=ev_msf(1,:,:)

copy_VarCoords(ev_msf(0,:,:), ev_sum)
printVarSummary(ev_sum)

  do n=0,neval-1
      ev_timef(n,:) = ev_timef(n,:) + ev_timef@ts_mean(n)  ; add time series mean
   end do      
 

;**********************************************************************
;画空间模态
;**********************************************************************
wks1=gsn_open_wks("png", ""+diri_pi+"msf_eof")
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString="era"
    res@gsnRightString="EOF-1(43.28%)"
    gsn_define_colormap(wks1,"hotcolr_19lev")
  
    res@cnFillColors=(/2,4,5,6,7,10,10,14,16,17,18,20/)
    ; res@lbLabelStrings=ispan(-6,6,1)
    ; res@cnFillPalette  = "GMT_polar"
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -0.05
    res@cnMaxLevelValF = 0.05
    res@cnLevelSpacingF = 0.01
    res@cnLineLabelFontHeightF=0.01
    res@cnFillOn = True 
    res@gsnContourNegLineDashPattern=1
    res@cnLinePalette="BlueDarkRed18"
    res@cnLineLabelsOn=True
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False
    res@lbLabelBarOn=False
    res@vpWidthF=0.7
    
     ; res@cnExplicitLineLabelsOn=True
     ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
       ; res@cnLabelScalingMode="ScaleFactor"
       ; res@cnLabelScaleValueF=10^9
      
    ; res@lbOrientation="Vertical"
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(100,290,40)
    res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","70~S~o~N~W"/)
    res@tmXBMinorOn=False
    ; res@vcRefMagnitudeF = 2.0                ; define vector ref mag
    ; res@vcRefLengthF    = 0.04        ; define length of vec ref
    ; res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    ; res@vcMinDistanceF  = 0.02       ; thin out vectors
    ; res@vcMapDirection  = False
    
    ; ; res@vcFillArrowHeadXF=0.1
    ; res@vcFillArrowMinFracWidthF=1.0
    ; res@vcFillArrowHeadMinFracYF=1.0
    ; ; res@vcFillArrowHeadXF=1.0
    ; res@vcFillArrowWidthF=0.05
    ; res@vcFillArrowHeadYF=0.15

    ; res@vcRefAnnoOrthogonalPosF  = -1.05
    ; res@vcRefAnnoPerimOn=False
    res@tmYRMode             = "Automatic"          ; turn off special labels on right axis

   plot=new(4,graphic)
   plot(0)=gsn_csm_pres_hgt(wks1,ev_msf(0,:,:),res)
    res@gsnRightString="EOF-2(26.29%)"
   plot(1)=gsn_csm_pres_hgt(wks1,ev_msf(1,:,:),res)

  
   ;---Resources for the polyline  添加下垫面的陆地
  ployres                  = True          ; polyline mods desired
  ployres@gsLineThicknessF = 15.0            ; line thickness   
  ployres@gsLineColor      = "black"         ; color of lines
  ployres@gsLineDashPattern=15
  ;---Attach the polyline
  ploylevel=1000
  ploylon1=fspan(100, 149, 25)
  ploylon2=fspan(280, 290, 11)
  ploylevel1=conform_dims(dimsizes(ploylon1), ploylevel, 0)
  ploylevel2=conform_dims(dimsizes(ploylon2), ploylevel, 0)
  
  dum1=new(2,graphic)
  dum2=new(2,graphic) 

   do i=0,1
  dum1(i)= gsn_add_polyline(wks1,plot(i),  ploylon1,ploylevel1,  ployres)
  end do
 ployres@gsLineDashPattern=0
  do i=0,1
 dum2(i)= gsn_add_polyline(wks1,plot(i),  ploylon2,ploylevel1,  ployres)
  end do

  ;************************************************************************
  ;画时间变化序列
  ;*************************************************************************

  res2=True
  res2@gsnFrame=False
  res2@gsnDraw=False 

  res2@gsnLeftString="era"
  res2@gsnRightString="PC-1"
  res2@gsnLeftStringFontHeightF=0.01
  res2@gsnRightStringFontHeightF=0.01


  res2@gsnAboveYRefLineColor = "gray25"              ; above ref line fill red
  res2@gsnBelowYRefLineColor = "gray75"             ; below ref line fill blue
  res2@gsnYRefLine           = 0.0             ; create a reference line
  res2@tmXBMode="Explicit"
  res2@tmXBValues=ispan(1,408,60)
  res2@tmXBLabels=ispan(1979,2012,5)
  res2@tmXBMinorValues=ispan(1,408,12)
  res2@tmXBMinorOn=True
  ; res2@tmXBMinorValues=ispan(1,38,1)
  res2@tmXBLabelFontAspectF=1.0
  res2@tmXBLabelFontHeightF=0.008
  res2@tmXBLabelFont=30     ;黑体
  res2@tmXBLabelFontThicknessF=2.0
  res2@tmXBMajorLengthF=0.005        ;調整刻度的长度
  res2@tmXBMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
  ; res2@tmYLLabelFontAspectF=2.0      ;字体的高除以长的比率

  res2@tmXBMinorLengthF=0.003
  res2@tmXBMinorOutwardLengthF=0.003
 
  res2@tmYLMajorLengthF=0.005        ;調整刻度的长度
  res2@tmYLMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
  res2@tmYLLabelFont=30
  res2@tmYLLabelFontHeightF=0.008
  res2@tmYLMinorOn=False
  res2@tmXTOn=False
  res2@tmYROn=False

  ; res2@vpHeightF=0.13
  ; res2@vpWidthF=0.7
  ; res2@vpXF=0.15
  ; res2@vpYF=0.95

plot(2)=gsn_csm_xy(wks1, ev_timef&time, ev_timef_std(0,:), res2)
plot(3)=gsn_csm_xy(wks1, ev_timef&time, ev_timef_std(1,:), res2)


pres=True
pres@gsnMaximize=True
pres@gsnPanelLabelBar=True
pres@pmLabelBarWidthF=0.4

gsn_panel(wks1, (/plot(0),plot(1)/), (/2,1/), pres)

; delete(  res@cnLevelSelectionMode)
;     delete( res@cnMinLevelValF )

; delete(   res@cnMaxLevelValF )
;     delete( res@cnLevelSpacingF )

;----------------------------------------------
;画正负累加的图
;-----------------------------------------------
wks2=gsn_open_wks("png", ""+diri_pi+"正负相位累加")
delete(res@cnFillPalette)
gsn_define_colormap(wks2,"hotcolr_19lev")
   res@cnMinLevelValF = -0.05
   res@cnMaxLevelValF = 0.05
   res@cnLevelSpacingF = 0.01
   res@cnFillColors=(/2,4,5,6,7,10,10,14,16,17,18,20/)
plot1=gsn_csm_pres_hgt(wks2,ev_sum,res)

dum3=new(1,graphic)
 ployres@gsLineDashPattern=15
dum3= gsn_add_polyline(wks2,plot1,  ploylon1,ploylevel1,  ployres)
 ployres@gsLineDashPattern=0
dum4= gsn_add_polyline(wks2,plot1,  ploylon2,ploylevel1,  ployres)

draw(plot1)
frame(wks2)


  timeres2=True

  timeres2@gsnFrame=False
  timeres2@gsnDraw=False 

  timeres2@gsnLeftString="era"
  timeres2@gsnRightString="EOF time series"
  timeres2@gsnLeftStringFontHeightF=0.01
  timeres2@gsnRightStringFontHeightF=0.01

  timeres2@xyLineColors=(/"blue","red"/)
  timeres2@xyMonoLineColor = False
 
  timeres2@xyMonoLineThickness = False
  timeres2@xyLineThicknesses = (/5, 5/)
  timeres2@xyMonoDashPattern = False
  timeres2@xyDashPatterns = 0
   


  timeres2@gsnYRefLine           = 0.0             ; create a reference line
  timeres2@tmXBMode="Explicit"
  timeres2@tmXBValues=ispan(1,408,60)
  timeres2@tmXBLabels=ispan(1979,2012,5)
  timeres2@tmXBMinorValues=ispan(1,408,12)
  timeres2@tmXBMinorOn=True
  ; timeres2@tmXBMinorValues=ispan(1,38,1)
  timeres2@tmXBLabelFontAspectF=1.0
  timeres2@tmXBLabelFontHeightF=0.008
  timeres2@tmXBLabelFont=30     ;黑体
  timeres2@tmXBLabelFontThicknessF=2.0
  timeres2@tmXBMajorLengthF=0.005        ;調整刻度的长度
  timeres2@tmXBMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
  ; timeres2@tmYLLabelFontAspectF=2.0      ;字体的高除以长的比率

  timeres2@tmXBMinorLengthF=0.003
  timeres2@tmXBMinorOutwardLengthF=0.003
 
  timeres2@tmYLMajorLengthF=0.005        ;調整刻度的长度
  timeres2@tmYLMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
  timeres2@tmYLLabelFont=30
  timeres2@tmYLLabelFontHeightF=0.008
  timeres2@tmYLMinorOn=False
  timeres2@tmXTOn=False
  timeres2@tmYROn=False

  timeres2@vpHeightF=0.3
  timeres2@vpWidthF=0.7
  timeres2@vpXF=0.15
  timeres2@vpYF=0.95

  timeres2@pmLegendDisplayMode    = "Always"            ; turn on legend
  timeres2@pmLegendSide           = "Top"               ; Change location of 
  timeres2@pmLegendParallelPosF   = 0.9                  ; move units right
  timeres2@pmLegendOrthogonalPosF = -0.35                ; move units down
  timeres2@pmLegendWidthF         = 0.10                ; Change width and
  timeres2@pmLegendHeightF        = 0.10                ; height of legend.
  timeres2@lgPerimOn              = False               ; turn off box around
  timeres2@lgLabelFontHeightF     = .01                 ; label font height
  timeres2@xyExplicitLegendLabels = (/"PC-1","PC-2"/)         ; create explicit labels
  timeres2@lgItemOrder = (/ 1, 0/) 
  
  wks3=gsn_open_wks("png", ""+diri_pi+"时间序列")
  wks4=gsn_open_wks("png", ""+diri_pi+"第三模态的时间序列")

  plot2=gsn_csm_xy(wks3, ev_timef&time, (/runave(ev_timef_std(0,:),3,0),runave(ev_timef_std(1,:),3,0)/), timeres2)
  ; plot2=gsn_csm_xy(wks3, ev_timef&time, (/band_pass(ev_timef_std(0,:)),band_pass(ev_timef_std(1,:))/), timeres2)
  delete(timeres2@gsnYRefLine)
  delete(timeres2@tmXBValues)
  delete(timeres2@tmXBLabels)
  delete(timeres2@tmXBMinorValues)
  delete(timeres2@gsnXRefLine)

  timeres2@gsnYRefLine=(/0,1,-1/)
  time_enino=(/1987,1990,1991,1992,1994,2002,2004,2009/)
  timef3=month_to_annual(ev_timef_std(2,:), 1)
  timef3!0="time"
  timef3&time=ispan(1979,2012,1)

  printVarSummary(timef3)
  timeres2@tmXBValues=ispan(1979,2012,5)
  timeres2@tmXBLabels=ispan(1979,2012,5)
  timeres2@tmXBMinorValues=ispan(1979,2012,1)
  timef3_new=dim_standardize(timef3,0)
  copy_VarMeta(timef3, timef3_new)
  plot3=gsn_csm_xy(wks4, timef3&time, dim_standardize(timef3,0), timeres2)
  scares=True
  sca=gsn_add_polymarker(wks4, plot3, time_enino ,timef3_new({time_enino}),  scares)

draw(plot2)
draw(plot3)
frame(wks3)
frame(wks4)



;-----------------------------------------------------------------
;pc1和pc2 scatter图：
;-------------------------------------------------------------------
 wks4=gsn_open_wks("png", ""+diri_pi+"pc1和pc2散点图")
   
    ress=True
    ress@gsnDraw=False
    ress@gsnFrame=False
    ress@xyLineThicknesses = (/ 1.0,2.0/)          ; make second line thicker
    ress@xyLineColors      = (/"black"/)          ; change line color
    ress@xyMarkerColor = "white"
    ress@xyMarkLineModes   = (/"Markers","lines"/)
    ; res@xyMonoLineColor     = False
    ; ress@xyExplicitLegendLabels = (/"Anual change", "mon_change", "Linear trend"/) 
    ress@xyMarkers         =  16                      ; choose type of marker  
    ress@xyMarkerSizeF     = 0.0001                 ; Marker size (default 0.01)
    ress@xyDashPatterns      = 1                      ; solid line 
    ress@gsnYRefLine           = 0.0             ; create a reference line
    ress@gsnXRefLine            =0.0
    ress@gsnYRefLineThicknesses=2.0
    ress@gsnXRefLineThicknessF=2.0
    ress@tiYAxisString="PC-2"
    ress@tiXAxisString="PC-1"
    ; ress@gsnLeftString="(a)dSLP"
    ; ress@gsnRightString="msf_zeroline"
    ; ress@gsnRightString=""+slp_data(3,1)+""
    ress@gsnLeftStringFontHeightF=0.01
    ress@gsnRightStringFontHeightF=0.01
   
    ; ress@trXMaxF=dimsizes(intesion&time)
    ; ress@tmXBValues=ispan(13,dimsizes(intesion&time),60)
    ; ress@tmXBLabels=ispan(1980,2015,5)
    
    ress@cnLevelSelectionMode ="ManualLevels"
    ress@cnMinLevelValF = -0.1
    ress@cnMaxLevelValF = 0.1
    ress@cnLevelSpacingF = 0.02






    ress@vpHeightF=0.5
    ress@vpWidthF=0.7
    ress@vpXF=0.15
    ress@vpYF=0.95
    
    ress@tmXBLabelFontAspectF=1.0
    ress@tmXBLabelFontHeightF=0.008
    ress@tmXBLabelFont=30     ;黑体
    ress@tmXBLabelFontThicknessF=2.0
    ress@tmXBMajorLengthF=0.005        ;調整刻度的长度
    ress@tmXBMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
    ; ress@tmYLLabelFontAspectF=2.0      ;字体的高除以长的比率

    ress@tmXBMinorLengthF=0.003
    ress@tmXBMinorOutwardLengthF=0.003
    
    ress@tmYLMajorLengthF=0.005        ;調整刻度的长度
    ress@tmYLMajorOutwardLengthF=0.005  ;调整刻度的长度，注意2个命令要同时使用
    ress@tmYLLabelFont=30
    ress@tmYLLabelFontHeightF=0.008
    ress@tmYLMinorOn=False
    ress@tmXTOn=False
    ress@tmYROn=False
    

 
   ; marker_str_enino=marker(ev_timef_std,str_enino)
   ; marker_str_lanina=marker(ev_timef_std,str_lanina)
   ; marker_weak_enino=marker(ev_timef_std,weak_enino)
   ; marker_weak_lanina=marker(ev_timef_std,weak_lanina)

 ;-------------------------------------------------------------
 ;散点拟合:对总体进行拟合
 ; ; ------------------------------------------------------------
   x=ev_timef(0,:)
   y=ev_timef(1,:)

   mono = 1                            ; ascending=1 , descending=-1
   ii   =  dim_pqsort_n(x,mono,0)      ; ascending order 
   xx   = x(ii)                        
   yy   = y(ii)
   
   n = 3
   c = lspoly_n(x, y, 1, n, 0)    ; all weights are set to one 
   Y = c(0) + c(1)*xx + c(2)*xx^2 
   

;-------------------------------------------------------------------
; 散点拟合：对极端时间进行拟合
;-------------------------------------------------------------
  ;  str=new((/dimsizes(str_enino)+dimsizes(str_lanina)/),typeof(str_enino))
  ;  str(0:dimsizes(str_enino)-1)=str_enino
  ;  str(dimsizes(str_enino):dimsizes(str_lanina)+dimsizes(str_enino)-1)=str_lanina
  ;  mono=1
  ;  ;                             ; ascending=1 , descending=-1
  ;  ii_str   =  dim_pqsort_n(str,mono,0)      ; ascending order 
  ;  str_new=str(ii_str)
 
   
  ;  x1=ev_timef(0,str_new)
  ;  y1=ev_timef(1,str_new)

  ;  mono = 1                            ; ascending=1 , descending=-1
  ;  ii   =  dim_pqsort_n(x1,mono,0)      ; ascending order 
  ;  xx   = x1(ii)                        
  ;  yy   = y1(ii)
  ; print("ii="+ii)
  ; print("xx="+xx)
  ;  n = 3
  ;  c = lspoly_n(x1, y1, 1, n, 0)    ; all weights are set to one 
  ;  Y = c(0) + c(1)*xx + c(2)*xx^2 
   
;-------------------------------------------------------
 
   ; x=ev_timef(0,:)
   ; y=ev_timef(1,:)

   ; mono = 1                            ; ascending=1 , descending=-1
   ; ii   =  dim_pqsort_n(x,mono,0)      ; ascending order 
   ; xx   = x(ii)                        
   ; yy   = y(ii)

   ; pc1_inter=ispan(-95,60,5)
   
   ; yy_avg=new(dimsizes(pc1_inter)-1,typeof(yy))
   ; xx_avg=new(dimsizes(pc1_inter)-1,typeof(xx))
   ; printVarSummary(yy_avg)

   ; do i=-95,55,5
   ;  m=ind(xx.ge.i.and.xx.le.(i+5))
   ;  yy_avg((i+95)/5)=dim_avg_n(yy(m),0)
   ;  xx_avg((i+95)/5)=dim_avg_n(xx(m),0)
   ;  delete(m)
   ; end do

;--------------------------------------------------------------------


   marker_str_enino=marker(ev_timef,str_enino)
   marker_str_lanina=marker(ev_timef,str_lanina)
   marker_weak_enino=marker(ev_timef,weak_enino)
   marker_weak_lanina=marker(ev_timef,weak_lanina)
   marker_neutral=marker(ev_timef,Neutral)

   polyres                   = True
   polyres@gsMarkerIndex     = 16          ; polymarker style
   polyres@gsMarkerSizeF     = 10         ; polymarker size
   polyres@gsMarkerColor     ="red"
    
   

   ; plot3=gsn_csm_xy(wks4,xx,(/yy/),ress)

   plot3=gsn_csm_xy(wks4,xx,(/yy,Y/),ress)

   dum1=gsn_add_polymarker(wks4, plot3, marker_str_enino(0,:), marker_str_enino(1,:), polyres)
   polyres@gsMarkerColor     ="pink"
   dum2=gsn_add_polymarker(wks4, plot3, marker_weak_enino(0,:), marker_weak_enino(1,:), polyres)
   polyres@gsMarkerColor     ="purple"
   dum3=gsn_add_polymarker(wks4, plot3, marker_str_lanina(0,:), marker_str_lanina(1,:), polyres)
   polyres@gsMarkerColor     ="blue"
   dum4=gsn_add_polymarker(wks4, plot3, marker_weak_lanina(0,:), marker_weak_lanina(1,:), polyres)
   polyres@gsMarkerSizeF     = 8         ; polymarker size
   polyres@gsMarkerColor     ="black"
   dum5=gsn_add_polymarker(wks4, plot3, marker_neutral(0,:), marker_neutral(1,:), polyres)
   ;  polyres@gsMarkerColor     ="green"
   ; dum6=gsn_add_polymarker(wks4, plot3, xx_avg,yy_avg , polyres)
     ; ploylineres=True
    ; dum7=gsn_add_polyline(wks4, plot3, xx_avg, Y, ploylineres)

   
  txid_tr = new(5, graphic)
  amid_tr = new(5, graphic)
  txres                       = True
  txres@txPerimOn             = True
  txres@txFontHeightF         = 0.015

  amres_tr                  = True
  amres_tr@amParallelPosF   =  0.38    ; This is the right edge of the plot.
 
  k=(/-0.47,-0.47,-0.42,-0.42,-0.37/)
  l=(/0.20,0.38,0.20,0.38,0.20/)
  letter = (/"str enino","str lanina","weak enino","weak lanina","neutral year"/)
  color=(/"red","purple","pink","blue","black"/)
  do i = 0, 4
     amres_tr@amOrthogonalPosF=k(i)
     amres_tr@amParallelPosF   =  l(i)    ; This is the right edge of the plot.
  ;---Create text strings
      tr_label = letter(i)
     
      txres@txBackgroundFillColor = "transparent"
      txres@txPerimOn = False
      txres@txFontColor=color(i)
      txid_tr(i) = gsn_create_text(wks4, tr_label, txres)

  ;---Attach text strings to plot
      amid_tr(i) = gsn_add_annotation(plot3, txid_tr(i), amres_tr)

  end do


   draw(plot3)
   frame(wks4)
  
 ;----------------------------------------------------------------
  ;重构eof
 ;------------------------------------------------------------------

 print("c(2)="+c(2))
 print("c(1)="+c(1))
 print("c(0)="+c(0))

 pc2_re= y - c(2)*x^2
 
 ev_timef(1,:)=(/pc2_re/)

 z20  = eof2data_n(ev_msf,ev_timef,2)             ; reconstruct array using the 8 EOFs
                                                  ; z8(ntim,nlat,mlon
copy_VarMeta(msf_time_de, z20)
printVarSummary(z20)
; print("z20="+z20(:,:,1))
; print("msf_time_de="+msf_time_de(:,:,1))
nlev=dimsizes(z20&level)
nlon=dimsizes(z20&lon)
ntime=dimsizes(z20&time)



;*********************************************写成nc文件,注意坐标变量之间的dimension也要命名。

   diro=diri
    filo1="msf_time_de_re.nc"
    system("/bin/rm -f "+diro + filo1+"")
    fout=addfile(diro + filo1,"c")
    setfileoption(fout,"DefineMode",True)
   ;======================================文件名
    fAtt               = True            ; assign file attributes
    fAtt@title         = "msf重构结果"    
    fAtt@creation_date = "2017.8.30"
    fileattdef( fout, fAtt)            ; copy file attributes    
    ;----------------------------坐标变量的定义
    dimnames=(/"level","lon","time"/)
    dimsize=(/nlev,nlon,ntime/)
    dimunlim=(/False,False,False/)
    filedimdef(fout, dimnames, dimsize, dimunlim)
 ;-----------------------------定义变量名，类型和维数
    filevardef(fout, "level",typeof(z20&level), (/"level"/))
    filevardef(fout, "lon", typeof(z20&lon),(/"lon"/))
    filevardef(fout, "time", typeof(z20&time),(/"time"/))
    filevardef(fout, "msf_time_de_re", typeof(z20),(/"level","lon","time"/))


    filevarattdef(fout, "level", z20&level)
    filevarattdef(fout, "lon", z20&lon)

    filevarattdef(fout, "msf_time_de_re", z20)
    
   
    setfileoption(fout,"DefineMode",False)

    fout->level=(/z20&level/)
    fout->lon=(/z20&lon/)
    fout->time=(/z20&time/)
   
    fout->msf_time_de_re=(/z20/)
   

end

