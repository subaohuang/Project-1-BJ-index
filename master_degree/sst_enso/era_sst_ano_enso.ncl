function one_test(pv_ano,east_time,pv_clm)
begin

    pv_east=dim_avg_n_Wrap( pv_ano(:,:,east_time),2)
    n=dimsizes(east_time)
    var_pv_ano_east    =  dim_variance_n(pv_ano(:,:,east_time), 2)    ; 13.477
    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  =where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)

    t_pv_ano_east = (pv_east-pv_clm)/var_pv_ano_east

    confi_pv_east=student_t(t_pv_ano_east , n-1)

    copy_VarCoords(pv_east, confi_pv_east)

    ; confi_pv_east_se=confi_pv_east(::4,::4)

 return confi_pv_east
 end

function normoalized(data,mean)
begin

msf_enino_tr=data/mean
copy_VarCoords(data, msf_enino_tr)
return msf_enino_tr
end 

begin
diri0="/mnt/f/share/era/data/"
diri_pi="/mnt/f/share/era/"
fils0=systemfunc("ls "+diri0+"era_nino3_4.ascii")
nino3_4=asciiread(fils0, 456, "float")
nino3_4_de=nino3_4(0:407)
; printVarSummary(nino3_4_de)
nino3_4_de!0="time"
nino3_4_de&time=ispan(1,408,1)

nino3_4_nor=dim_standardize_n(nino3_4_de, 1, 0)

enino=ind(nino3_4_de.gt.1.0)
lanila=ind(nino3_4_de.lt.-1.0)
Neutral=ind(nino3_4_de.le.1.0.and.nino3_4_de.ge.-1.0)

 str_enino=ind(nino3_4_de.gt.1.0)
 str_lanila=ind(nino3_4_de.lt.-1.0)
 weak_enino=ind(nino3_4_de.gt.0.5.and.nino3_4_de.lt.1.0)
 weak_lanila=ind(nino3_4_de.gt.-1.0.and.nino3_4_de.lt.-0.5)


 diri1="/mnt/f/share/era/data/"
 ; fils1=systemfunc("ls "+diri1+"era_u_sst_slp.nc")
 fils1=systemfunc("ls "+diri1+"era_sur_global.nc")
 f1=addfile(fils1, "r")
 time=f1->time
 lon=f1->longitude
 lat=f1->latitude
 date=cd_calendar(time, 1)

 tind1=ind(date.eq.197901)
 tind2=ind(date.eq.201612)
 latind1=ind(lat.eq.30)
 latind2=ind(lat.eq.-30)
 lonind1=ind(lon.eq.100)
 lonind2=ind(lon.eq.290)
 
  ; sst_r=f1->sst(tind1:tind2,latind1:latind2,lonind1:lonind2)  ;;;就画热带太平洋的情况
 
 sst_r=f1->sst(tind1:tind2,latind1:latind2,:)

 sst=sst_r*sst_r@scale_factor+sst_r@add_offset-273.15
 copy_VarMeta(sst_r, sst)


 sst_enino=sst(enino,:,:)
 sst_lanina=sst(lanila,:,:)
 sst_neutral=sst(Neutral,:,:)


 sst_de=dtrend_n(rmMonAnnCycLLT(sst(latitude|:,longitude|:,time|:)),True,2)
 
 copy_VarCoords(sst(latitude|:,longitude|:,time|:), sst_de)
 
 sst_clm_de=dim_avg_n_Wrap(sst_de, 2)

 sst_enino_de=dim_avg_n_Wrap(sst_de(:,:,enino),2)
 sst_lanina_de=dim_avg_n_Wrap(sst_de(:,:,lanila),2)
 sst_neutral_de=dim_avg_n_Wrap(sst_de(:,:,Neutral),2)
 sst_sum_de=sst_enino_de+sst_lanina_de

 
 printVarSummary(sst_enino_de)
 copy_VarCoords(sst_enino_de,sst_sum_de )

  sst_str_enino_de=dim_avg_n_Wrap(sst_de(:,:,str_enino),2)
  sst_str_lanina_de=dim_avg_n_Wrap(sst_de(:,:,str_lanila),2)

  sst_weak_enino_de=dim_avg_n_Wrap(sst_de(:,:,weak_enino),2)
  sst_weak_lanina_de=dim_avg_n_Wrap(sst_de(:,:,weak_lanila),2)

  sst_sum_str_de=sst_str_enino_de+sst_str_lanina_de
  sst_sum_weak_de=sst_weak_enino_de+sst_weak_lanina_de

copy_VarCoords(sst_str_enino_de,sst_sum_str_de )
copy_VarCoords(sst_weak_lanina_de,sst_sum_weak_de )

;;;; t检验
 confi_enino=one_test(sst_de,enino,sst_clm_de)
 confi_lanina=one_test(sst_de,lanila,sst_clm_de)
 confi_neutral=one_test(sst_de,Neutral,sst_clm_de)

;;;






;;;;




;;;标准化
   dim=dimsizes(sst_str_enino_de)
 
   mean1=new((/dim(0),dim(1)/),typeof(sst_str_enino_de))
   mean1=dim_avg_n_Wrap(nino3_4_de(str_enino), 0)
 
   mean2=new((/dim(0),dim(1)/),typeof(sst_str_enino_de))
   mean2=dim_avg_n_Wrap(nino3_4_de(str_lanila), 0)

  
   mean3=new((/dim(0),dim(1)/),typeof(sst_str_enino_de))
   mean3=dim_avg_n_Wrap(nino3_4_de(weak_enino), 0) 

   mean4=new((/dim(0),dim(1)/),typeof(sst_str_enino_de))
   mean4=dim_avg_n_Wrap(nino3_4_de(weak_lanila), 0)

   sst_str_enino_de_nor=normoalized(sst_str_enino_de,mean1)
   sst_str_lanina_de_nor=normoalized(sst_str_lanina_de,mean2)
   sst_weak_enino_de_nor=normoalized(sst_weak_enino_de,mean3)
   sst_weak_lanina_de_nor=normoalized(sst_weak_lanina_de,mean4)

   sst_sum_str_de_nor=sst_str_enino_de_nor-sst_str_lanina_de_nor
   sst_sum_weak_de_nor=sst_weak_enino_de_nor-sst_weak_lanina_de_nor

  copy_VarCoords(sst_str_enino_de, sst_sum_str_de_nor)
  copy_VarCoords(sst_weak_lanina_de, sst_sum_weak_de_nor)
;;;
;-------------------------------------------------------------------------------
 wks=gsn_open_wks("pdf", diri_pi+"era_sst_enso")
   res = True
    res@gsnAddCyclic=False
    res@gsnDraw=False
    res@gsnFrame=False

    ; res@cnFillPalette  = "temp_diff_18lev"
    ; res@cnFillPalette="BlueDarkRed18"
    res@cnFillPalette="NCV_blu_red"
    res@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
    res@cnFillOn = True 
    res@cnLinesOn = True
    res@cnLineLabelsOn=False
    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    res@cnMinLevelValF       = -1.6
    res@cnMaxLevelValF       =  1.6
    res@cnLevelSpacingF      =  0.2
    res@cnFillDrawOrder="PreDraw"
    res@gsnLeftString="era"
    res@lbLabelBarOn=True

    res@lbOrientation = "Vertical"
    res@lbBoxEndCapStyle="TriangleBothEnds"

    ; res@lbLabelStrings=ispan(-4,4,1)

    res@pmLabelBarWidthF=0.1     ;colorbar的宽度
 

    res@tmXBMode="Explicit"

    res@tmXBValues=fspan(0,357.5,5)

    res@tmYLTickSpacingF=10
    res@tmYLTickStartF=-30
    res@tmYLTickEndF=30
    res@tmYROn=False
    res@tmXTOn=False
    res@tmYLLabelFontHeightF=0.01
    res@tmXBLabelFontHeightF=0.01
    res@tmXBLabelFontThicknessF=0.02
    res@tmXBMinorOn=True

    ; res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W"/)

    res@tmXBLabels=(/"0~S~o~","90~S~o~N~E","180~S~o~N~E","90~S~o~N~W","0~S~o~"/)
    

    res@tmYLLabels=(/"30~S~o~N~S","20~S~o~N~S","10~S~o~N~S","0~S~o~N","10~S~o~N~N","20~S~o~N~N","30~S~o~N~N"/)

    ; res@mpDataSetName = "Earth..4"
    ; res@mpDataBaseVersion = "mediumRes"
    res@mpFillAreaSpecifiers=(/"Land","water"/)      ;确定需要填充的具体位置
    res@mpSpecifiedFillColors=(/"white","transparent"/) ;填充位置相应的具体颜色

    res@mpOutlineOn = True; turn off outline 
    ; res@mpGeophysicalLineThicknessF= 4.          ; double the thickness of geophysical boundaries 
    res@mpNationalLineThicknessF= 2.               ; double the thickness of national boundaries 
    res@mpOutlineSpecifiers = (/"Land"/)
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = min(sst_de&latitude)
    res@mpMaxLatF             = max(sst_de&latitude)
    res@mpMinLonF             = min(sst_de&longitude)
    res@mpMaxLonF             = max(sst_de&longitude)
    res@mpCenterLonF          = avg(sst_de&longitude)


    plot=new(4,graphic)
    res@gsnRightString=" Ei Nino"
    res@gsnLeftString="era"
    plot(0)= gsn_csm_contour_map(wks, sst_enino_de, res)
    res@gsnRightString="La Nina"
    res@gsnLeftString="era"
    plot(1)= gsn_csm_contour_map(wks, sst_lanina_de, res)
    res@cnMinLevelValF       = -0.32
    res@cnMaxLevelValF       =  0.32
    res@cnLevelSpacingF      =  0.08
    res@gsnRightString="Neutral"
    res@gsnLeftString="era"
    plot(2)= gsn_csm_contour_map(wks, sst_neutral_de, res)
    res@cnMinLevelValF       = -1.0
    res@cnMaxLevelValF       =  1.0
    res@cnLevelSpacingF      =  0.2
    res@gsnRightString="Ei Nino + La Nina"
    res@gsnLeftString="era"
    plot(3)= gsn_csm_contour_map(wks, sst_sum_de, res)
;---------------------------------------------------

;;;t检验显著性检验打点

  res_pro=True
  res_pro@gsnDraw=False
  res_pro@gsnFrame=False
  res_pro@cnLevelSelectionMode = "ExplicitLevels"
  res_pro@cnLinesOn=False
  res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res_pro@cnMinLevelValF      = 0.00        ; set min contour level
  res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
  res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing
  ; res_pro@cnLevels = ispan(-6,6,1)
  res_pro@cnLineLabelsOn=False
  res_pro@tmYRMode = "Automatic"
  res_pro@cnInfoLabelOn=False
  ; res_pro@cnLineLabelBackgroundColor="Transparent"
  ; res_pro@cnLineThicknessF = 3.
  ; res_pro@cnFillOn        = True

  iplot_enino=gsn_csm_contour(wks,confi_enino,res_pro)
  iplot_lanina=gsn_csm_contour(wks,confi_lanina,res_pro)
  iplot_neutral=gsn_csm_contour(wks,confi_neutral,res_pro)

  opt                      = True
  opt@gsnShadeFillType     = "pattern"
  opt@gsnShadeMid     = 17         ; stipple
  ; opt@gsnShadeHigh   = "green"
  opt@gsnShadeFillScaleF   = 0.2         ; add extra density
  opt@gsnShadeFillDotSizeF = 0.0      ; make dots larger

  iplot_enino = gsn_contour_shade(iplot_enino, 0,0.05, opt)
  iplot_lanina = gsn_contour_shade(iplot_lanina, 0,0.05, opt)
  iplot_neutral = gsn_contour_shade(iplot_neutral, 0,0.05, opt)

 overlay(plot(0), iplot_enino)
 overlay(plot(1), iplot_lanina)
 overlay(plot(2), iplot_neutral)



pres=True
pres@gsnPanelFigureStrings=(/"a","b","c","d"/)
; pres@gsnPanelFigureStringsBackgroundFillColor="transparent"
gsn_panel(wks, plot, (/2,2/), pres)




















;;;画6张图的情况----------------------------------------------------
wks1=gsn_open_wks("pdf", "/mnt/f/share/era/era_sst_ano_enso")
res2=True
res2 = True
res2@gsnAddCyclic=False
res2@gsnDraw=False
res2@gsnFrame=False

; res2@cnFillPalette  = "temp_diff_18lev"
 ; res2@cnFillPalette="BlueDarkRed18"
 ; res2@cnFillPalette="BlAqGrWh2YeOrReVi22"
res2@cnFillPalette="NCV_blu_red"
res2@cnLineColors="Transparent"  ;让填充是渐变的还是分隔的，colorbar是渐变还是间隔
res2@cnFillOn = True 
res2@cnLinesOn = False
; res2@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
; res2@cnMinLevelValF       =   -1.0
; res2@cnMaxLevelValF       =   1.0
; res2@cnLevelSpacingF      =  0.2
res2@lbOrientation = "Vertical"

res2@mpLimitMode           = "LatLon"
    res@mpMinLatF             = min(sst_de&latitude)
    res@mpMaxLatF             = max(sst_de&latitude)
    res@mpMinLonF             = min(sst_de&longitude)
    res@mpMaxLonF             = max(sst_de&longitude)
    res@mpCenterLonF          =avg(sst_de&longitude)

res2@vpHeightF=0.4
res2@vpXF=0.1

plot1=new(6,graphic)
 
 res@gsnRightString="strong Ei Nino"
 plot1(0)= gsn_csm_contour_map(wks1, sst_str_enino_de_nor, res)
 res@gsnRightString="strong La Nina"
 plot1(1)= gsn_csm_contour_map(wks1, sst_str_lanina_de_nor, res)
 res@gsnRightString="strong Ei Nino-La Nina"
 plot1(2)= gsn_csm_contour_map(wks1, sst_sum_str_de_nor, res)
  res@gsnRightString="weak Ei Nino"
 plot1(3)= gsn_csm_contour_map(wks1, sst_weak_enino_de_nor, res)
res@gsnRightString="weak La Nina"
 plot1(4)= gsn_csm_contour_map(wks1, sst_weak_lanina_de_nor, res)

; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
; res@cnMinLevelValF       =   -0.5
; res@cnMaxLevelValF       =   0.5
; res@cnLevelSpacingF      =  0.2

; delete(res@cnLevelSelectionMode)
; delete(res@cnMaxLevelValF)
; delete(res@cnMinLevelValF)
; delete(res@cnLevelSpacingF)
 res@gsnRightString="weak Ei Nino-La Nina"
 plot1(5)= gsn_csm_contour_map(wks1, sst_sum_weak_de_nor, res)

pres1=True
pres1@gsnMaximize=True
pres1@gsnPanelLabelBar=True
pres1@lbTitleString    = "k/k"
pres1@lbTitlePosition  = "Right"              ; title position
pres1@lbTitleFontHeightF= .015                ; make title smaller
pres1@lbTitleDirection = "Across"             ; title direction
pres1@gsnPanelFigureStrings=(/"a","b","c","d","e","f"/)
gsn_panel(wks1, plot1, (/2,3/), pres1)



end 