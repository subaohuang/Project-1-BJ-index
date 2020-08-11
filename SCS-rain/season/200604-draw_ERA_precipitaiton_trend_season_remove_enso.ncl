;;********************************************************
;;画GPCP的降水趋势
;;
;;lnshuheng 2019.12.24
;**********************************************************
; load "/home/ys17-19/lsh/Project/other/function_lin/function_lin.ncl"
; function cal_trend(data,season)
;   begin
;   dim = dimsizes(data)
;     p_rc=regCoef_n(ispan(1, dim(0), 1), data, 0, 0)
;     copy_VarCoords(data(1,:,:), p_rc)
;     p_rc=(/p_rc*dim(0)/)
;   return p_rc
; end 
function cal_trend(data)
  begin
    dim = dimsizes(data)
    p_mk=trend_manken(data, False, 0)
    p_rc = p_mk(1,:,:)
    copy_VarCoords(data(1,:,:), p_rc)
    p_rc=(/p_rc*dim(0)/)
    p_rc@prob = p_mk(0,:,:)
  return p_rc
end 
function remove_enso_signal(data,nino3_4)
  begin  
    data_de=(/dtrend_msg_n(ispan(1,dimsizes(data&time),1),data,False,True,0)/) 
    data_rc = regCoef_n(nino3_4,data_de, 0, 0)
    copy_VarCoords(data(1,:,:), data_rc)
    data_rc_3d  = conform_dims(dimsizes(data), data_rc, (/1,2/))
    m = reshape(data_rc@yintercept,dimsizes(data_rc_3d)) 
    nino3_4_3d = conform(data_rc_3d,nino3_4,(/0/))

    data_relate_enso  = nino3_4_3d*data_rc_3d + m 
    data_new = data 
    data_new = (/data - data_relate_enso/)
  return data_new
end 
begin
;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.197901)

    sst =  lonFlip(f1->sst(timeind,:,:))
 
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ssta_de = ssta
    ssta_de=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ssta_index = dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))
    printVarSummary(ssta_index)

    nino3_4 = month_to_season(ssta_index, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
    printVarSummary(nino3_4)
    
    ; nino3_4_new = nino3_4                ;;;;前期ENSO
    ; nino3_4_new(0:33) = (/nino3_4(1:)/)
    ; nino3_4_new(34) =ssta_index(35*12-1) 

    nino3_4 = (/dim_standardize_n(nino3_4, 1, 0)/)
;;;;;


    ;;;;;
    f2=addfile("/home/yangsong3/data-observation/ERA-interim/monthly/precipitation/ERAI-pre-1degree.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201512)
    lsp=f2->lsp(timeind,:,:)  ;;;units are m
    cp=f2->cp(timeind,:,:)   

    ;;convert units to mm/day
    p = cp
    p = (/lsp + cp/)
    p = p*10^3


    p_mam = month_to_season(p, "MAM")
    p_jja = month_to_season(p, "JJA")
    p_son = month_to_season(p, "SON")
    p_djf = month_to_season(p, "DJF")

    
    p_mam = (/remove_enso_signal(p_mam,nino3_4)/)
    p_jja = (/remove_enso_signal(p_jja,nino3_4)/)
    p_son  = (/remove_enso_signal(p_son,nino3_4)/)
    p_djf = (/remove_enso_signal(p_djf,nino3_4)/)
    
  
    ; rc_mam = cal_trend(p_mam,"MAM")
    ; rc_jja = cal_trend(p_jja,"JJA")
    ; rc_son = cal_trend(p_son,"SON")
    ; rc_djf = cal_trend(p_djf,"DJF")

    ; prob_mam=rc_test(rc_mam)
    ; prob_jja=rc_test(rc_jja)
    ; prob_son=rc_test(rc_son)
    ; prob_djf=rc_test(rc_djf)

    rc_mam = cal_trend(p_mam)
    rc_jja = cal_trend(p_jja)
    rc_son = cal_trend(p_son)
    rc_djf = cal_trend(p_djf)

    ; prob_mam=rc_test(rc_mam)
    ; prob_jja=rc_test(rc_jja)
    ; prob_son=rc_test(rc_son)
    ; prob_djf=rc_test(rc_djf)
    
    prob_mam = 1- rc_mam@prob
    prob_jja = 1- rc_jja@prob
    prob_son = 1- rc_son@prob
    prob_djf = 1- rc_djf@prob


    copy_VarCoords(p_mam(1,:,:), prob_mam)
    copy_VarCoords(p_mam(1,:,:), prob_jja)
    copy_VarCoords(p_mam(1,:,:), prob_son)
    copy_VarCoords(p_mam(1,:,:), prob_djf)
    
    
wks   = gsn_open_wks("pdf", "不同季节ERA降水的趋势-scs-1979-2013-去掉ENSO影响")
;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="NCV_blu_red"
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2

  res@mpLimitMode           = "LatLon"
  
  res@mpMinLatF             = 0
  res@mpMaxLatF             = 45
  res@mpMinLonF             =90
  res@mpMaxLonF             = 155



  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"


  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0

  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.7

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -2.5
  res@cnMaxLevelValF       =   2.5
  res@cnLevelSpacingF      =  0.25

  res@cnFillPalette="MPL_BrBG"
  res@gsnLeftString="Spring"

  plot = new(4,graphic)
  plot(0) = gsn_csm_contour_map(wks, rc_mam, res)
  res@gsnLeftString="Summer"
  plot(1) = gsn_csm_contour_map(wks, rc_jja, res)
  res@gsnLeftString="Autumn"
  plot(2) = gsn_csm_contour_map(wks, rc_son, res)
  res@gsnLeftString="Winter"
  plot(3) = gsn_csm_contour_map(wks, rc_djf, res)


  

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    
 
  ;;;添加赤道直线

   pgx=ispan(30,390,10)
   pgy=new(dimsizes(pgx),typeof(pgx))
   pgy=0

   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=1.0
   pgres@gsLineDashPattern=1
   
   str = ispan(1, 4, 1) + ""
   do i =0,3 
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
   end do 



  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

     str2 = ispan(21, 24, 1) +""
   do i =0,3 
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
  end do 


    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot=new(4,graphic)

    
 
    iplot(0) = gsn_csm_contour(wks,prob_mam,res_pro)    
    iplot(1) = gsn_csm_contour(wks,prob_jja,res_pro)
    iplot(2) = gsn_csm_contour(wks,prob_son,res_pro)
    iplot(3) = gsn_csm_contour(wks,prob_djf,res_pro)

    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

   
     do i=0,3
    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.2,opt)
     overlay(plot(i), iplot(i))
  end do

    ;;;  
   
pres = True
pres@gsnPanelLabelBar  = True
gsn_panel(wks, plot, (/2,2/), pres)

end