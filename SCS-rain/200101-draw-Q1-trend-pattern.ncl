;;********************************************************
;;draw diabatic heatspacial pattern 
;;q1 derived from NCEP_DOE calculate by eqution
;;                               lnshuheng 2020.1.1 
;**********************************************************
function cal_trend(x,y)
begin
  dim = x 
  q1_annual = y 
  q1_annual_rc = regCoef_n(ispan(1,dim,1),q1_annual,0,0)
  copy_VarCoords(q1_annual(1,:,:), q1_annual_rc)  
  q1_annual_rc = q1_annual_rc*dim(0)
return q1_annual_rc
end 
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)
return prob
end
begin
; ==============================================================
; Mission 1. read data
; ==============================================================
;;;read Q1 monthly data
  fils = systemfunc("ls /home/yangsong3/data-observation/NCEP_DOE/monthly/q1/q1."+"*"+".mon.mean.nc")
  f3 = addfiles(fils, "r")
  date = cd_calendar(f3[:]->time, 1)

  timeind = ind(date.le.201312.and.date.ge.197901)

  q1 = f3[:]->q1(timeind,:,:,:)

  printVarSummary(q1)
  pver  = q1&lev
 
  dp=dpres_plevel(pver,1013,min(pver),0)              ; caculate dp
  dp_cfm=conform_dims(dimsizes(q1),dp,1)          ; confirm dp to caculate sum(var*dp)
  wgt_q1 = q1(:,1,:,:)
  wgt_q1=(/dim_sum_n(dp_cfm*q1,1)/)               ; sum(var*dp)

  q1_annual = month_to_annual(wgt_q1, 1)
  q1_djf = month_to_season(wgt_q1, "DJF")
  q1_mam = month_to_season(wgt_q1, "MAM")
  q1_son = month_to_season(wgt_q1, "SON")
  q1_jja = month_to_season(wgt_q1, "JJA")
 

  dim = dimsizes(q1_annual)



  q1_djf_rc = cal_trend(dim(0),q1_djf)
  q1_mam_rc = cal_trend(dim(0),q1_mam)
  q1_son_rc = cal_trend(dim(0),q1_son)
  q1_jja_rc = cal_trend(dim(0),q1_jja)    

  prob_djf = rc_test(q1_djf_rc)
  prob_mam = rc_test(q1_mam_rc)
  prob_jja = rc_test(q1_jja_rc)
  prob_son = rc_test(q1_son_rc)
  
;;**************************************************************
;画图
;;***************************************************************
    wks   = gsn_open_wks("pdf", "Q1的趋势平面图")
    plot=new(4,graphic)
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnLeftStringFontHeightF=20.0
    res@cnFillPalette="BlueWhiteOrangeRed"
    res@cnFillOn=True
    res@cnLinesOn=False

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=True
    res@lbLabelStride=2
    res@lbOrientation="Vertical"
    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -5
    res@mpMaxLatF             = 45
    res@mpMinLonF             =100
    res@mpMaxLonF             = 160
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

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   -0.04
    res@cnMaxLevelValF       =   0.04
    res@cnLevelSpacingF      =  0.005


    plot(0) = gsn_csm_contour_map(wks, q1_djf_rc, res)
   
    plot(1) = gsn_csm_contour_map(wks, q1_mam_rc, res)
    
    plot(2) = gsn_csm_contour_map(wks, q1_jja_rc, res)

    plot(3) = gsn_csm_contour_map(wks, q1_son_rc, res)



  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0

   str = (/"1","2","3","4"/)
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   do i=0,3
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   end do

   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False

    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.0045
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(4,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_djf,res_pro)
    iplot(1)=gsn_csm_contour(wks,prob_mam,res_pro)
    iplot(2)=gsn_csm_contour(wks,prob_jja,res_pro)
    iplot(3)=gsn_csm_contour(wks,prob_son,res_pro)
    

    do i=0,3
    overlay(plot(i), iplot(i))
    end do

 pres = True
 gsn_panel(wks, plot, (/2,2/), pres)


end
  
  ;;;


  

;;;;

  