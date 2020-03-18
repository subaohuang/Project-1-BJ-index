;;********************************************************
;;回归ERA的u,v,slp,z到水汽收支方程的不同的变化上
;;
;;lnshuheng 2019.12.24
;**********************************************************
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
;;;读入Q-term
    f3 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Q-term-5.nc", "r")
    series_Q = f3->series_Q(:)
    series_thermaldyna_q = f3->series_thermaldyna_q
    series_dyna_q = f3->series_dyna_q
    series_quad_q = f3->series_quad_q
    series_eddy_q = f3->series_eddy_q
;;;

;==========================================================
;;;回归u，v,slp，500hpa Geophyoptensial heigh
;;==========================================================
    ;;;read data
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"wwnd","air"/)
    filein = path_wind + var_name +".mon.mean.nc"
        
    ; ;;;region 
    ;   lat1 = -10
    ;   lat2 = 90
    ;   lon1 = 0
    ;   lon2 = 360
    ; ;;;
    ; ;;;region 
    ;   lat1 = 0
    ;   lat2 = 15
    ;   lon1 = 110
    ;   lon2 = 155
    ; ;;;
    
    ;;;region 
      lat1 = -5
      lat2 = 45
      lon1 = 100
      lon2 = 160
    ;;;


    ;;;time span
    startmon =197901
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    dim = yearnum
    ;;


    f1 = addfile(filein(0),"r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.startmon.and.date.le.endmon)

    w = short2flt(f1->w(timeind,:,{lat1:lat2},{lon1:lon2}))
    w = (/w*-1/)

    f2 = addfile(filein(1),"r")
    t = short2flt(f2->t(timeind,:,{lat1:lat2},{lon1:lon2}))
    

    level = w&level

    ;;将温度处理成位势温度
     pt = pot_temp(level*100, t, 1, False)
     copy_VarCoords(t, pt)
     
  
    p=level*100
    dtdp = pt
    dtdp = center_finite_diff_n(pt, p,False , 0, 1)  
    dtdp = dtdp*-1
   ;;;处理成年异常数据


     ; q_annual  = month_to_annual(q, 1)
     ; t_annual  = month_to_annual(t, 1)
     MONTH = "JJA"

     ; w_annual = month_to_season(w, MONTH)
     ; t_annual = month_to_season(t, MONTH)
     ; dtdp_annual = month_to_season(dtdp, MONTH)
     
      w_annual = month_to_annual(w, 1) 
      t_annual = month_to_annual(pt, 1)
      dtdp_annual = month_to_annual(dtdp, 1)


     t_annual = dtdp_annual*100

    ;;;垂直积分
    dp=dpres_plevel(level,1013,min(level),0)              ; caculate dp
    dp_cfm=conform_dims(dimsizes(w_annual),dp,1)          ; confirm dp to caculate sum(var*dp)
    wgt_w = w_annual(:,1,:,:)
    wgt_w=(/dim_sum_n(dp_cfm*w_annual,1)/)               ; sum(var*dp)

    wgt_w = w_annual(:,{500},:,:)

    wgt_t = t_annual(:,1,:,:)
    wgt_t=(/dim_sum_n(dp_cfm*t_annual,1)/)               ; sum(var*dp)

    ;;;




     ;;;区域平均
     w_annual_scs  = dim_avg_n_Wrap(w_annual(:,:,{0:15},{110:155}), (/2,3/))
     t_annual_scs = dim_avg_n_Wrap(t_annual(:,:,{0:15},{110:155}), (/2,3/))
     ;;;



   ;;;计算趋势
    w_rc = regCoef_n(ispan(1,yearnum,1), w_annual_scs, 0, 0)
    t_rc = regCoef_n(ispan(1,yearnum,1), t_annual_scs, 0, 0)
    copy_VarCoords(w_annual_scs(1,:), w_rc)
    copy_VarCoords(t_annual_scs(1,:), t_rc)

    wgt_t_rc = regCoef_n(ispan(1,yearnum,1), wgt_t, 0, 0)
    wgt_w_rc = regCoef_n(ispan(1,yearnum,1), wgt_w, 0, 0)
    copy_VarCoords(wgt_t(1,:,:), wgt_t_rc)
    copy_VarCoords(wgt_w(1,:,:), wgt_w_rc)

    prob_t = rc_test(wgt_t_rc)
    prob_w = rc_test(wgt_w_rc)

    ; prob_t  = prob_t*10^6
    


   ;;;
   


   ;;;气候态
    w_clm_scs = dim_avg_n_Wrap(w_annual_scs, 0)
    t_clm_scs = dim_avg_n_Wrap(t_annual_scs, 0)
   ;;;
  
    w_rc = (/w_rc*yearnum/)
    w_change = w_clm_scs
    w_change = (/w_clm_scs + w_rc/)


    t_rc = (/t_rc*yearnum/)

    
    t_change = t_clm_scs
    t_change = (/t_clm_scs + t_rc/)




;;;=================================================
;;;画垂直坡面图
;;===================================================
 wks = gsn_open_wks("pdf","ERA-In的w和位温的气候态和趋势的profile"+MONTH)                  ; send graphics to PNG file
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       =  (/"black","blue"/)          ; change line color
    resxy@xyMarkLineModes    =  (/"lines","lines"/)
    resxy@xyExplicitLegendLabels  =  (/"main","dyn"/)     ;加线的标注
    resxy@trYReverse        = True               
    resxy@trYMinF = 100        ; reverse Y-axis

    plot1 = gsn_csm_xy(wks, (/w_clm_scs,w_change/),level ,resxy)
    resxy@xyLineThicknessF  =  (/ 1.0/)          ; make second line thicker
    delete(resxy@trYMinF)
    
    ; resxy@trXMinF = 0
    ; resxy@trXMaxF = 20
    plot2 = gsn_csm_xy(wks, (/t_clm_scs,t_change/),level ,resxy)

    pres = True

    gsn_panel(wks, (/plot1,plot2/), (/2,1/), pres)
 
;;;=================================================
;;;平面分布图
;;===================================================
 wks2   = gsn_open_wks("pdf", "ERA-In中的W和静力稳定度的趋势平面图-"+MONTH)
    plot=new(2,graphic)
    res=True
    res@gsnFrame=False
    res@gsnDraw=False
    res@gsnAddCyclic = False
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
    res@cnMinLevelValF       =   -0.001
    res@cnMaxLevelValF       =   0.001
    res@cnLevelSpacingF      =  0.0001


    plot(0) = gsn_csm_contour_map(wks2, wgt_w_rc, res)

    res@cnMinLevelValF       =   -0.6
    res@cnMaxLevelValF       =   0.6
    res@cnLevelSpacingF      =  0.1

   
    plot(1) = gsn_csm_contour_map(wks2, wgt_t_rc, res)

   

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
   do i=0,1
   plot@$str(i)$ = gsn_add_polyline(wks2, plot(i), pgx_pool, pgy_pool, pgres)
   end do

   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0.0,0.1/)
    res_pro@lbLabelBarOn = False

    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.0045
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(2,graphic)

    iplot(0)=gsn_csm_contour(wks2,prob_w,res_pro)
    iplot(1)=gsn_csm_contour(wks2,prob_t,res_pro)

    

    do i=0,1
    overlay(plot(i), iplot(i))
    end do

 pres = True
 gsn_panel(wks2, plot, (/2,1/), pres)

    
end 














