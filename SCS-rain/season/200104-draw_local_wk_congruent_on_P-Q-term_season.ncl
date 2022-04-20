;;********************************************************
;;计算ERA的u,v,slp的场到P,  Q-term以及其动力项，还有观测的U,V,SLP的趋势
;;
;;lnshuheng 2020.1.4
;**********************************************************
function  rc_test(m_zonala_rc)
begin
  tval = onedtond(m_zonala_rc@tval , dimsizes(m_zonala_rc))
   df   = onedtond(m_zonala_rc@nptxy, dimsizes(m_zonala_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(m_zonala_rc, prob)
return prob
end

begin
;;;读入Q-term
    f3 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Q-term-5.nc", "r")
    series_Q = f3->series_Q(:33)
    series_thermaldyna_q = f3->series_thermaldyna_q(:33)
    series_dyna_q = f3->series_dyna_q(:33)
    series_quad_q = f3->series_quad_q(:33)
    series_eddy_q = f3->series_eddy_q(:33)
;;;

;;;读入GPCP降水场
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.197901.and.date.le.201212)
    p_t=f2->precip(timeind,:,:)

;;;

;==========================================================
;;;回归local wk and had 
;;==========================================================


   ;;;读入Hadim_zonal海温数据
    diri="/home/ys17-19/lsh/data/era/"
    fils1=systemfunc("ls "+diri+"era_m_zonal_meri_new.nc")
    f1=addfile(fils1, "r")
   ;;;


    ;;;region 
      lat1 = -40
      lat2 = 40
      lon1 = 0
      lon2 = 360
    ;;;

    ;;;time span
    startmon =197901
    endmon = 201212

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
    dim = yearnum
    ;;

    m_zonal=f1->m_zonal
  
    
    m_zonal_clm = dim_avg_n_Wrap(m_zonal(:,{500},:,:), 0)
    m_zonal_500=rmMonAnnCycTLL(m_zonal(:,{500},{-60:60},:))
   
;;;;season 

    season = "JJA"


    m_zonal_annual=month_to_season(m_zonal_500, season)
  
    dims=dimsizes(m_zonal_annual)
;;
   
   ;;;处理成年异常数据
    
     p_t = (/rmMonAnnCycTLL(p_t)/)
     p = dim_avg_n_Wrap(p_t(:,{0:15},{110:155}), (/1,2/))
     
     p_annual = month_to_season(p, season)
     p_annual = (/dim_standardize(p_annual, 1)/)
   

     ;;;计算趋势,计算线性一致性之前先把趋势计算了
     
     m_zonal_rc  = regCoef_n(ispan(1, yearnum, 1), m_zonal_annual, 0, 0)
     ; m_zonal_rc = m_zonal_rc*10

      m_zonal_rc = m_zonal_rc*yearnum*10^4

      copy_VarCoords(m_zonal_annual(1,:,:), m_zonal_rc)



    
     p_rc  = regCoef_n(ispan(1, yearnum, 1), p_annual, 0, 0)
     p_rc = p_rc*yearnum
  

    
     m_zonal_annual = (/dtrend_msg_n(ispan(1,dim(0),1),m_zonal_annual,False,True,0)/)
     p_annual = (/dtrend_msg_n(ispan(1,dim(0),1),p_annual,False,True,0)/)


   ;;;回归前指数的趋势以及处理
      series_Q = (/dim_standardize(series_Q, 1)/)
      series_Q_rc  = regCoef_n(ispan(1,dim(0),1), series_Q, 0, 0)
      series_Q_rc = series_Q_rc*yearnum
      
      series_Q = (/dtrend_msg_n(ispan(1,dim(0),1),series_Q,False,True,0)/)

      
      series_dyna_q = (/dim_standardize(series_dyna_q, 1)/)
      series_dyna_q_rc  = regCoef_n(ispan(1,dim(0),1), series_dyna_q, 0, 0)
      series_dyna_q_rc = series_dyna_q_rc*yearnum
    
  
      series_dyna_q = (/dtrend_msg_n(ispan(1,dim(0),1),series_dyna_q,False,True,0)/)
    
   
   ;;回归
    ;;1: 回归到P
    m_zonal_p_rc = regCoef_n(p_annual, m_zonal_annual, 0, 0)
    copy_VarCoords(m_zonal_annual(1,:,:), m_zonal_p_rc)
    
    ;;2: 回归到Q
    m_zonal_Q_rc = regCoef_n(series_Q, m_zonal_annual, 0, 0)
    copy_VarCoords(m_zonal_annual(1,:,:), m_zonal_Q_rc)
    
    m_zonal_dyna_rc = regCoef_n(series_dyna_q, m_zonal_annual, 0, 0)
    copy_VarCoords(m_zonal_annual(1,:,:), m_zonal_dyna_rc)
    
   ;;;


  
   ;;;显著性检验
    
    prob_m_zonal = rc_test(m_zonal_rc)
  
    prob_m_zonal_p = rc_test(m_zonal_p_rc)

    prob_m_zonal_Q = rc_test(m_zonal_Q_rc)
  
    prob_m_zonal_dyna = rc_test(m_zonal_dyna_rc)
  

   ;;;;
   

    

   
    ;;;
    ; m_zonal_rc = where(prob_m_zonal.le.0.1, m_zonal_rc, m_zonal_rc@_FillValue)
    ; m_zonal_p_rc = where(prob_m_zonal_p.le.0.1, m_zonal_p_rc, m_zonal_p_rc@_FillValue)
    ; m_zonal_Q_rc = where(prob_m_zonal_Q.le.0.1, m_zonal_Q_rc, m_zonal_Q_rc@_FillValue)
    ; m_zonal_dyna_rc = where(prob_m_zonal_dyna.le.0.1, m_zonal_dyna_rc, m_zonal_dyna_rc@_FillValue)
    

   ; u_rc = where(prob_u.le.0.1.or.prob_v.le.0.1, u_rc, u_rc@_FillValue)

  ;;;趋势一致性

    m_zonal_p_rc = (/m_zonal_p_rc*p_rc*10^4/)
    m_zonal_Q_rc = (/m_zonal_Q_rc*series_Q_rc*10^4/)
    m_zonal_dyna_rc = (/m_zonal_dyna_rc*series_dyna_q_rc*10^4/)
  ;;;;
    



;==============================================================
; Mission 4. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "m_zonal_congruent_on_P,Q-dynaterm-"+season)
  plot=new(4,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="BlueRed"

  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -6.0
  res@cnMaxLevelValF       =   6.0
  res@cnLevelSpacingF      =  1.0
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             = lon1
  res@mpMaxLonF             = lon2 

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0



  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0
  
 
  plot(0) = gsn_csm_contour_map(wks,m_zonal_rc ,res)
    ;  res@cnMinLevelValF       =   -0.5
    ; res@cnMaxLevelValF       =   0.5
    ; res@cnLevelSpacingF      =  0.1

  plot(1) = gsn_csm_contour_map(wks,m_zonal_p_rc ,res)

  plot(2) = gsn_csm_contour_map(wks,m_zonal_Q_rc ,res)

  plot(3) = gsn_csm_contour_map(wks,m_zonal_dyna_rc ,res)
  
  
   
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
   do i = 0,3
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
    res_pro@cnFillDotSizeF = 0.002
  
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(4,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_m_zonal,res_pro)
    iplot(1)=gsn_csm_contour(wks,prob_m_zonal_p,res_pro)
    iplot(2)=gsn_csm_contour(wks,prob_m_zonal_Q,res_pro)
    iplot(3)=gsn_csm_contour(wks,prob_m_zonal_dyna,res_pro)
    do i = 0,3
    overlay(plot(i), iplot(i))
    end do

    pres = True
    pres@gsnPanelLabelBar = False
    gsn_panel(wks, plot, (/2,2/), pres)
end
    
  







    
   