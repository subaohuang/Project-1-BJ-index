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
   f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")
    date=cd_calendar(f2->time, 1)
    timeind=ind(date.ge.197901.and.date.le.201312)
    p=f2->precip(timeind,:,:)
    
    p = (/rmMonAnnCycTLL(p)/)
    p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)

    p_annual = month_to_annual(p, 1)
    dim = dimsizes(p_annual)
    
    ;;;计算区域的降水平均，得到一个index,先看降水由哪里因素决定
     pre_index =dim_avg_n_Wrap(p_annual(:,{-5:5},{210:270}), (/1,2/))
    ;;;

    pre_index = (/dim_standardize(pre_index, 1)/)
    
    delete(date)
    delete(timeind)
;==========================================================
;;;回归u，v,slp，500hpa Geophyoptensial heigh
;;==========================================================
    ;;;read data
    path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    var_name = (/"shum","air"/)
    filein = path_wind + var_name +".mon.mean.nc"
        
    ; ;;;region 
    ;   lat1 = -10
    ;   lat2 = 90
    ;   lon1 = 0
    ;   lon2 = 360
    ; ;;;
        ;;;region 
      lat1 = -40
      lat2 = 40
      lon1 = 0
      lon2 = 360
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

    q = short2flt(f1->q(timeind,:,{lat1:lat2},{lon1:lon2}))

    f2 = addfile(filein(1),"r")
    t = short2flt(f2->t(timeind,:,{lat1:lat2},{lon1:lon2}))
    t = (/t-273.16/)

    level = q&level
   ;;;处理成年异常数据


     ;;;区域平均
     q_scs = dim_avg_n_Wrap(q(:,:,{-5:5},{210:270}), (/2,3/))
     t_scs = dim_avg_n_Wrap(t(:,:,{-5:5},{210:270}), (/2,3/))
     ;;;
   
     q_annual_scs = month_to_annual(q_scs, 1)
     t_annual_scs = month_to_annual(t_scs, 1)
   ;;;

     q_annual_scs = (/dtrend_msg_n(ispan(1,yearnum,1),q_annual_scs,False,True,0)/)
     t_annual_scs = (/dtrend_msg_n(ispan(1,yearnum,1),t_annual_scs,False,True,0)/)


;;;读入Q-term
    f3 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Q-term-5.nc", "r")
    series_Q = f3->series_Q(:)
    series_thermaldyna_q = f3->series_thermaldyna_q
    series_dyna_q = f3->series_dyna_q
    series_quad_q = f3->series_quad_q
    series_eddy_q = f3->series_eddy_q
;;;
   

   series_dyna_q = (/dim_standardize(series_dyna_q, 1)/)
   series_thermaldyna_q = (/dim_standardize(series_thermaldyna_q, 1)/)
   

;;;
 

   ;;;计算回归
    q_rc = regCoef_n(pre_index, q_annual_scs, 0, 0)
    t_rc = regCoef_n(pre_index, t_annual_scs, 0, 0)
    copy_VarCoords(q_annual_scs(1,:), q_rc)
    copy_VarCoords(t_annual_scs(1,:), t_rc)
   ;;;
   


   ;;;气候态
    q_clm_scs = dim_avg_n_Wrap(q_annual_scs, 0)
    t_clm_scs = dim_avg_n_Wrap(t_annual_scs, 0)
   ;;;
  
    q_rc = (/q_rc*yearnum/)
    q_change = q_clm_scs
    q_change = (/q_clm_scs + q_rc/)


    t_rc = (/t_rc*yearnum/)

    
    t_change = t_clm_scs
    t_change = (/t_clm_scs + t_rc/)




;;;=================================================
;;;画垂直坡面图
;;===================================================
 wks = gsn_open_wks("pdf","q和ts的profile回归到降水")                  ; send graphics to PNG file
    resxy = True
    resxy@gsnDraw = False
    resxy@gsnFrame  =False
    resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
    resxy@xyMonoLineThickness = True
    resxy@xyLineColors       =  (/"black","blue"/)          ; change line color
    resxy@xyMarkLineModes    =  (/"lines","lines"/)
    resxy@xyExplicitLegendLabels  =  (/"main","dyn"/)     ;加线的标注
    resxy@trYReverse        = True                       ; reverse Y-axis

    plot1 = gsn_csm_xy(wks, (/q_rc/),level ,resxy)
    resxy@xyLineThicknessF  =  (/ 1.0/)          ; make second line thicker
    plot2 = gsn_csm_xy(wks, (/t_rc/),level ,resxy)

    pres = True

    gsn_panel(wks, (/plot1,plot2/), (/2,1/), pres)
    
end 














