;;********************************************************
;;计算水汽通量算度，并且分解水汽通量散度
;; ERA-interim:u,v，sh:specific humidity
;; used ncl function: uv2dv_cfd,
;;grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2019.12.24
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

season = "MAM"
  ;;;read era data
        path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
        var_name = (/"uwnd","vwnd","shum.total"/)
        filein = path_wind + var_name +".mon.mean.nc"
        
        ;;;region 
          lat1 = -10
          lat2 = 30
          lon1 = 100
          lon2 = 180
        ;;;
        ;   ;;;region 
        ;   lat1 = -10
        ;   lat2 = 10
        ;   lon1 = 100
        ;   lon2 = 160
        ; ;;;

        ;;;time span
        startmon =197901
        endmon = 201512

        startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
        endyear_r = str_split_by_length(endmon, 4 )

        startyear = stringtoint(startyear_r(0))
        endyear = stringtoint(endyear_r(0))

        yearnum = endyear - startyear +1
        ;;


        f1 = addfile(filein(0),"r")

        date=cd_calendar(f1->time, 1)
        timeind=ind(date.ge.startmon.and.date.le.endmon)

        u = short2flt(f1->u(timeind,:,{lat1:lat2},{lon1:lon2}))

        f2 = addfile(filein(1),"r")
        v = short2flt(f2->v(timeind,:,{lat1:lat2},{lon1:lon2}))

        
        f3 = addfile(filein(2),"r")
        date :=cd_calendar(f3->time, 1)
        timeind :=ind(date.ge.startmon.and.date.le.endmon)
        q = short2flt(f3->q(timeind,:,{lat1:lat2},{lon1:lon2}))        
   ;;;
   print("=====finish read data==========")
;==============================================================
; Mission 1. calculate Q
; ==============================================================
    era_qu  = q
    era_qv  = q
    era_qu  = u*q
    printVarSummary(era_qu)
    era_qv  = v*q
    era_wgt_qu = era_qu(:,0,:,:)
    era_wgt_qv = era_qv(:,0,:,:)


   ;;;------Cal DivQ-------
    ;-----------Follow is the vertical integral-----------
    pver=era_qu&level
    pver=pver*100
    pver@units="Pa"


    dp = pver              ; caculate dp
    dp(0:dimsizes(dp)-2) = abs(pver(1:) - pver(0:dimsizes(pver)-2))
    dp(dimsizes(dp)-1) = dp(dimsizes(dp)-2)
    ; dp=dpres_plevel(pver,100000,min(pver),0)              ; caculate dp

    dp_cfm=conform_dims(dimsizes(era_qu),dp,1)          ; confirm dp to caculate sum(var*dp)



    era_wgt_qu=(/dim_sum_n(dp_cfm*era_qu,1)/)               ; sum(var*dp)
    era_wgt_qv=(/dim_sum_n(dp_cfm*era_qv,1)/)               ; sum(var*dp)

    era_wgt_qu = (/era_wgt_qu/9.8/)
    era_wgt_qv = (/era_wgt_qv/9.8/)

    g_lon   =   era_wgt_qu&longitude
    g_lat   =   era_wgt_qu&latitude
    g_nlat  =   dimsizes(g_lat)

    DivQ = uv2dv_cfd(era_wgt_qu,era_wgt_qv,g_lat,g_lon,2)
    ; DivQ = uv2dvF_Wrap(era_wgt_qu(:,::-1,:),era_wgt_qv(:,::-1,:)) 
    DivQ   =   (/-DivQ*86400/) 
    copy_VarCoords(era_wgt_qu, DivQ)
    printVarSummary(DivQ)
    ; print(DivQ(0,:,:))
    print("=====finish DivQ calculate==========")
    
    era_divQ_bar = month_to_season(DivQ, season)
    era_divQ_delta = dim_rmvmean_n_Wrap(era_divQ_bar, 0)
    
    printVarSummary(era_divQ_delta)



;==============================================================
; Mission 2. calculate dynamic and thermaldynamic term of Q
; ==============================================================
   
   ;;;prepare variable

   
   q_bar = month_to_season(q, season)
   u_bar = month_to_season(u, season)
   v_bar = month_to_season(v, season)

   q_delta = dim_rmvmean_n_Wrap(q_bar, 0)
   u_delta = dim_rmvmean_n_Wrap(u_bar, 0)
   v_delta = dim_rmvmean_n_Wrap(v_bar, 0)
   
   q_clm = dim_avg_n_Wrap(q_bar, 0)
   u_clm = dim_avg_n_Wrap(u_bar, 0)   
   v_clm = dim_avg_n_Wrap(v_bar, 0)
  
   q_clm_cfm = conform_dims(dimsizes(q_bar),q_clm,(/1,2,3/))
   u_clm_cfm = conform_dims(dimsizes(u_bar),u_clm,(/1,2,3/))          
   v_clm_cfm = conform_dims(dimsizes(v_bar),v_clm,(/1,2,3/))

   
   dim = dimsizes(q_bar(1,:,:,:))
   q_mon = q
   u_mon = u
   v_mon = v
   do i = 0,yearnum-1,1
   q_mon(12*i:12*i+11,:,:,:) = conform_dims((/12,dim(0),dim(1),dim(2)/), q_bar(i,:,:,:), (/1,2,3/))
   u_mon(12*i:12*i+11,:,:,:) = conform_dims((/12,dim(0),dim(1),dim(2)/), u_bar(i,:,:,:), (/1,2,3/))
   v_mon(12*i:12*i+11,:,:,:) = conform_dims((/12,dim(0),dim(1),dim(2)/), v_bar(i,:,:,:), (/1,2,3/))
   end do 
   
    q_eddy = q 
    q_eddy = (/q - q_mon/)     ;;;;eddy term
    u_eddy = u 
    u_eddy = (/u - u_mon/)     ;;;;eddy term
    v_eddy = v 
    v_eddy = (/v - v_mon/)     ;;;;eddy term

   
    
    uq_eddy = q_eddy*u_eddy
    copy_VarCoords(q_eddy, uq_eddy)

   
    
    vq_eddy = q_eddy*v_eddy
    copy_VarCoords(q_eddy, vq_eddy)
    
    

    uq_eddy_bar = month_to_season(uq_eddy, season)
    vq_eddy_bar = month_to_season(vq_eddy, season)
    


    uq_eddy_delta = dim_rmvmean_n_Wrap(uq_eddy_bar, 0)
    vq_eddy_delta = dim_rmvmean_n_Wrap(vq_eddy_bar, 0)

    print("=====end variable set==========")
   


  
  ;;;;;
    print("=====now start componte calculate==========")
   ;;; Calculate the MFC dynamic term
   ;;
   ;; MFC_dynamic_adv = δvd(q)  
   
      dq_clm = grad_latlon_cfd(q_clm,g_lat,g_lon,False,False)
      era_dynamic_adv_q = u_delta*conform_dims(dimsizes(u_delta), dq_clm[1], (/1,2,3/)) \
      + v_delta*conform_dims(dimsizes(v_delta), dq_clm[0], (/1,2,3/))

      dp_cfm_new = conform_dims(dimsizes(u_bar),dp,1)          ; confirm dp to caculate sum(var*dp)
      era_wgt_dyna_adv_q = (/dim_sum_n(dp_cfm_new*era_dynamic_adv_q,1)/)               ; sum(var*dp)
      era_wgt_dyna_adv_q = (/era_wgt_dyna_adv_q*86400/(-9.8)/)       ;;convert to day
      copy_VarCoords(q_bar(:,1,:,:), era_wgt_dyna_adv_q)
    ;;

    ;; MFC_dynamic_con = qd(δv)

      duv_delta = uv2dv_cfd(u_delta,v_delta,g_lat,g_lon,2)
      era_dynamic_con_q = duv_delta*conform_dims(dimsizes(duv_delta), q_clm, (/1,2,3/))
      
      era_wgt_dyna_con_q = (/dim_sum_n(dp_cfm_new*era_dynamic_con_q,1)/)               ; sum(var*dp)
      era_wgt_dyna_con_q = (/era_wgt_dyna_con_q*86400/(-9.8)/)       ;;convert to day
      copy_VarCoords(q_bar(:,1,:,:), era_wgt_dyna_con_q)
    ;;;;

    print("=====end dynamic componte calculate==========")

    ; Calculate the MFC thermaldynamic term
    ;   MFC_thermaldynamic_adv = vd(δq)
       
      dq_delta = grad_latlon_cfd(q_delta,g_lat,g_lon,False,False)
      era_thermaldynamic_adv_q = u_clm_cfm*dq_delta[1] + v_clm_cfm*dq_delta[0]

      era_wgt_thermaldyna_adv_q = (/dim_sum_n(dp_cfm_new*era_thermaldynamic_adv_q,1)/)               ; sum(var*dp)
      era_wgt_thermaldyna_adv_q = (/era_wgt_thermaldyna_adv_q*86400/(-9.8)/)       ;;convert to day
      copy_VarCoords(q_bar(:,1,:,:), era_wgt_thermaldyna_adv_q)
    ;;;

    ;;  MFC_thermaldynamic_con = qd(δv)

      duv_clm = uv2dv_cfd(u_clm,v_clm,g_lat,g_lon,2)
      era_thermaldynamic_con_q = q_delta*conform_dims(dimsizes(q_delta), duv_clm, (/1,2,3/))

      era_wgt_thermaldyna_con_q = (/dim_sum_n(dp_cfm_new*era_thermaldynamic_con_q,1)/)               ; sum(var*dp)
      era_wgt_thermaldyna_con_q = (/era_wgt_thermaldyna_con_q*86400/(-9.8)/)       ;;convert to day
      copy_VarCoords(q_bar(:,1,:,:), era_wgt_thermaldyna_con_q)
      print("=====end thermaldynamic componte calculate==========")
    ;;;


  ;;;quadratic adv term δv*d(δq)
      era_quad_adv_q = u_delta*dq_delta[1] + v_delta*dq_delta[0]    

      era_wgt_quad_adv_q = (/dim_sum_n(dp_cfm_new*era_quad_adv_q,1)/)               ; sum(var*dp)
      era_wgt_quad_adv_q = (/era_wgt_quad_adv_q*86400/(-9.8)/)       ;;convert to day
      copy_VarCoords(q_bar(:,1,:,:), era_wgt_quad_adv_q)
  ;;;

  ;;;;;;quadratic con term δq*d(δv)
      era_quad_con_q = q_delta*uv2dv_cfd(u_delta,v_delta,g_lat,g_lon,2)

      era_wgt_quad_con_q = (/dim_sum_n(dp_cfm_new*era_quad_con_q,1)/)               ; sum(var*dp)
      era_wgt_quad_con_q = (/era_wgt_quad_con_q*86400/(-9.8)/)       ;;convert to day
      copy_VarCoords(q_bar(:,1,:,:), era_wgt_quad_con_q)
     print("=====end quadratic componte calculate==========")
  ;;;


  ; ;;;quandratic climate term d((δv*δq)clm)

  ;   era_quad_clm_q = uv2dv_cfd(era_quad_qu_clm ,era_quad_qv_clm,g_lat,g_lon,2)

  ;   dp_cfm_new2 = conform_dims(dimsizes(u_bar(1,:,:,:)),dp,0)          ; confirm dp to caculate sum(var*dp)
  ;   era_wgt_quad_q_clm = (/dim_sum_n(dp_cfm_new2*era_quad_clm_q,0)/)               ; sum(var*dp)
  ;   era_wgt_quad_q_clm = (/era_wgt_quad_q_clm*86400/(-9.8)/)       ;;convert to day
  ;   copy_VarCoords(q_bar(1,1,:,:), era_wgt_quad_q_clm)
  ;   printVarSummary(era_wgt_quad_q_clm)

  ;   print("=====end quandratic clmate componte calculate==========")

  ;;eddy term δd((bar(v'q'))

    ; era_eddy_q_annual = uv2dv_cfd(uq_eddy_bar,vq_eddy_bar,g_lat,g_lon,2)
    ; era_eddy_q = dim_rmvmean_n_Wrap(era_eddy_q_annual, 0)  
    era_eddy_q = uv2dv_cfd(uq_eddy_delta,vq_eddy_delta,g_lat,g_lon,2)

    era_wgt_eddy_q = (/dim_sum_n(dp_cfm_new*era_eddy_q,1)/)               ; sum(var*dp)
    era_wgt_eddy_q = (/era_wgt_eddy_q*86400/(-9.8)/)       ;;convert to day
    copy_VarCoords(q_bar(:,1,:,:), era_wgt_eddy_q)

    print("=====end eddy componte calculate==========")

  ;;;

; ;==============================================================
; ; Mission 3.  area avg(x)
; ; ==============================================================
   ;;era_wgt_dyna_q ; era_wgt_thermaldyna_q ;;era_wgt_quad_q;;era_wgt_eddy_q 
  ;  lat1 = 0
   ;lat2 = 15
   ;lon1 = 110
   ;lon2 = 155
    series_era_divQ_delta  = dim_avg_n_Wrap(era_divQ_delta(:,{0:20},{125:160}), (/1,2/))
    
    series_dyna_adv_q = dim_avg_n_Wrap(era_wgt_dyna_adv_q(:,{0:20},{125:160}),(/1,2/))
    series_dyna_con_q = dim_avg_n_Wrap(era_wgt_dyna_con_q(:,{0:20},{125:160}),(/1,2/))
    
    series_thermaldyna_adv_q = dim_avg_n_Wrap(era_wgt_thermaldyna_adv_q(:,{0:20},{125:160}),(/1,2/))
    series_thermaldyna_con_q = dim_avg_n_Wrap(era_wgt_thermaldyna_con_q(:,{0:20},{125:160}),(/1,2/))
    
    series_quad_adv_q = dim_avg_n_Wrap(era_wgt_quad_adv_q(:,{0:20},{125:160}),(/1,2/))
    series_quad_con_q = dim_avg_n_Wrap(era_wgt_quad_con_q(:,{0:20},{125:160}),(/1,2/))

    series_eddy_q = dim_avg_n_Wrap(era_wgt_eddy_q(:,{0:20},{125:160}), (/1,2/))

    
    ; series_eddy_q = (/series_thermaldyna_q + series_dyna_q + series_quad_q + series_eddy_q - avg(dim_avg_n_Wrap(era_wgt_quad_q_clm({0:20},{125:160}),1))/)


    ;;;写成nc 文件输出

    outpath = "/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/review/new/Q-term-8.nc"
       ;;-----------------------写入nc文件---------------------------
    system("rm " + outpath)   ; remove any pre-existing file
    ncdf = addfile(outpath ,"c")  ; open output netCDF file

    filedimdef(ncdf,"time",-1,True)  ; make time an UNLIMITED dimension; recommended  for most applications
   ;===================================================================
    ; create global attributes of the file (optional)
    ;===================================================================
    fAtt               = True            ; assign file attributes
    fAtt@title         = "Q-term,热力平流，热力幅合，动力平流，动力幅合，平方以及eddy项"
    fAtt@Conventions   = "CHAZHI"
    fAtt@creation_date = systemfunc ("date")
    fileattdef( ncdf, fAtt )            ; copy file attributes==
    ;=================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
    ncdf->series_Q  =  series_era_divQ_delta
    ncdf->series_dyna_adv_q = series_dyna_adv_q
    ncdf->series_dyna_con_q = series_dyna_con_q
    ncdf->series_thermaldyna_adv_q = series_thermaldyna_adv_q
    ncdf->series_thermaldyna_con_q = series_thermaldyna_con_q
    ncdf->series_quad_adv_q = series_quad_adv_q
    ncdf->series_quad_con_q = series_quad_con_q
    ncdf->series_eddy_q = series_eddy_q

    ncdf->time = ispan(1979,2015,1)
    ;;;

; ;==============================================================
; ; Mission 4.  trend
; ; ==============================================================




    q_delta_rc = regline(ispan(1,yearnum,1), series_era_divQ_delta)
    q_delta_rc = (/q_delta_rc*yearnum/)

    print("q_delta_rc="+q_delta_rc)

    q_dyna_adv_rc = regline(ispan(1,yearnum,1), series_dyna_adv_q)
    q_dyna_adv_rc = (/q_dyna_adv_rc*yearnum/)
    
    print("q_dyna_adv_rc="+q_dyna_adv_rc)

    q_dyna_con_rc = regline(ispan(1,yearnum,1), series_dyna_con_q)
    q_dyna_con_rc = (/q_dyna_con_rc*yearnum/)
    
    print("q_dyna_con_rc="+q_dyna_con_rc)

    q_thermaldyna_adv_rc = regline(ispan(1,yearnum,1), series_thermaldyna_adv_q)
    q_thermaldyna_adv_rc = (/q_thermaldyna_adv_rc*yearnum/)

    print("q_thermaldyna_adv_rc="+q_thermaldyna_adv_rc)
    
    q_thermaldyna_con_rc = regline(ispan(1,yearnum,1), series_thermaldyna_con_q)
    q_thermaldyna_con_rc = (/q_thermaldyna_con_rc*yearnum/)
    
    print("q_thermaldyna_con_rc= "+ q_thermaldyna_con_rc)
    
    q_quad_adv_rc = regline(ispan(1,yearnum,1), series_quad_adv_q)
    q_quad_adv_rc = (/q_quad_adv_rc*yearnum/)
    
    q_quad_con_rc = regline(ispan(1,yearnum,1), series_quad_con_q)
    q_quad_con_rc = (/q_quad_con_rc*yearnum/)

    
    q_eddy_rc = regline(ispan(1,yearnum,1), series_eddy_q)
    q_eddy_rc = (/q_eddy_rc*yearnum/)



;==============================================================
; Mission 5. bar plot
; ==============================================================
  
  x = (/1,2,3,4,5,6,7,8/)
  y = (/q_delta_rc,q_dyna_adv_rc,q_dyna_con_rc,q_thermaldyna_adv_rc,\
       q_thermaldyna_con_rc,q_quad_adv_rc,q_quad_con_rc,q_eddy_rc/)

  labels = (/"TOTAL","ADV(DY)","CON(DY)","ADV(TH)","CON(TH)",\
  "ADV(QD)","CON(QD)","EDDY"/)
;********************************************************
; create the plot
;********************************************************   
  wks1 = gsn_open_wks("pdf","bar-q-component-trend")                  ; send graphics to PNG file
  res                       = True                 ; plot mods desired         
  res@gsnFrame              = False                ; don't advance frame yet
  res@gsnXYBarChart         = True                 ; turn on bar chart
  res@gsnXYBarChartBarWidth = 0.75                 ; change bar widths
  ; res@gsnXYBarChartColors = (/"firebrick","red","orange","green", \
  ; "navy","blue","SkyBlue","SlateBlue"/)            ; choose colors
  
  res@tmXBOn                = True                ; turn off tickmarks at bot
  res@trYMinF               = -0.4                    ; bring bars down to zero
  res@trYMaxF               = 1.2                   ; bring bars down to zero

  res@trXMinF               = 0                    ; adds space on either end
  res@trXMaxF               = 9                    ; of the 1st and last bars

  res@gsnYRefLine = 0.      ; set the reference line equal to 0
  res@gsnAboveYRefLineColor = "Gray" 
  res@gsnBelowYRefLineColor = "white"
  
  ; res@tiMainString          = "Estimated Cancer Deaths for 2002"
;---Put city labels on Y axis
  res@tmXBMode                 = "Explicit"
  res@tmXBValues               = x
  res@tmXBLabels               = labels
  res@tmXBLabelFontHeightF = 0.01

  plot = gsn_csm_xy (wks1,x,y,res)                  ; create plot
  
 draw(plot)
 frame(wks1)

 ;==============================================================
; Mission 5. draw timeseries
; ==============================================================
  wks2   = gsn_open_wks("pdf", "Q的动力项和热力项的时间序列-precise")
  plot2 = new(1,graphic)
  
  resxy = True
  resxy@gsnDraw = False
  resxy@gsnFrame = False
    
  resxy@tmXBMode  =  "Explicit"
  resxy@trXMaxF = endyear
  resxy@trXMinF = startyear
  resxy@tmXBValues = ispan(startyear,endyear,5)
  resxy@tmXBLabels = ispan(startyear,endyear,5)
  resxy@tmXBMinorOn = True
  resxy@tmXBMinorValues = ispan(startyear,endyear,1)

  resxy@tmYLMinorOn = False

  resxy@vpHeightF = 0.4
  resxy@vpWidthF = 0.6
  ; resxy@vpXF = 0.15
  ; resxy@vpYF = 0.95

  resxy@xyLineThicknessF  =  (/ 3.0/)          ; make second line thicker
  resxy@xyMonoLineThickness = True
  resxy@xyLineColors       =  (/"black","blue","red","green","gray","pink","purple","orange"/)          ; change line color
  resxy@xyMarkLineModes    =  (/"lines"/)
  resxy@xyMonoMarkLineMode = True
  resxy@xyExplicitLegendLabels  =  (/"main","dyn","thermaldynamic","quandratic", "eddy"/)     ;加线的标注
  

  time = ispan(startyear, endyear, 1)
  plot2 = gsn_csm_xy(wks2,time,(/series_era_divQ_delta,series_dyna_adv_q,\
                              series_dyna_con_q,series_thermaldyna_adv_q,series_thermaldyna_con_q,series_quad_adv_q,series_quad_con_q,\
                          series_eddy_q/), resxy)
  draw(plot2)
  frame(wks2)


end 
    
  







