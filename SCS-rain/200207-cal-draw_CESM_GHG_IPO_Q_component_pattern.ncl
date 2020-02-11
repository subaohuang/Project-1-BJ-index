;;********************************************************
;;CESM模式试验
;;计算水汽通量算度，并且分解水汽通量散度
;; IPO敏感性试验:U,V，Q:specific humidity 
;; used ncl function: uv2dv_cfd,
;;grad_latlon_cfd:lon must be in increasing order
;;linshuheng 2020.02.07
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

function cal_trend(yearnum,era_divQ_delta,era_wgt_qu_delta,era_wgt_qv_delta)
begin
  ;;;;Q component
    Q_rc  = regCoef_n(ispan(1, yearnum, 1), era_divQ_delta, 0, 0)
    copy_VarCoords(era_divQ_delta(1,:,:), Q_rc)
    Q_rc = (/Q_rc*yearnum/)
    
    era_wgt_qu_delta_rc  = regCoef_n(ispan(1, yearnum, 1), era_wgt_qu_delta, 0, 0)
    copy_VarCoords(era_divQ_delta(1,:,:), era_wgt_qu_delta_rc)
    era_wgt_qu_delta_rc = (/era_wgt_qu_delta_rc*yearnum/)

    era_wgt_qv_delta_rc  = regCoef_n(ispan(1, yearnum, 1), era_wgt_qv_delta, 0, 0)
    copy_VarCoords(era_divQ_delta(1,:,:), era_wgt_qv_delta_rc)
    era_wgt_qv_delta_rc = (/era_wgt_qv_delta_rc*yearnum/)
 return [/Q_rc,era_wgt_qu_delta_rc,era_wgt_qv_delta_rc/]
 end

function cal_divQ(data_qu,data_qv)
begin 
    ;;;------Cal DivQ-------
    ;-----------Follow is the vertical integral-----------
    pver=data_qu&lev

    dp=dpres_plevel(pver,100000,min(pver),0)              ; caculate dp
    dp_cfm=conform_dims(dimsizes(data_qu),dp,0)          ; confirm dp to caculate sum(var*dp)

    era_wgt_qu = data_qu(0,:,:)
    era_wgt_qv = data_qv(0,:,:)

    era_wgt_qu=(/dim_sum_n(dp_cfm*data_qu,0)/)               ; sum(var*dp)
    era_wgt_qv=(/dim_sum_n(dp_cfm*data_qv,0)/)               ; sum(var*dp)

    era_wgt_qu = (/era_wgt_qu/9.8/)
    era_wgt_qv = (/era_wgt_qv/9.8/)

    g_lon   =   era_wgt_qu&lon
    g_lat   =   era_wgt_qu&lat
    g_nlat  =   dimsizes(g_lat)

    DivQ = uv2dv_cfd(era_wgt_qu,era_wgt_qv,g_lat,g_lon,2)
  
    DivQ   =   (/-DivQ*86400/) 
    copy_VarCoords(era_wgt_qu, DivQ)
return DivQ
end     

begin
  ;;;read CESM data 
		path_wind= (/"/home/ys17-19/lsh/CESM-data/F/F_2000_GHG/","/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/"/)
    ; path_wind= (/"/home/ys17-19/lsh/CESM-data/F/F_2000_IPO/","/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/"/)
		file_name = (/"F_2000_ghg.cam.h0.0101-4012_chazhi.nc","F_2000_CTRL.cam.h0.0101-4012_chazhi.nc"/)
		filein = path_wind + file_name 
        
    ; ;;;region 
    ;   lat1 = -10
    ;   lat2 = 50
    ;   lon1 = 60
    ;   lon2 = 160
    ; ;;;

    ;;;region 
      lat1 = -50
      lat2 = 50
      lon1 = 0
      lon2 = 360
    ;;;

  

  
        ;;;time span
         
        startmon =601
        endmon = 4012

        ;;


		f1 = addfile(filein(0),"r")

		date=cd_calendar(f1->time, 1)
		timeind=ind(date.ge.startmon.and.date.le.endmon)

		u_ghg = f1->U(timeind,:,{lat1:lat2},{lon1:lon2})
    v_ghg = f1->V(timeind,:,{lat1:lat2},{lon1:lon2})
    q_ghg = f1->Q(timeind,:,{lat1:lat2},{lon1:lon2})



		f2 = addfile(filein(1),"r")

    u_CT = f2->U(timeind,:,{lat1:lat2},{lon1:lon2})
    v_CT = f2->V(timeind,:,{lat1:lat2},{lon1:lon2})
    q_CT = f2->Q(timeind,:,{lat1:lat2},{lon1:lon2})
  
   ;;;
   print("=====finish read data==========")
;==============================================================
; Mission 1. calculate Q  σ(qv)  =v*σq + q*σv + σq*σv ;;chage is for long-term mean differece 
; ==============================================================
    ct_qu  = q_CT
    ct_qv  = q_CT
    ct_qu  = u_CT*q_CT
    printVarSummary(ct_qu)
    ct_qv  = v_CT*q_CT

    ghg_qu = q_ghg
    ghg_qv = q_ghg
    ghg_qu = u_ghg*q_ghg
    ghg_qv = v_ghg*q_ghg
    

    ct_qv_clm = dim_avg_n_Wrap(ct_qv, 0)
    ct_qu_clm = dim_avg_n_Wrap(ct_qu, 0)
    
    ghg_qv_clm = dim_avg_n_Wrap(ghg_qv, 0)
    ghg_qu_clm = dim_avg_n_Wrap(ghg_qu, 0)

    u_CT_clm = dim_avg_n_Wrap(u_CT, 0)
    v_CT_clm = dim_avg_n_Wrap(v_CT, 0)

    q_CT_clm = dim_avg_n_Wrap(q_CT, 0)

    diff_qu = ghg_qu_clm - ct_qu_clm
    diff_qv = ghg_qv_clm - ct_qv_clm
    copy_VarCoords(ghg_qu_clm, diff_qu)
    copy_VarCoords(ghg_qv_clm, diff_qv)

    diff_q = dim_avg_n_Wrap(q_ghg, 0) - dim_avg_n_Wrap(q_CT, 0)
    diff_u = dim_avg_n_Wrap(u_ghg, 0) - dim_avg_n_Wrap(u_CT, 0)
    diff_v = dim_avg_n_Wrap(v_ghg, 0) - dim_avg_n_Wrap(v_CT, 0)

    copy_VarCoords(q_ghg(0,:,:,:), diff_q)
    copy_VarCoords(u_ghg(0,:,:,:), diff_u)
    copy_VarCoords(v_ghg(0,:,:,:), diff_v)

    ;;;thermaldynamic term  v*σq 


    thermal_qu = u_CT_clm*diff_q 
    thermal_qv = v_CT_clm*diff_q 
    copy_VarCoords(u_CT_clm, thermal_qu)
    copy_VarCoords(v_CT_clm, thermal_qv)

    dyna_qu = q_CT_clm*diff_u
    dyna_qv = q_CT_clm*diff_v
    
    copy_VarCoords(q_CT(0,:,:,:), dyna_qu)
    copy_VarCoords(q_CT(0,:,:,:), dyna_qv)

    quad_qu = diff_u*diff_q
    quad_qv = diff_v*diff_q
    
    copy_VarCoords(diff_u, quad_qu)
    copy_VarCoords(diff_v, quad_qv)

    
    ;;;
;==============================================================
; Mission 2. calculate divergence and integral   ;;;计算水汽通量的散度以及积分
; ==============================================================

    diff_divQ  = cal_divQ(diff_qu,diff_qv)

    thermal_divQ = cal_divQ(thermal_qu , thermal_qv)

    dyna_divQ = cal_divQ(dyna_qu,dyna_qv)

    eddy_divQ  = cal_divQ(quad_qu,quad_qv)
    print("=====finish DivQ calculate==========")
  

   
;==============================================================
; Mission 6. draw plot
; ==============================================================
  wks   = gsn_open_wks("pdf", "CESM-diff2CTrun-ghg-global")
  plot=new(4,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
  res@cnFillPalette="MPL_RdBu"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -4.0
  res@cnMaxLevelValF       =   4.0
  res@cnLevelSpacingF      =  0.5
  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             = 90
  ; res@mpMaxLonF             = 160
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 30
  res@mpMaxLonF             = 280


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

  res@mpFillOn = False
  ; res@gsnScalarContour        = True
  ;   ;;;;;;;LineArrow
  ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  ; ; res@vcGlyphStyle="CurlyVector"

  ; res@vcMinDistanceF           = 0.02           ; thin vectors
  ; ; res@vcMinMagnitudeF          =0.04
  ; ; res@vcRefMagnitudeF          = 1.0
 
  ; res@vcRefLengthF             = 0.08     ; define length of v
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  ; res@vcLineArrowThicknessF=1.0
  ; res@vcVectorDrawOrder="PostDraw"
  ; res@vcRefMagnitudeF          = 1000.0


  res@gsnLeftString = "Moisture Convergence Flux"


  plot(0)=gsn_csm_contour_map(wks, diff_divQ, res)
  res@gsnLeftString = "Dynamic component"

  plot(1)=gsn_csm_contour_map(wks, dyna_divQ, res)
  res@gsnLeftString = "Thermaldynamic component"

  plot(2)=gsn_csm_contour_map(wks, thermal_divQ, res)
  res@gsnLeftString = "eddy component"

  plot(3)=gsn_csm_contour_map(wks, eddy_divQ, res)
  


  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = "" + ispan(1,5,1) + ""
   
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

  
 pres = True
 pres@gsnPanelLabelBar = True
 gsn_panel(wks, plot, (/2,2/), pres)
 

end 
    
   
   

  





    
   ;;;

   






