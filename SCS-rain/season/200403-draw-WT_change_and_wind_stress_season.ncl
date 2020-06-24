  ;;********************************************************
;;计算海洋平流项的时间变化的比较，分为动力项和热力项的作用,只关注
;;南北两侧的水平平流的动力项与热力项以及垂直扩散项
;;                                        
;;                                     linshuheng 2020.03.11
;;
;;
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
function  wgt_annual(total_month)
begin
    weight1 = ispan(1,12,1)*1.0
    weight1 = (2*weight1-1)/24

    weight2 = ispan(1,12,1)*1.0
    weight2 = (25-2*weight2)/24
    total_annual_test = new(34,"float")
    total_annual_test@_FillValue = default_fillvalue("float")
    printVarSummary(total_annual_test)

    do i =0,32,1
      total_annual_test(i+1) = sum(total_month(i*12:(i+1)*12-1)*weight1) + sum(total_month((i+1)*12:(i+2)*12-1)*weight2)
    end do    
return total_annual_test
end    
    


;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 201312

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region

; ;;;region
;    lat1 = -5
;    lat2 = 45
;    lon1 = 90
;    lon2 = 160
; ;;;

;;;region
   lat1 = -10
   lat2 = 20
   lon1 = 100
   lon2 = 160
;;;

   
;;;read soda data and then calculate ocean dynamic : bar(v*σT)， bar(w*σT)，bar denotes 
;;;averaged in mixed layer
    
    f2=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_u_mn_1980-2015_chazhi.nc", "r")   
    f3=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_v_mn_1980-2015_chazhi.nc", "r")   
    f4=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_wt_mn_1980-2015_chazhi.nc", "r")   
    f5=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_temp_mn_1980-2015_chazhi.nc", "r")
    ; f6=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_mlt_mn_1980-2015_chazhi.nc", "r") 
     f6=addfile("/home/ys17-19/lsh/Project/SCS-rain/SODA3.4.2_mld_density_criterion_0.8.nc", "r")  

   ;;; 选择季节
   season = "DJF"


    indyear = (2013-1980+1)*12-1
    u_r = f2->u(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    v_r = f3->v(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    wt_r = f4->wt(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    T_r= f5->temp(0:indyear,{0:330},{lat1:lat2},{lon1:lon2})
    ; mlt_r = f6->mlt(0:indyear,{lat1:lat2},{lon1:lon2})
    mlt = f6->mld(0:indyear,{lat1:lat2},{lon1:lon2})
    

    u = tofloat(u_r)
    v = tofloat(v_r)
    wt = tofloat(wt_r)
    T = tofloat(T_r)
    ; mlt = tofloat(mlt_r)
    copy_VarCoords(u_r, u)
    copy_VarCoords(v_r, v)
    copy_VarCoords(wt_r, wt)
    copy_VarCoords(T_r, T)
    ; copy_VarCoords(mlt_r, mlt)

    delete(u_r)
    delete(v_r)
    delete(wt_r)
    delete(T_r)
    ; delete(mlt_r)

   
    depth = T&depth
printVarSummary(wt)



;;;calculate lateral  and bottom boundaries of the domain Heat transport, bottom is 45m , lateral is 
  ;;; four horizatal boundaries

    ;;parameter 
     pi = get_pi("float")
     dx = u&lat
     dx = 6378388.*cos(0.0174533*u&lat)
     dx!0 = "lat"
     copy_VarCoords(u(1,1,:,1), dx)

     dy = (pi/180)*1*6378388.
   


      dz_3d_test = new(dimsizes(u),float)
      do i = 0,dimsizes(mlt&time)-1
        dz_3d_test(i,:,:,:) = dpres_plevel(depth, mlt(i,:,:), 0, 0)
      end do 
      copy_VarCoords(u, dz_3d_test)
      printVarSummary(dz_3d_test)
  
   
     dx_3d_test = conform_dims(dimsizes(T(:,:,{0:15},{110:155})) ,dx({0:15}) , 2)


     dx_3d_test  = where(ismissing(T(:,:,{0:15},{110:155})), 0.0, dx_3d_test)
     dz_3d_test  = where(ismissing(T), 0.0, dz_3d_test)
     
     volume = dim_sum_n(dx_3d_test*dz_3d_test(:,:,{0:15},{110:155})*dy, (/1,2,3/))


     Tbox =dim_sum_n_Wrap(T(:,:,{0:15},{110:155})*dy*dx_3d_test*dz_3d_test(:,:,{0:15},{110:155}), (/1,2,3/))/volume
     copy_VarCoords(T(:,1,1,1), Tbox)

     ; volume = sum(45*dy*45.58*dx({0:15})) 

     dT = center_finite_diff_n(Tbox, 1, False, 0, 0)
     copy_VarCoords(Tbox, dT)
     dT_annual  = month_to_season(dT, season) 
    
     Tbox_annual = month_to_season(Tbox, season)


  
     dT_test = center_finite_diff_n(Tbox_annual, 1, False, 0, 0)
     dT_wgt = dT_test
     dT_wgt(1:33) = Tbox_annual(1:33) - Tbox_annual(0:32)
   
     

    ;; lateral heat transport north and south v*(T -Tbox)dxdz,south: 0N ,110-155,north: 15N ,110-155
       x1 = 0
       x2 = 15
       y1  = 110.5
       y2 = 155
       z1  = 45.6
    
       dim =dimsizes(v(:,:,{x1},{y1:y2}))
    
       dz_3d1 = dz_3d_test(:,:,{x1},{y1:y2})
       dz_3d2 = dz_3d_test(:,:,{x2},{y1:y2})
       Tbox_3d = conform_dims(dim ,Tbox , 0)


    ;;;;;;;垂直扩散变化项
     
    ;;; lateral heat transport vertical w*(T -Tbox)dxdy 
       
       dim:= dimsizes(wt(:,4,{x1:x2},{y1:y2}))
       dx_3d :=conform_dims(dim ,dx({0:15}) , 1)
       Tbox_3d :=conform_dims(dim, Tbox, 0)       

       wtT = wt(:,4,{x1:x2},{y1:y2})
       mlt_scs = mlt(:,{x1:x2},{y1:y2})

       printVarSummary(mlt_scs)
      
       mlt_scs = where(ismissing(mlt_scs), avg(mlt_scs), mlt_scs)
      
       printMinMax(mlt_scs, 1)
 
  
   ;;δ(WT)
    ; f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/SODA3.4.2_wt_T_in_MLD.nc", "r")
    f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/SODA3.4.2_wt_T_in_MLD_density_defined_0.8.nc", "r")

    wt_mlt = f1->wt_mlt
    T_mlt = f1->T_mlt_diff

    wtT = wt_mlt
    wtT = (/wt_mlt*T_mlt/)
    


   wtT&lat@units = "degrees_north"
   wtT&lon@units = "degrees_east"


    dwt  = wtT(1,:,:)

    wtT_annual = month_to_season(wtT, season)



   dwt = (/dim_avg_n_Wrap(wtT_annual(17:,:,:), 0) - dim_avg_n_Wrap(wtT_annual(0:17-1,:,:), 0)/)
   dwt&lat@units = "degrees_north"
   dwt&lon@units = "degrees_east"

    printVarSummary(dwt)

   ;;

   ;;δ(W)T
    tdw = wtT(1,:,:)
    
    T_mlt_annual = month_to_season(T_mlt, season) 

    T_clm  = dim_avg_n_Wrap(T_mlt_annual(0:17-1,:,:), 0)


    wt_mlt_annual = month_to_season(wt_mlt, season)

    w_delta = dim_avg_n_Wrap(wt_mlt_annual(17:,:,:), 0) - dim_avg_n_Wrap(wt_mlt_annual(0:17-1,:,:), 0)
    
    tdw = (/T_clm*w_delta/)
   
    ;;;

    ;;δ(T)w
    wdt = wtT(1,:,:)



    w_clm  = dim_avg_n_Wrap(month_to_season(wt_mlt(0:12*17-1,:,:), season),0)
    T_delta = dim_avg_n_Wrap(T_mlt_annual(17:,:,:), 0) - dim_avg_n_Wrap(T_mlt_annual(0:17-1,:,:), 0)
    
    wdt = (/T_delta*w_clm/)


    ;;

   ;;δ(T)δ(w)
    dwdt = wtT(1,:,:)
    dwdt = (/w_delta*T_delta/)

   ;;

   ;;eddy 
   eddy = wtT(1,:,:)
   eddy = (/dwt - wdt - tdw - dwdt/)
   eddy&lat@units = "degrees_north"
   eddy&lon@units = "degrees_east"

   ;;
  
    dwt = (/dwt*10^5/)
    wdt = (/wdt*10^5/)
    tdw = (/tdw*10^5/)
    dwdt = (/dwdt*10^5/)
    eddy = (/eddy*10^5/)
  ;;



   T_clm&lat@units = "degrees_north"
   T_clm&lon@units = "degrees_east"
   
  copy_VarCoords(T_clm, T_delta )
  copy_VarCoords(T_clm, w_delta )
  copy_VarCoords(T_clm, w_clm )

  w_clm = (/w_clm*10^5/)
  w_delta = (/w_delta*10^5/)


;;;;wind stress

  
;;;read soda data and then calculate EKMAN transport and curl of wind stress
    
       
    f7=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_taux_mn_1980-2015_chazhi.nc", "r")   
    f8=addfile("/home/yangsong3/data-observation/SODA/3.4.2/soda3.4.2_tauy_mn_1980-2015_chazhi.nc", "r")   


    taux_r = f7->taux(0:indyear,{lat1:lat2},{lon1:lon2})
    tauy_r = f8->tauy(0:indyear,{lat1:lat2},{lon1:lon2})

    taux = tofloat(taux_r)
    tauy = tofloat(tauy_r)
    
    copy_VarCoords(taux_r, taux) 
    copy_VarCoords(tauy_r, tauy)
   
    delete(taux_r)
    delete(tauy_r)
    
    taux&lat@units = "degrees_north"
    taux&lon@units = "degrees_east"
    tauy&lat@units = "degrees_north"
    tauy&lon@units = "degrees_east"

    lat = taux&lat
    pi = get_pi("float")
    omega = 2*pi/86400
    f = 2*omega*sin(pi*lat/180)
    density = 1025
    f_3d = conform_dims(dimsizes(taux), f, 1)
    
    ue = tauy
    ve = taux  
    ue  = (/tauy/(density*f_3d)/)   ;;;EKMAN
    ve = (/-1*taux/(density*f_3d)/)   ;;;Ekman


    vr = uv2vr_cfd((/taux/(density*f_3d)/) ,(/tauy/(density*f_3d)/) ,taux&lat,taux&lon,2)
     ; vr = uv2vr_cfd(taux ,tauy ,taux&lat,taux&lon,2)

    copy_VarCoords(taux, vr)
    printVarSummary(vr)
    
    wg = vr 
    copy_VarCoords(taux, wg)


    ;;;年平均
       ue_annual = month_to_season(taux, season)
       ve_annual = month_to_season(tauy, season)
       wg_annual = month_to_season(wg, season)
    ;;;
    

    ;;;clm diff 
    
    ue_annual_1 = dim_avg_n_Wrap(ue_annual(0:12,:,:), 0)
    ue_annual_2 = dim_avg_n_Wrap(ue_annual(13:,:,:), 0)

    ue_annual_diff = ue_annual_1 
    ue_annual_diff = (/ue_annual_2 - ue_annual_1/)

    ve_annual_1 = dim_avg_n_Wrap(ve_annual(0:12,:,:), 0)
    ve_annual_2 = dim_avg_n_Wrap(ve_annual(13:,:,:), 0)

    ve_annual_diff = ve_annual_1 
    ve_annual_diff = (/ve_annual_2 - ve_annual_1/)


    wg_annual_1 = dim_avg_n_Wrap(wg_annual(0:12,:,:), 0)
    wg_annual_2 = dim_avg_n_Wrap(wg_annual(13:,:,:), 0)

    wg_annual_diff = wg_annual_1 
    wg_annual_diff = (/wg_annual_2 - wg_annual_1/)
    ;;;;


    ;;;;clm 
    ue_clm = dim_avg_n_Wrap(ue_annual, 0)
    ve_clm = dim_avg_n_Wrap(ve_annual, 0)
    wg_clm = dim_avg_n_Wrap(wg_annual, 0)
    wg_clm  = (/wg_clm*10^5/)
    ;;;


  
 
;;;calculate correlation coef

   rc_ue = ue_annual_diff
   rc_ve = ve_annual_diff
   rc_wg = wg_annual_diff
   rc_wg  = (/rc_wg*10^5/)
   ;;;;



    
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "WT的变化以及气候态以及wind-stress-mld-0.8-define"+season)

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="cmp_b2r"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160

  res@mpMinLatF             = 0
  res@mpMaxLatF             = 15
  res@mpMinLonF             = 110
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

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2
 

  plot = new(6,graphic)

  res@gsnLeftString = "Wclm"
  plot(0) = gsn_csm_contour_map(wks, w_clm, res)
  
  res@gsnLeftString = "dW"
  res@cnMinLevelValF       =   -0.2
  res@cnMaxLevelValF       =   0.2
  res@cnLevelSpacingF      =  0.04  
  plot(1) = gsn_csm_contour_map(wks, w_delta, res)
  
  res@gsnLeftString = "Tclm"
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -4
  res@cnMaxLevelValF       =   0
  res@cnLevelSpacingF      =  0.5
  
  res@cnFillPalette = "MPL_PuBu"
  res@gsnSpreadColorStart = 110
  res@gsnSpreadColorEnd = 3

  plot(2) = gsn_csm_contour_map(wks, T_clm, res)

  res@gsnLeftString = "dT"
  res@cnMinLevelValF       =   -2.0
  res@cnMaxLevelValF       =   0.0
  res@cnLevelSpacingF      =  0.25
  plot(3) = gsn_csm_contour_map(wks, T_delta, res)
   
  ;;;;;draw wind stress
  res2 = res
  res2@gsnLeftString = "Wind stress and curl climatology "


  res2@gsnScalarContour        = True
   ;;;;;;;LineArrow
  res2@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
  res2@vcGlyphStyle="CurlyVector"


  res2@vcRefMagnitudeF = 0.05
  res2@vcMinMagnitudeF = 0.005
  res2@vcMinDistanceF           = 0.018      ; thin vectors
  res2@vcRefLengthF             = 0.03 ; define length of v
 
  res2@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res2@vcLineArrowThicknessF=1.0


  res2@vcLineArrowHeadMaxSizeF=0.006
  res2@vcLineArrowHeadMinSizeF=0.006
  res2@vcVectorDrawOrder="PostDraw"
  res2@mpFillOn =False


  res2@cnMinLevelValF       =   -0.5
  res2@cnMaxLevelValF       =   0.5
  res2@cnLevelSpacingF      =  0.1
  delete(res2@gsnSpreadColorEnd)
  delete(res2@gsnSpreadColorStart)
  

  color = read_colormap_file("/home/ys17-19/lsh/data/colordata/bluewhiteorangered_11color.rgb")

  ; res2@cnFillPalette = "GMT_polar"
  delete(res2@cnFillPalette)
  res2@cnFillPalette = color
  plot(4) = gsn_csm_vector_scalar_map(wks, ue_clm ,ve_clm,wg_clm ,res2)

  res2@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res2@cnMinLevelValF       =   -0.05
  res2@cnMaxLevelValF       =   0.05
  res2@cnLevelSpacingF      =  0.01

  res2@vcRefMagnitudeF = 0.005
  res2@vcMinMagnitudeF          =0.0005
  
  
 
  res2@gsnLeftString = "Wind stress and curl change "
  plot(5) = gsn_csm_vector_scalar_map(wks,rc_ue ,rc_ve,rc_wg ,res2)

  
  


  ;;;添加赤道直线
    pgx=ispan(30,390,10)
    pgy=new(dimsizes(pgx),typeof(pgx))
    pgy=0

    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=1.0
    pgres@gsLineDashPattern=1
  
    str = "" + ispan(1,dimsizes(plot),1)

   do i = 0,dimsizes(str)-1
     plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx, pgy, pgres)
   end do


  ;;;添加indo-pacific warm pool
    ; ;;添加区域
   ;;warm pool
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   str2 = "0" + ispan(1,dimsizes(plot),1)
   do i = 0,dimsizes(str2)-1
     plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   end do


  pres = True
  pres@gsnPanelLabelBar  = False

  gsn_panel(wks, plot, (/3,2/), pres)
 


end







   ;;;

;;;;;





