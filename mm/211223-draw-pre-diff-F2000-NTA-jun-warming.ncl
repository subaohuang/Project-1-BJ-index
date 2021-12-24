;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
end 

function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 

begin
  modelname = "F2000_TPP_indNTA_warming"
  ; modelname = "F2000_TPP_indALLTA_junwarming"

  fils1 = "/home/yangsong3/CMIP6/linshh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h1.0101-4012_chazhi.nc"

  fils2 = "/home/yangsong3/CMIP6/linshh/CESM-data/F_mm/"+modelname+"/"+modelname+".cam.h1.0001-0030_chazhi.nc"
  
  wks=gsn_open_wks("pdf", modelname+"-Jun-PREDIFF-alt-Jun-warming.pdf")
  
  plot = new(1,graphic)
  iplot = new(1,graphic)
  ; str3 = ispan(1,4,1) + ""

  f1=addfile(fils1, "r")
  f2=addfile(fils2, "r")

  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.601.and.date.le.3012)

  PRECL_CT  =  f1->PRECL(timeind,:,:)
  PRECC_CT  =  f1->PRECC(timeind,:,:)

  PRE_CT = PRECC_CT
  PRE_CT = (/PRECC_CT+PRECL_CT/)
  ; PRE_CT = (/PRECL_CT/)  ;;;大尺度降水
  PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day

  PRECL_id  =  f2->PRECL(timeind,:,:)
  PRECC_id  =  f2->PRECC(timeind,:,:)   

  PRE_id = PRECC_id
  PRE_id = (/PRECC_id+PRECL_id/)
  ; PRE_id = (/PRECL_id/)
  PRE_id = (/PRE_id*1000*3600*24/)  ;; m/s conver to mm/day
   
 

  ;;;;U,V
  u_CT  =  f1->U(timeind,{850},:,:)
  v_CT  =  f1->V(timeind,{850},:,:)

  u_id  =  f2->U(timeind,{850},:,:)
  v_id  =  f2->V(timeind,{850},:,:)

  PRE_CT_12 = split_12(PRE_CT)
  PRE_id_12 = split_12(PRE_id)
  
  u_CT_12 = split_12(u_CT)
  u_id_12 = split_12(u_id)

  v_CT_12 = split_12(v_CT)
  v_id_12 = split_12(v_id)
  
 

  PRE_CT_jja =  dim_avg_n_Wrap(PRE_CT_12(5:6,:,:,:), 0)
  PRE_indian_jja = dim_avg_n_Wrap(PRE_id_12(5:6,:,:,:), 0)


  PRE_CT_jja =  PRE_CT_12(5,:,:,:)
  PRE_indian_jja = PRE_id_12(5,:,:,:)


  PRE_CT_clm =  dim_avg_n_Wrap(PRE_CT_jja(:,:,:), 0)
  PRE_indian_clm = dim_avg_n_Wrap(PRE_indian_jja(:,:,:), 0)
  

 
  PRE_diff = PRE_CT_clm
  PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)


  u_CT_jja = dim_avg_n_Wrap(u_CT_12(5:6,:,:,:), 0)
  u_indian_jja = dim_avg_n_Wrap(u_id_12(5:6,:,:,:), 0)

  u_CT_jja =  u_CT_12(5,:,:,:)
  u_indian_jja = u_id_12(5,:,:,:)


  u_CT_clm =  dim_avg_n_Wrap(u_CT_jja(:,:,:), 0)
  u_indian_clm = dim_avg_n_Wrap(u_indian_jja(:,:,:), 0)
  u_diff = u_CT_clm
  u_diff = (/u_indian_clm - u_CT_clm/)



  v_CT_jja = dim_avg_n_Wrap(v_CT_12(5:6,:,:,:), 0)
  v_indian_jja = dim_avg_n_Wrap(v_id_12(5:6,:,:,:), 0)
 
  v_CT_jja =  v_CT_12(5,:,:,:)
  v_indian_jja = v_id_12(5,:,:,:)


  v_CT_clm =  dim_avg_n_Wrap(v_CT_jja(:,:,:), 0)
  v_indian_clm = dim_avg_n_Wrap(v_indian_jja(:,:,:), 0)

  v_diff = v_CT_clm
  v_diff = (/v_indian_clm - v_CT_clm/)
 


  ;;;t_test 
    prob_p = t_test( PRE_CT_jja,PRE_indian_jja)
    prob_u = t_test( u_CT_jja ,u_indian_jja)
    prob_v = t_test( v_CT_jja ,v_indian_jja)
  ;;;
   PRE_diff@_FillValue = 10^9
   u_diff@_FillValue = 10^9
   v_diff@_FillValue = 10^9
   ; PRE_diff = where(prob_p.le.0.1, PRE_diff, PRE_diff@_FillValue)
    u_diff = where(prob_u.le.0.1.or.prob_v.le.0.1, u_diff, u_diff@_FillValue)
    v_diff = where(prob_v.le.0.1.or.prob_u.le.0.1, v_diff, v_diff@_FillValue)
;;;


;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=  True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineColor = "White"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "MPL_BrBG"


    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    ; res@cnMinLevelValF       =   -2
    ; res@cnMaxLevelValF       =   2
    ; res@cnLevelSpacingF      =  0.5
    
    res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
    res@cnLevels := (/-1.5,-1.0,-0.5,-0.2,0.2,0.5,1.0,1.5/)

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -20.0
    res@mpMaxLatF             = 80.0
    ; res@mpMinLonF             = 40.0
    res@mpMinLonF             = -180.0
    res@mpMaxLonF             = 180.0
    ; res@mpMinLatF             = -60.0
    ; res@mpMaxLatF             = 60.0
    ; res@mpMinLonF             = 0.0
    ; res@mpMaxLonF             = 360.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2.0
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3
    ; ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
    res@gsnLeftString = "Jun NTA warming induced J PRE&850U,V diff (SEN - CTRL)"


    ; res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.01           ; thin vectors
    res@vcMinMagnitudeF          =0.2
    ; res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 2.0
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    ; res@vcLineArrowThicknessF=1.5
    ; res@vcLineArrowHeadMaxSizeF=0.005
    ; res@vcLineArrowHeadMinSizeF=0.005
    res@vcVectorDrawOrder="PostDraw"
    
    
    res@gsnRightString = ""
    plot = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PRE_diff,res)

    ; plot(i) = gsn_csm_contour_map(wks,PRE_diff,res)


     ; pgres=True
     ; pgres@gsLineColor="black"
     ; pgres@gsLineThicknessF=2.0
     ; pgres@gsLineDashPattern=0

     
     ; ;;;画box
     ; lat1 = 0
     ; lat2 = 15
     ; lon1 = 125
     ; lon2 = 160

     ; pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
     ; pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
     ; ;;;
     ;  plot@$str3(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)

    ;;;打点
    ; res_pro=True
    ; res_pro@gsnAddCyclic = False
    ; res_pro@gsnDraw=False
    ; res_pro@gsnFrame=False
    ; res_pro@cnFillOn = True

    ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
    ; res_pro@cnLevels = (/0,0.01/)
    ; res_pro@lbLabelBarOn = False
    ; res_pro@cnLinesOn=False
    ; ; res_pro@cnMonoLineColor =False
    ; ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    ; res_pro@cnMonoFillPattern = True
    ; res_pro@cnFillPattern = 17
    ; res_pro@cnFillDotSizeF = 0.003
  
    ; res_pro@cnMonoFillColor = False
    ; res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLineLabelsOn=False
    ; res_pro@tmYRMode = "Automatic"
    ; res_pro@cnInfoLabelOn=False
   
    ; res_pro@gsnCenterString = ""
    ; ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; ; res_pro@cnLineThicknessF = 3.
    ; ; res_pro@cnFillOn        = True
    ; iplot(i)=gsn_csm_contour(wks,prob_p,res_pro)


     ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.01     ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    iplot = gsn_csm_contour(wks,prob_p,res_pro)

     
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 6         ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越
    iplot = gsn_contour_shade(iplot , 0,0.01,opt)
    ; overlay(plot, iplot)
    print("finish subset figure")

    draw(plot)
    frame(wks)

 
end
  

