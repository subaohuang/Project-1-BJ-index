;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
begin
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h0.0101-4012.nc"
  ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindianwarming/F_2000_indian.cam.h0.0101-4012.nc"  
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addatlanticwarming/F_2000_atlantic.cam.h0.0101-4012.nc"  
  

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.601.and.date.le.4012)
    


    PSL_CT  =  f1->Z3(timeind,{500},:,:)
 

    u_CT  =  f1->U(timeind,{200},:,:)
    v_CT  =  f1->V(timeind,{200},:,:)

    PSL_indian  =  f2->Z3(timeind,{500},:,:)

    u_indian  =  f2->U(timeind,{200},:,:)
    v_indian  =  f2->V(timeind,{200},:,:)
    
    PSL_CT_clm = dim_avg_n_Wrap(PSL_CT, 0)
    PSL_indian_clm = dim_avg_n_Wrap(PSL_indian, 0)
    
    u_CT_clm = dim_avg_n_Wrap(u_CT, 0)
    u_indian_clm = dim_avg_n_Wrap(u_indian, 0)

    v_CT_clm = dim_avg_n_Wrap(v_CT, 0)
    v_indian_clm = dim_avg_n_Wrap(v_indian, 0)


    PSL_diff = PSL_CT_clm
    PSL_diff = (/PSL_indian_clm - PSL_CT_clm/)

    PSL_diff_zonmean = dim_avg_n_Wrap(PSL_diff,0)
    PSL_diff_zonmean_com = conform_dims(dimsizes(PSL_diff), PSL_diff_zonmean, 1)

    ; PSL_diff = (/PSL_diff - PSL_diff_zonmean_com/)

    
    u_diff = u_CT_clm
    u_diff = (/u_indian_clm - u_CT_clm/)
    
    v_diff = v_CT_clm
    v_diff = (/v_indian_clm - v_CT_clm/)
    
   
     wks=gsn_open_wks("pdf", "CESM-F2000-rossby-atlantic-new")
    ; wks=gsn_open_wks("pdf", "CESM-slp_u_v_pattern-F-pday-rmindian-rmalt-rmallocean")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True
    res@cnFillPalette = "BlRe"

    ; res@cnFillPalette = "BlueWhiteOrangeRed"


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    ; res@cnMinLevelValF       =   -0.4
    ; res@cnMaxLevelValF       =   0.4
    ; res@cnLevelSpacingF      =  0.1
    
    res@cnMinLevelValF       =   -2
    res@cnMaxLevelValF       =   2
    res@cnLevelSpacingF      =  0.2

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -90.0
    res@mpMaxLatF             = 90.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF          = 180.0

    ; res@mpMinLatF             = -30.0
    ; res@mpMaxLatF             = 30.0
    ; res@mpMinLonF             = 0.0
    ; res@mpMaxLonF             = 360.0
    ; res@mpCenterLonF          = 180.0


    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    


    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.7
    res@vpHeightF = 0.3
;
    ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.03           ; thin vectors
    res@vcMinMagnitudeF          =0.04
    res@vcRefMagnitudeF          = 1.0
    res@vcRefLengthF             = 0.04     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01
    res@vcVectorDrawOrder="PostDraw"
    
    res@gsnLeftString = "Geopotential height at the level of 510(difference between indian ocean warming 1degree and control run)"
    res@gsnRightString = ""
    ; plot = gsn_csm_vector_scalar_map(wks,u_diff,v_diff, PSL_diff,res)
    plot = gsn_csm_contour_map(wks, PSL_diff, res)


    draw(plot)
    frame(wks)
    print("done script")
  end
  

