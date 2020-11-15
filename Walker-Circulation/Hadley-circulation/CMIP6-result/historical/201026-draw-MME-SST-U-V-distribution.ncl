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


plot = new(1,graphic)
wks = gsn_open_wks("pdf","201026-MME的CESM模式的SST,U,V趋势的水平分布")         ; send graphics to PNG file
f1=addfile("/home/ys17-19/lsh/Project/Walker-Circulation/Hadley-circulation/CMIP6-result/historical/sst_u_v_1000hpa_trend_1900-2014.nc", "r")

rc_u_all = lonFlip(f1->rc_u)
rc_v_all = lonFlip(f1->rc_v)
rc_sst_all = lonFlip(f1->rc_sst)


rc_u_mme = dim_avg_n_Wrap(rc_u_all, 0)
rc_v_mme = dim_avg_n_Wrap(rc_v_all, 0)
rc_sst_mme = dim_avg_n_Wrap(rc_sst_all, 0)



;;;
   
;;;;;;
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn = False

    res@gsnContourNegLineDashPattern = 1
    res@gsnContourPosLineDashPattern = 0
    res@gsnContourZeroLineThicknessF = 2.0



    res@cnLineColor = "black"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False

    ; res@cnFillPalette = "precip_diff_12lev"
    ; res@cnFillPalette = "BlueWhiteOrangeRed"
    res@cnFillPalette = "BlueDarkRed18"



    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    res@cnMinLevelValF       =   0.0
    res@cnMaxLevelValF       =   1.0
    res@cnLevelSpacingF      =  0.05
    ; res@cnFillPalette = "NCV_jaisnd"
    ; res@cnFillPalette = "sunshine_9lev"
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"

    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
   
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)/2


    res@mpFillOn = True
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3
    ;;


    ;;;LineArrow
    res@gsnScalarContour        = True
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.02           ; thin vectors
    ; res@vcMinMagnitudeF          =0.5
    ; res@vcMaxMagnitudeF = 7.0
    res@vcRefMagnitudeF          = 0.5
    res@vcRefLengthF             = 0.05     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.5
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005
    ; res@vcVectorDrawOrder="PostDraw"
    res@gsnRightString = ""
    plot = gsn_csm_vector_scalar_map(wks,rc_u_mme,rc_v_mme, rc_sst_mme,res)
    ; plot = gsn_csm_contour_map(wks, rc_mme,res)

draw(plot)
frame(wks)
end 


 


  