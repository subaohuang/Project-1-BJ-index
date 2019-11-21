begin
fils1  =  "/home/ys17-19/lsh/HEAT_MC_SON.cam.h0.25001-28912_chazhi.nc"

fils2 = "/home/ys17-19/lsh/B2000_f19g16_CP_CTRL.cam.h0.25001-28912_chazhi.nc"

f1 = addfile(fils1, "r")
f2 = addfile(fils2, "r")

 date  =  cd_calendar(f1->time, 1)

 timeind  =  ind(date.le.28912.and.date.ge.26001)



u_heat = f1->U(timeind,{850},:,:)
u_ct = f2->U(timeind,{850},:,:)

v_heat = f1->V(timeind,{850},:,:)
v_ct = f2->V(timeind,{850},:,:)

slp_heat = f1->PSL(timeind,:,:)
slp_ct = f2->PSL(timeind,:,:)


u_heat_mam = month_to_season(u_heat, "SON")
u_ct_mam = month_to_season(u_ct, "SON")

v_heat_mam = month_to_season(v_heat, "SON")
v_ct_mam = month_to_season(v_ct, "SON")

slp_heat_mam = month_to_season(slp_heat, "SON")
slp_ct_mam = month_to_season(slp_ct, "SON")


u_heat_mam_pa = dim_avg_n_Wrap(u_heat_mam(:,{-5:5},{150:210}), (/1,2/))

u_heat_mam_id = dim_avg_n_Wrap(u_heat_mam(:,{-5:5},{70:100}), (/1,2/))


u_ct_mam_pa = dim_avg_n_Wrap(u_ct_mam(:,{-5:5},{150:210}), (/1,2/))

u_ct_mam_id = dim_avg_n_Wrap(u_ct_mam(:,{-5:5},{70:100}), (/1,2/))

corr_heat = escorc_n(u_heat_mam_pa, u_heat_mam_id, 0, 0)

corr_ct = escorc_n(u_ct_mam_pa, u_ct_mam_id, 0, 0)

print("corr_heat="+corr_heat)
print("corr_ct="+corr_ct)




u_diff = dim_avg_n_Wrap(u_heat_mam, 0) - dim_avg_n_Wrap(u_ct_mam, 0)

v_diff = dim_avg_n_Wrap(v_heat_mam, 0) - dim_avg_n_Wrap(v_ct_mam, 0)

slp_diff = dim_avg_n_Wrap(slp_heat_mam, 0) - dim_avg_n_Wrap(slp_ct_mam, 0)

copy_VarCoords(u_heat(1,:,:), u_diff)
copy_VarCoords(v_heat(1,:,:), v_diff)
copy_VarCoords(slp_heat(1,:,:), slp_diff)


;;;画图

;;---------------画趋势图--------------------------------------
    wks=gsn_open_wks("pdf", "DIFF-HEAT-MC-U-V-SLP")
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

    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    ; res@cnMinLevelValF       =   -0.3
    ; res@cnMaxLevelValF       =   0.3
    ; res@cnLevelSpacingF      =  0.05

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF          = 180.0

    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.6
    res@vpHeightF = 0.3

    ;;
     ;;;LineArrow
    res@gsnScalarContour        = True
    
    res@vcGlyphStyle             = "LineArrow"      ; turn on curly vectors
    ; res@vcGlyphStyle="CurlyVector"

    res@vcMinDistanceF           = 0.03           ; thin vectors
    ; res@vcMinMagnitudeF          =0.05
    ; res@vcRefMagnitudeF          = 0.6
    res@vcRefLengthF             = 0.04     ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF=1.0
    res@vcLineArrowHeadMaxSizeF=0.01
    res@vcLineArrowHeadMinSizeF=0.01
    res@vcVectorDrawOrder="PostDraw"

    ;;;CurlyVector
    ; res@vcLineArrowThicknessF=1.5
    ; res@vcMinFracLengthF=0.3
    ; res@vcMinDistanceF           = 0.015           ; thin vectors
    ; res@vcLineArrowThicknessF=1.0
    ; res@vcMinDistanceF           = 0.02          ; thin vectors
    ; res@vcMinMagnitudeF          =3.0
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    ;;;
    res@gsnLeftString = "Observation"
    plot = gsn_csm_vector_scalar_map(wks,u_diff,v_diff,slp_diff,res)
   
   draw(plot)
   frame(wks)

end




;;;
