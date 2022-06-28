begin 
diri_in = "/home/ys17-23/lsh/LBM/ln_solver/data/output"
; f1 = addfile(diri_in+"/1k_monsoon_heat.nc", "r")
; f1 = addfile(diri_in+"/t21_outdata_monh_3.nc", "r")
; f1 = addfile(diri_in+"/t42_outdata_monh_3.lar.nc", "r")
; f1 = addfile(diri_in+"/t42_outdata_monh_3_wnp_1.5_new.nc", "r")
; f1 = addfile(diri_in+"/t42_outdata_monh_3.nc", "r")
f1 = addfile(diri_in+"/t21_outdata_monh_3.nc", "r")

u = f1->u
v = f1->v
; ta = f1->psi
ta = f1->p
; ta = ta/100.0
; ta  = ta/10^6
printMinMax(ta, 1)
printVarSummary(ta)


; f2 = addfile(diri_in+"/t42_outdata_monh_3_wnp_1.5.nc", "r")
; f2 = addfile(diri_in+"/t21_outdata_monh_3_wnp_1.5.nc", "r")
f2 = addfile(diri_in+"/t21_outdata_monh_3_wnp_1.5_125-155.nc", "r")

u2 = f2->u
v2 = f2->v
; ta = f1->psi
ta2 = f2->p
; ta2 = ta2/100.0
; ta2  = ta2/10^6



;;;forcing data 
    diri_in = "/home/ys17-23/lsh/LBM/draw/"
    ; f1 = addfile(diri_in+"/frc.test_1k.nc", "r")
    ; f1 = addfile(diri_in+"/t21_frc_5_25_3k.nc", "r")
    ; f1 := addfile("frc.t21l20.test.new2.nc", "r")
    ; f1 := addfile("t21_frc_indh_wnpc_1k.nc", "r")
    ; f1 := addfile("t42_frc_indh3_wnpc1.5.nc", "r")
    f1 := addfile(diri_in+"t21_frc_indh3.nc", "r")
    t1 = f1->t(0,:,:,:) 

    ;;;forcing
    t_850 = t1({0.45},:,:)
    t_850 = t_850*86400.0
    t_850 = t_850
    ;;;;
    diri_in = "/home/ys17-23/lsh/LBM/draw/"
    ; f1 = addfile(diri_in+"/frc.test_1k.nc", "r")
    ; f1 = addfile(diri_in+"/t21_frc_5_25_3k.nc", "r")
    ; f1 := addfile("frc.t21l20.test.new2.nc", "r")
    ; f1 := addfile("t21_frc_indh_wnpc_1k.nc", "r")
    ; f2 := addfile("t42_frc_indh3_wnpc1.5.nc", "r")
    ; f2 := addfile("t21_frc_indh3_wnpc1.5.nc", "r")
    f2 := addfile(diri_in+"t21_frc_indh3wnp1.5.125-155.nc", "r")
    t2 = f2->t(0,:,:,:) 

    ;;;forcing
    t_200 = t2({0.45},:,:)
    t_200 = t_200*86400.0
    t_200 = t_200
    ;;;;

; t = f1->t(:,:,::-1,:)
; v = f1->v(:,:,::-1,:)
; d = f1->d(:,:,::-1,:)
; p = f1->p(:,:,::-1,:)
  printMinMax(t2, 1)
;   printMinMax(v, 1)
;   printMinMax(d, 1)
;   printMinMax(p, 1)



; u_mean = dim_avg_n_Wrap(u(14:19,:,:,:),0)
; v_mean = dim_avg_n_Wrap(v(14:19,:,:,:),0)

u_mean = dim_avg_n_Wrap(u(19:29,:,:,:),0)
v_mean = dim_avg_n_Wrap(v(19:29,:,:,:),0)
ta_mean = dim_avg_n_Wrap(ta(19:29,0,:,:),0)


u2_mean = dim_avg_n_Wrap(u2(19:29,:,:,:),0)
v2_mean = dim_avg_n_Wrap(v2(19:29,:,:,:),0)
ta2_mean = dim_avg_n_Wrap(ta2(19:29,0,:,:),0)


u_850 = u_mean({850},:,:)
v_850 = v_mean({850},:,:)
ta_850 = ta_mean
; ta_850 = ta_mean({850},:,:)




printMinMax(t_850, 1)


u_200 = u2_mean({850},:,:)
v_200 = v2_mean({850},:,:)
ta_200 = ta2_mean


;;;draw 

plot = new(2,graphic)
wks = gsn_open_wks("pdf", "UVPSI-heating-over-5-25")
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  res@cnFillPalette = "BlueWhiteOrangeRed"

  res@cnFillOn=True
  res@cnLinesOn=True
  ; res@cnMonoLineColor = True
  ; res@cnLineColor = "red"

  res@gsnContourZeroLineThicknessF = 0.01

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -20
  res@mpMaxLatF             = 40
  res@mpMinLonF             = 30
  res@mpMaxLonF             = 290
  res@mpFillOn = False

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

   
  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.5
  res@tmYLMajorThicknessF=1.5
  res@tmXBMinorThicknessF = 1.5
  res@tmYLMinorThicknessF = 1.5
  res@tmBorderThicknessF=1.5
  res@vcRefAnnoOn = False

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -3.0
  ; res@cnMaxLevelValF       =   3.0
  ; res@cnLevelSpacingF      =  .25
  res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  res@cnLevels = (/-1.5,-1.,-0.5,-0.25,0.25,0.5,1,1.5,2.0,2.5,3.0/)
  res@cnFillColors= (/4,21,37,53,117,177,193,209,225,234,241,255/) - 2 


   ; res@cnMinLevelValF       =   -1.0
  ; res@cnMaxLevelValF       =   1.0
  ; res@cnLevelSpacingF      =  .1
  res@gsnLeftString= ""

  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; ; res@cnLevels := (/-12,-9,-6,-3,3,6,9,12/)*0.01
  ; res@cnLevels := (/-20,-15,-10,-5,5,10,15,20/)*0.01
  ; res@cnInfoLabelOn = False

  res@gsnLeftStringFontHeightF = 0.015

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

   res@gsnScalarContour        = True

    ;;;;;;;LineArrow
    ; res@vcGlyphStyle = "CurlyVector"

    res@vcMonoLineArrowColor = True
    ; res@vcLineArrowColor = "darkgreen"
    res@vcLineArrowColor = "black"
    res@vcMinDistanceF           = 0.016          ; thin vectors
    ; res@vcMinMagnitudeF          = 0.05
    ; res@vcMinMagnitudeF          = 0.1
    ; res@vcRefMagnitudeF          = 0.5
    res@vcMinMagnitudeF          = 0.10
    res@vcRefMagnitudeF          = 1.0

    ; res@vcMinMagnitudeF          = 0.05
    ; res@vcRefMagnitudeF          = 0.5

    res@vcRefLengthF             = 0.04     ; define length of v
    ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcLineArrowThicknessF = 1.5
    ; res@vcLineArrowHeadMaxSizeF=0.008
    ; res@vcLineArrowHeadMinSizeF=0.008
    ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    ; res@vcLineArrowHeadMinSizeF = 0.006

    res@vcLineArrowHeadMaxSizeF = 0.006           ;change the vector line arrows
    res@vcLineArrowHeadMinSizeF = 0.006
    

    res@gsnLeftStringFontHeightF = 0.015
    res@gsnRightStringFontHeightF = 0.015

   ; res@gsnRightString = "Strong SEASM"
    res@gsnLeftString = "(a) ISM Heating"
     res@tmXBLabelsOn = False
    plot(0) = gsn_csm_vector_scalar_map(wks,u_850, v_850,t_850, res)

     res@tmXBLabelsOn = True
    res@gsnLeftString = "(b) ISM Heaitng + WNP Cooling"

    res@vcRefAnnoOn = True
    res@vcRefAnnoString1 = "1 m/s"
    res@vcRefAnnoOrthogonalPosF = 0.3
    res@vcRefAnnoParallelPosF = 0.95
    res@vcRefAnnoFontHeightF = 0.01
    res@vcRefAnnoPerimOn = False

    plot(1) = gsn_csm_vector_scalar_map(wks,u_200, v_200,t_200, res) 

    rescn = res 
    rescn@cnFillOn =False
    rescn@cnLinesOn = True
    rescn@cnInfoLabelOn  = False

    ; rescn@cnMinLevelValF       =   -1.0
    ; rescn@cnMaxLevelValF       =   1.0
    ; rescn@cnLevelSpacingF      =  .2
     rescn@cnLevelSelectionMode = "ManualLevels" 
    rescn@cnMinLevelValF       =   -0.5
    rescn@cnMaxLevelValF       =   0.5
    rescn@cnLevelSpacingF      =  0.1
    rescn@gsnContourNegLineDashPattern= 1
    rescn@gsnLeftString= ""
    rescn@gsnRightString = ""


   plot_add = gsn_csm_contour(wks, ta_850, rescn)
   plot_add2 = gsn_csm_contour(wks, ta_200 , rescn)
  
   ; overlay(plot(0), plot_add)
   ; overlay(plot(1), plot_add2)

pres = True
pres@gsnPanelLabelBar  = True
pres@pmLabelBarWidthF = 0.6
pres@pmLabelBarHeightF = 0.05
pres@pmLabelBarOrthogonalPosF = 0.05

; pres@lbLabelStride=1
gsn_panel(wks, plot, (/2,1/), pres)

end

