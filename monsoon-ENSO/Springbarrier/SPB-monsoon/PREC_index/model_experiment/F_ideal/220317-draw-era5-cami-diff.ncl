
;;*************************************************************
;Aim: Output global forcing in certain period
;Mission:
;
;   1. Read data from datasets
;   2. Post-post process to deserved results
;   3. Output the data
;                     
;*************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1)/),float)
    do i=0,11   
    data_12(i,:,:) = data(i::12,:)
    end do 
return data_12
end 

function get_clm_data(diri,var,timestart,timeend)
begin

  fils = diri + var+ "_mon_1x1_1979_2020.nc"

  f1 = addfile(fils, "r")
  time = f1->time 
  date = cd_calendar(time, 1)
  timeind = ind(date.ge.timestart.and.date.le.timeend)

  data = short2flt(f1->$var$(timeind,:,{-90:90},{0:360}))
  data_clm = clmMonTLLL(data)

  ; data_clm_jjas = dim_avg_n_Wrap(data_clm(5:8,:,:,:), 0)
  data_clm_jjas = dim_avg_n_Wrap(data_clm(:,:,:,:), 0)

return data_clm_jjas
end 

begin
 
; ==============================================================
; Mission 1. read data
; ==============================================================

;;; read data from model intial file 

  fils = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703.nc") ;;;daily数据计算的q1

  f3 = addfile(fils, "r")

  US = f3->US(0,12,:,:)
  VS = f3->VS(0,12,:,:)
  T = f3->T(0,12,:,:)
  printVarSummary(US)
  printVarSummary(VS)
;;;;

;;; read data from model intial file 

  fils2 = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703_JJASMEAN.nc") ;;;daily数据计算的q1
  ; fils2 = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703_era5.nc") ;;;daily数据计算的q1

  f4 = addfile(fils2, "r")

  US_2 = f4->US(0,12,:,:)
  VS_2 = f4->VS(0,12,:,:)
  T_2 = f4->T(0,12,:,:)
  printVarSummary(US_2)
  printVarSummary(VS_2)
;;;;


  T_diff = T_2
  T_diff = (/T_2 - T/)

  US_diff = US_2
  US_diff = (/US_2 - US/)

  VS_diff = VS_2
  VS_diff = (/VS_2 - VS/)

; print(T_diff)
; exit()
printMinMax(US_diff, 1)
printMinMax(VS_diff, 1)
printMinMax(T_diff, 1)

; print(US_diff)
; exit()
plots = new(2,graphic)
  wks = gsn_open_wks("pdf", "cami-diff-test")
;;;;

    res=True 
    res@gsnDraw = False
    res@gsnFrame = False
    res@cnFillPalette = "BlueWhiteOrangeRed"
    ; res@cnFillPalette = "GMT_polar"
    res@cnFillOn=True
    res@cnLinesOn=True

    res@cnLineLabelsOn=False
    res@lbLabelBarOn=False
    res@lbOrientation="Vertical"

    res@mpLimitMode           = "LatLon"

    res@mpMinLatF             = -60
    res@mpMaxLatF             = 60
    res@mpMinLonF             = 0
    res@mpMaxLonF             = 360

    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpFillOn = False

    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBOn=True
    res@tmYROn=False

    res@tmXBMode="Explicit"
    res@tmYLMode="Explicit"

    res@tmXBLabelFontHeightF=0.01
    res@tmYLLabelFontHeightF=0.01
    res@tmXBMajorThicknessF=1.0
    res@tmYLMajorThicknessF=1.0
    res@tmBorderThicknessF=1.0

    ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
    ; ; res@cnMinLevelValF       =   -0.4
    ; ; res@cnMaxLevelValF       =   0.4
    ; ; res@cnLevelSpacingF      =  0.05

    ; res@cnMinLevelValF       =   -50.
    ; res@cnMaxLevelValF       =   50.
    ; res@cnLevelSpacingF      =  10.
    
    ; res@cnFillColors= (/3,4,5,6,8,8,10,11,13,14/) - 2 

    res@vcGlyphStyle="CurlyVector"
    ; res@vcGlyphStyle = "LineArrow"

    res@vcMonoLineArrowColor = True
    res@vcLineArrowColor = "black"
    res@vcMinDistanceF           = 0.018        ; thin vectors


    res@gsnScalarContour        = True

    res@vcRefLengthF             = 0.05    ; define length of v

    ; res@vcLineArrowThicknessF=1.2
    ; res@vcLineArrowHeadMaxSizeF=0.01
    ; res@vcLineArrowHeadMinSizeF=0.01
    ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    ; res@vcLineArrowHeadMinSizeF = 0.006

    ; res@vcRefMagnitudeF          =  0.2
    ; res@vcMinMagnitudeF          = 0.02
    res@vcRefMagnitudeF          =  0.5
    res@vcMinMagnitudeF          = 0.5

    ; res@gsnLeftStringFontHeightF = 0.02

    res@vcRefAnnoString1On = True
    res@vcRefAnnoString1 = "0.2 m/s"

    res@vcRefAnnoFontHeightF = 0.013
    res@vcRefAnnoOrthogonalPosF  = -1.2              ; move ref vector down -1.0在图的右上
    res@vcRefAnnoParallelPosF = 1.0

    res@vcVectorDrawOrder="PostDraw"

    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False


    res@gsnRightString = ""

    ; res@gsnLeftString = wksstring(i)
    ; res@gsnLeftString = timeday(i)
    res@gsnLeftStringFontHeightF = 0.0185
    res@gsnLeftStringOrthogonalPosF = -0.1

    ; res@gsnRightString = "OND mean"
    ; res@gsnRightStringFontHeightF = 0.0185

    ; plot(i) = gsn_csm_vector_scalar_map(wks, U_diff, V_diff,PS_diff, res)
    ; plot = gsn_csm_vector_scalar_map(wks, US_diff,VS_diff,T_diff, res)
    plots(0) = gsn_csm_contour_map(wks, US, res)
    plots(1) = gsn_csm_contour_map(wks, US_2, res)

pres =  True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plots, (/2,1/), pres)
  

end 