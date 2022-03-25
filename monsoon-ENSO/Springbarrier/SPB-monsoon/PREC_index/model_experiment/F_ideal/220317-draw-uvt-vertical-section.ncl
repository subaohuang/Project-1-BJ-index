
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

  data_clm_jjas = dim_avg_n_Wrap(data_clm(5:8,:,:,:), 0)

return data_clm_jjas
end 

begin
 
; ==============================================================
; Mission 1. read data
; ==============================================================

;;; read data from model intial file 

  fils = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703.nc") ;;;daily数据计算的q1

  f3 = addfile(fils, "r")

  US = f3->US(0,:,:,:)
  VS = f3->VS(0,:,:,:)
  T = f3->T(0,:,:,:)
  printVarSummary(US)
  printVarSummary(VS)
;;;;

;;; read data from model intial file 

  fils2 = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703_JJASMEAN.nc") ;;;daily数据计算的q1
  ; fils2 = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703_era5.nc") ;;;daily数据计算的q1

  f4 = addfile(fils2, "r")

  US_2 = f4->US(0,:,:,:)
  VS_2 = f4->VS(0,:,:,:)
  T_2 = f4->T(0,:,:,:)
  printVarSummary(US_2)
  printVarSummary(VS_2)
;;;;

; ;; read data from model intial file 

;   fils2 = systemfunc("ls F_dynamic_f19_ctrl_1year.cam.h0.0001-08.nc") ;;;daily数据计算的q1
;   ; fils2 = systemfunc("ls /home/ys17-23/lsh/Project/monsoon-ENSO/Springbarrier/SPB-monsoon/PREC_index/model_experiment/F_ideal/cami_0000-01-01_1.9x2.5_L30_c070703_era5.nc") ;;;daily数据计算的q1

;   f4 = addfile(fils2, "r")

;   US_2 = f4->U(0,:,:,:)
;   VS_2 = f4->V(0,:,:,:)
;   T_2 = f4->T(0,:,:,:)
;   printVarSummary(US_2)
;   printVarSummary(VS_2)
; ;;;;

;   T_diff = T_2
;   T_diff = (/T_2 - T/)

;   US_diff = US_2
;   US_diff = (/US_2 - US/)

;   VS_diff = VS_2
;   VS_diff = (/VS_2 - VS/)

; printMinMax(US_diff, 1)
; printMinMax(VS_diff, 1)
; printMinMax(T_diff, 1)

;;;;;;
US_mean = dim_avg_n_Wrap(US, 2)
VS_mean = dim_avg_n_Wrap(VS, 2)
T_mean = dim_avg_n_Wrap(T, 2)
  

US_2_mean = dim_avg_n_Wrap(US_2, 2)
VS_2_mean = dim_avg_n_Wrap(VS_2, 2)
T_2_mean = dim_avg_n_Wrap(T_2, 2)

; T_mean = (/T_2_mean - T_mean/)
; T_2_mean = T_mean

; US_mean = (/US_2_mean - US_mean/)
; US_2_mean = US_mean


plots = new(2,graphic)

  wks = gsn_open_wks("pdf", "cami-vertical-section")
;;;;

 ;set up plotting resources
 res     = True
 res@gsnDraw = False ; do not draw until the end
 res@gsnFrame = False ; do not automatically advance 'frame'

 res@trYReverse  = True ; reverse Y axis

 res@cnFillOn = True ; color fill
 res@lbLabelBarOn = False

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level


  ; res@cnMinLevelValF       =   -50.
  ; res@cnMaxLevelValF       =   50.
  ; res@cnLevelSpacingF      =  10.

 ; res@cnLevelSelectionMode = "ExplicitLevels"
 ; res@cnLevels = flevels
 ; contour every second line
 res@cnMonoLevelFlag=False
 ; res@cnLevelFlags = lconts
 res@gsnContourZeroLineThicknessF=0. ;don't plot zero line
 res@gsnContourNegLineDashPattern = 2
 res@cnLineThicknessF=2.

 ;contour labelling
 res@cnLineLabelsOn = True
 res@cnExplicitLineLabelsOn=False

;controlling label bar
  res@lbLabelAutoStride   = True         ; Clean up labelbar labels.
  res@lbBoxLinesOn        = True         ; labelbar box lines.
  res@lbTitlePosition     = "Bottom"
  res@lbTitleFontHeightF = 0.02
  res@lbLabelFontHeightF = 0.015
  res@pmLabelBarHeightF = 0.1
  res@lbBoxMinorExtentF = 0.2
  res@lbBoxSeparatorLinesOn=False

  plots(0) = gsn_csm_contour(wks,T_mean,res)
  plots(1) = gsn_csm_contour(wks,T_2_mean,res)

  ;   plots(0) = gsn_csm_contour(wks,US_mean,res)
  ; plots(1) = gsn_csm_contour(wks,US_2_mean,res)

pres =  True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plots, (/2,1/), pres)
  

end 