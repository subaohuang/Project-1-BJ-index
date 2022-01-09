
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"


begin 

modelname = "B2000_decoupling_indian_ocean"
timestart = "251"
timeend = "280"


; f1 = systemfunc("ls /WORK/sysu_hjkx_ys/linshh/B/decouple/"+modelname+"/run/"+modelname+".pop.h.0{251..280}"+"*"+".nc")
f1 = systemfunc("ls /WORK/sysu_hjkx_ys/linshh/B/decouple/"+modelname+"/run/"+modelname+".pop.h.0{"+timestart+".."+timeend+"}"+"*"+".nc")


; modelname = "B2000_f19g16_CTRL"
; f1 = systemfunc("ls /WORK/sysu_hjkx_ys/linshh/data/"+modelname+"/"+modelname+".pop.h.0{251..280}"+"*"+".nc")
TEMP_file = addfiles (f1, "r")


TLONG = TEMP_file[1]->TLONG
TLAT = TEMP_file[1]->TLAT
TEMP = TEMP_file[:]->TEMP(:,0,:,:)

TEMP@lat2d = TLAT
TEMP@lon2d = TLONG
;;;

;;;差值TEMP 
  sst_sen = PopLatLon(TEMP,"gx1v6","fv1.9x2.5","bilin","da","090206")
  printVarSummary(sst_sen)
;;;;

out = addfile(modelname+"_sst_0251-0280_fvgrid.nc", "c")
out->TEMP = sst_sen


modelname = "B2000_f19g16_CTRL"
f1 = systemfunc("ls /WORK/sysu_hjkx_ys/linshh/data/"+modelname+"/"+modelname+".pop.h.0{251..280}"+"*"+".nc")


TEMP_file = addfiles (f1, "r")


TLONG = TEMP_file[1]->TLONG
TLAT = TEMP_file[1]->TLAT
TEMP_ctrl = TEMP_file[:]->TEMP(:,0,:,:)

TEMP_ctrl@lat2d = TLAT
TEMP_ctrl@lon2d = TLONG


;;;差值TEMP 
    sst_ctrl = PopLatLon(TEMP_ctrl,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst_ctrl)
;;;;

out = addfile(modelname+"_sst_0251-0280_fvgrid.nc", "c")
out->TEMP = sst_ctrl

exit()

sst_ctrl = (/rmMonAnnCycTLL(sst_ctrl)/)


nino34_ctrl = dim_avg_n_Wrap(sst_ctrl(:,{-5:5},{190:240}), (/1,2/))


sst_ctrl_djf = month_to_season(sst_ctrl, "JJA")

nino34_ctrl_djf = month_to_season(nino34_ctrl, "DJF")




sst_sen = (/rmMonAnnCycTLL(sst_sen)/)


nino34_sen = dim_avg_n_Wrap(sst_sen(:,{-5:5},{190:240}), (/1,2/))


sst_sen_djf = month_to_season(sst_sen, "JJA")

nino34_sen_djf = month_to_season(nino34_sen, "DJF")


dim = dimsizes(nino34_sen_djf)


sst_ind_reg =  regCoef_n(nino34_sen_djf(0:dim-2), sst_sen_djf(1:dim-1,:,:) , 0, 0)
copy_VarCoords(sst_sen(1,:,:), sst_ind_reg)

sst_ctrl_reg =  regCoef_n(nino34_ctrl_djf(0:dim-2), sst_ctrl_djf(1:dim-1,:,:), 0, 0)
copy_VarCoords(sst_ctrl(1,:,:), sst_ctrl_reg)

;;;
wks = gsn_open_wks("pdf", "regpattern-indian_sst")
;;**************************************************************
;画图
;;***************************************************************
   

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
; 
  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "BlueRed"

  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 30
  res@mpMaxLonF             = 300

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

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -1.0
  res@cnMaxLevelValF       =   1.0
  res@cnLevelSpacingF      =  0.2


  res@vcGlyphStyle="CurlyVector"
  res@vcMonoLineArrowColor = True
  res@vcLineArrowColor = "black"
  res@vcMinDistanceF           = 0.02          ; thin vectors




  res@gsnLeftString = "Heating"

  plot = new(2,graphic)
  plot(0) =gsn_csm_contour_map(wks, sst_ind_reg, res)
  plot(1) =gsn_csm_contour_map(wks, sst_ctrl_reg, res)
  ; plot(0) = gsn_csm_vector_scalar_map(wks, uwnd_diff_jja, vwnd_diff_jja,p_diff_jja, res)
  
  ; res@gsnLeftString = "Cooling"
  ; plot(1) = gsn_csm_vector_scalar_map(wks, uwnd_diff, vwnd_diff,, res)

  pres = True
  gsn_panel(wks, plot, (/2,1/), pres)
end