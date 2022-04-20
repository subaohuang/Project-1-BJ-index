;--------------------------------------------------------
; DKRZ NCL example: contour_plot.ncl
;
; Description:      Create a contour plot in one frame, 
;                   title and labelbar in another frame each.
;
;                    This script is part of the Pseudo_3D_plots
;                    examples to generate 3D looking images.
;
; NCL Version:        6.3.0
;
; 22.07.16  meier-fleischer(at)dkrz.de
;--------------------------------------------------------
begin

  title = ""
  
  f   = addfile("/home/ys17-19/lsh/data/sst/HadISST_sst.nc","r")                             ;-- open input file
  date=cd_calendar(f->time, 1)
  timeind:=ind(date.le.201612.and.date.ge.197901)
  var = lonFlip(f->sst(timeind,:,:))                                       ;-- read variable first time step
  

  sst_clm = dim_avg_n_Wrap(var, 0)

  ;;;-- contour plot workstation
  wks_type          = "png"
  wks_type@wkWidth  =  2500
  wks_type@wkHeight =  2500
  wks = gsn_open_wks(wks_type, "plot_contour") 


;-- title plot workstation
  wkst_type          = "png"
  wkst_type@wkWidth  =  2500
  wkst_type@wkHeight =  2500
  wkst = gsn_open_wks(wkst_type, "plot_title") 
  
;-- labelbar plot workstation
  wksl_type          = "png"
  wksl_type@wkWidth  =  2500
  wksl_type@wkHeight =  2500
  wksl = gsn_open_wks(wksl_type, "plot_labelbar") 
  


;;**************************************************************
;画图
;;***************************************************************
  plot=new(1,graphic)
  res=True
 
  res@gsnLeftStringFontHeightF=20.0
  
  ; res@cnFillPalette = "NCV_jaisnd"
 ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillOn=True
  res@cnLinesOn= False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False

  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 70

  res@mpMaxLonF             = 300
  ; res@tmXBLabelFontThicknessF  = 0.05

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=False
  res@tmYROn=False
  res@tmYLOn =False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"
  
  res@tmYLLabelsOn = False
  res@tmXBLabelsOn = False

  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   12.0
  res@cnMaxLevelValF       =   30.0
  res@cnLevelSpacingF      =  2.0

  contour = gsn_csm_contour_map(wks, sst_clm, res)


  ;-------------------------------
;-- Title plot:
;-------------------------------
  tires                         =  True
  tires@txJust                  = "CenterCenter"            ;-- text justification
  tires@txFontHeightF           =  0.017                    ;-- text font size

  tix = 0.5
  tiy = 0.8
  
  gsn_text_ndc(wkst,title,tix,tiy,tires)                    ;-- draw title
  print("       - title plot done")

  frame(wkst)                                               ;-- advance title frame
  
;-------------------------------
;-- Labelbar plot:
;-------------------------------
;-- retrieve the contour resources for the labelbar
  getvalues contour@contour
     "cnLevels"             :   cnlevs                      ;-- contour levels
     "cnFillColors"         :   cncols                      ;-- contour colors
     "cnMonoFillPattern"    :   mpattern                    ;-- contour fill pattern
  end getvalues

  labels = ""+cnlevs                                        ;-- labelbar labels
  nboxes =  dimsizes(cnlevs)+1                              ;-- number of color boxes
 
  lbx = 0.3                                                 ;-- labelbar x-position
  lby = 0.9                                                 ;-- labelbar y-position
  lbw = 0.2                                                 ;-- labelbar width
  lbh = 0.8                                                 ;-- labelbar height

;-- labelbar resources
  lbres                         =  True
  lbres@lbAutoManage            =  False                    ;-- we want to control the labelbar
  lbres@lbOrientation           = "vertical"                ;-- labelbar orientation
  lbres@lbLabelFontHeightF      =  0.024                    ;-- increase label font size
  lbres@lbLabelFontThicknessF   =  4.                       ;-- increase label font thickness
  lbres@lbLabelStride           =  1                        ;-- label every value
  lbres@lbPerimOn               =  False                    ;-- no box around labelbar
  lbres@lbLabelFontColor        = "black"                   ;-- labe font color
  lbres@lbLabelAlignment        = "InteriorEdges"           ;-- where to draw the labels
  lbres@lbLabelOffsetF          =  0.05                     ;-- move labels to the right
  lbres@lbMonoFillPattern       =  mpattern                 ;-- use same fill pattern as contour plot
  lbres@lbMonoFillColor         =  False                    ;-- use multiple colors
  lbres@lbFillColors            =  cncols                   ;-- contour colors
  lbres@lbLabelStrings          =  labels                   ;-- contour labels

  lbres@vpWidthF                =  lbw                      ;-- labelbar viewport width
  lbres@vpHeightF               =  lbh                      ;-- labelbar viewport height
   
  gsn_labelbar_ndc(wksl,nboxes,labels,lbx,lby,lbres)        ;-- draw labelbar

;-- text resources labelbar annotations
  txres = True
  txres@txFontHeightF          = 0.018                      ;-- units text font size
  
  gsn_text_ndc(wksl, "[C~S~o~N~]", 0.440, 0.12, txres)      ;-- add units to labelbar

  print("       - labelbar plot done")

  frame(wksl)                                               ;-- advance labelbar frame

end

