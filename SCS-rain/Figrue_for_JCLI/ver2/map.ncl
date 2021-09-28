;;********************************************************
;;画GPCP的降水趋势
;;
;;lnshuheng 2019.12.24
;**********************************************************
begin
   
wks   = gsn_open_wks("pdf", "map")
;;**************************************************************
;画图
;;***************************************************************
      
  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -40
  res@mpMaxLatF             = 40
  res@mpMinLonF             = -80
  res@mpMaxLonF             = 280

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmXBValues = (/-60,0,60,120,180,240/) 
  res@tmXBLabels = (/"60W","0","60E","120E","180","120W"/) 

  res@tmXBLabelFontHeightF=0.01
  res@tmYLLabelFontHeightF=0.01
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmYLMinorOn = True
  
  ; res@tmXBValues = (/100,120,140,160/) 
  ; res@tmXBLabels = (/100,120,140,160/) +""+"E" 


  ; res@mpShapeMode="FreeAspect"
  ; res@vpHeightF=0.3
  ; res@vpWidthF=0.7
  plot = new(2,graphic)
  res@gsnLeftString = "(a) Boreal Spring"
  res@gsnLeftStringFontHeightF = 0.015
  plot(0) = gsn_csm_map(wks, res)
  res@gsnLeftString = "(b) Other seasons"
  plot(1) = gsn_csm_map(wks, res)
  

  pres = True
  gsn_panel(wks, plot, (/2,1/), pres)
  

end