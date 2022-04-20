;;********************************************************
;;画散点图,得到关系
;;
;;lsh 2020-11-04
;**********************************************************
; function chazhi(data)
;  begin
;    xi = data&longitude
;    yi = data&latitude
;    xo = fspan(0.5,359.5,360)
;    yo = fspan(-89.5,89.5,180)
 
;    data_chazhi = linint2_Wrap(xi, yi, data, True, xo, yo, 0)
  
;    data_chazhi!1= "lat"
;    data_chazhi!2= "lon"
   
;    data_chazhi&lat = yo
;    data_chazhi&lon = xo
;    copy_VarAtts(data, data_chazhi)
 
;    return data_chazhi
; end 
begin
  ;;;read data
  lat1 = -40
  lat2 =  40
  lon1 =  30
  lon2 = 290

  ; startmon = 198201
  ; endmon = 201912

  startmon = 197901
  endmon = 201512

  f1 = addfile("/home/yangsong3/data-observation/NOAA-OLR/olr.mon.mean.197406-201912.nc", "r")
  date:=cd_calendar(f1->time, 1)
  timeind :=ind(date.ge.startmon.and.date.le.endmon)
  olr = short2flt(f1->olr(timeind,{lat1:lat2},{lon1:lon2})) ;;;2.5 ×2.5 grid

  printVarSummary(olr)



  f2 = addfile("/home/yangsong3/data-observation/linshh/data/sst/HadISST_sst.nc", "r")
  date:=cd_calendar(f2->time, 1)

  timeind:=ind(date.le.201512.and.date.ge.197901)

  sst_r =  lonFlip(f2->sst(timeind,:,:))

  sst = sst_r(:,{lat1-1:lat2+1},{lon1-1:lon2+1})

  printVarSummary(sst)





;;;;chazhi 
   print("**********start interp************")
    ; xi = sst&lon
    ; yi = sst&lat

    xi = sst&longitude
    yi = sst&latitude

    xo = olr&lon
    yo = olr&lat



    sst_chazhi = linint2_Wrap(xi, yi, sst, True, xo, yo, 0)

    sst_chazhi!1= "lat"
    sst_chazhi!2= "lon"

    sst_chazhi&lat = yo
    sst_chazhi&lon = xo
;;;;
   
   
;;;;
    printVarSummary(sst_chazhi)
     olr = where(ismissing(sst_chazhi), olr@_FillValue, olr)
;;;;




;;;;Climatology bin
    ; sst_clm = dim_avg_n_Wrap(sst_chazhi(:,{0:20},{125:160}), 0)   
    sst_clm = dim_avg_n_Wrap(month_to_season(sst_chazhi(:,{0:20},{125:160}),"MAM"), 0)   
    olr_clm = dim_avg_n_Wrap(olr(:,{0:20},{125:160}), 0)


    sst_clm_one = ndtooned(sst_clm)
    olr_clm_one = ndtooned(olr_clm)

;;;


;;;;anomaly 

    olr = (/rmMonAnnCycTLL(olr)/)
    sst_chazhi = (/rmMonAnnCycTLL(sst_chazhi)/)
;;;

;;;;Monthly bin

   ; sst_chazhi = (/runave_n_Wrap(sst_chazhi, 3, 0, 0)/)
   ; olr = (/runave_n_Wrap(olr, 3, 0, 0)/)

   dim = dimsizes(sst_chazhi)
   
   sst_chazhi_june = sst_chazhi(3:dim(0)-1:12,:,:)
   sst_chazhi_july = sst_chazhi(4:dim(0)-1:12,:,:)
   sst_chazhi_aug = sst_chazhi(5:dim(0)-1:12,:,:)

   sst_chazhi_jj =  array_append_record(sst_chazhi_june,sst_chazhi_july,0)
   sst_chazhi_jja = array_append_record(sst_chazhi_jj , sst_chazhi_aug,0)


    olr_june = olr(3:dim(0)-1:12,:,:)
    olr_july = olr(4:dim(0)-1:12,:,:)
    olr_aug = olr(5:dim(0)-1:12,:,:)

    olr_jj =  array_append_record(olr_june,olr_july,0)
    olr_jja = array_append_record(olr_jj , olr_aug,0)
    
   ;  sst_one = ndtooned(sst_chazhi_jja(:,{0:15},{125:160}))
   ;  olr_one = ndtooned(olr_jja(:,{0:15},{125:160}))

    ; sst_one = ndtooned(sst_chazhi(:,{0:20},{125:160}))
;;;;


;;;;计算单点回归系数

    corr = escorc_n(sst_chazhi_jja, olr_jja, 0, 0)                    ;;计算相关性
    copy_VarCoords(sst_chazhi(1,:,:), corr)
    printVarSummary(corr)
    std1 = dim_stddev_n(olr_jja, 0)
    std2 = dim_stddev_n(sst_chazhi_jja, 0)


    rc = corr

    ; rc = (std1/std2)*corr
;;;
wks = gsn_open_wks("pdf", "单点回归-sst-OLR-MAM")
;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnAddCyclic = False
  
 
   res@cnFillPalette = "BlueWhiteOrangeRed"
   ; res@cnFillPalette = "BlueDarkRed18"

  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn= True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpFillOn  =False
  res@mpFillDrawOrder = "PostDraw"
  res@mpMinLatF             = -5
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 110
  res@mpMaxLonF             = 170

  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 30
  ; res@mpMinLonF             = 100
  ; res@mpMaxLonF             = 170

  ; res@mpMinLatF             = 0
  ; res@mpMaxLatF             = 15
  ; res@mpMinLonF             = 125
  ; res@mpMaxLonF             = 160

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"



  res@tmXBLabelFontHeightF=0.025
  res@tmYLLabelFontHeightF=0.025
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0

  ; res@tmXBValues = (/100,120,140,160/) 
  ; res@tmXBLabels = (/100,120,140,160/) +""+"E" 

  res@gsnLeftStringFontHeightF = 0.03


  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -10.0
  ; res@cnMaxLevelValF       =   10.0
  ; res@cnLevelSpacingF      =  2.0
  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.5
  res@cnMaxLevelValF       =   0.5
  res@cnLevelSpacingF      =  .1

  res@cnInfoLabelOn = False
  

  res@gsnLeftString="CORR between SST and OLR in each point"
  res@gsnRightString = ""


  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

  plot=new(1,graphic)


  ; plot(0) = gsn_csm_contour_map(wks, olr_rela_mam({-0.5:15.5},{124.5:160.5}), res)
  ; res@gsnLeftString="(b) JJA OLR change related to SST change"
  ; plot(1) = gsn_csm_contour_map(wks, olr_rela_jja({-0.5:15.5},{124.5:160.5}), res)
  ; res@gsnLeftString="(c) SON OLR change related to SST change"
  ; plot(2) = gsn_csm_contour_map(wks, olr_rela_son({-0.5:15.5},{124.5:160.5}), res)
  ; res@gsnLeftString="(d) DJF OLR change related to SST change"
  ; plot(3) = gsn_csm_contour_map(wks, olr_rela_djf({-0.5:15.5},{124.5:160.5}), res)

  plot = gsn_csm_contour_map(wks,rc, res)


draw(plot)
frame(wks)
end
