load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl" 
;**********************************************************
begin
 ;;;time span
    startmon =25601
    endmon = 28012

    startyear_r = str_split_by_length(startmon, 3)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 3 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
   lat1 = -60
   lat2 = 60
   lon1 = 0
   lon2 = 360

;;;;
 

;;;;model name 
   modelname = "B2000_alt_north_year_CTRL"
;;;;


;;;;读入控制实验结果
    diri="/home/ys17-19/lsh/CESM-data/B/"+ modelname+"/"
    fils1=systemfunc("ls "+diri+""+modelname+".TEMP.h0.0251-0280.nc")
    f1=addfile(fils1, "r")
    date=cd_calendar(f1->time, 1)

    timeind=ind(date.le.endmon.and.date.ge.startmon)

    sst_r = f1->TEMP(timeind,0,:,:)  
   ;;;;差值TEMP 
    sst = PopLatLon(sst_r,"gx1v6","fv1.9x2.5","bilin","da","090206")
    printVarSummary(sst)
    

    sst_ct = sst(:,{lat1:lat2},{lon1:lon2})
    
;;;read NCEP-NCAR heat flux data
    fils3 = diri + modelname+".cam.h0.0251-0280_OLR_heatflux.nc"
    f3 = addfile(fils3, "r")

    date:=cd_calendar(f3->time, 1)

    timeind:=ind(date.le.endmon.and.date.ge.startmon)
  
    
    lhtfl_ct = f3->LHFLX(timeind,{lat1:lat2},{lon1:lon2})
    shtfl_ct = f3->SHFLX(timeind,{lat1:lat2},{lon1:lon2})
    nswrs_ct = f3->FSNS(timeind,{lat1:lat2},{lon1:lon2})
    nlwrs_ct = f3->FLNS(timeind,{lat1:lat2},{lon1:lon2})
;;;

    lhtfl_ct_clm = dim_avg_n_Wrap(lhtfl_ct, 0)
    shtfl_ct_clm = dim_avg_n_Wrap(shtfl_ct, 0)
    nswrs_ct_clm = dim_avg_n_Wrap(nswrs_ct, 0)
    nlwrs_ct_clm = dim_avg_n_Wrap(nlwrs_ct, 0)

;;;;
plot = new(4,graphic)

  wks   = gsn_open_wks("pdf", "B2000_HEATFLUX_CLM")
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlueWhiteOrangeRed"
  ; res@cnFillPalette = "NCV_jaisnd"
  ; res@gsnSpreadColorStart = 130
  ; res@gsnSpreadColorEnd = 253

  res@cnFillOn=True
  res@cnLinesOn=False
  res@cnLineLabelInterval = 2
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpFillDrawOrder = "PostDraw"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160
   ; lat1 = 0
   ; lat2 = 15
   ; lon1 = 110
   ; lon2 = 155

   lat1 = -60
   lat2 = 60
   lon1 = 0
   lon2 = 360

  res@mpMinLatF             = lat1
  res@mpMaxLatF             = lat2
  res@mpMinLonF             =lon1
  res@mpMaxLonF             = lon2

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0
  res@gsnRightString = ""

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -300.0
  res@cnMaxLevelValF       =   300.0
  res@cnLevelSpacingF      =  30.0


  res@gsnLeftString = "SST trend "
  plot(0) = gsn_csm_contour_map(wks, lhtfl_ct_clm, res)
  res@cnMinLevelValF       =   -50.0
  res@cnMaxLevelValF       =   50.0
  res@cnLevelSpacingF      =  10.0

  plot(1) = gsn_csm_contour_map(wks, shtfl_ct_clm, res)
  res@cnMinLevelValF       =   -300.0
  res@cnMaxLevelValF       =   300.0
  res@cnLevelSpacingF      =  30.0

  plot(2) = gsn_csm_contour_map(wks, nswrs_ct_clm, res)
  res@cnMinLevelValF       =   -100.0
  res@cnMaxLevelValF       =   100.0
  res@cnLevelSpacingF      =  10.0

  plot(3) = gsn_csm_contour_map(wks, nlwrs_ct_clm, res)

 pres = True
 gsn_panel(wks, plot, (/2,2/), pres)
 end

  ;;;添加赤道直线
