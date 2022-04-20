;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
begin
    fils1 = "/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_ERA5_197901-202012.nc" 
    ; fils1 = "/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_NCEP1_197901-202012.nc" 
        ; fils1 = "/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_NCEP1_197901-202012-test.nc"
    ; fils1 = "/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_vdmethod_NCEP1_197901-202012-test.nc"
    fils2 = "/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_NCEP1_197901-202012-test.nc"
    ; fils1 = "/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_vdmethod_ERA5_197901-202012.nc"
    ; fils2 = "/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_ERA5_197901-202012.nc"
        ; fils2 = "/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_NCEP1_197901-202012-test.nc"
    ; fils1 = "/home/yangsong3/data-observation/linshh/obs-derived/msf/msf_HC_vdmethod_ERA5_197901-202012.nc"  

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f1->time, 1)

    timestart = 197901
    timeend = 202012
    timeind = ind(date.ge.timestart.and.date.le.timeend)

    msf1=dim_avg_n_Wrap(f1->msf(timeind,{1000:100},:,:),3)

    msf2 = dim_avg_n_Wrap(f2->msf(timeind,{1000:100},:,:),3)
    
    
    ; msf_time = lonFlip(f1->msf(timeind,{1000:100},:,:))
    ; msf_time_ctrl = dim_avg_n_Wrap(msf_time(:,:,:,{-80:20}),3)



    msf1 = (/msf1/10^10/)
    msf2 = (/msf2/10^10/)
    

    msf1_annual = month_to_annual(msf1, 1)
    msf2_annual = month_to_annual(msf2, 1)

    msf_clm = dim_avg_n_Wrap(msf1_annual, 0)

    
    size = dimsizes(msf1_annual)
    yearnum = size(0)

    msf1_rc = msf_clm
    msf1_rc = (/regCoef_n(ispan(1,yearnum,1), msf1_annual , 0, 0)/)
    msf2_rc = (/regCoef_n(ispan(1,yearnum,1), msf2_annual , 0, 0)/)

    msf1_rc  = msf1_rc*yearnum
    msf2_rc  = msf2_rc*yearnum

    ; msf_diff = msf1_rc
    ; msf_diff = (/msf2_rc - msf1_rc/)

    ; msf_diff = (/dim_avg_n_Wrap(msf1_annual, 0) - dim_avg_n_Wrap(msf2_annual, 0)/)
;


    
    ;;;;画图
     wks = gsn_open_wks("pdf", "F_2000-msf-change-ERA5-vdmethod-"+timestart + "-"+ timeend)
     ; wks = gsn_open_wks("pdf", "F_2000-msf-change-vd-method-ERA5-"+timestart + "-"+ timeend)

    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLinesOn=True

    res@cnFillOn =True

    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -1.
    res@cnMaxLevelValF = 1.
    res@cnLevelSpacingF = 0.2
    ;  res@cnMinLevelValF = -0.1
    ; res@cnMaxLevelValF = 0.1
    ; res@cnLevelSpacingF = 0.01

    res@lbLabelBarOn =True
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    res@cnLineThicknessF = 2.0

    res@tmXBMajorThicknessF=1.5
    res@tmYLMajorThicknessF=1.5
    res@tmBorderThicknessF=1.5

    res@tmXBLabelFontHeightF=0.02
    res@tmYLLabelFontHeightF=0.02

    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(-30,30,5)


    res@tmYRMode="Automatic"

    res@vpHeightF=0.3
    res@vpWidthF=0.7

    ;;;;
  
    res@gsnLeftStringFontHeightF = 0.02
    plot = gsn_csm_pres_hgt(wks,msf1_rc,res)

    ;;;画气候态的等值线
    res@gsnRightString=""
    res@gsnLeftString=""
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn=True
    res@cnLevels = fspan(-10,10,11)
    res@cnLineLabelsOn=True
    res@cnLineLabelBackgroundColor="Transparent"
    res@cnLineThicknessF = 2.
    res@cnFillOn = False
    res@gsnContourNegLineDashPattern=1
    res@gsnContourZeroLineThicknessF=2.2

 
    plot_clm = gsn_csm_pres_hgt(wks, msf_clm, res)

    overlay(plot, plot_clm)
  draw(plot)
  frame(wks)
  end 

