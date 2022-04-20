;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
begin
  fils1 = "/home/yangsong3/data-observation/linshh/data/wc-result/msf_HC_F_2000_CTRL_0101-4012.nc"  
  ; fils1 = "/home/yangsong3/data-observation/linshh/data/wc-result/msf_HC_FC5_CTRL_0101-4012.nc"  
  ; fils1 = "/home/yangsong3/data-observation/linshh/data/wc-result/msf_F_2000_allocean-tropical_0101-4012.nc"  
  ; fils2 = "/home/yangsong3/data-observation/linshh/data/wc-result/msf_HC_lsh_F2000_all_sstano_0101-4012.nc"  
  fils2 = "/home/yangsong3/data-observation/linshh/data/wc-result/msf_HC_lsh_F2000_tro_sstano_0101-4012.nc"
  ; fils2 = "/home/yangsong3/data-observation/linshh/data/wc-result/msf_HC_lsh_FC5_all_sstano_0101-4012.nc"  

  f1=addfile(fils1, "r")
  f2=addfile(fils2, "r")

  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.601.and.date.le.3012)
  
  

  ; msf_time_ctrl= dim_avg_n_Wrap(f1->msf(timeind,{1000:100},:,{120:270}),3)

  ; msf_time_addindian= dim_avg_n_Wrap(f2->msf(timeind,{1000:100},:,{120:270}),3)

   msf_time_ctrl= dim_avg_n_Wrap(f1->msf(timeind,{1000:100},:,:),3)

   msf_time_addindian= dim_avg_n_Wrap(f2->msf(timeind,{1000:100},:,:),3)


  ; msf_time1 = lonFlip(f1->msf(timeind,{1000:100},:,:))
  ; msf_time_ctrl = dim_avg_n_Wrap(msf_time1(:,:,:,{-80:20}),3)

  ; msf_time2 = lonFlip(f2->msf(timeind,{1000:100},:,:))
  ; msf_time_addindian = dim_avg_n_Wrap(msf_time2(:,:,:,{-80:20}),3)


  msf_time_ctrl&lev@units = "hpa"
  msf_time_ctrl&lev@long_name = "pressure level"

  msf_time_addindian&lev@units = "hpa"
  msf_time_addindian&lev@long_name = "pressure level"
 


  msf_time_ctrl=(/msf_time_ctrl/10^10/)
  msf_time_addindian=(/msf_time_addindian/10^10/)


  printMinMax(msf_time_ctrl, 1)

  
  msf_clm_ctrl = dim_avg_n_Wrap(msf_time_ctrl, 0)
  msf_clm_addindian = dim_avg_n_Wrap(msf_time_addindian, 0)
  
  msf_annual_1 = msf_clm_ctrl
  msf_annual_2 = msf_clm_addindian

  msf_diff = msf_annual_1
  msf_diff = (/msf_clm_addindian - msf_clm_ctrl/)

  printMinMax(msf_diff, 1)






    ;;;;画图
    wks = gsn_open_wks("pdf", "F2000-HC-msf-diff-tro-sstano")
    ;;;;;draw 
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
    ;  res@cnMinLevelValF = -0.8
    ; res@cnMaxLevelValF = 0.8
    ; res@cnLevelSpacingF = 0.1

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
    plot = gsn_csm_pres_hgt(wks,msf_diff,res)

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

 
    plot_clm = gsn_csm_pres_hgt(wks, msf_clm_ctrl, res)

    overlay(plot, plot_clm)


    draw(plot)
    frame(wks)
  end 

