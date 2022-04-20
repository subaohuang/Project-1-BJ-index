;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
begin
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL/F_2000_CTRL.cam.h1.0101-4012_local_wk_hc_500.nc"
  fils2 = "/home/ys17-19/lsh/CESM-data/F/lsh_F2000_ind_sstano/lsh_F2000_ind_sstano.cam.h1.0101-4012_local_wk_hc_500.nc"  


  f1=addfile(fils1, "r")
  f2=addfile(fils2, "r")

  date=cd_calendar(f1->time, 1)
  timeind=ind(date.ge.601.and.date.le.3012)
 
  wk_ctrl  =  f1->m_zonal(timeind,:,:)
  hc_ctrl  =  f1->m_meri(timeind,:,:)
  
  wk_ind  =  f2->m_zonal(timeind,:,:)
  hc_ind  =  f2->m_meri(timeind,:,:)

  
  wk_ctrl_clm = dim_avg_n_Wrap(wk_ctrl, 0)
  hc_ctrl_clm = dim_avg_n_Wrap(hc_ctrl, 0)

  wk_ind_clm = dim_avg_n_Wrap(wk_ind, 0)
  hc_ind_clm = dim_avg_n_Wrap(hc_ind, 0)


  wk_diff = wk_ctrl_clm
  wk_diff = (/wk_ind_clm - wk_ctrl_clm/)
  

  hc_diff = hc_ctrl_clm
  hc_diff = (/hc_ind_clm - hc_ctrl_clm/)

  

   wks=gsn_open_wks("pdf", "CESM-F2000-WK-HC-indian")

    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =True
    res@cnLineColor = "White"
    res@cnLineLabelInterval = 2
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette = "precip_diff_12lev"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -0.002
    res@cnMaxLevelValF       =   0.002
    res@cnLevelSpacingF      =  0.0005
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 30.0
    res@mpMaxLonF             = 160.0
    res@mpCenterLonF          = (/res@mpMinLonF + res@mpMaxLonF/)
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3

   
    res@gsnRightString = ""

    plot = new(2,graphic)
    plot(0) = gsn_csm_contour_map(wks, wk_diff,res)
    plot(1) = gsn_csm_contour_map(wks, hc_diff,res)
    
    
    pres = True
    pres@gsnPanelLabelBar =  True

    gsn_panel(wks, plot, (/2,1/), pres)

    print("done script")
  end
  

