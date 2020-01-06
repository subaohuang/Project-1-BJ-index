;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
begin
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h0.0101-4012_chazhi.nc"
  ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindianwarming/F_2000_indian.cam.h0.0101-4012_chazhi.nc"  
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addpacifictrend/F_2000_pacific.cam.h0.0101-4012_chazhi.nc"  
    ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addatlanticwarming/F_2000_atlantic.cam.h0.0101-4012_chazhi.nc" 
  

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.601.and.date.le.4012)
    


      PRECL_CT  =  f1->PRECL(timeind,:,:)
      PRECC_CT  =  f1->PRECC(timeind,:,:)
     
      PRE_CT = PRECC_CT
      PRE_CT = (/PRECC_CT+PRECL_CT/)
      PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day

   
      
      PRECL_id  =  f2->PRECL(timeind,:,:)
      PRECC_id  =  f2->PRECC(timeind,:,:)
     
      PRE_id = PRECC_id
      PRE_id = (/PRECC_id+PRECL_id/)
      PRE_id = (/PRE_id*1000*3600*24/)  ;; m/s conver to mm/day

     PRE_CT_clm = dim_avg_n_Wrap(PRE_CT, 0)
     PRE_indian_clm = dim_avg_n_Wrap(PRE_id, 0)

     PRE_diff = PRE_CT_clm
     PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)

   
   
     wks=gsn_open_wks("pdf", "CESM-F2000-p_pattern-pacific")
    ; wks=gsn_open_wks("pdf", "CESM-slp_u_v_pattern-F-pday-rmindian-rmalt-rmallocean")
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=True
    res@cnFillOn = True 
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    res@cnFillPalette = "BlueWhiteOrangeRed"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -3
    res@cnMaxLevelValF       =   3
    res@cnLevelSpacingF      =  0.5
    
    ; res@cnMinLevelValF       =   -16
    ; res@cnMaxLevelValF       =   16
    ; res@cnLevelSpacingF      =  2

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF          = 180.0
    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    


    res@mpShapeMode = "FreeAspect"
    res@vpWidthF = 0.7
    res@vpHeightF = 0.3
;
    res@gsnRightString = ""
    plot = gsn_csm_contour_map(wks, PRE_diff, res)
    
    draw(plot)
    frame(wks)
    print("done script")
  end
  

