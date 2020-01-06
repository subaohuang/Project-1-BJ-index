;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
begin
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h0.0101-4012_chazhi.nc"
  ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindianwarming/F_2000_indian.cam.h0.0101-4012_chazhi.nc"  
  ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addpacifictrend/F_2000_pacific.cam.h0.0101-4012_chazhi.nc"  
    ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addatlanticwarming/F_2000_atlantic.cam.h0.0101-4012_chazhi.nc" 
      fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addallocean_tropical/F_2000_allocean.cam.h0.0101-4012.nc"
  

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

      ;;;ttest 
      simple_size =  PRE_CT_clm
      simple_size = dimsizes(PRE_CT&time)
    
      simple_variance1  = dim_variance_n_Wrap(PRE_id, 0)
      simple_variance2  = dim_variance_n_Wrap(PRE_CT, 0)

      prob_p = ttest(PRE_indian_clm ,simple_variance1, simple_size, PRE_CT_clm ,simple_variance2 ,simple_size,False,False)
     
     copy_VarCoords(PRE_CT_clm, prob_p)
     
   
   
     wks=gsn_open_wks("pdf", "CESM-F2000-p_pattern-allocean-SCS-philipieng region")
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

    res@cnFillPalette="MPL_BrBG"

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
    res@mpMinLatF             = -5
    res@mpMaxLatF             = 45
    res@mpMinLonF             =90
    res@mpMaxLonF             = 160
    res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    res@gsnLeftString  = "AGCM forced with Lanina-like SST"
    


;     res@mpShapeMode = "FreeAspect"
;     res@vpWidthF = 0.7
;     res@vpHeightF = 0.3
; ;
   

    res@gsnRightString = ""
    plot = gsn_csm_contour_map(wks, PRE_diff, res)
    



    ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   str3=unique_string("dum3")
   str4=unique_string("dum4")
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
    plot@$str3$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres)
   

     ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.002
  
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
   
    res_pro@gsnCenterString = ""
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
    iplot = new(1,graphic)

    iplot(0)=gsn_csm_contour(wks,prob_p,res_pro)

    i = 0
    overlay(plot(i), iplot(i))


    draw(plot)
    frame(wks)
    print("done script")
  end
  

 