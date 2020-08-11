;;;************************************************************
;;
;; linshuheng 2019.10.22
;;;*********************************************************
function t_test(series1,series2)
begin
  ;;;ttest 
      avg1 = dim_avg_n_Wrap(series1, 0)
      avg2 = dim_avg_n_Wrap(series2, 0)
      
      dim = dimsizes(avg1)
      simple_size =  avg1
      simple_size = dim(0)
    
      simple_variance1  = dim_variance_n_Wrap(series1, 0)
      simple_variance2  = dim_variance_n_Wrap(series2, 0)

      prob_p = ttest(avg1 ,simple_variance1, simple_size, avg2 ,simple_variance2 ,simple_size,False,False)
     copy_VarCoords(avg1, prob_p)
return prob_p
end 
begin
  fils1 = "/home/ys17-19/lsh/CESM-data/F/F_2000_CTRL_new/F_2000_CTRL.cam.h0.0101-4012_chazhi.nc"
  ; fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addindianwarming/F_2000_indian.cam.h0.0101-4012_chazhi.nc"  
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_IPO/F_2000_ipo.cam.h0.0101-4012_chazhi.nc"
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

  ;;;

    ;;;显著性检验
    prob_pre = t_test(PRE_id,PRE_CT)
    copy_VarCoords(PRE_id(0,:,:), prob_pre)
   


  ;;;;

    ; wks=gsn_open_wks("pdf", "CESM-F2000-p_pattern-pacific-GHG")
    wks=gsn_open_wks("pdf", "CESM-F2000-p_pattern-pacific-IPO-SCS")
    
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

    res@cnFillPalette = "MPL_BrBG"

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -0.4
    res@cnMaxLevelValF       =   0.4 
    res@cnLevelSpacingF      =  0.1
    

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
    res@mpMinLatF             = -5
    res@mpMaxLatF             = 45
    res@mpMinLonF             =90
    res@mpMaxLonF             = 160

    ; res@mpMinLatF             = -60
    ; res@mpMaxLatF             = 60
    ; res@mpMinLonF             =0
    ; res@mpMaxLonF             = 360




    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    
;
    res@gsnRightString = ""
    plot = gsn_csm_contour_map(wks, PRE_diff, res)
    

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

    iplot(0)=gsn_csm_contour(wks,prob_pre,res_pro)

    i = 0
    overlay(plot(i), iplot(i))


  ; ;;添加区域
   ;;warm pool
   pgres = True
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
  
   str = (/"1","2","3"/)
   
   ;;;画box
   lat1 = 0
   lat2 = 15
   lon1 = 110
   lon2 = 155

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   i = 0
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   


    draw(plot)
    frame(wks)
    print("done script")
  end
  

