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
 
  fils2 = "/home/ys17-19/lsh/CESM-data/F/F_2000_addatlanticwarming/F_2000_atlantic.cam.h0.0101-4012_chazhi.nc"  
  
  

    f1=addfile(fils1, "r")
    f2=addfile(fils2, "r")
  
 

    date=cd_calendar(f1->time, 1)
    timeind=ind(date.ge.601.and.date.le.4012)


    PSL_CT  =  f1->PSL(timeind,:,:)
    PSL_CT = (/PSL_CT/100/)


    PSL_indian  =  f2->PSL(timeind,:,:)
    PSL_indian = (/PSL_indian/100/)


    PSL_CT_clm = dim_avg_n_Wrap(PSL_CT, 0)
    PSL_indian_clm = dim_avg_n_Wrap(PSL_indian, 0)
   

    PSL_diff_indian = PSL_CT_clm
    PSL_diff_indian = (/PSL_indian_clm - PSL_CT_clm/)

    

    ;;;t_test 
     prob_slp_indian = t_test(PSL_indian,PSL_CT)
     printVarSummary(prob_slp_indian)
    
     

    wks=gsn_open_wks("pdf", "CESM-F2000-slp-pattern-atlantic")
    
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True 
    res@cnLinesOn =True
    res@cnLineLabelsOn=False
    res@cnRasterSmoothingOn=True

    ; res@cnFillPalette = "hotcolr_19lev"


    gsn_define_colormap(wks, "hotcolr_19lev")
 


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res@cnMinLevelValF       =   -1
    ; res@cnMaxLevelValF       =   1
    ; res@cnLevelSpacingF      =  0.1
    
    res@cnMinLevelValF       =   -1.5
    res@cnMaxLevelValF       =   1.5
    res@cnLevelSpacingF      =  0.1
    
    ; res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11,0,0,16,18,20,21/)

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
     
    ; res@mpMinLatF             = -60.0
    ; res@mpMaxLatF             = 60.0
    ; res@mpMinLonF             = 0.0
    ; res@mpMaxLonF             = 360.0

    res@mpMinLatF             = -10.0
    res@mpMaxLatF             = 50.0
    res@mpMinLonF             = 60.0
    res@mpMaxLonF             = 160.0




    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    


    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3

    
    
    res@gsnRightString = ""
    res@gsnLeftString = "Indian SST forcing"
 
    plot = new(1,graphic)
    plot(0) = gsn_csm_contour_map(wks, PSL_diff_indian, res)
     

    ;;;画box
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=2.0
   pgres@gsLineDashPattern=0
   
    lat1 = 0
    lat2 = 15
    lon1 = 110
    lon2 = 155

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
    str3 = ispan(1, 3, 1) +""
    plot1 =new(1,graphic)
    ;;;
    
    plot1@$str3(0)$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres) 
   

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

    iplot(0)=gsn_csm_contour(wks,prob_slp_indian,res_pro)
    

    
    overlay(plot(0), iplot(0))
    
    
   getvalues plot@contour
   "cnFillColors" :cnFillColors
   end getvalues
   print(cnFillColors)

   draw(plot)
   frame(wks)
    
    print("done script")
  end
  

