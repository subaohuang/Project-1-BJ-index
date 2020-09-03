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

    fils1 = "/home/ys17-19/lsh/CESM-data/B/" + (/"B2000_alt_north_year_CTRL/B2000_alt_north_year_CTRL.cam.h0.0251-0280_chazhi.nc"/)

    fils2 = "/home/ys17-19/lsh/CESM-data/B/" + (/"B2000_alt_north_year/B2000_alt_north_year.cam.h0.0251-0280_chazhi.nc", \
      "B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer/B2000_alt_north_year_WNPCLM_nudge2NTACLM_2buffer.cam.h0.0251-0280_chazhi.nc"/)

    leftstring = (/"AW - CTRL", "WNPCLM - CTRL","AW - WNPCLM"/)



    wks = gsn_open_wks("pdf", "diff-pre-restore-altantic-year_WNPCLM__2buffer_panel")
    plot = new(3,graphic)

    str3 = ispan(1, 3, 1) +""
    plot1 =new(1,graphic)
    iplot = new(3,graphic)


do i=0,2
    
    if i.ne.2 then
      f1=addfile(fils1, "r")
      f2=addfile(fils2(i), "r")
    else 
      f1=addfile(fils2(1), "r")
      f2=addfile(fils2(0), "r")
    end if

    time1 = f1->time
    time2 = f2->time
    time1 = (/time1 -1/)
    time2 = (/time2 -1/)
    date1=cd_calendar(time1, 1)
    date2=cd_calendar(time2, 1)


    timeind1=ind(date1.ge.25601.and.date1.le.28012)
    timeind2=ind(date2.ge.25601.and.date2.le.28012)
    
    PRECL_CT  =  f1->PRECL(timeind1,:,:)
    PRECC_CT  =  f1->PRECC(timeind1,:,:)

    PRE_CT = PRECC_CT
    PRE_CT = (/PRECC_CT+PRECL_CT/)
    ; PRE_CT = (/PRECL_CT/)  ;;;大尺度降水
    PRE_CT = (/PRE_CT*1000*3600*24/)  ;; m/s conver to mm/day

    PRECL_id  =  f2->PRECL(timeind2,:,:)
    PRECC_id  =  f2->PRECC(timeind2,:,:)

    PRE_id = PRECC_id
    PRE_id = (/PRECC_id+PRECL_id/)
    ; PRE_id = (/PRECL_id/)
    PRE_id = (/PRE_id*1000*3600*24/)  ;; m/s conver to mm/day

    PRE_CT_clm = dim_avg_n_Wrap(PRE_CT, 0)
    PRE_indian_clm =dim_avg_n_Wrap(PRE_id, 0)

    PRE_diff = PRE_CT_clm
    PRE_diff = (/PRE_indian_clm - PRE_CT_clm/)

    ;;;t_test 
    prob_p = t_test( PRE_CT ,PRE_id)
    PRE_diff@_FillValue = 10^9
     

  
    res=True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnAddCyclic=True


    res@cnLineColors="Transparent"  ;
    res@lbLabelBarOn=False
    res@cnFillOn = True
    res@cnLinesOn =False
    res@cnLineLabelsOn=False
    ; res@cnFillMode  = 
    ; res@cnRasterSmoothingOn=True
    res@trGridType = "TriangularMesh"

    ; res@cnFillPalette = "hotcolr_19lev"


    gsn_define_colormap(wks, "MPL_BrBG")
 


    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels
    ; res@cnMinLevelValF       =   -0.5
    ; res@cnMaxLevelValF       =   0.5
    ; res@cnLevelSpacingF      =  0.1
    
    res@cnMinLevelValF       =   -2.5
    res@cnMaxLevelValF       =   2.5
    res@cnLevelSpacingF      =  0.5
    
    ; res@cnFillColors = (/2,3,4,5,6,7,8,9,10,11,0,0,16,18,20,21/)

    res@tmXTOn=False
    res@tmYROn=False

    res@mpLimitMode           = "LatLon"
     
    res@mpMinLatF             = -30.0
    res@mpMaxLatF             = 30.0
    res@mpMinLonF             = 0.0
    res@mpMaxLonF             = 360.0
    res@mpCenterLonF = 180.0

    ; res@mpMinLatF             = -10.0
    ; res@mpMaxLatF             = 50.0
    ; res@mpMinLonF             = 60.0
    ; res@mpMaxLonF             = 160.0




    res@mpFillOn = False
    res@mpGeophysicalLineColor      = "Grey18"
    res@mpGeophysicalLineThicknessF = 1.5
    


    ; res@mpShapeMode = "FreeAspect"
    ; res@vpWidthF = 0.7
    ; res@vpHeightF = 0.3

    
    
    res@gsnRightString = ""
    res@gsnLeftString = leftstring(i)
 
    plot(i) = gsn_csm_contour_map(wks, PRE_diff, res)
     

    ;;;画box
    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=2.0
    pgres@gsLineDashPattern=0

    lat1 = 0
    lat2 = 15
    lon1 = 125
    lon2 = 160

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   
    ;;;
    
    plot1@$str3(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres) 
   

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
 

    iplot(i)=gsn_csm_contour(wks,prob_p ,res_pro)
    overlay(plot(i), iplot(i))
end do

   pres = True
   pres@gsnPanelLabelBar = True
   pres@lbOrientation = "Vertical"
   pres@pmLabelBarWidthF = 0.15
   gsn_panel(wks, plot, (/3,1/), pres)
  end
  
  

