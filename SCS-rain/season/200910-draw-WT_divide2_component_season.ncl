  ;;********************************************************
;;计算海洋平流项的时间变化的比较，分为动力项和热力项的作用,只关注
;;南北两侧的水平平流的动力项与热力项以及垂直扩散项
;;                                        
;;                                     linshuheng 2020.03.11
;;
;;
;**********************************************************
function  rc_test(ssta_rc)
begin
  tval = onedtond(ssta_rc@tval , dimsizes(ssta_rc))
   df   = onedtond(ssta_rc@nptxy, dimsizes(ssta_rc)) - 2
   b = tval    ; b must be same size as tval (and df)
   b = 0.5
   prob = betainc(df/(df+tval^2),df/2.0,b)    ; prob(nlat,nlon)
   copy_VarCoords(ssta_rc, prob)

return prob
end
function  wgt_annual(total_month)
begin
    weight1 = ispan(1,12,1)*1.0
    weight1 = (2*weight1-1)/24

    weight2 = ispan(1,12,1)*1.0
    weight2 = (25-2*weight2)/24
    total_annual_test = new(34,"float")
    total_annual_test@_FillValue = default_fillvalue("float")
    printVarSummary(total_annual_test)

    do i =0,32,1
      total_annual_test(i+1) = sum(total_month(i*12:(i+1)*12-1)*weight1) + sum(total_month((i+1)*12:(i+2)*12-1)*weight2)
    end do    
return total_annual_test
end    
    


;;**************************************************************
;;;main program
;;***************************************************************
begin

    season = "MAM"
 ;;;time span
    startmon =198001
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
;;;region
;;;region
   lat1 = -5
   lat2 = 20
   lon1 = 100
   lon2 = 170
;;;

  
    

   ;;δ(WT)
    f1 = addfile("/home/ys17-19/lsh/Project/SCS-rain/extend_to_2015/SODA3.4.2_wt_T_in_MLD_density_defined_0.8_2015.nc", "r")

    wt_mlt = f1->wt_mlt(:,{lat1:lat2},{lon1:lon2})
    T_mlt = f1->T_mlt_diff(:,{lat1:lat2},{lon1:lon2})

    wtT = wt_mlt
    wtT = wt_mlt*T_mlt

    T_mlt_season = month_to_season(T_mlt, season)
    wt_mlt_season = month_to_season(wt_mlt, season)

    wtT&lat@units = "degrees_north"
    wtT&lon@units = "degrees_east"

    wtT_season = month_to_season(wtT, season)


    dwt  = wtT_season(1,:,:)
    dwt = (/dim_avg_n_Wrap(wtT_season(18:,:,:), 0) - dim_avg_n_Wrap(wtT_season(0:18-1,:,:), 0)/)
    dwt&lat@units = "degrees_north"
    dwt&lon@units = "degrees_east"

    printVarSummary(dwt)

   ;;

   ;;δ(W)T
    tdw = wtT(1,:,:)
    T_clm  = dim_avg_n_Wrap(T_mlt_season(0:18-1,:,:), 0)
    w_delta = dim_avg_n_Wrap(wt_mlt_season(18:,:,:), 0) - dim_avg_n_Wrap(wt_mlt_season(0:18-1,:,:), 0)
    
    tdw = (/T_clm*w_delta/)
   
    ;;;

    ;;δ(T)w
    wdt = wtT(1,:,:)
    w_clm  = dim_avg_n_Wrap(wt_mlt_season(0:18-1,:,:), 0)
    T_delta = dim_avg_n_Wrap(T_mlt_season(18:,:,:), 0) - dim_avg_n_Wrap(T_mlt_season(0:18-1,:,:), 0)
    
    wdt = (/T_delta*w_clm/)


    ;;

   ;;δ(T)δ(w)
    dwdt = wtT(1,:,:)
    dwdt = (/w_delta*T_delta/)

   ;;

   ;;eddy 
   eddy = wtT(1,:,:)
   eddy = (/dwt - wdt - tdw - dwdt/)
   eddy&lat@units = "degrees_north"
   eddy&lon@units = "degrees_east"

   ;;
  
  dwt = (/dwt*10^5/)
  wdt = (/wdt*10^5/)
  tdw = (/tdw*10^5/)
  dwdt = (/dwdt*10^5/)  ;;365*86400
  eddy = (/eddy*10^5/)
;;;;draw

 avg1 = avg(tdw)/avg(dwt)
 avg2 = avg(eddy)/avg(dwt)
 avg3 = avg(wdt)/avg(dwt)
 avg4 = avg(dwdt)/avg(dwt)
 print(avg1)
 print(avg2)
 print(avg3)
 print(avg4)
    exit()
;;**************************************************************
;画图
;;***************************************************************
  wks   = gsn_open_wks("pdf", "WT的分解-philipian-"+season)

  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  res@cnFillPalette="BlRe"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 45
  ; res@mpMinLonF             =90
  ; res@mpMaxLonF             = 160

  res@mpMinLatF             = 0
  res@mpMaxLatF             = 15
  res@mpMinLonF             = 125
  res@mpMaxLonF             = 160



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

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.3
  res@cnMaxLevelValF       =   0.3
  res@cnLevelSpacingF      =  0.05
 

  plot = new(5,graphic)

  res@gsnLeftString = "Total"
  plot(0) = gsn_csm_contour_map(wks, dwt, res)
  res@gsnLeftString = "dW.Tclm"
  plot(1) = gsn_csm_contour_map(wks, tdw, res)
  res@gsnLeftString = "dT.Wclm"
  plot(2) = gsn_csm_contour_map(wks, wdt, res)
  res@gsnLeftString = "dT.dW"
  plot(3) = gsn_csm_contour_map(wks, dwdt, res)
  res@gsnLeftString = "EDDY"
  plot(4) = gsn_csm_contour_map(wks, eddy, res)
  



  pres = True
  pres@gsnPanelLabelBar  = True

  gsn_panel(wks, plot, (/3,2/), pres)
 


end







   ;;;

;;;;;





