  
;;**************************************************************
;;;main program
;;***************************************************************
begin
;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201512.and.date.ge.197901)

    yearnum = 2015-1979+1 

    sst  = lonFlip(f1->sst(timeind,:,:))

    ; sst = (/rmMonAnnCycTLL(sst)/)
    
    sst_clm = clmMonTLL(sst)

    ; sst_clm_avg = dim_avg_n_Wrap(sst_clm(:,{0:15},{125:160}), (/1,2/))
    sst_clm_avg = dim_avg_n_Wrap(sst_clm(:,{0:20},{300:}), (/1,2/))
    
    mean = avg(sst_clm_avg)

    ; sst_clm_avg = sst_clm_avg - mean
  ;;;

  ;;;cal change 
    ssta = sst
    ssta = (/rmMonAnnCycTLL(ssta)/)
    ; ssta_wnp = dim_avg_n_Wrap(ssta(:,{0:15},{125:160}), (/1,2/))
    ssta_wnp = dim_avg_n_Wrap(ssta(:,{0:20},{300:}), (/1,2/))
    ssta_wnp_12 = new((/yearnum,12/),float)

    do i=0,11   
    ssta_wnp_12(:,i) = ssta_wnp(i::12)
    end do 
    
    rc_12 = regCoef_n(ispan(1,yearnum,1), ssta_wnp_12, 0, 0)
    rc_12 = rc_12*yearnum


   ;;;

;;;read Q1 monthly data from linshuheng
    fils = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/q1/q1"+"*")

    f3 = addfiles(fils, "r")
    date := cd_calendar(f3[:]->time, 1)

    timeind = ind(date.le.201512.and.date.ge.197901)

    q1 = f3[:]->q1(timeind,:,:,:)

    q1 = (/q1*86400.0/)
    printVarSummary(q1)

    plev = q1&lev
    dim = dimsizes(q1)

    psfc =1013.0
    ptop = min(plev)
    dp   = dpres_plevel_Wrap(plev, psfc, ptop, 0)  ; dp(30)
    printVarSummary(dp)
    
    dp_4d = conform_dims(dimsizes(q1), dp, 1)
    
    q1dp = q1(:,1,:,:)
    q1dp =  dim_sum_n_Wrap(q1*dp_4d, 1)/dim_sum_n(dp_4d, 1)  ;;;积分

    ; q1dp = q1dp/9.8

    q_clm = clmMonTLL(q1dp)
    q_clm!1 = "lat"
    q_clm!2 = "lon"

    q_clm&lat = q1&lat
    q_clm&lon = q1&lon
  
    q_clm_avg = dim_avg_n_Wrap(q_clm(:,{0:15},{125:160}), (/1,2/))

    q_mean = avg(q_clm_avg)
    ; q_clm_avg = q_clm_avg - q_mean

    ;;;cal change 
    q1a = q1dp
    q1a = (/rmMonAnnCycTLL(q1dp)/)
    q1a_mean = dim_avg_n_Wrap(q1a(:,{0:15},{125:160}), (/1,2/))
    

    q1a_12 = new((/yearnum,12/),float)

    do i=0,11   
    q1a_12(:,i) = q1a_mean(i::12)
    end do 
    
    rc_q1_12 = regCoef_n(ispan(1,yearnum,1), q1a_12, 0, 0)
    rc_q1_12 = rc_q1_12*yearnum
   ;;;
   
   sst_change = sst_clm_avg + rc_12
   q_change = q_clm_avg + rc_q1_12


;;;
plot = new(2,graphic)


  wks   = gsn_open_wks("pdf", "seasonal cycle for SST over TWNP")

;;**************************************************************
;画图
;;***************************************************************
  res = True
  res@tmBorderThicknessF = 2.0
  res@gsnDraw        = False
  res@gsnFrame       = False


  ; res@vpWidthF  = 0.9
  ; res@vpHeightF = 0.4
  ; res@vpXF      = 0.06
  ; res@vpYF      = 0.75

  res@tmXTLabelsOn   = False
  res@tmXTOn         = False
  
  res@trXMinF        = 0
  res@trXMaxF        = 11

  ; res@trYMinF        = 28.0
  res@trYMinF        = 27.0
  res@trYMaxF        = 30.0

  res@tmXBLabelFontHeightF = 0.015
  res@tmYLLabelFontHeightF = 0.015

  res@gsnLeftStringFontHeightF   = 0.015
  
  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, 11, 1)
  res@tmXBLabels = ispan(1,12,1)

  res@tmYLMinorOn = False


  res@tmYLMode = "Explicit"
  res@tmYLValues = fspan(28.0, 30.0, 11)
  res@tmYLLabels = fspan(28.0, 30.0, 11)

  res@tiYAxisString = "SST [~S~o~N~C]" 
  res@tiYAxisFontHeightF=0.015

  res@tiXAxisString = "Month"
  res@tiXAxisFontHeightF = 0.015


  res@xyMarkLineModes = (/"MarkLines","MarkLines"/)
  res@xyMarkers = (/1,1/)
  res@xyMarkerSizes = (/0.03,0.03/)
  res@xyMarkerColors= (/"black","black"/)
  res@xyMarkerThicknesses = (/4.0,4.0/)
  res@xyLineThicknesses  = (/2.0,2.0/)
  res@xyLineColors       = (/"black","black"/)
  res@xyDashPatterns     = (/0,1/)
  res@lgPerimOn = False
  res@lgLabelFontHeightF = 0.018

  res@gsnLeftString = "(a) Seasonal Cycle of area-averaged SST "


  res@xyLineColors  = (/"black","black"/)
  res@xyMarkerColors= (/"black","black"/)

  resR = res
  resR@gsnLeftString = ""
  ; resR@tmYLValues := fspan(0.3, 0.6, 7)
  ; resR@tmYLLabels := fspan(0.3, 0.6, 7)
  ;  resR@trYMinF        = 0.3
  ; resR@trYMaxF        = 0.6
  resR@tmYLValues := fspan(-0.3, 0.6, 10)
  resR@tmYLLabels := fspan(-0.3, 0.6, 10)
   resR@trYMinF        = -0.3
  resR@trYMaxF        = 0.6
  resR@xyMarkerColors= (/"red","red"/)
  resR@xyLineColors       = (/"red","red"/)
  resR@tmYRMajorLineColor = "red"
  resR@tmYRMinorLineColor = "red"
  resR@tmYRLabelFontColor = "red"

  resR@tiYAxisString = "SST change[~S~o~N~C]" 
  resR@tiYAxisFontColor = "red"


   plot(0)  = gsn_csm_xy2(wks, ispan(0,11,1), sst_clm_avg, rc_12, res,resR)

    res@tmYLValues := fspan(0.4, 2.4, 11)
    res@tmYLLabels := fspan(0.4, 2.4, 11)
    res@trYMinF        = 0.4
    res@trYMaxF        = 2.4

    res@tiYAxisString = "Q1 [K/day]" 
    res@tiXAxisString = "Month"
    
    res@gsnLeftString = "(b) Seasonal Cycle of area-averaged and vertical mean Q1"

    resR@tmYLValues := fspan(-0.5, 0.5, 11)
    resR@tmYLLabels := fspan(0.5, 0.5, 11)
    resR@trYMinF        = -0.5
    resR@trYMaxF        = 0.5


  resR@tiYAxisString = "Q1 change [K/day]" 

   plot(1)  = gsn_csm_xy2(wks, ispan(0,11,1),q_clm_avg,rc_q1_12, res,resR)

  

 

  ;;---------------- 参考线 --------------
  x=ispan(0, 11, 1)
  y1= x*0.1
  y2= x*0.1
  y1 = mean
  y2 = q_mean

  lnres   =    True
  lnres@gsLineColor="black"
  lnres@gsLineThicknessF=1.0
  lnres@gsLineDashPattern=0
 
  dum000 = gsn_add_polyline(wks,plot(0),x,y1,lnres)
  dum001 = gsn_add_polyline(wks,plot(1),x,y2,lnres)
  ; y= 0.5
  ; lnres@gsLineColor="red"
  ; lnres@gsLineDashPattern=11
  ; dum111 = gsn_add_polyline(wks,plot,x,y,lnres)
  ; y=-0.5
  ; lnres@gsLineColor="blue"
  ; dum222 = gsn_add_polyline(wks,plot,x,y,lnres)

   

  ;;---------------- 标题 --------------
  txres               = True
  txres@txFontHeightF = 0.019
  txres@txFont = "helvetica-bold"
  txres@txFontColor   = "black"
  ; gsn_text_ndc(wks, "Spring prediction barrier index",  0.5, 0.775, txres)
  
pres  = True
gsn_panel(wks, plot, (/1,2/), pres)

print("****done******")
;;;;

end 
