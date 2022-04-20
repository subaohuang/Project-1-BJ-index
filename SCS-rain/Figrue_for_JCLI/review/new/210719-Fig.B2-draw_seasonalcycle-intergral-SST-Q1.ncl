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

    sst_clm_avg = dim_avg_n_Wrap(sst_clm(:,{0:15},{125:160}), (/1,2/))
    ; sst_clm_avg = dim_avg_n_Wrap(sst_clm(:,{-10:10},{40:110}), (/1,2/))
    
    mean = avg(sst_clm_avg)

    ; sst_clm_avg = sst_clm_avg - mean
  ;;;

  ;;;cal change 
    ssta = sst
    ssta = (/rmMonAnnCycTLL(ssta)/)
    ssta_wnp = dim_avg_n_Wrap(ssta(:,{0:15},{125:160}), (/1,2/))
    ; ssta_wnp = dim_avg_n_Wrap(ssta(:,{-10:10},{40:110}), (/1,2/))
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

   T = cdft_t(0.975, yearnum-1)


    ; print(rc_12)
    ; print(rc_q1_12)

    std = dim_stddev(rc_12)
    print(T*rc_12@rstd)




; exit()

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

  res@gsnLeftString = ""
  res@gsnRightString = ""

  res@vpWidthF  = 0.6
  res@vpHeightF = 0.4
  ; res@vpXF      = 0.06
  ; res@vpYF      = 0.75

  res@tmXTLabelsOn   = False
  res@tmXTOn         = False
  res@tmYROn = False


  res@trXMinF        = 0
  res@trXMaxF        = 11

  res@trYMinF        = 28.0
  res@trYMaxF        = 30.0


  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02

  res@gsnLeftStringFontHeightF   = 0.02
  res@gsnRightStringFontHeightF  = 0.02

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, 11, 1)
  res@tmXBLabels = ispan(1,12,1)


  res@tmYLMode = "Explicit"
  res@tmYLValues = fspan(28.0, 30.0, 11)
  res@tmYLLabels = fspan(28.0, 30.0, 11)

  res@tiYAxisString = "SST [~S~o~N~C]" 
  res@tiYAxisFontHeightF=0.02

  res@tiXAxisString = "Month"
  res@tiXAxisFontHeightF = 0.02


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

  res@gsnLeftString = "(a) Seasonal cycle of SST "



  res@xyLineColors  = (/"black","black"/)
  res@xyMarkerColors= (/"black","black"/)

  plot(0)  = gsn_csm_xy (wks, ispan(0,11,1), (/sst_clm_avg,sst_change/), res)

    res@tmYLValues := fspan(0.4, 2.4, 11)
    res@tmYLLabels := fspan(0.4, 2.4, 11)
    res@trYMinF        = 0.4
    res@trYMaxF        = 2.4

    res@tiYAxisString = "Q1 [K/day]" 
    res@tiXAxisString = "Month"
    
   res@gsnLeftString = "(b) Seasonal cycle of vertically averaged Q1"
   plot(1)  = gsn_csm_xy (wks, ispan(0,11,1), (/q_clm_avg,q_change/), res)

 

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


    ; ;***********************************************
    ; ; legend resources
    ; ;***********************************************  
    res_lines = True
    res_lines@gsLineDashPattern = 0.                    ; solid line
    res_lines@gsMarkerIndex     = 1
    res_lines@gsMarkerSizeF     = 0.03
    res_lines@gsLineThicknessF = 2.0
    res_lines@gsLineColor       = "black"
    res_lines@gsMarkerColor = "black"
    res_lines@xyMarkerThicknessF  = 4.0

    
    x := ispan(40,50,10)*0.1
    y = x
    y = 28.4

    poly1=  gsn_add_polyline(wks, plot(0), x-0.3, y, res_lines)
    poly2=  gsn_add_polymarker(wks, plot(0), x-0.3, y, res_lines)

    y = 28.2
    res_lines@gsLineDashPattern = 1.
    poly3=  gsn_add_polyline(wks, plot(0), x-0.3, y, res_lines)
    poly4=  gsn_add_polymarker(wks, plot(0), x-0.3, y, res_lines)
    
    y =0.8
    res_lines@gsLineDashPattern = 0.
    poly5=  gsn_add_polyline(wks, plot(1), x, y, res_lines)
    poly6=  gsn_add_polymarker(wks, plot(1), x, y, res_lines)
   
    y =0.6
    res_lines@gsLineDashPattern = 1.
    poly7=  gsn_add_polyline(wks, plot(1), x, y, res_lines)
    poly8=  gsn_add_polymarker(wks, plot(1), x, y, res_lines)


    textres = True
    textres@txFontHeightF = 0.012
    textres@txJust = "CenterLeft"


 ; (/,"SST Climatology plus 37-yr linear change"/)

     text1 = gsn_add_text(wks, plot(0), "Climatological SST", x(dimsizes(x)-1)+0.2-0.3, 28.4 , textres)

     text2 = gsn_add_text(wks, plot(0), "Climatological SST plus 37-yr linear change", x(dimsizes(x)-1)+0.2-0.3, 28.2 , textres)


     text3 = gsn_add_text(wks, plot(1), "Climatological Q1", x(dimsizes(x)-1)+0.2, 0.8, textres)

     text4 = gsn_add_text(wks, plot(1), "Climatological Q1 plus 37-yr linear change", x(dimsizes(x)-1)+0.2, 0.6 , textres)

pres  = True
gsn_panel(wks, plot, (/1,2/), pres)

print("****done******")
;;;;

end 
