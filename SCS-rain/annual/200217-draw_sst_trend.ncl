;;********************************************************
;;画HadIsst的sst趋势
;;
;;lnshuheng 2020.02.17
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
;;**************************************************************
;;;main program
;;***************************************************************
begin
  
   ;;;读入HadiSST海温数据
    diri="/home/ys17-19/lsh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)

    timeind:=ind(date.le.201012.and.date.ge.199901)

    sst=lonFlip(f1->sst(timeind,:,:))
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
  


   ;;;先得到秋季和冬季eninio年的沃克环流情况SONDJF
    ssta_annual = month_to_annual(ssta, 1)

    delete(date)
    delete(timeind)
    ;;;
    
    dim = dimsizes(ssta_annual)
    ssta_rc=regCoef_n(ispan(1, dim(0), 1), ssta_annual, 0, 0)
    copy_VarCoords(ssta_annual(1,:,:), ssta_rc)
    prob_p=rc_test(ssta_rc)

    ; ssta_rc = where(prob_p.le.0.05, ssta_rc, ssta_rc@_FillValue)
    ssta_rc=(/ssta_rc*dim(0)/)
    printVarSummary(ssta_rc)
    
    ; ssta_rc_avg = avg(dim_avg_n_Wrap(ssta_rc({0:15},{125:160}), 0)) 
    ; print(ssta_rc_avg)

    ; ssta_rc = (/ssta_rc - ssta_rc_avg/)
    ; print(ssta_rc)
    
    
    wks   = gsn_open_wks("pdf", "年平均SST趋势-global-1999-2010")
;;**************************************************************
;画图
;;***************************************************************
  plot=new(1,graphic)
  res=True
  res@gsnFrame=False
  res@gsnDraw=False
  res@gsnLeftStringFontHeightF=20.0
  
  ; res@cnFillPalette = "NCV_jaisnd"
 res@cnFillPalette = "BlueWhiteOrangeRed"

  ; res@gsnSpreadColorStart = 130
  ; res@gsnSpreadColorEnd = 253
  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=True
  res@lbLabelStride=2
  res@lbOrientation="Vertical"
  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -60
  res@mpMaxLatF             = 60
  res@mpMinLonF             =0
  res@mpMaxLonF             = 360

  ; res@mpMinLatF             = 0
  ; res@mpMaxLatF             = 15
  ; res@mpMinLonF             =125
  ; res@mpMaxLonF             = 160



  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  ; res@tmXBValues=(/30,90,150,210,270,330,387.5/)
  ; res@tmXBLabels=(/"30~S~o~N~E","90~S~o~N~E","150~S~o~N~E","150~S~o~N~W"\
  ;             ,"90~S~o~N~W","30~S~o~N~W","30~S~o~N~E"/)

  ; res@tmYLValues=(/30,20,10,0,-10,-20,-30/)
  ; res@tmYLLabels=(/"30~S~o~N~N","20~S~o~N~N","10~S~o~N~N","0"\
  ;       ,"10~S~o~N~S","20~S~o~N~S","30~S~o~N~S"/)
  ; res@tmYLLabels=(/"30~S~o~C~~Z130~~H5~N","20~S~o","10~S~o","0"\
  ;               ,"10~S~o","20~S~o~","~Z70~~H5~~V2~S~C~~Z100~~V20~30~S~o"/)

  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=3.0
  res@tmYLMajorThicknessF=3.0
  res@tmBorderThicknessF=4.0

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  res@cnMinLevelValF       =   -0.6
  res@cnMaxLevelValF       =   0.6
  res@cnLevelSpacingF      =  0.05
  
  ; res@cnMinLevelValF       =   0.2
  ; res@cnMaxLevelValF       =   0.6
  ; res@cnLevelSpacingF      =  0.05

 
  ; res@cnMinLevelValF       =   -0.2
  ; res@cnMaxLevelValF       =   0.2
  ; res@cnLevelSpacingF      =  0.025

  ; res@gsnLeftString="HadISST sst trend"
  ; res@gsnLeftString="SST trend (remove regional average)"
  res@gsnLeftString="SST trend 1990-2016"
  plot = gsn_csm_contour_map(wks, ssta_rc, res)

  res2 = res
  res2@cnLevelSelectionMode = "ExplicitLevels"
  res2@cnFillOn = False
  res2@cnLinesOn = True
  res2@cnLevels = (/0.55/)
  res2@cnLineColors = "black"
  res2@cnMonoLineColor  = True
  res2@cnLineThicknessF = 5
  res2@cnMonoLineThickness  = True
  res2@gsnLeftString = ""
  res2@cnInfoLabelOn = False

  plot1 = gsn_csm_contour(wks, ssta_rc, res2)
  overlay(plot, plot1)  

   str1 = unique_string("dum1")  ; "unique_string" will return a unique    


  ;;;添加赤道直线
   pgx=ispan(30,390,10)
   pgy=new(dimsizes(pgx),typeof(pgx))
   pgy=0
   str1=unique_string("dum1")
   str2=unique_string("dum2")
    str3=unique_string("dum3")
   pgres=True
   pgres@gsLineColor="black"
   pgres@gsLineThicknessF=1.0
   pgres@gsLineDashPattern=1
   plot@$str1$ = gsn_add_polyline(wks, plot(0), pgx, pgy, pgres)
  ;;;添加indo-pacific warm pool
   ;  ; ;;添加区域
   ; ;;warm pool
   ; pgres@gsLineThicknessF=2.0
   ; pgres@gsLineDashPattern=0
   ; str3=unique_string("dum3")
   ; str4=unique_string("dum4")
   
   ; ;;;画box
   ; lat1 = 0
   ; lat2 = 15
   ; lon1 = 110
   ; lon2 = 155

   ; pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   ; pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ; ;;;
   ;  plot@$str3$ = gsn_add_polyline(wks, plot(0), pgx_pool, pgy_pool, pgres)





    
    
   ;;;打点
    res_pro=True
    res_pro@gsnAddCyclic = False
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.01/)
    res_pro@lbLabelBarOn = False
    res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.003
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","White","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False
    ; res_pro@cnLineLabelBackgroundColor="Transparent"
    ; res_pro@cnLineThicknessF = 3.
    ; res_pro@cnFillOn        = True
  
  
     ; iplot = gsn_csm_contour(wks,prob_p  ,res_pro)

    ;;;  
    ; overlay(plot, iplot)
draw(plot)
frame(wks)


end