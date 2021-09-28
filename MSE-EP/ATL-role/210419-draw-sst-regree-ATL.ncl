;;********************************************************
;;计算SST倾向以及SST值与GPCP降水的趋势
;;其中热通量的方向都是向上为正，向下为负
;;的相关系数场，
;;NCEP-NCAR的精度为1.87×1.9
;;
;;                                         lnshuheng 2020.02.18
;**********************************************************
function chazhi(data)
 begin
   xi = data&lon
   yi = data&lat
   xo = fspan(0.5,359.5,360)
   yo = fspan(-89.5,89.5,180)
 
   data_chazhi = linint2_Wrap(xi, yi, data, True, xo, yo, 0)
  
   data_chazhi!1 = "lat"
   data_chazhi!2= "lon"
   
   data_chazhi&lat = yo
   data_chazhi&lon = xo
   copy_VarAtts(data, data_chazhi)
 
   return data_chazhi
end 
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
 ;;;time span
    startmon =197901
    endmon = 201912

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
    order  = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"/)
 

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ;;;
    delete(date)
    delete(timeind)

    ssta_de = ssta
    ; ssta_de=(/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ssta_in = dim_avg_n_Wrap(ssta_de(:,{5:-5},{190:240}), (/1,2/))
    printVarSummary(ssta_in)

    nino3_4 = month_to_season(ssta_in, "DJF")   ;;;第一年是JF平均，第二年是DJF平均
    printVarSummary(nino3_4)


    ; nino3_4_new = nino3_4                ;;;;前期ENSO
    ; nino3_4_new(0:33) = (/nino3_4(1:)/)
    ; nino3_4_new(34) =ssta_index(35*12-1) 

    nino3_4 = (/dim_standardize_n(nino3_4, 1, 0)/)
;;;

;;;area mean 
      rad    = 4.0*atan(1.0)/180.0
      re     = 6371220.0
      rr     = re*rad

      dlon   = abs(lon(2)-lon(1))*rr
      dx     = dlon*cos(lat*rad)
      ;          
      jlat  = dimsizes( lat )
      dy     = new ( jlat, typeof(dx))      ; close enough
      
      dy(0)  = abs(lat(2)-lat(1))*rr
      dy(1:jlat-2)  = abs(lat(2:jlat-1)-lat(0:jlat-3))*rr*0.5   
      dy(jlat-1)    = abs(lat(jlat-1)-lat(jlat-2))*rr

      area   = dx*dy                                ; cell area function of latitude only
      copy_VarCoords(lat, area)
      printVarSummary(area)
;;;;

  

    ; wks   = gsn_open_wks("pdf", "CMAP&olr回归到印度洋SST的回归系数")
    wks   = gsn_open_wks("pdf", "SST 回归大西洋SST的回归系数")
  
    plot=new(4,graphic)
    plot1 = plot

    str = ispan(1,4,1)+ ""
    str2 = ispan(5,8,1)+ ""
 
    iplot = new(4,graphic)

    season =  (/"MAM","JJA","SON","DJF"/)
;;;;;



do i =0,3

  ;;;region
   lat1 = -50
   lat2 = 50
   lon1 = 0
   lon2 = 360
  ;;;

 
  ssta_annual = month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), season(i))
  printVarSummary(ssta_annual)

  ssta_mam =  month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), "MAM")


  ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_mam&time),1), wgt_areaave(ssta_mam(:,{0:30},{280:340}), area({0:30}), 1.0, 0) , 0, 0)
  ssta_index_rc  =ssta_index_rc*yearnum


  print("*******ssta_index_rc="+ssta_index_rc)


  ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)
 
  ssta_index = wgt_areaave_Wrap(ssta_mam(:,{0:30},{280:340}), area({0:30}), 1.0, 0) 


  ; ;;;计算各变量的实际趋势
  ;  p_rc := regCoef_n(ispan(1,yearnum,1),p_annual,0,0)
  ;  p_rc = (/p_rc*yearnum/)

  ; ;;;;去趋势
  ; p_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),p_annual,False,True,0)/)
  ; u_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),u_annual,False,True,0)/)
  ; v_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),v_annual,False,True,0)/)
  


  ;;;correlation with sst
    r_sst_p := regCoef_n(ssta_index, ssta_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_p)

    r_sst_p = (/r_sst_p*ssta_index_rc/)
  
   ;;;;
   

    print("*******显著性检验*****************")

    ;;;t test 
    prob_p1= rc_test(r_sst_p)
 

    ; prob_p1 = where(p_percent.ge.0.9, 0.01, 2)
    ; ; prob_p2 = where(u_percent.ge.0.5, 0.01, 2)
    ; ; prob_p3 = where(v_percent.ge.0.5, 0.01, 2)
    ; prob_p4 = where(olr_percent.ge.0.9, 0.01, 2)


    ; prob_p1 = where(p_percent.le.-0.9, 0.06, prob_p1)
    ; ; prob_p2 = where(u_percent.le.-0.5, 0.06, 2)
    ; ; prob_p3 = where(v_percent.le.-0.5, 0.06, 2)
    ; prob_p4 = where(olr_percent.le.-0.9, 0.06, prob_p4)



    ; prob_p1 = where(prob_p1.eq.0, 2, prob_p1)
    ;;;;


   



  print("*********开始画图********")

;;**************************************************************
;画图
;;***************************************************************


  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  ; res@cnFillPalette="MPL_BrBG"
  res@cnFillOn=True
  res@cnLinesOn=False
  res@gsnAddCyclic = False

  res@cnLineLabelsOn=False

  res@lbLabelBarOn= True

 
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbLabelStride=1
  res@lbOrientation="horizontal"
  res@mpLimitMode           = "LatLon"
  ; res@mpMinLatF             = lat1
  ; res@mpMaxLatF             = lat2
  ; res@mpMinLonF             = lon1
  ; res@mpMaxLonF             = lon2
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 30
  ; res@mpMinLonF             = 100
  ; res@mpMaxLonF             = 170
   
  res@mpMinLatF             = -50
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360
  


  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False
  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

  res@tmBorderThicknessF=1.0

  res@tmXBLabelFontHeightF=0.02
  res@tmYLLabelFontHeightF=0.02
  res@tmXBMajorThicknessF=1.5
  res@tmYLMajorThicknessF=1.5
  res@tmBorderThicknessF = 1.5

  ; res@tmXBValues = (/100,120,140,160/) 
  ; res@tmXBLabels = (/100,120,140,160/) +""+"E" 
  
  ; res@tmYLValues = (/-30,-20,-10,0,10,20,30/) 
  ; res@tmYLLabels = (/"30S","20S","10S","0","10N","20N","30N"/) 

  ; res@tmYLMinorValues = ispan(-30,30,5)




  res@gsnScalarContour        = True
  res@mpFillOn = False


  ; res@cnFillPalette = "precip_diff_12lev"
  ; res@cnFillPalette = "CBR_drywet"
  res@cnFillPalette  = "BlueWhiteOrangeRed"
  res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)
  ; res@cnLevels := (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)/2.0
  ; res@cnLevels := (/-1.0,-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8,1.0/)/2.0
  ; res@cnLevels := (/-0.8,-0.6,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.6,0.8/)/2.0

  ; res@cnLevels := (/-2.5,-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)


    ; res@cnLevels := (/-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0,2.5/)
    ; res@cnFillColors= (/3,4,6,8,9,10,11,12,13/) - 2 
      ;  res@cnLevels := (/-1.5,-1.0,-0.5,-0.25,0.25,0.5,1.0,1.5,2.0/)
  ; res@cnFillColors= (/2,3,4,6,8,9,10,11,12,13/) - 2 
  ; res@cnLevels := (/-1.5,-1.0,-0.5,-0.2,0.2,0.5,1.0,1.5,2.0/)

    ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
    res@cnLevels := (/-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4/)



  ; gsn_define_colormap(wks, "precip_diff_12lev")
  ; res@cnFillColors= (/2,3,4,5,6,8,9,10,11,12,13/) - 2
  ; res@cnFillColors= (/2,3,4,5,6,-1,9,10,11,12,13/) -2
  res@lbLabelStride=1


  res@cnExplicitLabelBarLabelsOn = True
  ; res@lbLabelStrings =  res@cnLevels
  res@vcExplicitLabelBarLabelsOn = True
  res@stExplicitLabelBarLabelsOn = True
  res@lbLabelFontHeightF = 0.016
  res@lbAutoManage = False

  
  ; res@cnFillPalette = "BlRe"
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

  ; res@cnMinLevelValF       =   -0.5
  ; res@cnMaxLevelValF       =   0.5
  ; res@cnLevelSpacingF      =  0.1



    res@gsnLeftString = "The components of "+ season(i)+ " SST trend linked to MAM TNAST"

    ; res@gsnLeftStringOrthogonalPosF  = -0.12

  
  plot(i) = gsn_csm_contour_map(wks,r_sst_p, res)

  
  ; delete(res@cnLevelSelectionMode)
  delete(res@cnFillPalette)
  delete(res@cnLevels)
  delete(res@cnFillColors)


   
  ;;;添加indo-pacific warm pool 
    ; ;;添加区域
   ;;warm pool
    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF= 2.0
    pgres@gsLineDashPattern= 1

    pgres@gsLineThicknessF=2.0
    pgres@gsLineDashPattern=1
    str3=unique_string("dum3")
    str4=unique_string("dum4")
 
   ;;;画box
   lat1 = 0
   lat2 = 20
   lon1 = 210
   lon2 = 280

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)

    plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   ;;;
   pgres@gsLineColor="red"

    lat1 = 0
   lat2 = 30
   lon1 = 290
   lon2 = 340

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;
   
    plot@$str2(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   
   


    ;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLinesOn=False
    res_pro@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
    res_pro@cnMinLevelValF      = 0.00        ; set min contour level
    res_pro@cnMaxLevelValF      = 1.00        ; set max contour level
    res_pro@cnLevelSpacingF     = 0.1       ; set contour spacing

    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False


    iplot(i) = gsn_csm_contour(wks,prob_p1,res_pro)

   
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17        ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    
         
        ; res_pro=True
        ; res_pro@gsnDraw=False
        ; res_pro@gsnFrame=False
        ; res_pro@cnFillOn = True

        ; res_pro@cnLevelSelectionMode = "ExplicitLevels"
        ; res_pro@cnLevels = (/0,0.05,1.0/)
        ; res_pro@lbLabelBarOn = False

        ; res_pro@cnLinesOn=False
        ; res_pro@cnFillDrawOrder = "PostDraw"

        ; res_pro@cnMonoFillPattern = True
        ; res_pro@cnFillPattern = 17
        ; res_pro@cnFillDotSizeF = 0.0028
        ; res_pro@cnMonoFillColor = False
        ; res_pro@cnFillColors = (/"Transparent","red","blue","Transparent"/)
        ; ; res_pro@cnLevels = ispan(-6,6,1)
        ; res_pro@cnLineLabelsOn=False
        ; res_pro@tmYRMode = "Automatic"
        ; res_pro@cnInfoLabelOn=False

        ; iplot(i*2) = gsn_csm_contour(wks,prob_p1,res_pro)

        ; iplot(i*2+1) = gsn_csm_contour(wks,prob_p4,res_pro)
    
    

    ;;;  
    
      overlay(plot(i), iplot(i))
      
   
      delete(lat1)
      delete(lat2)
      delete(lon1)
      delete(lon2)
end do

  pres = True
  pres@gsnMaximize = True
  
  pres@gsnPanelBottom = 0.05

  ; gsn_panel(wks, plot, (/2,2/), pres)
  draw(plot(0))
  frame(wks)
  
; system("convert -geometry 1800x1800 -density 600x600 -trim  GPCP&olr回归到印度洋SST的回归系数.png  GPCP&olr回归到印度洋SST的回归系数_high.png")     

 
 print("done")

end