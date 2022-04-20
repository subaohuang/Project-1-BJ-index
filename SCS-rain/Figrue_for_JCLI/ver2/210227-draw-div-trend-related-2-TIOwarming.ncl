;;********************************************************
;;;;计算1000hpa-850hpa平均的wind-divergence 的趋势与
;;;;TIO SST 增温有关的分量为多少，
;;;;TIO 增温能够解释多少 观测的divergence的趋势
;**********************************************************
function chazhi(data)
 begin
   xi = data&longitude
   yi = data&latitude
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

undef("remove_enso_signal")
function remove_enso_signal(nino3_4,  data)
  begin  
   ;;;;;
      ; data_de=(/dtrend_msg_n(ispan(1,dimsizes(data&time),1),data,False,True,0)/) 
      data_de  = data
      data_rc = regCoef_n(nino3_4,data_de, 0, 0)
      copy_VarCoords(data(1,:,:), data_rc)
      data_rc_3d  = conform_dims(dimsizes(data), data_rc, (/1,2/))
      m = reshape(data_rc@yintercept,dimsizes(data_rc_3d)) 
      nino3_4_3d = conform(data_rc_3d,nino3_4,(/0/))
  
      data_relate_enso  = nino3_4_3d*data_rc_3d + m 
      data_new = data 
      data_new = (/data - data_relate_enso/)
  return data_new
end 


;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =197901
    endmon = 201512

    startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
    endyear_r = str_split_by_length(endmon, 4 )

    startyear = stringtoint(startyear_r(0))
    endyear = stringtoint(endyear_r(0))

    yearnum = endyear - startyear +1
  ;;
   order  = (/"(a)","(b)","(c)","(d)"/)
   season = (/"MAM","JJA","SON","DJF"/)
   str = ispan(1, 4, 1) + ""

   plot = new(4,graphic)
   iplot = plot
 ;;;
 

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
;==============================================================
; Mission 2 div data
; ==============================================================

  ;   fils3 = systemfunc("ls /home/yangsong3/data-observation/ERA-interim/monthly/divergence/DIV_monthly_*"+".nc")
  ;   f3 = addfiles(fils3, "r")

  ;   date:=cd_calendar(f3[:]->time, 1)
  ;   timeind:=ind(date.ge.startmon.and.date.le.endmon)
  ;   time := f3[:]->time(timeind)
  ;   lat := f3[1]->latitude
  ;   lon := f3[1]->longitude

  ;   div_r = new((/dimsizes(time),7, dimsizes(lat),dimsizes(lon)/),float)
  ;   do i = 0,dimsizes(timeind)/12-1
  ;      div_r(i*12:(i+1)*12-1,:,:,:) = short2flt(f3[i]->d(:,{1000:850},::-1,:))
  ;   end do

   
  ;   div = dim_avg_n_Wrap(div_r,1)

  ;   div = (/div*10^6/)
  ;   div = (/-1*div/)
  ;   div_chazhi = chazhi(div)
  ;   printVarSummary(div_chazhi)
    
  ;   system("/bin/rm -f era_div_chazhi_1979-2015.nc")
  ;   ncout1 = addfile("era_div_chazhi_1979-2015.nc", "c")
  ;   ncout1->div_chazhi = div_chazhi

  ; exit()


    f3=addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/ver2/era_div_chazhi_1979-2015.nc", "r")   
    date:=cd_calendar(f3->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    
    div_chazhi = f3->div_chazhi(timeind,:,:)

    printVarSummary(div_chazhi)
 
    div_chazhi=(/rmMonAnnCycTLL(div_chazhi)/)

    wks   = gsn_open_wks("pdf", "div-related-TIOwarming.pdf")

do i =0,3

  ; ;;;region
  ;  lat1 = -30
  ;  lat2 = 30
  ;  lon1 = 30
  ;  lon2 = 180
  ; ;;;
   ;;;region
   lat1 = -30
   lat2 = 30
   lon1 = 0
   lon2 = 360
  ;;;

 
  ssta_annual = month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), season(i))
  printVarSummary(ssta_annual)




  ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{-20:20},{40:110}), area({-20:20}), 1.0, 0) , 0, 0)
  ; ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{0:30},{290:360}), area({0:30}), 1.0, 0) , 0, 0)
  ssta_index_rc  =ssta_index_rc*yearnum

  if i.eq.3 then 
    ssta_index_rc = 0.42
  end if

  if i.eq.0 then 
      ssta_index_rc = 0.23
    end if
     if i.eq.1 then 
      ssta_index_rc = 0.41
    end if
  print("*******ssta_index_rc="+ssta_index_rc)

  ;;;remove_enso_signal 
  ssta_annual = (/remove_enso_signal(nino3_4, ssta_annual)/)
  ;;;

  ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)
 
 
  ;;;;;要分析的变量
  div_annual = month_to_season(div_chazhi(:,{lat1:lat2},{lon1:lon2}), season(i))


  ;;;remove_enso_signal 
  div_annual = (/remove_enso_signal(nino3_4, div_annual)/)
  ;;;;;;;




  ssta_index = wgt_areaave_Wrap(ssta_annual(:,{-20:20},{40:100}), area({-20:20}), 1.0, 0) ;;;indian 
  ; ssta_index = wgt_areaave_Wrap(ssta_annual(:,{0:30},{290:360}), area({0:30}), 1.0, 0) ;;;atlantic


  ;;;计算各变量的实际趋势

   p_rc = regCoef_n(ispan(1,yearnum,1),div_annual,0,0)
   p_rc = (/p_rc*yearnum/)
   
   copy_VarCoords(div_annual(1,:,:), p_rc)


  ;;;;去趋势
    div_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),div_annual,False,True,0)/)


  ;;;correlation with sst
    r_sst_p := regCoef_n(ssta_index, div_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_p)

    r_sst_p = (/r_sst_p*ssta_index_rc/)
   ;;;;;

  

   
    ;;;percent
     p_percent  = r_sst_p
     p_percent = r_sst_p/p_rc 
    ;;;


    ;;;;计算正趋势的平均
    p_rc_neg = where(p_rc({0:20},{115:155}).le.0, p_rc({0:20},{115:155}), 0)
    p_rc_avg = avg(p_rc_neg)

    r_sst_p_neg = where(r_sst_p({0:20},{115:155}).le.0, r_sst_p({0:20},{115:155}), 0)
    
    r_sst_p_avg = avg(r_sst_p_neg)

    percent_avg = r_sst_p_avg / p_rc_avg

  print("******解释多少负的趋势*****************"+""+season(i)+""+percent_avg)
 

 
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
  if i.eq.3 then
  res@lbLabelBarOn= True
  else
  res@lbLabelBarOn= False
  end if
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbLabelStride=1
  res@lbOrientation="horizontal"
  res@mpLimitMode           = "LatLon"

  res@mpMinLatF             = -5
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 100
  res@mpMaxLonF             = 170
  
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


    res@lbLabelStride=1

    res@cnExplicitLabelBarLabelsOn = True
    ; res@lbLabelStrings =  res@cnLevels
    res@vcExplicitLabelBarLabelsOn = True
    res@stExplicitLabelBarLabelsOn = True
    res@lbLabelFontHeightF = 0.016
    res@lbAutoManage = False



    res@cnLevelSelectionMode = "ExplicitLevels"      ; set manual contour level
    res@cnFillPalette = "MPL_PuOr"
    res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)

    res@gsnLeftString = order(i)+" "+ season(i)+ " linear congruency of div trend with TIOST"

    res@gsnLeftStringOrthogonalPosF  = -0.12

  
    plot(i) = gsn_csm_contour_map(wks, r_sst_p, res)

  
  ; plot(i) = gsn_csm_vector_scalar_map(wks,r_sst_u,r_sst_v , r_sst_w, res)
  

   ; ;;;画气候态的等值线
   ;  res@gsnRightString=""
   ;  res@gsnLeftString=""
   ;  res@cnLevelSelectionMode = "ExplicitLevels"

   ;  ; res@cnLevelSelectionMode = "AutomaticLevels"
   ;  res@cnLinesOn=True
   ;  res@cnLevels := fspan(-0.6,0.2,13)
   ;  res@cnLineLabelsOn=True
   ;  res@cnLineLabelBackgroundColor="White"
   ;  res@cnLineLabelInterval = 1
   ;  res@cnLineThicknessF = 2.
   ;  res@cnFillOn = False
   ;  res@gsnContourNegLineDashPattern=1
   ;  res@gsnContourZeroLineThicknessF=1.0
   ;  res@gsnLeftString  = ""
   ;  res@gsnRightString =""
   ;  res@cnInfoLabelOn = False    

   ;  plot1(i) = gsn_csm_contour(wks, r_sst_w, res)
   ;  ;;;  
   ;  overlay(plot(i), plot1(i))
    
 

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
   lon1 = 125
   lon2 = 160

   pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

   
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
 


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

  gsn_panel(wks, plot, (/2,2/), pres)
  
; system("convert -geometry 1800x1800 -density 600x600 -trim  GPCP&olr回归到印度洋SST的回归系数.png  GPCP&olr回归到印度洋SST的回归系数_high.png")     

 
 print("done")

end