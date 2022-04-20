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

  

;;;read GPCP precipitation data
    
    ; f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")   
    ; date:=cd_calendar(f2->time, 1)
    ; timeind:=ind(date.ge.startmon.and.date.le.endmon)
    ; p=f2->precip(timeind,:,:)
    ; printVarSummary(p)

    ; p_chazhi = chazhi(p)
  
    ; system("/bin/rm -f gpcp_pre_chazhi_1979-2015.nc")
    ; ncout1 = addfile("gpcp_pre_chazhi_1979-2015.nc", "c")
    ; ncout1->p_chazhi = p_chazhi




    f2=addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/ver1/gpcp_pre_chazhi_1979-2015.nc", "r")   
    ; f2=addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/ver1/cmap_pre_chazhi_1979-2015.nc", "r")   

    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    
    p_chazhi = f2->p_chazhi(timeind,:,:)

    printVarSummary(p_chazhi)
 
    p_chazhi=(/rmMonAnnCycTLL(p_chazhi)/)
;;;; 


;;;read era data
    ; path_wind = "/home/yangsong3/data-observation/ERA-interim/monthly/"
    ; var_name = (/"uwnd","vwnd","wwnd"/)
    ; filein = path_wind + var_name +".mon.mean.nc"

    ; ;;;region 
    ; lat1 = -90
    ; lat2 = 90
    ; lon1 = 0
    ; lon2 = 360
    ; ;;;


    ; f1 = addfile(filein(0),"r")

    ; date :=cd_calendar(f1->time, 1)

    ; timeind :=ind(date.ge.197901.and.date.le.201512)

    ; u = short2flt(f1->u(timeind,{850},{lat1:lat2},{lon1:lon2}))

    ; f2 = addfile(filein(1),"r")
    ; v = short2flt(f2->v(timeind,{850},{lat1:lat2},{lon1:lon2}))

    ; u!1 = "lat"
    ; u!2 = "lon"

    ; v!1 = "lat"
    ; v!2 = "lon"
    
    ; printVarSummary(u)

    ; f3 = addfile(filein(2),"r")

    ; date := cd_calendar(f3->time, 1)

    ; timeind:=ind(date.ge.197901.and.date.le.201512)

    ; w = short2flt(f3->w(timeind,{500},{lat1:lat2},{lon1:lon2}))
    ; w = (/w*-100.0/)

    ; w!1 = "lat"
    ; w!2 = "lon"

    ; u_chazhi = chazhi(u)
    ; printVarSummary(u_chazhi)
    ; u_chazhi=(/rmMonAnnCycTLL(u_chazhi)/)

    ; ; system("/bin/rm -f era_u_chazhi_1979-2015.nc")
    ; ; ncout1 = addfile("era_u_chazhi_1979-2015.nc", "c")
    ; ; ncout1->u_chazhi = u_chazhi

    ; v_chazhi = chazhi(v)
    ; printVarSummary(v_chazhi)
    ; v_chazhi=(/rmMonAnnCycTLL(v_chazhi)/)

    ; ; system("/bin/rm -f era_v_chazhi_1979-2015.nc")
    ; ; ncout1 = addfile("era_v_chazhi_1979-2015.nc", "c")
    ; ; ncout1->v_chazhi = v_chazhi

    ; w_chazhi = chazhi(w)
    ; printVarSummary(w_chazhi)
    ; w_chazhi=(/rmMonAnnCycTLL(w_chazhi)/)

;      system("/bin/rm -f era_w_chazhi_1979-2015.nc")
;     ncout1 = addfile("era_w_chazhi_1979-2015.nc", "c")
;     ncout1->w_chazhi = w_chazhi
; ;;;;

    ; print("=====finish read data==========")

    f1=addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/ver1/era_u_chazhi_1979-2015.nc", "r")   
    date:=cd_calendar(f1->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    
    u_chazhi = f1->u_chazhi(timeind,:,:)

    printVarSummary(u_chazhi)
 
    u_chazhi=(/rmMonAnnCycTLL(u_chazhi)/)

    f2=addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/ver1/era_v_chazhi_1979-2015.nc", "r")   
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    
    v_chazhi = f2->v_chazhi(timeind,:,:)

    printVarSummary(v_chazhi)
 
    v_chazhi=(/rmMonAnnCycTLL(v_chazhi)/)


    f3=addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/ver1/era_w_chazhi_1979-2015.nc", "r")   
    date:=cd_calendar(f3->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    
    w_chazhi = f3->w_chazhi(timeind,:,:)

    printVarSummary(w_chazhi)
 
    w_chazhi=(/rmMonAnnCycTLL(w_chazhi)/)
;;;;;


    f6 = addfile("/home/ys17-19/lsh/Project/SCS-rain/Figrue_for_JCLI/ver1/olr_chazhi_1979-2015.nc", "r")
    ;;;region
    lat1 = -90
    lat2 = 90
    lon1 = 0
    lon2 = 360

    date:=cd_calendar(f6->time, 1)
    timeind :=ind(date.ge.startmon.and.date.le.endmon)
    olr_chazhi = f6->olr_chazhi

;;;;;


;;;;
    ; wks   = gsn_open_wks("pdf", "CMAP&olr回归到印度洋SST的回归系数")
    wks   = gsn_open_wks("pdf", "GPCP&olr回归到印度洋SST的回归系数")
  
    plot=new(8,graphic)
    plot1 = plot

    str = ispan(1,8,1)+ ""
 
    iplot = new(8,graphic)

    season =  (/"MAM","JJA","SON","DJF"/)
;;;;;



do i =0,3

  ;;;region
   lat1 = -30
   lat2 = 30
   lon1 = 30
   lon2 = 180
  ;;;

 
  ssta_annual = month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), season(i))
  printVarSummary(ssta_annual)


  ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{-20:20},{40:110}), area({-20:20}), 1.0, 0) , 0, 0)
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
  p_annual = month_to_season(p_chazhi(:,{lat1:lat2},{lon1:lon2}), season(i))


  ;;;remove_enso_signal 
  p_annual = (/remove_enso_signal(nino3_4, p_annual)/)
  ;;;;;;;


 
  ;;;;;要分析的变量
  u_annual = month_to_season(u_chazhi(:,{lat1:lat2},{lon1:lon2}), season(i))
  ;;;remove_enso_signal 
  u_annual = (/remove_enso_signal(nino3_4, u_annual)/)
  ;;;
  


  ;;;;;要分析的变量
  v_annual = month_to_season(v_chazhi(:,{lat1:lat2},{lon1:lon2}), season(i))

  ;;;remove_enso_signal 
  v_annual = (/remove_enso_signal(nino3_4, v_annual)/)
  ;;;
  
  
  ;;;;;要分析的变量
  w_annual = month_to_season(w_chazhi(:,{lat1:lat2},{lon1:lon2}), season(i))



  ;;;;;要分析的变量
  olr_annual = month_to_season(olr_chazhi(:,{lat1:lat2},{lon1:lon2}), season(i))

    


  ssta_index = wgt_areaave_Wrap(ssta_annual(:,{-20:20},{40:100}), area({-20:20}), 1.0, 0) 


  ;;;计算各变量的实际趋势
   p_rc = regCoef_n(ispan(1,yearnum,1),p_annual,0,0)
   u_rc = regCoef_n(ispan(1,yearnum,1),u_annual,0,0)
   v_rc = regCoef_n(ispan(1,yearnum,1),v_annual,0,0)
   olr_rc = regCoef_n(ispan(1,yearnum,1),olr_annual,0,0)

   p_rc = (/p_rc*yearnum/)
   u_rc = (/u_rc*yearnum/)
   v_rc = (/v_rc*yearnum/)
   olr_rc = (/olr_rc*yearnum/)

   



  ;;;;去趋势
  p_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),p_annual,False,True,0)/)
  u_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),u_annual,False,True,0)/)
  v_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),v_annual,False,True,0)/)
  olr_annual =(/dtrend_msg_n(ispan(1,yearnum,1),olr_annual,False,True,0)/)


  ;;;correlation with sst
    r_sst_p := regCoef_n(ssta_index, p_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_p)

    r_sst_p = (/r_sst_p*ssta_index_rc/)


   ;;;;;

    ;;;correlation with sst
    r_sst_olr := regCoef_n(ssta_index, olr_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_olr)

    r_sst_olr = (/r_sst_olr*ssta_index_rc/)

    region_mean = avg(r_sst_olr({0:15},{125:160}))
    ;;;;;

     print("*********"+season(i)+"**********r_sst_olr="+region_mean)

    ; ;;;normalized before regression 
    ; ssta_index = (/dim_standardize_n(ssta_index, 1, 0)/)
    ; ;;;;

    ;;;correlation with sst
    r_sst_u := regCoef_n(ssta_index, u_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_u)

    r_sst_u = (/r_sst_u*ssta_index_rc/)
    ;;;;;


    ;;;correlation with sst
    r_sst_v := regCoef_n(ssta_index, v_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_v)

    r_sst_v = (/r_sst_v*ssta_index_rc/)
    ;;;;;


    ;;;correlation with sst
    r_sst_w := regCoef_n(ssta_index, w_annual, 0, 0)
    copy_VarCoords(ssta_annual(0,:,:), r_sst_w)

    r_sst_w = (/r_sst_w*ssta_index_rc/)
    ;;;;;
   
   
    ;;;percent
    
     u_percent  = r_sst_u
     u_percent = r_sst_u/u_rc  

     v_percent  = r_sst_v
     v_percent = r_sst_v/v_rc  

     p_percent  = r_sst_p
     p_percent = r_sst_p/p_rc  

   

     olr_percent  = r_sst_olr
     olr_percent = r_sst_olr/olr_rc  
     

    ; p_percent = (/abs(p_percent)/)
    u_percent = (/abs(u_percent)/)
    v_percent = (/abs(v_percent)/)
    ;       olr_percent = (/abs(olr_percent)/)



    print("*******显著性检验*****************")

    ;;;t test 
    prob_p1= rc_test(r_sst_p)
    prob_p2= rc_test(r_sst_u)
    prob_p3= rc_test(r_sst_v)
    prob_p4= rc_test(r_sst_olr)


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



    if i.ne.0 then 
      prob_p2({0:30},{100:140}) =  where(prob_p2({0:30},{100:140}).le.0.6, 0.1,prob_p2({0:30},{100:140}))
    end if
    
    if i.eq.1 then 
      prob_p1({0:15},{125:160}) =  where(r_sst_p({0:15},{125:160}).le.-0.5, 0.005, prob_p1({0:15},{125:160}))
    end if

    r_sst_u = where(prob_p2.le.0.1.or.prob_p3.le.0.1, r_sst_u, r_sst_u@_FillValue)





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
  ; res@mpMinLatF             = lat1
  ; res@mpMaxLatF             = lat2
  ; res@mpMinLonF             = lon1
  ; res@mpMaxLonF             = lon2
  ; res@mpMinLatF             = -5
  ; res@mpMaxLatF             = 30
  ; res@mpMinLonF             = 100
  ; res@mpMaxLonF             = 170
   
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 30
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


  res@cnFillPalette = "precip_diff_12lev"
  ; res@cnFillPalette = "CBR_drywet"
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

    res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
    res@cnLevels := (/-1.5,-1.0,-0.5,-0.2,0.2,0.5,1.0,1.5/)



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

  res@vcGlyphStyle="CurlyVector"
  res@vcMonoLineArrowColor = True
  res@vcLineArrowColor = "black"
  res@vcMinDistanceF           = 0.02          ; thin vectors

  res@vcRefLengthF             = 0.05    ; define length of v
  res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcVectorDrawOrder="PostDraw"
  
  ; resv@vcMonoLineArrowColor = True
  ; resv@vcLineArrowColor = "blue"
  res@vcLineArrowThicknessF=1.2
  ; res@vcLineArrowHeadMaxSizeF=0.01
  ; res@vcLineArrowHeadMinSizeF=0.01
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006

    res@vcRefMagnitudeF          = 0.2*10
    ; res@vcMinMagnitudeF          =.3
    res@gsnLeftStringFontHeightF = 0.0185
    ; res@gsnLeftStringFontHeightF = 0.02

    res@vcRefAnnoString1On = True
    res@vcRefAnnoString1 = "2 m/s"
    res@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
    res@vcRefAnnoParallelPosF = 1.02
    res@vcRefAnnoFontHeightF = 0.013
   
    res@vcRefAnnoBackgroundColor = "Transparent"
    res@vcRefAnnoPerimOn = False

    res@gsnLeftString = order(i*2)+" "+ season(i)+ " linear congruency of Precip & 850UV trend with TIOST"

    res@gsnLeftStringOrthogonalPosF  = -0.12

  
  plot(i*2) = gsn_csm_vector_scalar_map(wks,r_sst_u,r_sst_v , r_sst_p, res)

  
  ; delete(res@cnLevelSelectionMode)
  delete(res@cnFillPalette)
  delete(res@cnLevels)
  delete(res@cnFillColors)
  ; res@cnFillPalette = "BlueWhiteOrangeRed"
  res@cnFillPalette = "MPL_bwr"

  res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    res@lbLabelStride=1

  ; res@cnLevels := (/-2.0,-1.5,-1.0,-0.5,0.5,1.0,1.5,2.0/)*10

  res@cnLevels := (/-16,-12,-8,-4,-2,2,4,8,12,16/)

  ; res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  ; res@cnMinLevelValF      = -2.5       ; set min contour level
  ; res@cnMaxLevelValF      = 2.5      ; set max contour level
  ; res@cnLevelSpacingF     = 0.5       ; set contour spacing

  res@gsnLeftStringOrthogonalPosF  = 0.06


  res@gsnLeftString = order(i*2+1)+" "+ season(i)+ " linear congruency of OLR trend with TIOST"

  plot(i*2+1) = gsn_csm_contour_map(wks, r_sst_olr, res)

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

   
   plot@$str(i*2)$ = gsn_add_polyline(wks, plot(i*2), pgx_pool, pgy_pool, pgres)
   plot@$str(i*2+1)$ = gsn_add_polyline(wks, plot(i*2+1), pgx_pool, pgy_pool, pgres)
   


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


    iplot(i*2) = gsn_csm_contour(wks,prob_p1,res_pro)

    iplot(i*2+1) = gsn_csm_contour(wks,prob_p4,res_pro)
    

    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17        ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

    iplot(i*2) = gsn_contour_shade(iplot(i*2) , 0,0.1,opt)
    iplot(i*2+1) = gsn_contour_shade(iplot(i*2+1) , 0,0.1,opt)
         
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
    
      overlay(plot(i*2), iplot(i*2))
      overlay(plot(i*2+1), iplot(i*2+1))
   
      delete(lat1)
      delete(lat2)
      delete(lon1)
      delete(lon2)
end do

  pres = True
  pres@gsnMaximize = True
  
  pres@gsnPanelBottom = 0.05

  gsn_panel(wks, plot, (/4,2/), pres)
  
; system("convert -geometry 1800x1800 -density 600x600 -trim  GPCP&olr回归到印度洋SST的回归系数.png  GPCP&olr回归到印度洋SST的回归系数_high.png")     

 
 print("done")

end