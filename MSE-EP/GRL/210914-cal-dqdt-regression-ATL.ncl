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


function partial_reg(data1,data2,data3) ;;;data1:index data2: varible data3:effect want to remove
begin 
    ;;data1 :x1 data2 :y ; data3: x2
      r12 = escorc_n(data1, data2, 0, 0) ;;;二维 ry1
      r13 = escorc_n(data1, data3, 0, 0) ;;;二维 r12
      r23 = escorc_n(data3, data2, 0, 0) ;;;二维 ry2
      std2 = dim_stddev_n_Wrap(data2, 0) ;;;二维 sy
      std1 = dim_stddev_n_Wrap(data1, 0) ;;;一个数字 s1
      std1_2d = std2
      std1_2d = (/std1/)

      reg12 = (/r12 - r23*r13/)/(1-r13^2)*(std2/std1_2d)
      corr = (/r12 - r23*r13/)/sqrt((1-r23^2)*(1-r13^2))
      copy_VarCoords(data2(1,:,:), reg12)
      ;偏相关的有效自由度
        ; r1_auto = esacr_n(data1,1,0)
        ; r2_auto = esacr_n(data2,1,0)
        
        ; one = r2_auto(0,:,:)
        ; one =1
        ; N = dimsizes(data1&time)*(one- r1_auto(0)*r2_auto(0,:,:))/(one+r1_auto(0)*r2_auto(0,:,:))
      N = dimsizes(data1)
      ; N = equiv_sample_size(data1,0.05,0)

      prob_mam = rtest(corr, N-1, 0)
      copy_VarCoords(data2(0,:,:), prob_mam)
      reg12@corr = corr
      reg12@prob = prob_mam

return reg12
end 

function partial_reg2(data1,data2,data3) ;;;data1:index data2: varible data3:effect want to remove
begin 
    ;;data1 :x1 1d data2 :y 3d; data3: x2 1d
      r13 = regCoef_n(data3, data1, 0, 0)

      data1 = (/data1 - (r13*data3 + r13@yintercept)/)


      r23 = regCoef_n(data3, data2, 0, 0)
      r23_3d = conform_dims(dimsizes(data2), r23, (/1,2/))
      
      m = reshape(r23@yintercept,dimsizes(data2)) 
      data3_3d = conform(data2,data3,(/0/))

      data2 = (/data2 - (r23_3d*data3_3d + m)/)
      
      reg12 = regCoef_n(data1, data2,0, 0)
      copy_VarCoords(data2(1,:,:),reg12)
      prob_mam = rc_test(reg12)
      copy_VarCoords(data2(0,:,:), prob_mam)
      reg12@prob = prob_mam

return reg12
end 


function split_12(data)
begin
    dim = dimsizes(data) 
    yearnum = dim(0)/12

    data_12 = new((/12,yearnum,dim(1),dim(2)/),float)
    do i=0,11   
    data_12(i,:,:,:) = data(i::12,:,:)
    end do 
return data_12
end 

;;**************************************************************
;;;main program
;;***************************************************************
begin
 ;;;time span
    startmon =198001
    endmon = 202012

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
    
    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")   
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p=f2->precip(timeind,:,:)
    printVarSummary(p)
    p= (/rmMonAnnCycTLL(p)/)
;;;; 


    ; f2=addfile("/home/yangsong3/data-observation/NCEP_NCAR/monthly/prate.sfc.mon.mean.nc", "r")
    ; date:=cd_calendar(f2->time, 1)
    ; timeind:=ind(date.ge.startmon.and.date.le.endmon)
    ; p_r=f2->prate(timeind,{-90:90},{0:360})
    ; p_r = (/p_r*86400/)

    ; ;;;处理数据
    ; ;;; 水平插值 lat:-88.75：88.75 ； lon:1.25 2.5×2.5度

    ; yo = fspan(-88.75, 88.75, 72)
    ; xo = fspan(1.25, 358.75, 144)
    ; p = linint2_Wrap(p_r&lon, p_r&lat, p_r, True, xo, yo, 0)

    ; delete(p&LAT)
    ; delete(p&LON)
    ; p!1 = "lat"
    ; p!2 = "lon"
    ; p&lat = yo
    ; p&lon = xo
    ; p&lat@units = "degree_north"
    ; p&lon@units = "degree_east"
  
    ; p= (/rmMonAnnCycTLL(p)/)
    ; printVarSummary(p)
; ;;;; 

;;;
    f3=addfile("/home/ys17-19/lsh/Project/MSE-EP/GRL/evapr_oaflux_chazhi.nc", "r")
    date:=cd_calendar(f3->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    eva=f3->evapr(timeind,:,:)
    eva = (/eva*10/365/) ;;convert to mm/day
    ; eva = (/rmMonAnnCycTLL(eva)/)

    p_e  = p
    p_e = (/p - eva/)
;;;;

;;;;;choose season 

 season = "MAM"



  ;;;read era data
   
        path_wind = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
        var_name = (/"shum.mon.mean.nc"/)

        filein = path_wind + var_name 
    
        ;;;region 
          lat1 = -50
          lat2 = 50
          lon1 = 0
          lon2 = 360
        ;;;

        ;;;time span
        startmon =198001
        endmon = 202012

        startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
        endyear_r = str_split_by_length(endmon, 4 )

        startyear = stringtoint(startyear_r(0))
        endyear = stringtoint(endyear_r(0))

        yearnum = endyear - startyear +1
        ;;
        f3 = addfile(filein(0),"r")
        date :=cd_calendar(f3->time, 1)
        timeind :=ind(date.ge.startmon.and.date.le.endmon)
        q = f3->shum(timeind,:,{lat1:lat2},{lon1:lon2}) 
        q= q/1000.0  



      ;;;MSE 
      plev = q&level
      ptop = min(plev)
      psfc = 1013.0
      dz = plev
      dz = dpres_plevel(plev*100, psfc*100, ptop*100, 0)
      dz&level = plev
      printVarSummary(dz)

    ;;;;;;;;;;;;;;;;;;;
    q_col =  dim_sum_wgt_n_Wrap(q(:,{1000:500},:,:), dz({1000:500}), 1, 1)
    q_col = (/rmMonAnnCycTLL(q_col)/)

    q_col = (/center_finite_diff_n(q_col, 1, False, 0, 0)/) ;;;mm/month

    q_col = q_col/(9.8*30)

    q_12 = split_12(q_col)

    q_annual = month_to_season(q_col(:,{lat1:lat2},{lon1:lon2}), season)
    q_annual = (/dim_avg_n(q_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)

    q_annual_rc = regCoef_n(ispan(1,yearnum,1), q_annual, 0,0)
    q_annual_rc = (/q_annual_rc*yearnum/) 
    copy_VarCoords(q_annual(1,:,:), q_annual_rc)
    q_annual_rc_mean  =  avg(q_annual_rc({0:20},{250:280}))

    q_annual = (/dtrend_msg_n(ispan(1,yearnum,1),q_annual,False,True,0)/)
;;;;;


;;;;;
    diri="/home/ys17-19/lsh/Project/MSE-EP/GRL/"
    fils1=systemfunc("ls "+diri+"3D-JFMAM-Q-term-8-NEP-2020.nc")
    f1=addfile(fils1, "r")

    series_Q = f1->series_Q
    series_dyna_adv_q = f1->series_dyna_adv_q
    series_dyna_con_q = f1->series_dyna_con_q

    series_thermaldyna_adv_q = f1->series_thermaldyna_adv_q
    series_thermaldyna_con_q = f1->series_thermaldyna_con_q

    series_quad_adv_q = f1->series_quad_adv_q
    series_quad_con_q = f1->series_quad_con_q
    series_eddy_q = f1->series_eddy_q
  
print("=====finish read data==========")



;;;;;;split to 12 month  is a array size (12,yearnum)
   ssta_12 = split_12(ssta)
   p_12 = split_12(p_e)
   ; eva_12 = split_12(eva_12)

   
;;;;;
   season = "MAM"

  ;;;region
   lat1 = -50
   lat2 = 50
   lon1 = 0
   lon2 = 360
  ;;;

  p_annual = month_to_season(p(:,{lat1:lat2},{lon1:lon2}), season)
  
  p_annual = (/dim_avg_n(p_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)

  ssta_annual = month_to_season(ssta(:,{lat1:lat2},{lon1:lon2}), season)
  
  ssta_annual = (/dim_avg_n(ssta_12(0:4,:,{lat1:lat2},{lon1:lon2}),0)/)

  ssta_index_rc = regCoef_n(ispan(1,dimsizes(ssta_annual&time),1), wgt_areaave(ssta_annual(:,{0:30},{280:340}), area({0:30}), 1.0, 0) , 0, 0)
  ssta_index_rc  =ssta_index_rc*yearnum


  print("*******ssta_index_rc="+ssta_index_rc)


  ssta_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),ssta_annual,False,True,0)/)
 
  ssta_index = wgt_areaave_Wrap(ssta_annual(:,{10:30},{280:330}), area({10:30}), 1.0, 0)
  ; ssta_index = wgt_areaave_Wrap(ssta_annual(:,{-10:20},{40:100}), area({-10:20}), 1.0, 0)
  ; ssta_index = wgt_areaave_Wrap(ssta_annual(:,{0:20},{250:280}), area({0:20}), 1.0, 0)  
  ; ssta_local_index = dim_avg_n_Wrap(ssta_annual(:,{0:20},{250:280}), (/1,2/))
  ssta_local_index = dim_avg_n_Wrap(ssta_annual(:,{-5:5},{190:240}), (/1,2/))
  ; ssta_local_index = wgt_areaave_Wrap(ssta_annual(:,{10:30},{280:330}), area({10:30}), 1.0, 0)
  ; ssta_local_index = ssta_index
  


  ;;;trend 
    series_Q_trend = regCoef_n(ispan(1,yearnum,1),series_Q,0,0)
    series_Q_trend  = (/series_Q_trend*yearnum/)
    copy_VarCoords(series_Q(1,:,:), series_Q_trend)

    series_dyna_con_q_rc = regCoef_n(ispan(1,yearnum,1),series_dyna_con_q,0,0)
    series_dyna_con_q_rc  = (/series_dyna_con_q_rc*yearnum/)
    copy_VarCoords(series_dyna_con_q(1,:,:), series_dyna_con_q_rc)


    p_rc =  regCoef_n(ispan(1,yearnum,1),p_annual,0,0)
    p_rc = (/p_rc*yearnum/)
    copy_VarCoords(p_annual(1,:,:), p_rc)

    p_rc_mean = avg(p_rc({0:20},{250:280}))


  ;;;



  ;;;;去趋势
  p_annual =  (/dtrend_msg_n(ispan(1,yearnum,1),p_annual,False,True,0)/)

  series_Q =  (/dtrend_msg_n(ispan(1,yearnum,1),series_Q,False,True,0)/)
  series_dyna_adv_q =  (/dtrend_msg_n(ispan(1,yearnum,1),series_dyna_adv_q,False,True,0)/)
  series_dyna_con_q =  (/dtrend_msg_n(ispan(1,yearnum,1),series_dyna_con_q,False,True,0)/)
  series_thermaldyna_adv_q =  (/dtrend_msg_n(ispan(1,yearnum,1),series_thermaldyna_adv_q,False,True,0)/)
  series_thermaldyna_con_q =  (/dtrend_msg_n(ispan(1,yearnum,1),series_thermaldyna_con_q,False,True,0)/)

  series_quad_adv_q =  (/dtrend_msg_n(ispan(1,yearnum,1),series_quad_adv_q,False,True,0)/)
  series_quad_con_q =  (/dtrend_msg_n(ispan(1,yearnum,1),series_quad_con_q,False,True,0)/)

  series_eddy_q =  (/dtrend_msg_n(ispan(1,yearnum,1),series_eddy_q,False,True,0)/)
  ;;;;;

  dim = dimsizes(series_Q)

  data1 = new((/2,dim(0),dim(1),dim(2)/),typeof(series_Q))
  
  ; data = (/series_Q, series_dyna_con_q ,series_dyna_adv_q,series_thermaldyna_con_q,series_thermaldyna_adv_q,series_quad_con_q,\
  ;   series_quad_adv_q, series_eddy_q/)

  data1 = (/series_Q, series_dyna_con_q/)

  data2 = (/series_Q_trend,series_dyna_con_q_rc/)
  

  prob_rc_Q = rc_test(series_Q_trend)
  prob_rc_dyna_con_q_rc = rc_test(series_dyna_con_q_rc)

  data3 = (/prob_rc_Q,prob_rc_dyna_con_q_rc/)
  


  
    wks   = gsn_open_wks("pdf", "Q-p-e-observerd trend and linear congruent")
    plot=new(6,graphic)
    plot1 = plot

    str = ispan(1,6,1)+ ""
    str2 = ispan(1,6,1)+ "a"
 
    iplot = new(6,graphic)
    leftsting = (/"(a) P-E","(b) P-E","(c) DIVQ","(d) DIVQ","(e) DYN(CON)","(f) DYN(CON)"/)
    rightString = (/"Observed trend","Trends linearly congruent with ATL warming",\
                   "Observed trend","Trends linearly congruent with ATL warming",\
                   "Observed trend","Trends linearly congruent with ATL warming"/)

do i=0,5

  if i.eq.0 then 
  r_sst_Q := p_rc 
  prob_p1 := rc_test(p_rc)
  end if 

  if i.eq.1 then 

    r_sst_Q := partial_reg(ssta_index,q_annual,ssta_local_index)
    r_sst_Q = (/r_sst_Q*ssta_index_rc/)
    copy_VarCoords(q_annual(1,:,:), r_sst_Q) 

    prob_p1 := r_sst_Q@prob
    copy_VarCoords(q_annual(1,:,:), prob_p1) 

    print("percent="+avg(r_sst_Q({0:20},{250:280}))/q_annual_rc_mean)
    print(avg(r_sst_Q({0:20},{250:280})))
    print(q_annual_rc_mean)
    exit()
  end if

  if mod(i,2).eq.0.and.i.gt.1 then 
    r_sst_Q := data2(i/2-1,:,:)
    copy_VarCoords(series_Q(1,:,:), r_sst_Q) 
    prob_p1 := data3(i/2-1,:,:)
    copy_VarCoords(series_Q(1,:,:),prob_p1) 

  end if 
   
  if mod(i,2).ne.0.and.i.gt.1 then 
    data_select :=data1(i/3-1,:,:,:)
    copy_VarCoords(series_Q, data_select) 
    ;;;correlation with sst
    ; r_sst_p := regCoef_n(ssta_index, p_annual, 0, 0)
    ; copy_VarCoords(p_annual(0,:,:), r_sst_p)
    ; r_sst_p = (/r_sst_p*ssta_index_rc/)
    r_sst_Q := partial_reg(ssta_index,data_select,ssta_local_index)
    r_sst_Q = (/r_sst_Q*ssta_index_rc/)
    prob_p1 := r_sst_Q@prob
    copy_VarCoords(series_Q(1,:,:), prob_p1) 
     ;;;;
   end if
 
if i.eq.3 then 
    percent = avg(r_sst_Q({0:20},{250:280}))/avg(series_Q_trend({0:20},{250:280}))
    print(percent)
    exit()
end if
  

;;;;
   

;;;;;
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

  res@lbLabelBarOn= False

 
  res@pmLabelBarOrthogonalPosF = 0.2

  res@lbLabelStride=1
  res@lbOrientation="horizontal"
  res@mpLimitMode           = "LatLon"

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

  ; res@gsnScalarContour        = True
  res@mpFillOn = True
  res@mpFillDrawOrder = "PostDraw"


  res@cnFillPalette = "precip_diff_12lev"
  ; res@cnFillPalette = "CBR_drywet"
  res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels


  res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; res@cnLevels := (/-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4/)*0.5
  ; if i.gt.1 then 
  res@cnLevels := (/-0.25,-0.2,-0.15,-0.1,0.1,0.15,0.2,0.25/)*8
  ; else 
  ; res@cnLevels := (/-0.25,-0.2,-0.15,-0.1,0.1,0.15,0.2,0.25/)*4
  ; end if

  res@lbLabelStride=1


  res@cnExplicitLabelBarLabelsOn = True
  ; res@lbLabelStrings =  res@cnLevels

  res@lbLabelFontHeightF = 0.016
  res@lbAutoManage = False

 
  ; res@cnFillPalette = "BlRe"
  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

  ; res@cnMinLevelValF       =   -0.5
  ; res@cnMaxLevelValF       =   0.5
  ; res@cnLevelSpacingF      =  0.1
  ;   res@vcExplicitLabelBarLabelsOn = True
  ; res@stExplicitLabelBarLabelsOn = True
  ; res@vcGlyphStyle="CurlyVector"
  ; res@vcGlyphStyle="LineArrow"
  ; res@vcMonoLineArrowColor = True
  ; res@vcLineArrowColor = "black"
  ; res@vcMinDistanceF           = 0.02          ; thin vectors

  ; res@vcRefLengthF             = 0.05    ; define length of v
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  ; res@vcVectorDrawOrder="PostDraw"

  ; ; resv@vcMonoLineArrowColor = True
  ; ; resv@vcLineArrowColor = "blue"
  ; res@vcLineArrowThicknessF=1.2
  ; ; res@vcLineArrowHeadMaxSizeF=0.01
  ; ; res@vcLineArrowHeadMinSizeF=0.01
  ; res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  ; res@vcLineArrowHeadMinSizeF = 0.006

  ; res@vcRefMagnitudeF          = 0.2*10
  ; ; res@vcMinMagnitudeF          =.3
  ; res@gsnLeftStringFontHeightF = 0.0185
  ; ; res@gsnLeftStringFontHeightF = 0.02

  ; res@vcRefAnnoString1On = True
  ; res@vcRefAnnoString1 = "2.0 m/s"
  ; res@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
  ; res@vcRefAnnoParallelPosF = 1.02
  ; res@vcRefAnnoFontHeightF = 0.013

  ; res@vcRefAnnoBackgroundColor = "Transparent"
  ; res@vcRefAnnoPerimOn = False

  res@gsnLeftString = leftsting(i)
  res@gsnRightString = rightString(i)
  res@gsnLeftStringFontHeightF = 0.015
  res@gsnRightStringFontHeightF = 0.015
  ; res@gsnLeftString = "The components of JFMAM shum & 850UV trend linked to TNA"

  ; res@gsnLeftStringOrthogonalPosF  = -0.12

  plot(i) = gsn_csm_contour_map(wks,r_sst_Q, res)

  
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
    lon1 = 250
    lon2 = 280

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)

    plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
    ;;;
    pgres@gsLineColor="red"

    lat1 = 10
    lat2 = 30
    lon1 = 280
    lon2 = 330

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
    res_pro@cnLevelSpacingF     = 0.05       ; set contour spacing

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

    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.05,opt)
    

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
  pres@gsnPanelLabelBar = True

  gsn_panel(wks, plot, (/4,2/), pres)
  ; draw(plot(0))
  ; frame(wks)

 print("done")

end