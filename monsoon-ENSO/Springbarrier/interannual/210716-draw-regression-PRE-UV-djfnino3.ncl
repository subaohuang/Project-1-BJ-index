;;*************************************************
;; draw circulation pattern during monsoon year 
;;*******************************************************************

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


begin
    startmon = 197901
    endmon = 201912
    startyear = 1979
    endyear = 2019

;;;读入HadiSST海温数据
    diri="/home/yangsong3/data-observation/linshh/data/sst/"
    fils1=systemfunc("ls "+diri+"HadISST_sst.nc")
    f1=addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)  ;;; 1900-2019

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)
    
    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)

    ;;detrend 
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))

    nino34_djf = month_to_season(nino34, "DJF")

    nino34_djf(0:dimsizes(nino34_djf)-2) = nino34_djf(1:) 
    nino34_djf(dimsizes(nino34_djf)-1) = nino34(dimsizes(nino34)-1) 

    std_nino = dim_stddev(nino34_djf)
    sst_total_year = ispan(startyear, endyear, 1)
    nino_year = sst_total_year(ind(nino34_djf.ge.std_nino))
    lanina_year = sst_total_year(ind(nino34_djf.le.-1*std_nino))

;;;


;;;读入ERA wind data

    fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
    ; fils1 = "/home/yangsong3/data-observation/20thCRV2/uwnd.mon.mean.nc"
    f1 = addfile(fils1, "r")
    time1 = f1->time
    date1=cd_calendar(time1, 1)

    endmon = endmon                 ;;; 1979-2019
    startmon = startmon 
    timeind1=ind(date1.ge.startmon.and.date1.le.endmon)  

 
;;;;U,V

    u_200  = short2flt(f1->u(timeind1,{200},:,:))
    u_850  = short2flt(f1->u(timeind1,{850},:,:))

    ; u_200  = f1->uwnd(timeind1,{200},:,:)
    ; u_850  = f1->uwnd(timeind1,{850},:,:)
    u_850 = (/rmMonAnnCycTLL(u_850)/)
    u_200 = (/rmMonAnnCycTLL(u_200)/)
;;;

   ; fils2 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/v_mon_1x1_1979_2020.nc"
   ; f2 = addfile(fils2, "r")
   ; time2 = f2->time
   ; date2 = cd_calendar(time2, 1)
   ; timeind2 = ind(date2.ge.197901.and.date2.le.201912)
;;;;



;;;
    ; fils2 = "/home/yangsong3/data-observation/20thCRV2/vwnd.mon.mean.nc"
    fils2  = "/home/yangsong3/data-observation/ERA5-monthly/pressure/v_mon_1x1_1979_2020.nc"
    f2 = addfile(fils2, "r")
    time2 = f2->time
    date2 = cd_calendar(time2, 1)

    timeind2 = ind(date2.ge.startmon.and.date2.le.endmon)  ;;;风场数据从1979年开始

;;;;U,V

    ; v_200  = f2->vwnd(timeind1,{200},:,:)
    ; v_850  = f2->vwnd(timeind1,{850},:,:)

    ; v_850 = (/rmMonAnnCycTLL(v_850)/)
    ; v_200 = (/rmMonAnnCycTLL(v_200)/)
    v_200  = short2flt(f2->v(timeind2,{200},:,:))
    v_850  = short2flt(f2->v(timeind2,{850},:,:))

    v_850 = (/rmMonAnnCycTLL(v_850)/)
    v_200 = (/rmMonAnnCycTLL(v_200)/)

    u_850 = (/dtrend_msg_n(ispan(1,dimsizes(u_850&time),1),u_850,False,True,0)/)
    v_850 = (/dtrend_msg_n(ispan(1,dimsizes(v_850&time),1),v_850,False,True,0)/)
;;;



;;;;

    f2=addfile("/home/yangsong3/data-observation/GPCP/precip.mon.mean.nc", "r")  
    date:=cd_calendar(f2->time, 1)
    timeind:=ind(date.ge.startmon.and.date.le.endmon)
    p=f2->precip(timeind,:,:)
    printVarSummary(p)
  
    p = (/rmMonAnnCycTLL(p)/)

    ;;detrend 
    p = (/dtrend_msg_n(ispan(1,dimsizes(p&time),1),p,False,True,0)/)



; ;;; JJ mean 
;     u_850_june =  dim_avg_n_Wrap(u_850(5::12,:,:), 0)   ;;1979-2019
;     u_850_july =  dim_avg_n_Wrap(u_850(6::12,:,:), 0)

;     u_850_jj = u_850_june
;     u_850_jj = (/(u_850_june+u_850_july)*0.5/)

;     v_850_june =  dim_avg_n_Wrap(v_850(5::12,:,:), 0)
;     v_850_july =  dim_avg_n_Wrap(v_850(6::12,:,:), 0)

;     v_850_jj = v_850_june
;     v_850_jj = (/(v_850_june+v_850_july)*0.5/)

    
;     ssta_june =  dim_avg_n_Wrap(ssta(5::12,:,:), 0)
;     ssta_july =  dim_avg_n_Wrap(ssta(6::12,:,:), 0)

;     ssta_jj = ssta_june
;     ssta_jj = (/(ssta_june+ssta_july)*0.5/)
; ;;;;;



;;;;read SPB index data 
   ; opt1="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB_index_interannual.ascii"
    opt1="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB_index_interannual_1980-2019_method1.ascii" 
    spb_year = asciiread(opt1, -1, "float")  ;;1980-2019 ;;;



;;;strong and weak SPB year 
    total_year = ispan(startyear+1,endyear,1) ;;1980-2019
    inds1 = ind(spb_year.ge.0.5)
    indw1= ind(spb_year.le.-0.5)
    indnor = ind(spb_year.gt.-0.5.and.spb_year.lt.0.5)

    ; inds = ind(spb_year.ge.1.0)
    ; indw = ind(spb_year.le.-1.0)

    ; indnor = ind(spb_year.gt.-1.0.and.spb_year.lt.1.0)

    monstryear = total_year(inds1)
    monwekyear = total_year(indw1)
    noryear = total_year(indnor)
;;;;


 ;;; 
    sst_year = ispan(startyear,endyear,1)
    sst_year@_FillValue = 99
  
    ind_sst_s = get1Dindex_Exclude(sst_year, monstryear)    ;;;去掉某一年数据,数据从1979年开始-2019年结束
    ind_sst_w = get1Dindex_Exclude(sst_year, monwekyear)    ;;;去掉某一年数据


    nino34_djf_s = nino34_djf(ind_sst_s)  ;;指数从1979-2019，期间去掉几个异常年份
    nino34_djf_w = nino34_djf(ind_sst_w)


;;;


    plot = new(9,graphic)
    plot1 = new(9,graphic)
    plot2 = new(9,graphic)

    iplot = plot
    iplot1 = plot
    iplot2 = plot
    ; iplot1 = plot1
    str = ispan(1,9,1) + ""

    ; wks = gsn_open_wks("pdf", "UVSST_regression-SPB-before")
    wks = gsn_open_wks("pdf", "UVPRE_regression-DJFnino34-removeSTRSPByear")
    wks1 = gsn_open_wks("pdf", "UVPRE_regression-DJFnino34-removeWEAKSPByear")
    wks2 = gsn_open_wks("pdf", "UVPRE_regression-DJFnino34-allyear")
;;;


; ;;;计算回归系数 JJ UV regressed onto DJF(-1) nino34
  

;   rc_u_s =  regCoef_n(nino34_djf_s(0:yearnum-2), u_850_jj_s(1:yearnum-1,:,:) , 0, 0)
;   rc_u_w =  regCoef_n(nino34_djf_w(0:yearnum-2), u_850_jj_w(1:yearnum-1,:,:) , 0, 0)

;   rc_v_s =  regCoef_n(nino34_djf_s(0:yearnum-2), v_850_jj_s(1:yearnum-1,:,:) , 0, 0)
;   rc_v_w =  regCoef_n(nino34_djf_w(0:yearnum-2), v_850_jj_w(1:yearnum-1,:,:) , 0, 0)
  
;   rc_sst_s =  regCoef_n(nino34_djf_s(0:yearnum-2), ssta_jj_s(1:yearnum-1,:,:) , 0, 0)
;   rc_sst_w =  regCoef_n(nino34_djf_w(0:yearnum-2), ssta_jj_w(1:yearnum-1,:,:) , 0, 0)
  
 month = (/"Dec(-1)","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug"/)
do i=0,8
;; exclude anomaous SPB year ssta , u ,v 
   if i.le.0 then 
      ssta_month =  p(i+11:dimsizes(p&time)-1-12:12,:,:)  ;;; 1979 D 1980D .....2018 D
      u_850_month = u_850(i+11:dimsizes(u_850&time)-1-12:12,:,:)
      v_850_month = v_850(i+11:dimsizes(v_850&time)-1-12:12,:,:)
    else
      ssta_month =  p(i+11:dimsizes(p&time)-1:12,:,:)     ;;  1980 J 1981 J ......2019 J
      u_850_month = u_850(i+11:dimsizes(u_850&time)-1:12,:,:)
      v_850_month = v_850(i+11:dimsizes(v_850&time)-1:12,:,:)
    end if
;;;;
   
;;;;去掉某一年
  if i.eq.0 then 
    or_year = ispan(startyear,endyear-1,1)
    or_year@_FillValue = 99

    yearnum = endyear - startyear 
   
    inds = get1Dindex_Exclude(or_year, monstryear)    ;;;去掉某一年数据
    indw = get1Dindex_Exclude(or_year, monwekyear)    ;;;去掉某一年数据

  else 
    or_year = ispan(startyear+1,endyear,1)
    or_year@_FillValue = 99

    yearnum = endyear - startyear
   
    inds = get1Dindex_Exclude(or_year, monstryear)    ;;;去掉某一年数据
    indw = get1Dindex_Exclude(or_year, monwekyear)    ;;;去掉某一年数据
  end if 



;;;;;;;;
     u_850_jj_s = u_850_month(inds,:,:)    
     v_850_jj_s = v_850_month(inds,:,:)    

     u_850_jj_w = u_850_month(indw,:,:)    
     v_850_jj_w = v_850_month(indw,:,:)

     ssta_jj_s = ssta_month(inds,:,:)    
     ssta_jj_w = ssta_month(indw,:,:)
;;;;





    ssta_month_str = regCoef_n(nino34_djf_s(0:dimsizes(nino34_djf_s)-2), ssta_jj_s, 0, 0)
    copy_VarCoords(ssta_month(0,:,:), ssta_month_str)

    u_month_str = regCoef_n(nino34_djf_s(0:dimsizes(nino34_djf_s)-2), u_850_jj_s, 0, 0)
    copy_VarCoords(u_850_month(0,:,:), u_month_str)
    
    v_month_str = regCoef_n(nino34_djf_s(0:dimsizes(nino34_djf_s)-2), v_850_jj_s, 0, 0)
    copy_VarCoords(v_850_month(0,:,:), v_month_str)
    
   
    ssta_month_weak = regCoef_n(nino34_djf_w(0:dimsizes(nino34_djf_w)-2), ssta_jj_w, 0, 0)
    copy_VarCoords(ssta_month(0,:,:), ssta_month_weak)

    u_month_weak = regCoef_n(nino34_djf_w(0:dimsizes(nino34_djf_w)-2), u_850_jj_w, 0, 0)
    copy_VarCoords(u_850_month(0,:,:), u_month_weak)
    
    v_month_weak = regCoef_n(nino34_djf_w(0:dimsizes(nino34_djf_w)-2), v_850_jj_w, 0, 0)
    copy_VarCoords(v_850_month(0,:,:), v_month_weak)
   

    ssta_month_all = regCoef_n(nino34_djf(0:dimsizes(nino34_djf)-2), ssta_month, 0, 0)
    copy_VarCoords(ssta_month(0,:,:), ssta_month_all)

    u_month_all = regCoef_n(nino34_djf(0:dimsizes(nino34_djf)-2), u_850_month, 0, 0)
    copy_VarCoords(u_850_month(0,:,:), u_month_all)

    v_month_all = regCoef_n(nino34_djf(0:dimsizes(nino34_djf)-2), v_850_month, 0, 0)
    copy_VarCoords(v_850_month(0,:,:), v_month_all)


   ; ;;;;后期
   ;   if i.le.6 then 
   ;    ssta_month =  ssta(i+5+12:dimsizes(ssta&time)-1-12:12,:,:)
   ;    u_850_month = u_850(i+5+12:dimsizes(u_850&time)-1-12:12,:,:)
   ;    v_850_month = v_850(i+5+12:dimsizes(v_850&time)-1-12:12,:,:)
   ;  else
   ;    ssta_month =  ssta(i+5+12:dimsizes(ssta&time)-1:12,:,:)
   ;    u_850_month = u_850(i+5+12:dimsizes(u_850&time)-1:12,:,:)
   ;    v_850_month = v_850(i+5+12:dimsizes(v_850&time)-1:12,:,:)
   ;  end if
   
   ;  ssta_month_str = regCoef_n(spb_year(0:dimsizes(spb_year)-2), ssta_month, 0, 0)
   ;  copy_VarCoords(ssta_month(0,:,:), ssta_month_str)

   ;  u_month_str = regCoef_n(spb_year(0:dimsizes(spb_year)-2), u_850_month, 0, 0)
   ;  copy_VarCoords(u_850_month(0,:,:), u_month_str)
    
   ;  v_month_str = regCoef_n(spb_year(0:dimsizes(spb_year)-2), v_850_month, 0, 0)
   ;  copy_VarCoords(v_850_month(0,:,:), v_month_str)
    

    prob_sst_str = rc_test(ssta_month_str)
    prob_u_str = rc_test(u_month_str)
    prob_v_str = rc_test(v_month_str)


    prob_sst_weak = rc_test(ssta_month_weak)
    prob_u_weak = rc_test(u_month_weak)
    prob_v_weak = rc_test(v_month_weak)

    
    prob_sst_all = rc_test(ssta_month_all)
    prob_u_all = rc_test(u_month_all)
    prob_v_all = rc_test(v_month_all)


  
   ; ssta_month_str = where(prob_sst_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)

   ; ssta_month_weak = where(prob_sst_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)


    u_month_str = where(prob_u_str.le.0.1.or.prob_v_str.le.0.1, u_month_str, u_month_str@_FillValue)

    v_month_str = where(prob_v_str.le.0.1.or.prob_u_str.le.0.1, v_month_str, v_month_str@_FillValue)

    u_month_weak = where(prob_u_weak.le.0.1.or.prob_v_weak.le.0.1, u_month_weak, u_month_weak@_FillValue)

    v_month_weak = where(prob_v_weak.le.0.1.or.prob_u_weak.le.0.1, v_month_weak, v_month_weak@_FillValue)

    u_month_all = where(prob_u_all.le.0.1.or.prob_v_all.le.0.1, u_month_all, u_month_all@_FillValue)

    v_month_all = where(prob_v_all.le.0.1.or.prob_u_all.le.0.1, v_month_all, v_month_all@_FillValue)
    

  
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  res@cnFillPalette = "precip_diff_12lev"
  ; res@cnFillPalette = "temp_diff_18lev"

  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -20
  res@mpMaxLatF             = 50
  res@mpMinLonF             = 30
  res@mpMaxLonF             = 250

  res@mpFillOn = False

  res@mpCenterLonF          = (res@mpMinLonF+res@mpMaxLonF)/2.0

  res@tmXTOn=False
  res@tmYROn=False
  res@tmXBOn=True
  res@tmYROn=False

  res@tmXBMode="Explicit"
  res@tmYLMode="Explicit"

   
  res@tmXBLabelFontHeightF=0.015
  res@tmYLLabelFontHeightF=0.015
  res@tmXBMajorThicknessF=1.0
  res@tmYLMajorThicknessF=1.0
  res@tmBorderThicknessF=1.0

  ; res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -1.5 
  ; res@cnMaxLevelValF       =   1.5
  ; res@cnLevelSpacingF      =  0.5

  ; res@cnMinLevelValF       =   -0.8
  ; res@cnMaxLevelValF       =   0.8
  ; res@cnLevelSpacingF      =  0.1

  res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  res@cnLevels := (/-1.5,-1.0,-0.5,-0.1,0.1,0.5,1.0,1.5/)
  ; res@cnLevels := (/-20,-15,-10,-5,5,10,15,20/)*0.01
  ; res@cnInfoLabelOn = False
  
  res@gsnLeftString = "PRE&850UV " + month(i)
  res@gsnLeftStringFontHeightF = 0.02

  ; res@mpShapeMode = "FreeAspect"
  ; res@vpHeightF = 0.2
  ; res@vpWidthF = 0.7

  res@gsnScalarContour        = True

  ;;;;;;;LineArrow
  res@vcMinDistanceF           = 0.02           ; thin vectors
  res@vcMinMagnitudeF          =0.4
  res@vcRefMagnitudeF          = 4.0
  res@vcRefLengthF             = 0.1     ; define length of v
  ; res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
  res@vcLineArrowThicknessF = 1.0
  ; res@vcLineArrowHeadMaxSizeF=0.008
  ; res@vcLineArrowHeadMinSizeF=0.008
  res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
  res@vcLineArrowHeadMinSizeF = 0.006
    
    
    if i.eq.8 then 
      res@vcRefAnnoOn = True
      res@vcRefAnnoOrthogonalPosF = 0.4
    else 
      res@vcRefAnnoOn = False
    end if

  
  res@gsnRightString = "Persistent ENSO (n=41-8)"
  plot(i) = gsn_csm_vector_scalar_map(wks,u_month_str, v_month_str, ssta_month_str, res)
  res@gsnRightString = " Transitive ENSO (n=41-8)"
  plot1(i) = gsn_csm_vector_scalar_map(wks1,u_month_weak, v_month_weak, ssta_month_weak, res)
  res@gsnRightString = "ALL ENSO"
  plot2(i) = gsn_csm_vector_scalar_map(wks2,u_month_all, v_month_all, ssta_month_all, res)



   ; ;;;draw box 
   ;  pgres=True
   ;  pgres@gsLineColor="black"
   ;  pgres@gsLineThicknessF=0.8
   ;  pgres@gsLineDashPattern=0

   ;  lat1 = -5
   ;  lat2 = 5
   ;  lon1 = 190
   ;  lon2 = 240

   ;  pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
   ;  pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ; ;;;

   
   ;  plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
 

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
    res_pro@cnFillDotSizeF = 0.001
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)
    ; res_pro@cnLevels = ispan(-6,6,1)
    res_pro@cnLineLabelsOn=False
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot(i) = gsn_csm_contour(wks,prob_sst_str,res_pro)
     iplot1(i) = gsn_csm_contour(wks1,prob_sst_weak,res_pro)
      iplot2(i) = gsn_csm_contour(wks2,prob_sst_all,res_pro)
   
  
    overlay(plot(i), iplot(i))
    overlay(plot1(i), iplot1(i))
    overlay(plot2(i), iplot2(i))

delete(ssta_month_str)
delete(u_month_str)
delete(v_month_str)
delete(ssta_month_weak)
delete(u_month_weak)
delete(v_month_weak)
delete([/ssta_month_all,u_month_all,v_month_all/])
end do
pres = True
pres@gsnPanelLabelBar  = True

; pres@pmLabelBarOrthogonalPosF = -0.008

; pres@lbLabelStride=1
gsn_panel(wks, plot, (/3,3/), pres)
gsn_panel(wks1, plot1, (/3,3/), pres)
gsn_panel(wks2, plot2, (/3,3/), pres)


end