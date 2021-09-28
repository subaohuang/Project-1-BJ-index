;;*************************************************
;; draw circulation pattern during monsoon year 
;;*******************************************************************

function one_test(pv_ano,east_time)
begin
    pv_east=dim_avg_n_Wrap(pv_ano(east_time,:,:),0)
    n=dimsizes(east_time)
    ; var_pv_ano_east    =  dim_variance_n(pv_ano(east_time,:,:), 0)    ; 13.477
    var_pv_ano_east    =  dim_variance_n(pv_ano(:,:,:), 0)    ; 13.477

    var_pv_ano_east  =  sqrt(var_pv_ano_east/n) 
    var_pv_ano_east  = where(var_pv_ano_east.eq.0, var_pv_ano_east@_FillValue, var_pv_ano_east)
    
    t_pv_ano_east = (pv_east)/var_pv_ano_east
    confi_pv_east = student_t(t_pv_ano_east , n-1)
    copy_VarCoords(pv_east, confi_pv_east)

 return confi_pv_east
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


    ; opt1="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB_index_interannual.ascii"
    opt1="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/interannual/SPB_index_interannual_1980-2019_method1.ascii" 
    spb_year = asciiread(opt1, -1, "float")  ;;1980-2019



;;;strong and weak monsoon year 
    total_year = ispan(startyear+1,endyear,1)
   

    ; inds = ind(spb_year.ge.0.5)
    ; indw = ind(spb_year.le.-0.5)

    ; indnor = ind(spb_year.gt.-0.5.and.spb_year.lt.0.5)
    inds = ind(spb_year.ge.1.0)
    indw = ind(spb_year.le.-1.0)

    indnor = ind(spb_year.gt.-1.0.and.spb_year.lt.1.0)

    monstryear = total_year(inds)
    monwekyear = total_year(indw)
    noryear = total_year(indnor)
    
   
    yearnum = dimsizes(nino34)/12 -1
    dim_sst = dimsizes(ssta)
    dim_u = dimsizes(u_850)
; exit()


; ;;;split to 12 month nino34_12 is a array size (12,yearnum)

   sst_year = ispan(startyear,endyear,1) 

   ind_ts1 =  get1Dindex(sst_year,monstryear-1)
   ind_tw1 =  get1Dindex(sst_year,monwekyear-1)
   ind_tnor1 =  get1Dindex(sst_year,noryear-1)
  
   u_year = ispan(startyear,endyear,1) 

   ind_us1 =  get1Dindex(u_year,monstryear-1)
   ind_uw1 =  get1Dindex(u_year,monwekyear-1)
   ind_unor1 =  get1Dindex(u_year,noryear-1)

  plot = new(12,graphic)
  plot1 = new(12,graphic)
  iplot = plot
  iplot1 = plot1
  str = ispan(1,12,1) + ""
 
  wks = gsn_open_wks("pdf", "UVSST_strong-SPB-before-1.0std")
  wks1 = gsn_open_wks("pdf", "UVSST_weak-SPB-before-1.0std")

 month = (/"Jun(-1)","Jul(-1)","Aug(-1)","Sep(-1)","Oct(-1)","Nov(-1)","Dec(-1)","Jan","Feb","Mar","Apr","May"/) 


do i=0,11
;;; strong and weak monsoon year ssta , u ,v 
   if i.le.6 then 
      ssta_month =  ssta(i+5:dimsizes(ssta&time)-1-12:12,:,:)
      u_850_month = u_850(i+5:dimsizes(u_850&time)-1-12:12,:,:)
      v_850_month = v_850(i+5:dimsizes(v_850&time)-1-12:12,:,:)
    else
      ssta_month =  ssta(i+5:dimsizes(ssta&time)-1:12,:,:)
      u_850_month = u_850(i+5:dimsizes(u_850&time)-1:12,:,:)
      v_850_month = v_850(i+5:dimsizes(v_850&time)-1:12,:,:)
    end if

    ssta_month_str = dim_avg_n_Wrap(ssta_month(ind_ts1,:,:),0)

    u_month_str = dim_avg_n_Wrap(u_850_month(ind_us1,:,:),0)

    v_month_str = dim_avg_n_Wrap(v_850_month(ind_us1,:,:),0)

    
    ssta_month_weak = dim_avg_n_Wrap(ssta_month(ind_tw1,:,:),0)

    u_month_weak = dim_avg_n_Wrap(u_850_month(ind_uw1,:,:),0)

    v_month_weak = dim_avg_n_Wrap(v_850_month(ind_uw1,:,:),0)

    prob_sst_str = one_test(ssta_month,ind_ts1)
    prob_u_str = one_test(u_850_month,ind_us1)
    prob_v_str = one_test(v_850_month,ind_us1)

    prob_sst_weak = one_test(ssta_month,ind_tw1)
    prob_u_weak = one_test(u_850_month,ind_uw1)
    prob_v_weak = one_test(v_850_month,ind_uw1)
   

   ; ssta_month_str = where(prob_sst_str.le.0.1, ssta_month_str, ssta_month_str@_FillValue)

   ; ssta_month_weak = where(prob_sst_weak.le.0.1, ssta_month_weak, ssta_month_weak@_FillValue)

    u_month_str = where(prob_u_str.le.0.1.or.prob_v_str.le.0.1, u_month_str, u_month_str@_FillValue)

    v_month_str = where(prob_v_str.le.0.1.or.prob_u_str.le.0.1, v_month_str, v_month_str@_FillValue)


    u_month_weak = where(prob_u_weak.le.0.1.or.prob_v_weak.le.0.1, u_month_weak, u_month_weak@_FillValue)

    v_month_weak = where(prob_v_weak.le.0.1.or.prob_u_weak.le.0.1, v_month_weak, v_month_weak@_FillValue)


  
;;**************************************************************
;画图
;;***************************************************************

  res=True
  res@gsnFrame=False
  res@gsnDraw=False

  res@cnFillPalette = "BlueWhiteOrangeRed"
  ; res@cnFillPalette = "temp_diff_18lev"

  res@cnFillOn=True
  res@cnLinesOn=False

  res@cnLineLabelsOn=False
  res@lbLabelBarOn=False
  res@lbOrientation="Vertical"

  res@mpLimitMode           = "LatLon"
  res@mpMinLatF             = -30
  res@mpMaxLatF             = 30
  res@mpMinLonF             = 0
  res@mpMaxLonF             = 360

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

  res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour level
  ; res@cnMinLevelValF       =   -0.4
  ; res@cnMaxLevelValF       =   0.4
  ; res@cnLevelSpacingF      =  0.1

  res@cnMinLevelValF       =   -0.8
  res@cnMaxLevelValF       =   0.8
  res@cnLevelSpacingF      =  0.1

  ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
  ; ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
  ; ; res@cnLevels := (/-12,-9,-6,-3,3,6,9,12/)*0.01
  ; res@cnLevels := (/-20,-15,-10,-5,5,10,15,20/)*0.01
  ; res@cnInfoLabelOn = False
  
  res@gsnLeftString = "SST&850UV " + month(i)
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
    
    
    if i.eq.11 then 
      res@vcRefAnnoOn = True
      res@vcRefAnnoOrthogonalPosF = 0.4
    else 
      res@vcRefAnnoOn = False
    end if

  
  res@gsnRightString = "Strong SPB"
  plot(i) = gsn_csm_vector_scalar_map(wks,u_month_str, v_month_str, ssta_month_str, res)
  res@gsnRightString = "Weak SPB"
  res@cnMinLevelValF       =   -0.6
  res@cnMaxLevelValF       =   0.6
  res@cnLevelSpacingF      =  0.1

  plot1(i) = gsn_csm_vector_scalar_map(wks1,u_month_weak, v_month_weak, ssta_month_weak, res)
  
  
   ;;;draw box 
    pgres=True
    pgres@gsLineColor="black"
    pgres@gsLineThicknessF=0.8
    pgres@gsLineDashPattern=0

    lat1 = -5
    lat2 = 5
    lon1 = 190
    lon2 = 240

    pgx_pool=(/lon1,lon2,lon2,lon1,lon1/)
    pgy_pool=(/lat1,lat1,lat2,lat2,lat1/)
   ;;;

   
   plot@$str(i)$ = gsn_add_polyline(wks, plot(i), pgx_pool, pgy_pool, pgres)
   plot1@$str(i)$ = gsn_add_polyline(wks, plot1(i), pgx_pool, pgy_pool, pgres)
   

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
  
    overlay(plot(i), iplot(i))
    overlay(plot1(i), iplot1(i))

end do
pres = True
pres@gsnPanelLabelBar  = True

; pres@pmLabelBarOrthogonalPosF = -0.008

; pres@lbLabelStride=1
gsn_panel(wks, plot, (/4,3/), pres)
gsn_panel(wks1, plot1, (/4,3/), pres)

end