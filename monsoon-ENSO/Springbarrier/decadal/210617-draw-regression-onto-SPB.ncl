;;*************************************************
;;draw SPB 10 years windows time series 
;;
;;*************************************************************
function cal_mean(length,nave,data,season)
begin
  dim = dimsizes(data)

  if season.eq."annual" then 
    data_sea = month_to_annual(data, 1) 
  else
    data_sea = month_to_season(data, season)
  end if 

  EAESM_series = new((/length-nave+1,dim(1),dim(2)/),float)
  do t = nave,length  
    EAESM  = dim_avg_n_Wrap(data_sea(t-nave:t-1,:,:),0)
    EAESM_series(t-nave,:,:) =  EAESM
    delete(EAESM)
  end do 
return EAESM_series
end 


begin
;;;read SPB data
    diri="/home/ys17-19/lsh/Project/monsoon-ENSO/Springbarrier/"
    ; fils1=systemfunc("ls "+diri+"SPB_series.nc")
    fils1=systemfunc("ls "+diri+"SPB_series_1899-2019-rmtrend.nc")
    ; fils1=systemfunc("ls "+diri+"SPB_series_1899-2019-rmtrend-fang-method.nc")

    f1=addfile(fils1, "r")

    SPB_series = f1 ->SPB_series
;;;


; ;;;读入ERA wind data

;    fils1 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/u_mon_1x1_1979_2020.nc"
;    f1 = addfile(fils1, "r")
;    time1 = f1->time
;    date1=cd_calendar(time1, 1)
;    timeind1=ind(date1.ge.197901.and.date1.le.201912)
 
  
; ;;;;U,V
;     u_200  = short2flt(f1->u(timeind1,{200},:,:))
;     u_850  = short2flt(f1->u(timeind1,{850},:,:))

;     u_850 = (/rmMonAnnCycTLL(u_850)/)
;     u_200 = (/rmMonAnnCycTLL(u_200)/)
; ;;;

;    fils2 = "/home/yangsong3/data-observation/ERA5-monthly/pressure/v_mon_1x1_1979_2020.nc"
;    f2 = addfile(fils2, "r")
;    time2 = f2->time
;    date2 = cd_calendar(time2, 1)
;    timeind2 = ind(date2.ge.197901.and.date2.le.201912)
 
  
; ;;;;U,V
;     v_200  = short2flt(f2->v(timeind2,{200},:,:))
;     v_850  = short2flt(f2->v(timeind2,{850},:,:))

;     v_850 = (/rmMonAnnCycTLL(v_850)/)
;     v_200 = (/rmMonAnnCycTLL(v_200)/)
; ;;;



;;;; read  HadISST data
    startmon = 189901
    endmon = 201912

    diri = "/home/yangsong3/data-observation/linshh/data/sst/"
    fils1 = systemfunc("ls "+diri+"HadISST_sst.nc")
    f1 = addfile(fils1, "r")
    date:=cd_calendar(f1->time, 1)
    lat = f1->latitude
    lon = f1->longitude

    timeind:=ind(date.le.endmon.and.date.ge.startmon)

    sst=lonFlip(f1->sst(timeind,:,:))
    printVarSummary(sst)

    ssta=(/rmMonAnnCycTLL(sst)/)
    copy_VarCoords(sst, ssta)
    ssta = (/dtrend_msg_n(ispan(1,dimsizes(ssta&time),1),ssta,False,True,0)/)
    ;;;
    ;;;
    delete(date)
    delete(timeind)
    
    nino34 =  dim_avg_n_Wrap(ssta(:,{5:-5},{190:240}), (/1,2/))
;;;

    season =(/"MAM","JJA","SON","DJF","annual"/)
    order = (/"(a)","(b)","(c)","(d)","(e)"/)
    wks = gsn_open_wks("pdf", "Corr-U,V,SST&SPB_intensity-fang-method")
    plot = new(5,graphic)
    iplot = plot
;;;;


do i=0,4
    ;;;get 10years running windows u,v,sst data 
    length = 2019-1899+1
    nave = 10 
    ;;;
    ; u_run = cal_mean(length,nave,u_850,season(i))
    ; v_run = cal_mean(length,nave,v_850,season(i))
    sst_run = cal_mean(length,nave,ssta,season(i))
    
;;;;

;;;corr 
    ; u_corr = escorc_n(SPB_series, u_run,0,0)
    ; v_corr = escorc_n(SPB_series, v_run,0,0)
    sst_corr = escorc_n(SPB_series,sst_run,0,0)

    ; copy_VarCoords(u_run(0,:,:), u_corr)
    ; copy_VarCoords(v_run(0,:,:), v_corr)
    copy_VarCoords(sst_run(0,:,:), sst_corr)
  

    yearnum = dimsizes(sst_run(:,1,1))

    ; prob_u = rtest(u_corr, yearnum, 0)
    ; prob_v = rtest(v_corr, yearnum, 0)
    prob_sst = rtest(sst_corr, yearnum, 0)

    ; copy_VarCoords(u_run(0,:,:), prob_u)
    ; copy_VarCoords(v_run(0,:,:), prob_v)
    copy_VarCoords(sst_run(0,:,:), prob_sst)

    ; u_corr = where(prob_u.le.0.05.or.prob_v.le.0.05, u_corr, u_corr@_FillValue)
    ; v_corr = where(prob_u.le.0.05.or.prob_v.le.0.05, v_corr, v_corr@_FillValue)
    ; sst_corr = where(prob_sst.le.0.01, sst_corr, sst_corr@_FillValue)
    ; delete(u_run)
    ; delete(v_run)
    delete(sst_run)
    
;;;


;;;draw plot 
    res=True
    res@gsnFrame=False
    res@gsnDraw=False

    res@cnFillOn=True
    res@cnLinesOn=False
    res@gsnAddCyclic = True

    res@cnLineLabelsOn=False
    res@lbLabelBarOn  = False

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

    ; res@tmXBLabelFontHeightF=0.02
    ; res@tmYLLabelFontHeightF=0.02
    ; res@tmXBMajorThicknessF=1.5
    ; res@tmYLMajorThicknessF=1.5
    ; res@tmBorderThicknessF = 1.5


    res@gsnScalarContour        = True
    res@mpFillOn = False

    ; res@cnFillPalette = "precip_diff_12lev"
    res@cnFillPalette = "BlueWhiteOrangeRed"
    ; res@cnLevelSelectionMode  = "ExplicitLevels" ; set manual contour levels
    

    ; res@cnFillColors= (/2,3,4,6,8,9,10,11,13/) - 2 
    ; res@cnLevels := (/-1.5,-1.0,-0.5,-0.2,0.2,0.5,1.0,1.5/)

    res@lbLabelStride=1

    ; res@cnExplicitLabelBarLabelsOn = True
    ; res@vcExplicitLabelBarLabelsOn = True
    ; res@stExplicitLabelBarLabelsOn = True
    ; res@lbLabelFontHeightF = 0.016
    ; res@lbAutoManage = False

    res@cnLevelSelectionMode = "ManualLevels"      ; set manual contour levels

    res@cnMinLevelValF       =   -1.0
    res@cnMaxLevelValF       =   1.0
    res@cnLevelSpacingF      =  0.2

    ; res@vcGlyphStyle="CurlyVector"
    res@vcGlyphStyle = "LineArrow"
    res@vcMonoLineArrowColor = True
    res@vcLineArrowColor = "black"
    res@vcMinDistanceF           = 0.02          ; thin vectors

    res@vcRefLengthF             = 0.02    ; define length of v
    res@vcRefAnnoOrthogonalPosF  = -1.0              ; move ref vector down -1.0在图的右上
    res@vcVectorDrawOrder="PostDraw"

    ; resv@vcMonoLineArrowColor = True
    ; resv@vcLineArrowColor = "blue"
    res@vcLineArrowThicknessF=1.2
    ; res@vcLineArrowHeadMaxSizeF=0.01
    ; res@vcLineArrowHeadMinSizeF=0.01
    res@vcLineArrowHeadMaxSizeF = 0.013           ;change the vector line arrows
    res@vcLineArrowHeadMinSizeF = 0.006

    res@vcRefMagnitudeF          = 1.0
 
    ; res@vcRefAnnoString1On = True
    ; res@vcRefAnnoString1 = "2 m/s"
    ; res@vcRefAnnoOrthogonalPosF  = -1.18              ; move ref vector down -1.0在图的右上
    ; res@vcRefAnnoParallelPosF = 1.02
    ; res@vcRefAnnoFontHeightF = 0.013

    ; res@vcRefAnnoBackgroundColor = "Transparent"
    ; res@vcRefAnnoPerimOn = False

    res@gsnLeftString = order(i)+" "+ season(i)+ "SPB&U,V,SST CORR"

    ; res@gsnLeftStringOrthogonalPosF  = -0.12

  
    ; plot(i) = gsn_csm_vector_scalar_map(wks,u_corr,v_corr,sst_corr, res)
    plot(i) = gsn_csm_contour_map(wks,sst_corr, res)

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


    iplot(i) = gsn_csm_contour(wks,prob_sst,res_pro)

    
    opt                      = True
    opt@gsnShadeFillType     = "pattern"
    opt@gsnShadeMid     = 17        ; stipple
    ; opt@gsnShadeFillDotSizeF =  0.0025
    ; opt@gsnShadeHigh   = "green"
    opt@gsnShadeFillScaleF   = 0.5         ; add extra density ,越大越係數

    iplot(i) = gsn_contour_shade(iplot(i) , 0,0.1,opt)
    
end do 

pres = True
pres@gsnPanelLabelBar = True
gsn_panel(wks, plot,(/3,2/), pres)
end 