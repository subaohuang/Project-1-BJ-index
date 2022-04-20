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

function mse_trend(data,year1,year2)
begin
  rc1 =regCoef_n(ispan(1, (year2-year1+1),1), data,0,0)
  copy_VarCoords(data(1,:,:), rc1)
  rc1 = rc1*(year2-year1+1)
  return rc1
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
begin


;;;;
  path_wind = "/home/yangsong3/data-observation/NCEP_NCAR/monthly/"
  var_name = (/"uwnd.mon.mean.nc","vwnd.mon.mean.nc","omega.mon.mean.nc"/)

  filein = path_wind + var_name 


  startmon = 198001
  endmon = 202012

  startyear = 1980
  endyear = 2020

  ;;;region 
  lat1 = -30
  lat2 = 30
  lon1 = 210
  lon2 = 280
  ;;;

  ; ;;;global region 
  ; lat1 = -50
  ; lat2 = 50
  ; lon1 = 0
  ; lon2 = 360
  ; ;;;
  startyear_r = str_split_by_length(startmon, 4)  ;; 将时间string分解成年和月
  endyear_r = str_split_by_length(endmon, 4 )

  startyear = stringtoint(startyear_r(0))
  endyear = stringtoint(endyear_r(0))

  yearnum = endyear - startyear +1
;;;;

;;;
  f3 = addfile(filein(2),"r")
  date:=cd_calendar(f3->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  w = dim_avg_n_Wrap(f3->omega(timeind,{1000:100},{lat1:lat2},{lon1:lon2}),3)

  w = (/rmMonAnnCycTLL(w)/)
  w = (/dtrend_msg_n(ispan(1,dimsizes(w&time),1),w,False,True,0)/)
;;;;

;;;;
  level = w&level
  

  f2 = addfile(filein(1),"r")
  date:=cd_calendar(f2->time, 1)
  timeind:=ind(date.ge.startmon.and.date.le.endmon)

  level_v = f2->level

  levind = get1Dindex(level_v,level)

  v = dim_avg_n_Wrap(f2->vwnd(timeind,levind,{lat1:lat2},{lon1:lon2}),3)

  v = (/rmMonAnnCycTLL(v)/)
  v = (/dtrend_msg_n(ispan(1,dimsizes(v&time),1),v,False,True,0)/)
;;;;


;;;;;;split to 12 month  is a array size (12,yearnum)
   w_12 = split_12(w)
   v_12 = split_12(v)
;;;;   

   
   
;;;等到年平均的数据
    w_annual = month_to_annual(w, 1)
    v_annual = month_to_annual(v, 1)

    w_annual = (/dim_avg_n(w_12(0:4,:,:,:),0)/)
    v_annual = (/dim_avg_n(v_12(0:4,:,:,:),0)/)
;;;;;

;;;;
    w_rc = mse_trend(w_annual,startyear,endyear)

    w_rc = w_rc*-50.0

    v_rc = mse_trend(v_annual,startyear,endyear)
;;;;


;;;prob 

  prob_w = rc_test(w_rc)
  ; prob_v = rc_test(v_rc)
;;;
wks = gsn_open_wks("pdf", "w-v-trend-latitude-section-150w-80w")

;;;;draw plot 
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@gsnLeftString= "Trends of JFMAM vertical and meridional wind averaged over 150w-80w"
    res@gsnLeftStringFontHeightF  = 0.015

    res@lbLabelStride = 1
    res@lbLabelFontHeightF = 0.015
   
    res@gsnMaximize=False
    ; res@lbLabelStrings=ispan(-6,6,1)
    ; res@cnFillPalette  = "BlueDarkRed18"
    res@cnFillPalette  = "BlueWhiteOrangeRed"
    res@cnLinesOn=False

    ; res@cnLevelSelectionMode ="ManualLevels"
    ; res@cnMinLevelValF = -1
    ; res@cnMaxLevelValF = 1
    ; res@cnLevelSpacingF = 0.2

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-1.0,-0.8,-0.6,-0.4,-0.2,0.2,0.4,0.6,0.8,1.0/)

    res@cnFillOn = True 
    ; res@gsnYRefLine=500
    ; res@gsnYRefLineColor="black"

    res@gsnContourNegLineDashPattern=1
    res@gsnContourZeroLineThicknessF=2.5

    res@cnLineLabelsOn=False
    res@cnMonoFillScale=True
    res@cnInfoLabelOn=False

    ; res@cnExplicitLineLabelsOn=True
    ; res@cnLineLabelStrings=""+ispan(-2,6,1)+""
    ; res@cnLabelScalingMode="ScaleFactor"
    ; res@cnLabelScaleValueF=10^9

    res@lbLabelBarOn=True
    ; res@cnLineColors=(/0,0,0,0,0,0,18,18,18,18,18,18/)
    res@tmYRMode="Automatic"
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
    res@tmXBValues=ispan(-30,30,10)
    res@tmXBLabels=(/"30S","20S","10S","0","10N","20N","30N"/) ;;;~S~o~N~W

    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(-30,30,5)
     
    res@tiYAxisString=""

    res@vcRefMagnitudeF = 0.6                ; define vector ref mag
    res@vcRefLengthF    = 0.03        ; define length of vec ref
    res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    res@vcMinDistanceF  = 0.02       ; thin out vectors
    res@vcMapDirection  = False

    ; res@vcFillArrowHeadXF=0.1
    res@vcFillArrowMinFracWidthF=1.0
    res@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadXF=1.0
    res@vcFillArrowWidthF=0.05
    res@vcFillArrowHeadYF=0.15

    res@vcRefAnnoOrthogonalPosF  = -1.14

    res@vcRefAnnoPerimOn=False
    res@vcRefAnnoString2On=False


    res@vpHeightF=0.3
    res@vpWidthF=0.7

    plot=gsn_csm_pres_hgt_vector(wks,w_rc,v_rc,w_rc,res)


    ;;;打点
    res_pro=True
    res_pro@gsnDraw=False
    res_pro@gsnFrame=False
    res_pro@cnFillOn = True
    res_pro@cnLinesOn = False
     res_pro@cnLineLabelsOn=False

    res_pro@cnLevelSelectionMode = "ExplicitLevels"
    res_pro@cnLevels = (/0,0.1/)
    res_pro@lbLabelBarOn = False

    ; res_pro@cnLinesOn=False
    ; res_pro@cnMonoLineColor =False
    ; res_pro@cnLineColors = ((/"Transparent","White","Transparent"/))

    res_pro@cnMonoFillPattern = True
    res_pro@cnFillPattern = 17
    res_pro@cnFillDotSizeF = 0.003
    res_pro@cnMonoFillColor = False
    res_pro@cnFillColors = (/"Transparent","black","Transparent"/)

   
    res_pro@tmYRMode = "Automatic"
    res_pro@cnInfoLabelOn=False

    iplot  = gsn_csm_contour(wks, prob_w, res_pro)
    ;;;  
    overlay(plot, iplot)

    ; ;;;画气候态的等值线
    ; res@gsnRightString=""
    ; res@gsnLeftString=""
    ; res@cnLevelSelectionMode = "ExplicitLevels"
    ; res@cnLinesOn=True
    ; res@cnLevels = ispan(-8,8,2)
    ; res@cnLineLabelsOn=True
    ; res@cnLineLabelBackgroundColor="Transparent"
    ; res@cnLineThicknessF = 2.
    ; res@cnFillOn = False
    ; res@gsnContourNegLineDashPattern=1
    ; res@gsnContourZeroLineThicknessF=2.2

    ; plot1(k) = gsn_csm_pres_hgt(wks, msf_mean_clm, res)
  
 
    ;;;  
    ; overlay(plot(k), plot1(k))
    draw(plot)
    frame(wks)

; pres = True
; pres@gsnMaximize = True
; pres@gsnPanelLabelBar = True
; pres@pmLabelBarHeightF = 0.05
; pres@pmLabelBarWidthF  = 0.6
; ; pres@pmLabelBarOrthogonalPosF = 0.5
; pres@pmLabelBarParallelPosF = 0.02
; ; pres@lbOrientation = "Vertical"
; ; pres@lbLabelFontHeightF = 0.01
; gsn_panel(wks, plot, (/6,5/), pres)
end 



