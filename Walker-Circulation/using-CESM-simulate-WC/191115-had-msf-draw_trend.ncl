 ;;;;;
function get_had(v_a)
begin
    a = 6371000
    g = 9.8
    pi = get_pi("double")
    dfi = pi
    cosfi = cos(v_a&latitude/180.0)
    print(cosfi)
    
    cosfi_new = conform_dims(dimsizes(vint), cosfi, 2)


   ;;; vibeta(p, x, linlog, psfc, pbot, ptop)积分纬向平均的结果，得到mass
   vint_beta = v_a
   psfc = v_a(:,1,:)
   psfc = 1013.0
   
    do i = 0,dimsizes(level) -1
      vint_beta(:,i,:)  = (/vibeta( level, v_a(time|:,latitude|:,level|:), 1, psfc, level(i), 0)/)
    end do

    vint_beta_fi = v_a
    vint_beta_fi := (/((2*a*dfi*cosfi_new)/g)*vint_beta*100/)
    vint_beta_fi = (/vint_beta_fi/10^10/)
    copy_VarCoords(v_a, vint_beta_fi)

return vint_beta_fi
end 

function get_trend(fi_annual)
begin
;;;trend 
  rc1 = regCoef_n(ispan(1, dimsizes(fi_annual(:,1,1)), 1), fi_annual, 0, 0)
   
  rc1 = (/rc1*10/)  ;;;convert to decadal trend mean

  copy_VarCoords(fi(0,{1000:100},{-30:30}), rc1)
;;;
return rc1
end 
;;;;;;
;;*****main program 
;;
;;;;;;
  begin 
	; fils1  =  "/home/ys17-19/lsh/data/era/era_pre_global.nc"
    fils1  =  "/home/ys17-19/lsh/data/wc-result/1850-PDAY-rmpacific_197912-201212.nc"
    
	f1  =  addfile(fils1, "r")
	timestart = 197901
	timeend = 201612
	date  =  cd_calendar(f1->time, 1)
	timeind  =  ind(date.le.timeend.and.date.ge.timestart)

	startyear_r = str_split_by_length(timestart, 4)  ;; 将时间string分解成年和月
	endyear_r = str_split_by_length(timeend, 4 )

	startyear = stringtoint(startyear_r(0))
	endyear = stringtoint(endyear_r(0))

	v  =  f1->V(timeind,::-1,:,:)

	level = f1->level(::-1)
    print(level)
   
    ;; v avg 
    v_a = dim_avg_n_Wrap(v , 3) ;;;time,level,lat
    ;;
   
    a = 6371000
    g = 9.8
    pi = get_pi("double")
    dfi = pi
    cosfi = cos(v_a&latitude/180.0)
    print(cosfi)
    
    cosfi_new = conform_dims(dimsizes(vint), cosfi, 2)
;;;

;;; vibeta(p, x, linlog, psfc, pbot, ptop)积分
   
   vint_beta = v_a
   psfc = v_a(:,1,:)
   psfc = 1013.0
   
   do i = 0,dimsizes(level) -1
    vint_beta(:,i,:)  = (/vibeta( level, v_a(time|:,latitude|:,level|:), 1, psfc, level(i), 0)/)
  end do

    vint_beta_fi = v_a
    vint_beta_fi := (/((2*a*dfi*cosfi_new)/g)*vint_beta*100/)
    vint_beta_fi = (/vint_beta_fi/10^10/)
   copy_VarCoords(v_a, vint_beta_fi)
   vint_beta_clm = dim_avg_n_Wrap(vint_beta_fi(:,{1000:100},{-30:30}), 0)


    vint_beta_fi_win = month_to_season(vint_beta_fi(:,{1000:100},{-30:30}), "DJF")
    vint_beta_fi_win_clm = dim_avg_n_Wrap(vint_beta_fi_win, 0)
;;;;;;



;;;




;;;;;;
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlRe"

    res@cnLinesOn=True
 
    res@cnFillOn =True
   
    res@cnLinesOn=False
    ; res@cnLevelSelectionMode ="ManualLevels"
    ; res@cnMinLevelValF = -10
    ; res@cnMaxLevelValF = 10
    ; res@cnLevelSpacingF = 2
    res@lbLabelBarOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    res@lbLabelBarOn = True
    
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
  
    ; res@tmXBValues=(/0,90,180,270,357.5/)
    ; res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(-90,90,10)

    ; res@tmXBValues=ispan(100,310,40)
    ; res@tmXBLabels=(/"100~S~o~N~E","140~S~o~N~E","180~S~o~N~E","140~S~o~N~W","100~S~o~N~W","60~S~o~N~W"/)
    ; res@tmXBMinorOn=True
    ; res@tmXBMinorValues=ispan(100,310,10)
    res@tmYRMode="Automatic"
     
    res@vpHeightF=0.3
    res@vpWidthF=0.7


    ; res@vcRefMagnitudeF = 0.2               ; define vector ref mag
    ; res@vcRefLengthF    = 0.04       ; define length of vec ref
    ; res@vcGlyphStyle    = "FillArrow"     ; turn on curly vectors
    ; res@vcMinDistanceF  = 0.03      ; thin out vectors
    ; res@vcMapDirection  = False
    ; res@vcLineArrowThicknessF=2.0
    ; ; res@vcFillArrowHeadXF=0.1
    ; res@vcFillArrowMinFracWidthF=1.0
    ; res@vcFillArrowHeadMinFracYF=1.0
    ; ; res@vcFillArrowHeadXF=1.0
    ; res@vcFillArrowWidthF=0.05
    ; res@vcFillArrowHeadYF=0.15
    
    ; res@vcLineArrowHeadMaxSizeF=0.005
    ; res@vcLineArrowHeadMinSizeF=0.005
   

    ;;;;


   
    ; res@cnFillPalette = "MPL_seismic"
     wks2 = gsn_open_wks("pdf", "积分得到的哈德莱环流趋势")
    
    plot_diff = new(2,graphic)
    res@gsnLeftString = "Hadley circulation"
    plot_diff(0) = gsn_csm_pres_hgt(wks2,rc1 ,res)
    plot_diff(1) = gsn_csm_pres_hgt(wks2,rc2,res)
     ; draw(plot_diff(0))
     ; frame(wks2)
    res2 =res
    res2@gsnRightString=""
    res2@gsnLeftString=""
    res2@cnLevelSelectionMode = "ExplicitLevels"
    res2@cnLinesOn=True
    res2@cnLevels = ispan(-10,10,2)
    res2@cnLineLabelsOn=True
    res2@cnLineLabelBackgroundColor="Transparent"
    res2@cnLineThicknessF = 2.
    res2@cnFillOn = False
    res2@gsnContourNegLineDashPattern=1
    res2@gsnContourZeroLineThicknessF=4.0
    
    plot_clm = new(2,graphic)
    plot_clm(0) = gsn_csm_pres_hgt(wks2, fi_clm, res2)
    plot_clm(1) = gsn_csm_pres_hgt(wks2, vint_beta_clm, res2)


    do i = 0,1
    overlay(plot_diff(i), plot_clm(i))
    end do

   pres = True
   gsn_panel(wks2, plot_diff, (/2,1/), pres)


    end 