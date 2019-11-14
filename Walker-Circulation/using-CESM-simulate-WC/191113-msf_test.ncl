 ;;;;;
  begin 
	fils1  =  "/home/ys17-19/lsh/data/rea_msf/msf_time_era_global.nc"
	f1  =  addfile(fils1, "r")
	timestart = 197901
	timeend = 201212
	date  =  cd_calendar(f1->time, 1)
	timeind  =  ind(date.le.timeend.and.date.ge.timestart)

	startyear_r = str_split_by_length(timestart, 4)  ;; 将时间string分解成年和月
	endyear_r = str_split_by_length(timeend, 4 )

	startyear = stringtoint(startyear_r(0))
	endyear = stringtoint(endyear_r(0))

	ud  =  f1->ud_time(:,:,timeind)
    ud_new = ud(time|:,level|:,lon|:)


    
	level = f1->level
    print(level)
    msf_time = f1->msf_time
    msf_clm = dim_avg_n_Wrap(msf_time({1000:100},:,:), 2)
	printVarSummary(ud_new)
    vint = ud_new
    dp_or = dpres_plevel(level, 1013, 0, 0)  ; dp(30)
   
   
    dp = level

    level2 = level
    level2(0:dimsizes(level)-2) = level(1:)
    level2(dimsizes(level)-1) = 0
    dp = level - level2

    ud_2 = ud_new

    ud_2(:,0:dimsizes(level)-2,:) = ud_new(:,1:,:)
    ud_2(:,dimsizes(ud_new&level)-1,:) = 0
    ud_a = (/ud_new + ud_2/)*0.5

    ddp  = dpres_plevel(level(15:), 1000, 0, 0)  ; dp(30)
 
 ;;; 积分
      
	do i=0,dimsizes(level)-2
        ; dp  = dpres_plevel(level(i:), max(level(i)), 0, 0)  ; dp(30)
        ; dp = level(i:)
        ; level2 = level(i:)
        ; level2(0:dimsizes(level(i:))-2) = level(i+1:)
        ; level2(dimsizes(level(i:))-1) = 0
        ; dp = level(i:) - level2
        ; printVarSummary(ud_new(:,i+1:,:))
        ; printVarSummary(dp)
        dp_new = conform_dims(dimsizes(ud_a(:,i:,:)), dp(i:dimsizes(level2)-1), 1)
        print(dp)
        ; exit()
        printVarSummary(dp_new)
        vint(:,i,:)= dim_sum_n(ud_new(:,i:,:)*dp_new,1)
        ; delete(dp)
        delete(dp_new)
        ; delete(level2)
	end do 
    ; vint(:,dimsizes(level)-1,:) = 100*ud_new(:,dimsizes(level)-1,:)
    
    a = 6371000
    g = 9.8
    pi = get_pi("double")
    dfi = (pi*10)/180

    fi = ud_new
	fi := (/((a*dfi)/g)*vint*100/)
    copy_VarCoords(ud_new, fi)
	delete(vint)
    printVarSummary(fi)
    fi_clm = dim_avg_n_Wrap(fi(:,{1000:100},:), 0)
    printVarSummary(fi_clm)
    fi_clm = (/fi_clm/10^9/)
    msf_clm = (/msf_clm/10^9*0.98/)


;;;;;;/


;;;;;;
    res = True
    res@gsnDraw=False
    res@gsnFrame=False
    res@cnFillPalette = "BlueDarkRed18"

    res@cnLinesOn=True
 
    res@cnFillOn =True
   
    res@cnLinesOn=False
    res@cnLevelSelectionMode ="ManualLevels"
    res@cnMinLevelValF = -6
    res@cnMaxLevelValF = 6
    res@cnLevelSpacingF = 1
    res@lbLabelBarOn =False
    res@cnLineLabelsOn = False
    res@cnInfoLabelOn = False
    res@lbLabelBarOn = True
    
    res@tiYAxisString=""
    res@tmXTOn=False
    res@tmYROn=False
    res@tmXBMode="Explicit"
  
    res@tmXBValues=(/0,90,180,270,357.5/)
    res@tmXBLabels=(/"0","90~S~o~N~E","180","90~S~o~N~W","0"/)
    res@tmXBMinorOn=True
    res@tmXBMinorValues=ispan(0,360,10)

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
    res@vcRefMagnitudeF = 0.2              ; define vector ref mag
    res@vcRefLengthF    = 0.05        ; define length of vec ref

    res@vcGlyphStyle    = "CurlyVector"     ; turn on curly vectors
    res@vcMinDistanceF  = 0.03       ; thin out vectors
    res@vcMapDirection  = False
    res@vcLineArrowThicknessF=2.0
    ; res@vcFillArrowHeadXF=0.1
    res@vcFillArrowMinFracWidthF=1.0
    res@vcFillArrowHeadMinFracYF=1.0
    ; res@vcFillArrowHeadXF=1.0
    res@vcFillArrowWidthF=0.05
    res@vcFillArrowHeadYF=0.15
    
    res@vcLineArrowHeadMaxSizeF=0.005
    res@vcLineArrowHeadMinSizeF=0.005

    res@vcVectorDrawOrder="PostDraw"

    res@vcRefAnnoOrthogonalPosF  = -1.14

    res@vcRefAnnoPerimOn=False
    res@vcRefAnnoString2On=False

    ;;;;


   
    ; res@cnFillPalette = "MPL_seismic"
     wks2 = gsn_open_wks("pdf", "积分得到的msf")
    
    plot_diff = new(2,graphic)
    res@gsnLeftString = "Control run"
    plot_diff(0) = gsn_csm_pres_hgt(wks2,fi_clm,res)
    plot_diff(1) = gsn_csm_pres_hgt(wks2,msf_clm,res)
     
   pres = True
   gsn_panel(wks2, plot_diff, (/2,1/), pres)


    end 